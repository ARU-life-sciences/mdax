//! `irx`: fast, windowed inverted-repeat (IR) candidate discovery for assemblies.
//!
//! This binary wraps `mdax`'s foldback (palindrome / inverted-repeat) detector and
//! applies it to *assemblies* by scanning each contig in overlapping windows.
//!
//! Conceptually, for each window we:
//! 1) run a coarse minimizer-geometry foldback detector (`detect_foldbacks_topk`)
//! 2) refine breakpoint position with a lightweight alignment-based step
//! 3) derive left/right arm bounds from supporting matchpoints (`best_pts` / `hit.pts`)
//! 4) deduplicate candidates (best-first) to avoid window overlap inflation
//! 5) emit BED-like TSV rows describing candidate IR intervals
//!
//! ## Why windows?
//! Assemblies can be huge; scanning the whole contig at once can blur multiple
//! signals and can be expensive. Windowing lets you keep detection local and
//! makes the runtime roughly linear in contig length.
//!
//! ## Why two-stage dedup?
//! Adjacent/overlapping windows often rediscover the same IR. We use:
//! - Stage A: very cheap bin-level dedup on `(break_bin, diag_bin)`
//! - Stage B: greedy interval overlap dedup on quantised intervals
//!   (best-first so we keep the strongest representative)
//!
//! ## Parallelism strategy
//! We parallelise *windows within a contig* using Rayon. This keeps the contig’s
//! sequence data shared (borrowed slice) and avoids extra IO complexity.
//!
//! Important: we use a **local Rayon thread pool** (not the global pool) when
//! `--threads > 0`, so the thread count is deterministic and not affected by
//! other crates initialising Rayon.
//!
//! ## Output format
//! We write a TSV with a BED6-like prefix plus additional columns:
//! - contig, start, end, name, score, strand, break_pos, identity_est, matches,
//!   span, la0, la1, ra0, ra1, win_start, win_end, kept_pts, bin
//!
//! Coordinates are 0-based, half-open (BED convention).
//!
//! ## Notes / limitations
//! - This is a *candidate finder*: it will return imperfect hits, especially in
//!   repeat-rich or low-complexity regions.
//! - The coarse detector uses minimizer buckets with a cap to avoid quadratic
//!   blowups; highly repetitive k-mers may be ignored.
//! - The arm bounds are derived from minimizer matchpoints and are approximate;
//!   you can tighten/validate later with a dedicated alignment pass if needed.

use anyhow::Result;
use clap::{Arg, Command, value_parser};
use rayon::ThreadPoolBuilder;
use std::cell::RefCell;
use std::collections::HashSet;
use std::{io::Write, path::PathBuf};
use thread_local::ThreadLocal;

use mdax::{
    cfg::{FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg, SharedCfg, SigCfg},
    foldback, io,
    scratch::FoldScratch,
    utils::RefineMode,
};

/// Build the CLI argument parser.
///
/// We keep knobs intentionally minimal:
/// - window geometry (`--window-len`, `--step`)
/// - number of candidates per window (`--hits-per-window`)
/// - a simple identity filter (`--min-ident`)
/// - dedup quantisation (`--dedup-bp`)
/// - thread count (`--threads`)
///
/// If you later want `irx` to mirror `mdax` options, you can expose fields from
/// `SharedCfg` / `FoldOnlyCfg` / `RefineCfg` as additional flags.
fn build_cli() -> clap::ArgMatches {
    Command::new("irx")
        .version(clap::crate_version!())
        .about("Detect candidate inverted repeats in assemblies using mdax foldback logic")
        .long_about(
    "irx scans an assembly FASTA (optionally gzipped) for inverted-repeat (IR) \
candidates using mdax's foldback detector.\n\n\
Workflow per contig:\n\
  - generate overlapping windows\n\
  - in each window, find up to K candidate foldback bins (top-k anti-diagonals)\n\
  - refine breakpoint position for each candidate\n\
  - derive left/right arm bounds from minimizer matchpoints\n\
  - best-first deduplicate across windows\n\
  - emit BED-like TSV rows\n\n\
Output is 0-based, half-open coordinates (BED convention).\n\n\
TSV columns (header line begins with '#'):\n\
  #contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\tmatches\tspan\tla0\tla1\tra0\tra1\twin_start\twin_end\tkept_pts\tbin\n\n\
Column meanings:\n\
  contig        First token of FASTA record id\n\
  start,end     Candidate interval covering both arms (contig coords)\n\
  name         Constant 'IR'\n\
  score        identity_est scaled to 0..1000 (rounded)\n\
  strand       '.' (IR is not strand-specific here)\n\
  break_pos    Refined split position in contig coords\n\
  identity_est Refinement identity estimate (0..1)\n\
  matches      Coarse supporting minimizer matches (fb.matches)\n\
  span         Coarse span between arms (fb.span, bp)\n\
  la0,la1      Left arm bounds (contig coords)\n\
  ra0,ra1      Right arm bounds (contig coords)\n\
  win_start,end Window bounds that produced this hit (contig coords)\n\
  kept_pts     #matchpoints retained near refined anti-diagonal (0 if unfiltered fallback)\n\
  bin          Anti-diagonal bin id (debug / dedup key component)\n\n\
",
)
        .arg(
            Arg::new("input")
                .help("Input assembly FASTA(.gz)")
                .long_help(
                    "Input assembly in FASTA format. If the filename ends with .gz, \
irx will transparently decompress it.\n\n\
Note: Needletail may include the full FASTA header in the record id; irx will \
trim to the first whitespace-delimited token for output.",
                )
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("bed")
                .help("Output BED/TSV file ('-' for stdout)")
                .long_help(
                    "Output path for the BED-like TSV. Use '-' to write to stdout.\n\n\
The first line is a header prefixed with '#'.",
                )
                .long("bed")
                .short('o')
                .default_value("-")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("window_len")
                .help("Window length for synthetic reads (0 = whole contig)")
                .long_help(
                    "Length of the sliding window used to scan contigs.\n\n\
- If 0, the entire contig is processed as a single window.\n\
- Otherwise, windows are [start, start+window_len) clamped to contig length.\n\n\
Smaller windows can increase sensitivity to local IRs but may increase runtime.",
                )
                .long("window-len")
                .default_value("250000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("step")
                .help("Step size between windows")
                .long_help(
                    "Step size between consecutive window starts.\n\n\
Smaller step => more overlap => higher sensitivity but more duplicate discoveries.\n\
Deduplication is applied downstream, but runtime still increases with overlap.",
                )
                .long("step")
                .default_value("50000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("hits_per_window")
                .help("How many foldback candidates to emit per window (top-k bins)")
                .long_help(
                    "Maximum number of coarse foldback candidates per window.\n\n\
Internally, mdax bins minimizer matchpoints by quantised anti-diagonal \
(d = p1 + p2). `--hits-per-window` selects the top-K bins by a cheap proxy rank \
(count + span).\n\n\
Higher values may find multiple IRs within the same window, at increased cost.",
                )
                .long("hits-per-window")
                .default_value("4")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_ident")
                .help("Minimum identity of the left/right arms of the putative inverted repeat.")
                .long_help(
                    "Minimum identity estimate (0..1) from the refinement step.\n\n\
This is a quick, approximate similarity measure between the two arms near the \
refined split. It is *not* a full alignment of the entire arm intervals.\n\n\
Typical values:\n\
  - 0.20: permissive, good for candidate discovery\n\
  - 0.50+: stricter, fewer false positives but may miss diverged IRs",
                )
                .long("min-ident")
                .default_value("0.2")
                .value_parser(value_parser!(f32)),
        )
        .arg(
            Arg::new("dedup_bp")
                .help("Quantise starts/ends for faster dedup (bp). Used by the interval-overlap dedupper.")
                .long_help(
                    "Quantisation size (in bp) used in deduplication.\n\n\
We do two-stage dedup:\n\
  A) breakpoint/bin dedup on (break_pos/q, diag_bin)\n\
  B) overlap dedup on intervals quantised to (start/q, end/q)\n\n\
Larger values merge more aggressively (fewer outputs) and speed up overlap checks.",
                )
                .long("dedup-bp")
                .default_value("200")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("threads")
                .help("Number of worker threads (0 = Rayon default)")
                .long_help(
                    "Number of worker threads.\n\n\
- 0: use Rayon default (typically = number of logical CPUs)\n\
- N>0: use exactly N threads via a *local* Rayon thread pool\n\n\
Using a local pool ensures the requested thread count is respected even if \
another crate initialises the Rayon global pool.",
                )
                .long("threads")
                .default_value("0")
                .value_parser(value_parser!(usize)),
        )
        .get_matches()
}

/// Generate (start,end) windows over a contig.
///
/// If `window_len == 0`, yields the whole contig once.
/// Otherwise yields overlapping windows advancing by `step`
/// until the end of the contig is reached.
fn windows(len: usize, window_len: usize, step: usize) -> impl Iterator<Item = (usize, usize)> {
    let mut start = 0usize;
    std::iter::from_fn(move || {
        if start >= len {
            return None;
        }
        if window_len == 0 || window_len >= len {
            start = len;
            return Some((0, len));
        }
        let end = (start + window_len).min(len);
        let out = (start, end);
        if end == len {
            start = len;
        } else {
            start += step.max(1);
        }
        Some(out)
    })
}

/// Needletail sometimes returns the whole FASTA header in `id`.
///
/// For BED-like output we typically want the first token up to whitespace.
fn first_token_utf8(bytes: &[u8]) -> &str {
    let end = bytes
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(bytes.len());
    std::str::from_utf8(&bytes[..end]).unwrap_or("<nonutf8>")
}

/// Derive arm bounds from coarse matchpoints (no filtering).
///
/// `best_pts` are (p1, p2) minimizer start positions in forward coordinates.
/// We map each point to:
/// - left  = min(p1,p2)
/// - right = max(p1,p2)
///
/// Then take min/max across all points, expand end coords by +k to cover k-mer
/// length, and clamp to `[0, win_len]`.
///
/// Returns (la0, la1, ra0, ra1) in *window-local* coords.
fn arm_bounds_from_best_pts(
    best_pts: &[(i32, i32)],
    k: usize,
    win_len: usize,
) -> Option<(usize, usize, usize, usize)> {
    if best_pts.is_empty() {
        return None;
    }

    let mut lmin = i32::MAX;
    let mut lmax = i32::MIN;
    let mut rmin = i32::MAX;
    let mut rmax = i32::MIN;

    for &(p1, p2) in best_pts {
        let left = p1.min(p2);
        let right = p1.max(p2);

        lmin = lmin.min(left);
        lmax = lmax.max(left);
        rmin = rmin.min(right);
        rmax = rmax.max(right);
    }

    let la0 = lmin.max(0) as usize;
    let la1 = (lmax.max(0) as usize).saturating_add(k).min(win_len);
    let ra0 = rmin.max(0) as usize;
    let ra1 = (rmax.max(0) as usize).saturating_add(k).min(win_len);

    if la0 < la1 && ra0 < ra1 {
        Some((la0, la1, ra0, ra1))
    } else {
        None
    }
}

/// Derive arm bounds but keep only points consistent with the refined split.
///
/// Key idea:
/// For a foldback, matchpoints satisfy roughly:
///     p1 + p2 ≈ 2 * split
///
/// We keep points where `| (p1+p2) - 2*split | <= tol`,
/// where `tol = fold_diag_tol * widen`.
///
/// Returns (la0, la1, ra0, ra1, kept_pts) in *window-local* coords.
fn arm_bounds_near_split_diag(
    best_pts: &[(i32, i32)],
    split: usize,
    k: usize,
    win_len: usize,
    fold_diag_tol: i32,
    widen: i32,
) -> Option<(usize, usize, usize, usize, usize)> {
    if best_pts.is_empty() {
        return None;
    }

    let split2 = 2 * (split as i32);
    let tol = fold_diag_tol.max(1) * widen.max(1);

    let mut lmin = i32::MAX;
    let mut lmax = i32::MIN;
    let mut rmin = i32::MAX;
    let mut rmax = i32::MIN;
    let mut kept = 0usize;

    for &(p1, p2) in best_pts {
        let d = p1 + p2;
        if (d - split2).abs() > tol {
            continue;
        }

        kept += 1;
        let left = p1.min(p2);
        let right = p1.max(p2);

        lmin = lmin.min(left);
        lmax = lmax.max(left);
        rmin = rmin.min(right);
        rmax = rmax.max(right);
    }

    if kept == 0 {
        return None;
    }

    let la0 = lmin.max(0) as usize;
    let la1 = (lmax.max(0) as usize).saturating_add(k).min(win_len);
    let ra0 = rmin.max(0) as usize;
    let ra1 = (rmax.max(0) as usize).saturating_add(k).min(win_len);

    if la0 < la1 && ra0 < ra1 {
        Some((la0, la1, ra0, ra1, kept))
    } else {
        None
    }
}

/// Pack two u32 values into a u64 (cheap HashSet key).
///
/// Used for Stage-A dedup keys: `(break_bin, diag_bin)`.
#[inline]
fn pack_u32_pair(a: u32, b: u32) -> u64 {
    ((a as u64) << 32) | (b as u64)
}

/// One candidate IR call, computed in parallel per window and then deduplicated
/// serially per contig.
///
/// The struct is designed to carry:
/// - coordinate fields needed for output
/// - inexpensive precomputed dedup keys
/// - a single “rank” for deterministic best-first selection
#[derive(Clone, Debug)]
struct Cand {
    // ---- global coordinates (contig space) ----
    start: usize,
    end: usize,
    break_pos: usize,
    la0: usize,
    la1: usize,
    ra0: usize,
    ra1: usize,
    win_start: usize,
    win_end: usize,

    // ---- quality / evidence ----
    ident: f32,
    score: i32,
    matches: usize,
    span: usize,
    kept_pts: usize,
    bin: i32,

    // ---- best-first dedup / stability ----
    rank: i64,
    key_a: u64, // Stage A key: (break_bin, diag_bin)
    qs: usize,  // Stage B quantised interval start
    qe: usize,  // Stage B quantised interval end
}

fn main() -> Result<()> {
    let args = build_cli();

    let input = args.get_one::<PathBuf>("input").unwrap();
    let bed = args.get_one::<PathBuf>("bed").unwrap();
    let window_len = *args.get_one::<usize>("window_len").unwrap();
    let step = *args.get_one::<usize>("step").unwrap();
    let hits_per_window = *args.get_one::<usize>("hits_per_window").unwrap();
    let dedup_bp = *args.get_one::<usize>("dedup_bp").unwrap();
    let min_ident = *args.get_one::<f32>("min_ident").unwrap();
    let threads = *args.get_one::<usize>("threads").unwrap();

    // Build a *local* Rayon pool if requested; otherwise, use the global/default.
    // This makes `--threads N` deterministic and avoids global-pool initialisation issues.
    let pool = if threads > 0 {
        Some(ThreadPoolBuilder::new().num_threads(threads).build()?)
    } else {
        None
    };

    // Minimal cfg: reuse mdax detector + refinement.
    //
    // Notes:
    // - `forward_only=true` is important for mdax's foldback geometry.
    // - `min_matches`, `min_arm`, and `fold_diag_tol` control coarse hit density.
    // - `refine` controls the refinement band/window and the maximum allowed edit-rate.
    let cfg = MdaxCfg {
        shared: SharedCfg {
            minimizer: MinimizerCfg {
                k: 17,
                w: 21,
                forward_only: true,
            },
            min_matches: 20,
            end_guard: 0,
            refine: RefineCfg {
                window: 200,
                arm: 500,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.25,
            },
            fold_diag_tol: 120,
        },
        fold: FoldOnlyCfg { min_arm: 2000 },
        fold2: FoldSecondPassCfg {
            min_support: 0,
            min_identity: 0.0,
            min_support_ident: 0.0,
        },
        sig: SigCfg {
            flank_bp: 600,
            take: 8,
            value_shift: 0,
        },
    };

    // Main run function. We place this in a closure so we can optionally
    // `pool.install(run)` when a local pool is used.
    let run = || -> Result<()> {
        let mut out: Box<dyn Write> = if bed.to_string_lossy() == "-" {
            Box::new(std::io::BufWriter::new(std::io::stdout()))
        } else {
            Box::new(std::io::BufWriter::new(std::fs::File::create(bed)?))
        };

        // TSV header
        match writeln!(
            out,
            "#contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\tmatches\tspan\tla0\tla1\tra0\tra1\twin_start\twin_end\tkept_pts\tbin"
        ) {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
            Err(e) => return Err(e.into()),
        }

        let mut r = io::open_fasta_reader(input)?;

        // Thread-local scratch: each Rayon worker reuses a FoldScratch to avoid
        // repeated allocations when scanning many windows.
        //
        // IMPORTANT: construct TLS inside `run` so it is associated with the
        // threads used by the chosen pool.
        let tls: ThreadLocal<RefCell<FoldScratch>> = ThreadLocal::new();

        // Process each contig sequentially (IO + per-contig dedup is simpler).
        // Windows *within a contig* are processed in parallel.
        while let Some(rec) = r.next() {
            let rec = rec?;
            let id = first_token_utf8(rec.id());
            let seq = rec.seq();
            let len = seq.len();

            let wins: Vec<(usize, usize)> = windows(len, window_len, step).collect();
            if wins.is_empty() {
                continue;
            }

            let q = dedup_bp.max(1);
            let mm_k = cfg.shared.minimizer.k;
            let fold_diag_tol = cfg.shared.fold_diag_tol;

            // For each window:
            // - detect top-k coarse foldback bins
            // - refine split position per hit
            // - compute approximate arm interval bounds
            // - package into `Cand` records
            //
            // We do *no* dedup here (that happens serially per contig).
            use rayon::prelude::*;

            let mut cands: Vec<Cand> = wins
                .par_iter()
                .flat_map_iter(|&(wstart, wend)| {
                    let sub = &seq[wstart..wend];
                    let win_len = sub.len();

                    // Thread-local scratch reuse.
                    let scratch_cell = tls.get_or(|| RefCell::new(FoldScratch::new()));
                    let mut scratch = scratch_cell.borrow_mut();

                    // Find up to K foldback candidates in this window.
                    let hits = foldback::detect_foldbacks_topk(
                        sub,
                        &cfg.shared,
                        &cfg.fold,
                        &mut scratch,
                        hits_per_window,
                    );

                    if hits.is_empty() {
                        return Vec::new().into_iter();
                    }

                    // Collect per-window candidates.
                    let mut out_local: Vec<Cand> = Vec::new();

                    for hit in hits {
                        if hit.bin_count < cfg.shared.min_matches {
                            continue;
                        }
                        if hit.bin_span < cfg.fold.min_arm {
                            continue;
                        }

                        let fb = hit.fb;
                        let bin = hit.bin;

                        // Refine split position around coarse split.
                        // We swallow errors per-hit to avoid aborting the whole contig.
                        let Some(rf) = foldback::refine_breakpoint(
                            sub,
                            fb.split_pos,
                            &cfg.shared,
                            &mut scratch.refine,
                        )
                        .ok()
                        .flatten() else {
                            continue;
                        };

                        // Filter by identity estimate.
                        if rf.identity_est < min_ident {
                            continue;
                        }

                        // Global break position (contig coords).
                        let break_pos = wstart + rf.split_pos;

                        // Determine arm bounds:
                        //   1) prefer points near refined split (anti-diagonal proximity)
                        //   2) fall back to unfiltered bounds
                        //   3) fall back to break±span/2
                        let (start, end, la0g, la1g, ra0g, ra1g, kept_pts) =
                            if let Some((la0, la1, ra0, ra1, kept)) = arm_bounds_near_split_diag(
                                &hit.pts,
                                rf.split_pos,
                                mm_k,
                                win_len,
                                fold_diag_tol,
                                20,
                            ) {
                                let la0g = wstart + la0;
                                let la1g = wstart + la1;
                                let ra0g = wstart + ra0;
                                let ra1g = wstart + ra1;
                                let start = la0g.min(ra0g);
                                let end = la1g.max(ra1g);
                                (start, end, la0g, la1g, ra0g, ra1g, kept)
                            } else if let Some((la0, la1, ra0, ra1)) =
                                arm_bounds_from_best_pts(&hit.pts, mm_k, win_len)
                            {
                                let la0g = wstart + la0;
                                let la1g = wstart + la1;
                                let ra0g = wstart + ra0;
                                let ra1g = wstart + ra1;
                                let start = la0g.min(ra0g);
                                let end = la1g.max(ra1g);
                                (start, end, la0g, la1g, ra0g, ra1g, 0)
                            } else {
                                let half = fb.span / 2;
                                let start = break_pos.saturating_sub(half);
                                let end = (break_pos + half).min(len);
                                (start, end, 0, 0, 0, 0, 0)
                            };

                        // Score is 0..1000 scale of identity estimate (BED-ish).
                        let score = (rf.identity_est.clamp(0.0, 1.0) * 1000.0).round() as i32;

                        // Stage A key: (break_bin, diag_bin) packed into u64.
                        // Including `diag_bin` makes this less aggressive (keeps distinct anti-diagonals).
                        let break_bin = (break_pos / q) as u32;
                        let diag_u32 = bin as u32;
                        let key_a = pack_u32_pair(break_bin, diag_u32);

                        // Quantised interval for Stage B overlap dedup.
                        let qs = (start / q) * q;
                        let qe = (end / q) * q;

                        // Deterministic rank for best-first dedup.
                        //
                        // You want rank to reflect what we trust most. A reasonable ordering is:
                        //  - refined identity (score)
                        //  - bin_rank from coarse detector (support+span proxy)
                        //  - raw minimizer matches as a final tiebreaker
                        //
                        // Packing into i64 keeps sort cheap.
                        let rank = ((score as i64) << 48)
                            | (((hit.bin_rank as i64) & 0xFFFF_FFFF) << 16)
                            | ((fb.matches.min(65535) as i64) << 0);

                        out_local.push(Cand {
                            start,
                            end,
                            break_pos,
                            la0: la0g,
                            la1: la1g,
                            ra0: ra0g,
                            ra1: ra1g,
                            win_start: wstart,
                            win_end: wend,
                            ident: rf.identity_est,
                            score,
                            matches: fb.matches,
                            span: fb.span,
                            kept_pts,
                            bin,
                            rank,
                            key_a,
                            qs,
                            qe,
                        });
                    }

                    out_local.into_iter()
                })
                .collect();

            if cands.is_empty() {
                continue;
            }

            // back to the main thread
            // Sort best-first so we keep the strongest representative when windows overlap.
            cands.sort_by(|a, b| {
                b.rank
                    .cmp(&a.rank)
                    .then_with(|| a.break_pos.cmp(&b.break_pos))
            });

            // Stage A: cheap dedup on break/bin.
            let mut seen_break_bins: HashSet<u64> = HashSet::new();
            // Stage B: greedy overlap dedup on quantised intervals.
            let mut accepted_intervals: Vec<(usize, usize)> = Vec::new();

            for c in cands {
                // A: check for duplicate break/bin key. This is very fast and removes exact duplicates and many near-duplicates.
                if !seen_break_bins.insert(c.key_a) {
                    continue;
                }

                // B: check for overlap with accepted intervals. This is O(N^2) in the worst case but should be fine for small N.
                let mut overlap = false;
                for &(s, e) in &accepted_intervals {
                    if !(c.qe <= s || c.qs >= e) {
                        overlap = true;
                        break;
                    }
                }
                if overlap {
                    continue;
                }
                accepted_intervals.push((c.qs, c.qe));

                // emit TSV
                match writeln!(
                    out,
                    "{id}\t{}\t{}\tIR\t{}\t.\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    c.start,
                    c.end,
                    c.score,
                    c.break_pos,
                    c.ident,
                    c.matches,
                    c.span,
                    c.la0,
                    c.la1,
                    c.ra0,
                    c.ra1,
                    c.win_start,
                    c.win_end,
                    c.kept_pts,
                    c.bin
                ) {
                    Ok(_) => {}
                    Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
                    Err(e) => return Err(e.into()),
                }
            }
        }

        Ok(())
    };

    // Run either inside the local pool or directly.
    if let Some(pool) = pool {
        pool.install(run)
    } else {
        run()
    }
}
