//! `irx`: fast, windowed inverted-repeat (IR) candidate discovery for assemblies.
//!
//! This binary wraps `mdax` foldback (palindrome / inverted-repeat) detector and
//! applies it to genome *assemblies* by scanning each contig in overlapping windows.
//!
//! For each window we:
//! 1) run a coarse minimizer-geometry foldback detector (`mdax::foldback::detect_foldbacks_topk`)
//! 2) refine breakpoint position with a lightweight alignment-based step
//! 3) derive left/right arm bounds from supporting matchpoints (`best_pts` / `hit.pts`)
//! 4) deduplicate candidates (best-first) to avoid window overlap inflation
//! 5) emit BED-like TSV rows describing candidate IR intervals
//!
//! ## Why windows?
//! Assemblies can be huge; scanning the whole contig at once can blur multiple
//! signals and can be expensive. Windowing lets us keep detection local and
//! makes the runtime roughly linear in contig length.
//!
//! ## Why two-stage dedup?
//! Adjacent/overlapping windows often rediscover the same IR. We use:
//! - Stage A: very cheap bin-level dedup on `(break_bin, diag_bin)`
//! - Stage B: drop a candidate whose two arms both overlap an accepted call's
//!   arms (the same IR found from another window); nested and overlapping but
//!   distinct IRs are kept
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
//!   user can tighten/validate later with a dedicated alignment pass if needed/wanted.
//! - Currently only IR's > 2kb are retained. Though working on a TODO to have a mode to
//!   return shorter repeats should the user desire. I suspect this will be more
//!   computationally expensive.

use anyhow::Result;
use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use clap::{Arg, Command, value_parser};
use mdax::elog;
use mdax::utils::revcomp_into;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use serde::Serialize;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
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
/// - dedup quantisation (`--dedup-bp`)
/// - thread count (`--threads`)
///
/// We may later want `irx` to mirror `mdax` options and expose fields from
/// `SharedCfg` / `FoldOnlyCfg` / `RefineCfg` as additional flags.
fn build_cli() -> clap::ArgMatches {
    Command::new("irx")
        .version(clap::crate_version!())
        .about("Detect candidate inverted repeats in assemblies using mdax foldback logic")
        .long_about(
    "`irx` scans an assembly FASTA(.gz) for inverted-repeat (IR) \
candidates using mdax's foldback detector.\n\n\
Workflow per contig:\n\
  - generate overlapping windows\n\
  - in each window, find up to K candidate foldback bins (top-k anti-diagonals)\n\
  - refine breakpoint position for each candidate\n\
  - derive left/right arm bounds from minimizer matchpoints\n\
  - best-first deduplicate across windows\n\
  - emit BED-like TSV rows (and optionally GFF3 via --gff)\n\n\
Output is 0-based, half-open coordinates (BED convention).\n\n\
TSV columns (header line begins with '#'):\n\n\
Column meanings:\n\
  contig         First token of FASTA record id\n\
  start,end      Candidate interval covering both arms (contig coords)\n\
  name           Constant 'IR'\n\
  score          identity_est scaled to 0..1000 (rounded)\n\
  strand         '.' (IR is not strand-specific here)\n\
  break_pos      Refined split position in contig coords\n\
  identity_est   Refinement identity estimate (0..1). Computed locally around breakpoint during refinement.\n\
  tir_ident      Estimated sequence identity of the IR arms (la0..la1 vs revcomp(ra0..ra1))\n\
  matches        Coarse supporting minimizer matches (fb.matches)\n\
  span           Coarse span between arms (fb.span, bp)\n\
  la0,la1        Left arm bounds (contig coords)\n\
  ra0,ra1        Right arm bounds (contig coords)\n\
  contig_len     Length of the contig/chromosome (bp)\n\
  win_start,end  Window bounds that produced this hit (contig coords)\n\
  kept_pts       Matchpoints supporting the reported arm bounds (the chosen arm cluster; 0 if unfiltered fallback)\n\
  bin            Anti-diagonal bin id (debug / dedup key component)\n\
  arm_len, spacer, ir_class  Shorter arm (bp), right-arm start minus left-arm end, class\n\
  dir_ident      Identity of the arms in DIRECT orientation (left vs right as-is)\n\
  direct_better  1 if dir_ident > tir_ident (arms match better as a direct repeat)\n\
  locus_id       IRs whose arms share sequence (overlap, transitively) form one locus;\n\
                 copies of a repeat paired in inverted orientation give many IRs per locus\n\
  locus_n_irs    Number of IRs in that locus\n\
  pass           local (--window-len) or large (--large-window-len) pass that found it\n\n\
",
)
        .arg(
            Arg::new("input")
                .help("Input assembly FASTA(.gz)")
                .long_help(
                    "Input assembly in FASTA format. If the filename ends with .gz, \
`irx` will transparently decompress it.\n\n\
Note: `irx` will trim to the first whitespace-delimited token for output.",
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
                .short('b')
                .default_value("-")
                .value_parser(value_parser!(PathBuf)),
        )
                .arg(
            Arg::new("fasta")
                .help("Optional FASTA output for accepted IR intervals ('-' for stdout).")
                .long_help("Write the sequences corresponding to accepted inverted-repeat (IR) \
intervals to a FASTA file.

Each emitted IR interval is extracted from the input assembly using the \
reported `start`–`end` coordinates and written as a FASTA record. The FASTA \
header contains the contig name and interval coordinates so sequences can be \
easily traced back to their source region.

If '-' is supplied as the filename, FASTA output is written to stdout. \
This is useful for piping IR sequences directly into downstream tools \
(e.g. alignment, repeat analysis, or clustering).

If this option is not provided, no FASTA sequences are written and only the \
BED-like interval table is produced.")
                .long("fasta")
                .short('f')
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("fasta_wrap")
                .help("Wrap FASTA sequence lines to this width (0 = no wrap).")
                .long_help("Set the line width used when writing FASTA sequences.

Sequences are wrapped to the specified number of characters per line in the \
output FASTA file. This improves readability and compatibility with tools \
that expect fixed-width FASTA formatting.

A value of 0 disables wrapping and writes each sequence on a single line.

The default value of 60 follows a common FASTA formatting convention.")
                .long("fasta-wrap")
                .default_value("60")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max_interval_bp")
                .help("Maximum length of an emitted IR interval (bp). Helps prevent window-scale inflation.")
                .long("max-interval-bp")
                .default_value("250000")
                .value_parser(value_parser!(usize)),
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
Internally, `mdax` bins minimizer matchpoints by quantised anti-diagonal \
(d = p1 + p2). `--hits-per-window` selects the top-K bins by a cheap proxy rank \
(count + span).\n\n\
Higher values find more IRs in IR-dense regions, at modest cost. The default (32) was \
raised from 4 after a stress test: with 6-12 IRs within one window, a cap of 4 lost \
up to half of them (and ~20% of genuine IRs in Castanea_sativa), for ~10% runtime.",
                )
                .long("hits-per-window")
                .default_value("32")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_arm")
                .help("Minimum length (bp) of each IR arm.")
                .long_help(
                    "Minimum length (bp) of *each* IR arm.\n\n\
Every emitted row satisfies arm_len >= --min-arm, where arm_len is the shorter of the \
two reported arms (full extent of the supporting matchpoints; see --arm-gap-bp). \
Minimizer matches cannot reach the last few bases of an arm, so an arm exactly at the \
threshold may measure a few bp short.\n\n\
The coarse detector's span gate measures the whole IR footprint (arm + spacer + arm), \
so on its own it does not bound arm length; the per-arm check is applied after arm \
bounds are derived and before deduplication.",
                )
                .long("min-arm")
                .default_value("2000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_matches")
                .help("Minimum minimizer matchpoints supporting a candidate.")
                .long_help(
                    "Minimum number of minimizer matchpoints in an anti-diagonal bin for it to \
be considered.\n\n\
Each matching k-mer pair contributes two mirrored points (left->right and right->left), \
so this is roughly twice the number of distinct matching minimizers.",
                )
                .long("min-matches")
                .default_value("20")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("dedup_bp")
                .help("Breakpoint quantisation (bp) for Stage-A dedup.")
                .long_help(
                    "Quantisation size (in bp) used in deduplication.\n\n\
We do two-stage dedup:\n\
  A) breakpoint/bin dedup on (break_pos/q, diag_bin)\n\
  B) drop a candidate whose two arms both overlap an accepted call's arms (the same IR\n\
     rediscovered from another window); nested or merely overlapping IRs are kept\n\n\
Larger values merge more aggressively in Stage A.",
                )
                .long("dedup-bp")
                .default_value("200")
                .value_parser(value_parser!(usize)),
        )
        .arg(
    Arg::new("diag_widen")
        .help("Anti-diagonal tolerance widen factor used in split-consistent arm bounds.")
        .long_help(
            "When deriving arm bounds near the refined split, we keep matchpoints where:\n\
  |(p1+p2) - 2*split| <= fold_diag_tol * diag_widen\n\n\
Smaller values are stricter (tighter arms, fewer cap hits).\n\
Larger values are more permissive (can inflate intervals).",
        )
        .long("diag-widen")
        .default_value("4")
        .value_parser(value_parser!(i32)),
        )
        .arg(
            Arg::new("large_window_len")
                .help("Window length (bp) of the second, large-IR pass (0 disables it).")
                .long_help(
                    "Second pass for IRs too large for the local windows.\n\n\
The local pass (--window-len/--step) only sees an IR whole if its footprint (left-arm \
start to right-arm end) is <= window-len - step (200 kb by default); recall falls to 0 at \
window-len. The large pass scans --large-window-len windows every --large-step bp, \
keeps only IRs with footprint > window-len - step, and merges them through the same \
dedup; the TSV/GFF 'pass' column says which pass found each IR. Its coarse footprint \
gate is raised so small IRs do not use its top-K slots. Complete for footprints up to \
--large-window-len - --large-step (1 Mb by default).",
                )
                .long("large-window-len")
                .default_value("1200000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("large_step")
                .help("Step (bp) between large-pass windows.")
                .long("large-step")
                .default_value("200000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("large_max_interval_bp")
                .help("Maximum emitted IR interval (bp) in the large pass.")
                .long("large-max-interval-bp")
                .default_value("1000000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_tir_ident")
                .help("Minimum arm identity (tir_ident) for a call to be reported.")
                .long_help(
                    "Minimum arm identity (tir_ident: left arm vs reverse-complemented right arm).\n\n\
Calls below it are rejected before deduplication, so they cannot suppress a genuine \
overlapping IR. Unrelated sequence scores ~0.37 on this measure, and the smallest \
identity irx can detect at --min-arm 2000 is ~0.8; on real genomes calls are bimodal with \
nothing between 0.45 and 0.60, the low mode being matches between scattered repeat \
copies that share an anti-diagonal bin without being homologous arms. 0 disables.",
                )
                .long("min-tir-ident")
                .default_value("0.6")
                .value_parser(value_parser!(f32)),
        )
        .arg(
            Arg::new("arm_gap_bp")
                .help("Minimum gap (bp) between matchpoints that splits an arm.")
                .long_help(
                    "Arm bounds come from gap-based clustering of the matchpoints along the arm:\n\
points are split into clusters wherever consecutive ones are more than\n\
max(--arm-gap-bp, 10 x median gap) apart; the cluster spanning the most arm bp is kept\n\
and reported at its full extent. This keeps whole arms of clean IRs while excluding\n\
stray matches from other repeat copies on the same anti-diagonal.",
                )
                .long("arm-gap-bp")
                .default_value("1000")
                .value_parser(value_parser!(i32)),
        )
        .arg(
            Arg::new("classify_ir")
                .help("Disable IR class columns (arm_len, spacer, ir_class) in TSV.")
                .long("classify-ir")
                // on by default
                .action(clap::ArgAction::SetFalse),
        )
        .arg(
            Arg::new("immediate_bp")
                .help("Max spacer (bp) to call an IR 'immediate'.")
                .long("immediate-bp")
                .default_value("2000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("wide_bp")
                .help("Min spacer (bp) to call an IR 'wide'.")
                .long("wide-bp")
                .default_value("100000")
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
        .arg(
    Arg::new("cap_policy")
        .help("What to do when an inferred interval exceeds --max-interval-bp.")
        .long_help(
            "When the inferred IR interval length exceeds --max-interval-bp:\n\
  clamp     Clamp to the cap\
  drop      Discard the candidate entirely\n\
  penalize  Clamp, but strongly reduce its rank so non-capped hits win dedup\n\n\
Recommended: penalize (prevents huge piles of exactly-cap intervals).",
        )
        .long("cap-policy")
        .default_value("penalize")
        .value_parser(["clamp", "drop", "penalize"]),
    )
    .arg(
    Arg::new("html")
        .help("Optional self-contained HTML dashboard report.")
        .long_help(
            "Write an interactive HTML dashboard containing embedded IR data \
and diagnostic plots.\n\n\
The report is self-contained except for any external JS libraries referenced \
by the template.",
        )
        .long("html")
        .value_parser(value_parser!(PathBuf)),
    )
    .arg(
    Arg::new("gff")
        .help("Optional GFF3 output of IR calls (use '-' for stdout).")
        .long_help(
            "Write IR calls as GFF3 (1-based, closed coordinates), sorted by contig \
(input order) then start, so the file can be bgzip'd and tabix-indexed directly.\n\n\
Each IR is three features:\n\
  inverted_repeat (SO:0000294)   left-arm start .. right-arm end, strand '.'\n\
  repeat_component (SO:0000840)  left arm, strand '+', Parent=<IR id>\n\
  repeat_component (SO:0000840)  right arm, strand '-', Parent=<IR id>\n\n\
The right arm's '-' strand records that it is the reverse complement of the left. \
The spacer is not a feature (it is negative for overlapping arms); it is the \
parent's `spacer` attribute, with `ir_class`. Score column = tir_ident ('.' if \
unavailable). `break_pos` is the 1-based position of the first base after the split.\n\n\
The header lists a ##sequence-region for every contig scanned, including contigs \
with no calls.",
        )
        .long("gff")
        .value_parser(value_parser!(PathBuf)),
    )
    .get_matches()
}

#[inline]
fn classify_ir_str(spacer: isize, immediate_bp: usize, wide_bp: usize) -> &'static str {
    if spacer <= 0 {
        "overlap"
    } else if spacer as usize <= immediate_bp {
        "immediate"
    } else if spacer as usize >= wide_bp {
        "wide"
    } else {
        "spaced"
    }
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

/// Identity of two arm sequences by block alignment: a fast proxy first, and
/// the traceback only for promising pairs (as for tir_ident since its
/// introduction). Used for both orientations (tir_ident and dir_ident).
fn arm_identity(a: &[u8], b: &[u8]) -> f32 {
    let proxy = arm_identity_proxy_blockalign(a, b).unwrap_or(f32::NAN);
    if proxy.is_finite() && proxy >= 0.65 {
        arm_identity_trace_blockalign(a, b).unwrap_or(proxy)
    } else {
        proxy
    }
}

/// Connected components of IRs linked by overlapping arm sequence: each IR
/// contributes two segments (its arms); overlapping segments join their IRs,
/// transitively. Returns (locus index per IR, size of each locus).
fn assign_loci(arms: &[(usize, usize, usize, usize)]) -> (Vec<usize>, Vec<usize>) {
    let n = arms.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut segs: Vec<(usize, usize, usize)> = Vec::with_capacity(2 * n);
    for (i, a) in arms.iter().enumerate() {
        segs.push((a.0, a.1, i));
        segs.push((a.2, a.3, i));
    }
    segs.sort_unstable();
    // Sweep: a segment starting before the running max end of the current run
    // overlaps that run (interval-graph components are exactly these runs).
    let mut run: Option<(usize, usize)> = None; // (max end, representative IR)
    for &(s, e, i) in &segs {
        match run {
            Some((end, rep)) if s < end => {
                let (ri, rr) = (find(&mut parent, i), find(&mut parent, rep));
                parent[ri] = rr;
                run = Some((end.max(e), rep));
            }
            _ => run = Some((e, i)),
        }
    }
    let mut id_of_root = vec![usize::MAX; n];
    let mut locus_of = vec![0usize; n];
    let mut sizes: Vec<usize> = Vec::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        if id_of_root[r] == usize::MAX {
            id_of_root[r] = sizes.len();
            sizes.push(0);
        }
        locus_of[i] = id_of_root[r];
        sizes[locus_of[i]] += 1;
    }
    (locus_of, sizes)
}

/// Bin size (bp) of the left-arm index used by Stage-B dedup.
const DEDUP_BIN: usize = 50_000;

/// Gap threshold multiplier for `--arm-bounds cluster`: a gap larger than this
/// many median inter-point gaps splits the arm. Gaps between minimizer matches
/// are roughly exponential, so P(gap > 10 x median ≈ 7 x mean) ≈ 0.1% per gap.
const ARM_GAP_FACTOR: f32 = 10.0;

/// Arm bounds by gap-based clustering of split-consistent matchpoints.
///
/// Each kept point is a pair (left, right) of mirrored positions; their
/// separation `right - left` grows by 2 bp per bp moved outward from the IR
/// centre, so sorting by it orders points from the arms' inner edges outward
/// without needing an estimate of the centre. Consecutive points further apart
/// (in arm bp) than `max(min_gap_bp, gap_factor * median gap)` start a new
/// cluster; the cluster spanning the most arm bp is the IR, and its min/max
/// give the arm bounds in full (no quantile trimming).
///
/// Why: quantile trimming cut ~10% from each end of every genuine arm
/// (plastid IR 25.7 kb -> 20.6 kb), while untrimmed min/max let a single
/// stray match from another repeat copy on the same anti-diagonal stretch an
/// arm by tens of kb. A dense homologous core has small gaps at any identity
/// the detector can see (median ~70 bp at 90%, ~125 bp at 85%), whereas stray
/// copies sit beyond a gap of kb.
///
/// Returns (la0, la1, ra0, ra1, kept_pts) in window-local coords.
#[allow(clippy::too_many_arguments)]
fn arm_bounds_near_split_diag_clustered(
    best_pts: &[(i32, i32)],
    split: usize,
    k: usize,
    win_len: usize,
    fold_diag_tol: i32,
    widen: i32,
    min_gap_bp: i32,
    gap_factor: f32,
) -> Option<(usize, usize, usize, usize, usize)> {
    let split2 = 2 * (split as i32);
    let tol = fold_diag_tol.max(1) * widen.max(1);

    // (separation, left, right) for split-consistent points.
    let mut pts: Vec<(i32, i32, i32)> = best_pts
        .iter()
        .filter(|&&(p1, p2)| (p1 + p2 - split2).abs() <= tol)
        .map(|&(p1, p2)| {
            let (l, r) = (p1.min(p2), p1.max(p2));
            (r - l, l, r)
        })
        .collect();
    if pts.is_empty() {
        return None;
    }
    pts.sort_unstable();

    // Gaps along the arm in bp (separation changes 2 bp per arm bp).
    let gaps: Vec<i32> = pts.windows(2).map(|w| (w[1].0 - w[0].0) / 2).collect();
    let mut positive: Vec<i32> = gaps.iter().copied().filter(|&g| g > 0).collect();
    let median_gap = if positive.is_empty() {
        0
    } else {
        let mid = positive.len() / 2;
        *positive.select_nth_unstable(mid).1
    };
    let threshold = (min_gap_bp.max(1) as f32).max(gap_factor * median_gap as f32);

    // Clusters as [start, end) index ranges. Keep the one spanning the most arm
    // bp (then most points, then innermost): extent is what defines an arm, while
    // point density falls with divergence -- by point count, a perfect 300 bp
    // stray copy can outvote a genuine 2 kb arm at 90% identity.
    let extent = |s: usize, e: usize| (pts[e - 1].0 - pts[s].0) / 2;
    let (mut best_s, mut best_e) = (0usize, 0usize);
    let mut best_key = (-1i32, 0usize);
    let mut s = 0usize;
    for i in 0..=gaps.len() {
        let boundary = i == gaps.len() || gaps[i] as f32 > threshold;
        if boundary {
            let e = i + 1;
            let key = (extent(s, e), e - s);
            if key > best_key {
                best_key = key;
                best_s = s;
                best_e = e;
            }
            s = i + 1;
        }
    }
    let cluster = &pts[best_s..best_e];

    let l0 = cluster.iter().map(|p| p.1).min()?;
    let l1 = cluster.iter().map(|p| p.1).max()?;
    let r0 = cluster.iter().map(|p| p.2).min()?;
    let r1 = cluster.iter().map(|p| p.2).max()?;

    let la0 = l0.max(0) as usize;
    let la1 = (l1.max(0) as usize).saturating_add(k).min(win_len);
    let ra0 = r0.max(0) as usize;
    let ra1 = (r1.max(0) as usize).saturating_add(k).min(win_len);

    if la0 < la1 && ra0 < ra1 {
        Some((la0, la1, ra0, ra1, cluster.len()))
    } else {
        None
    }
}

/// Split position implied by a bin's matchpoints: median of (p1 + p2) / 2.
///
/// Every point of a foldback/IR with centre `c` satisfies p1 + p2 ≈ 2c, so the
/// median anti-diagonal is a robust estimate of the centre that does not depend
/// on the spacer sequence (unlike breakpoint refinement).
fn median_diag_split(pts: &[(i32, i32)]) -> Option<usize> {
    if pts.is_empty() {
        return None;
    }
    let mut d: Vec<i64> = pts.iter().map(|&(p1, p2)| p1 as i64 + p2 as i64).collect();
    let mid = d.len() / 2;
    let (_, &mut m, _) = d.select_nth_unstable(mid);
    Some((m.max(0) / 2) as usize)
}

#[inline]
fn normalize_arms(
    mut la0: usize,
    mut la1: usize,
    mut ra0: usize,
    mut ra1: usize,
) -> (usize, usize, usize, usize) {
    // Ensure la0<=la1 and ra0<=ra1
    if la1 < la0 {
        std::mem::swap(&mut la0, &mut la1);
    }
    if ra1 < ra0 {
        std::mem::swap(&mut ra0, &mut ra1);
    }

    // Ensure "left arm" really is left of "right arm" (by start)
    if ra0 < la0 {
        std::mem::swap(&mut la0, &mut ra0);
        std::mem::swap(&mut la1, &mut ra1);
    }
    (la0, la1, ra0, ra1)
}

/// Build a conservative single-interval span from arm bounds.
/// Uses the shorter arm length to avoid inflation by a long/poorly-bounded arm.
///
/// Returns (start,end) in contig coords.
#[inline]
fn interval_from_arms_conservative(
    la0: usize,
    la1: usize,
    ra0: usize,
    ra1: usize,
    contig_len: usize,
) -> Option<(usize, usize)> {
    let (la0, la1, ra0, ra1) = normalize_arms(la0, la1, ra0, ra1);

    if la0 >= la1 || ra0 >= ra1 {
        return None;
    }

    let left_len = la1 - la0;
    let right_len = ra1 - ra0;
    let arm_len = left_len.min(right_len);

    // right arm "end" = ra0 + arm_len
    let end = (ra0 + arm_len).min(contig_len);

    if la0 < end { Some((la0, end)) } else { None }
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

/// Upper bound on the arm length that either arm-bounds function can report
/// from `pts`: the shorter of the left-arm and right-arm extents (plus k).
///
/// Matchpoints come in mirrored pairs, so the raw `p1` range spans the whole
/// IR footprint; splitting each point into `min(p1,p2)` / `max(p1,p2)` gives
/// the per-arm extents instead.
fn arm_extent_upper_bound(pts: &[(i32, i32)], k: usize) -> usize {
    if pts.is_empty() {
        return 0;
    }
    let (mut lmin, mut lmax, mut rmin, mut rmax) = (i32::MAX, i32::MIN, i32::MAX, i32::MIN);
    for &(p1, p2) in pts {
        let (l, r) = (p1.min(p2), p1.max(p2));
        lmin = lmin.min(l);
        lmax = lmax.max(l);
        rmin = rmin.min(r);
        rmax = rmax.max(r);
    }
    let left = (lmax - lmin).max(0) as usize + k;
    let right = (rmax - rmin).max(0) as usize + k;
    left.min(right)
}

/// Pack two u32 values into a u64 (cheap HashSet key).
///
/// Used for Stage-A dedup keys: `(break_bin, diag_bin)`.
#[inline]
fn pack_u32_pair(a: u32, b: u32) -> u64 {
    ((a as u64) << 32) | (b as u64)
}

/// Integer writer
#[inline]
fn write_usize<W: Write>(w: &mut W, mut x: usize) -> std::io::Result<()> {
    // small, allocation-free integer formatting
    let mut buf = [0u8; 32];
    let mut i = buf.len();
    if x == 0 {
        i -= 1;
        buf[i] = b'0';
    } else {
        while x > 0 {
            let d = (x % 10) as u8;
            i -= 1;
            buf[i] = b'0' + d;
            x /= 10;
        }
    }
    w.write_all(&buf[i..])
}

/// FASTA writer that avoids any UTF-8 conversion and avoids `format!`.
/// `seq` is written as raw bytes; line wrapping is done via slicing.
fn write_fasta_record_bytes<W: Write>(
    w: &mut W,
    contig: &str,
    start: usize,
    end: usize,
    break_pos: usize,
    ident_milli: u16, // identity*1000 (rounded)
    tir_ident: f32,
    matches: usize,
    span: usize,
    bin: i32,
    class: &str,
    spacer_bp: i32,
    seq: &[u8],
    wrap: usize,
) -> std::io::Result<()> {
    // Header:
    // >contig:start-end|break=...|ident=...|matches=...|span=...|bin=...
    w.write_all(b">")?;
    w.write_all(contig.as_bytes())?;
    w.write_all(b":")?;
    write_usize(w, start)?;
    w.write_all(b"-")?;
    write_usize(w, end)?;

    w.write_all(b"|break=")?;
    write_usize(w, break_pos)?;

    w.write_all(b"|ident=")?;
    // ident_milli as X.XXX without float formatting
    write_usize(w, (ident_milli / 1000) as usize)?;
    w.write_all(b".")?;
    let frac = ident_milli % 1000;
    w.write_all(&[
        b'0' + (frac / 100) as u8,
        b'0' + ((frac / 10) % 10) as u8,
        b'0' + (frac % 10) as u8,
    ])?;

    w.write_all(b"|tir_ident=")?;
    // tir_ident as X.XXX without float formatting
    if tir_ident.is_finite() && tir_ident >= 0.0 {
        let tir_milli = (tir_ident.clamp(0.0, 1.0) * 1000.0).round() as u16;
        write_usize(w, (tir_milli / 1000) as usize)?;
        w.write_all(b".")?;
        let frac = tir_milli % 1000;
        w.write_all(&[
            b'0' + (frac / 100) as u8,
            b'0' + ((frac / 10) % 10) as u8,
            b'0' + (frac % 10) as u8,
        ])?;
    } else {
        w.write_all(b"NaN")?;
    }

    w.write_all(b"|matches=")?;
    write_usize(w, matches)?;

    w.write_all(b"|span=")?;
    write_usize(w, span)?;

    w.write_all(b"|bin=")?;
    // bin i32 (handle negative without allocation)
    if bin < 0 {
        w.write_all(b"-")?;
        write_usize(w, (-bin) as usize)?;
    } else {
        write_usize(w, bin as usize)?;
    }

    w.write_all(b"|class=")?;
    w.write_all(class.as_bytes())?;

    w.write_all(b"|spacer=")?;
    // spacer_bp i32 without allocation
    if spacer_bp < 0 {
        w.write_all(b"-")?;
        write_usize(w, (-spacer_bp) as usize)?;
    } else {
        write_usize(w, spacer_bp as usize)?;
    }

    w.write_all(b"\n")?;

    // Sequence
    if wrap == 0 {
        w.write_all(seq)?;
        w.write_all(b"\n")?;
        return Ok(());
    }

    let mut i = 0usize;
    while i < seq.len() {
        let j = (i + wrap).min(seq.len());
        w.write_all(&seq[i..j])?;
        w.write_all(b"\n")?;
        i = j;
    }
    Ok(())
}

/// Very fast proxy: block-aligner global score -> identity-like quantity.
/// Good for filtering/ranking; not exact.
///
/// `right_rc` must already be reverse-complemented.
fn arm_identity_proxy_blockalign(left: &[u8], right_rc: &[u8]) -> Option<f32> {
    if left.is_empty() || right_rc.is_empty() {
        return None;
    }

    // choose query/ref
    let (q_bytes, r_bytes) = if left.len() <= right_rc.len() {
        (left, right_rc)
    } else {
        (right_rc, left)
    };

    // block sizes: powers of two; start small, allow growth
    let min_block = 64usize;
    let max_len = q_bytes.len().max(r_bytes.len()).max(64);
    let max_block = 1024usize.min(max_len.next_power_of_two()).max(min_block);

    let gaps = Gaps {
        open: -5,
        extend: -1,
    };

    let q = PaddedBytes::from_bytes::<NucMatrix>(q_bytes, max_block);
    let r = PaddedBytes::from_bytes::<NucMatrix>(r_bytes, max_block);

    // No TRACE, no X-drop
    let mut a = Block::<false, false>::new(q.len(), r.len(), max_block);
    a.align(&q, &r, &NW1, gaps, min_block..=max_block, 0);
    let res = a.res();

    // crude normalisation: score per aligned base mapped to [0,1]
    let aln_len = (q.len().min(r.len()).max(1)) as f32;
    let s = (res.score as f32) / aln_len;

    let proxy = ((s + 1.0) * 0.5).clamp(0.0, 1.0);
    Some(proxy)
}

/// Slower, but more accurate than `arm_identity_proxy_blockalign` because it
/// uses the traceback CIGAR to count matches/mismatches/indels.
///
/// `right_rc` must already be reverse-complemented.
fn arm_identity_trace_blockalign(left: &[u8], right_rc: &[u8]) -> Option<f32> {
    if left.is_empty() || right_rc.is_empty() {
        return None;
    }

    let (q_bytes, r_bytes) = if left.len() <= right_rc.len() {
        (left, right_rc)
    } else {
        (right_rc, left)
    };

    let min_block = 64usize;
    let max_len = q_bytes.len().max(r_bytes.len()).max(64);
    let max_block = 1024usize.min(max_len.next_power_of_two()).max(min_block);

    let gaps = Gaps {
        open: -5,
        extend: -1,
    };

    let q = PaddedBytes::from_bytes::<NucMatrix>(q_bytes, max_block);
    let r = PaddedBytes::from_bytes::<NucMatrix>(r_bytes, max_block);

    // TRACE=true, X_DROP=false
    let mut a = Block::<true, false>::new(q.len(), r.len(), max_block);
    a.align(&q, &r, &NW1, gaps, min_block..=max_block, 0);
    let res = a.res();

    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    a.trace()
        .cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);

    let mut m: u64 = 0;
    let mut x: u64 = 0;
    let mut ins: u64 = 0;
    let mut del: u64 = 0;

    for op in cigar.to_vec().iter() {
        match op.op {
            Operation::Eq => m += op.len as u64,
            Operation::X => x += op.len as u64,
            Operation::I => ins += op.len as u64,
            Operation::D => del += op.len as u64,
            _ => {}
        }
    }

    let denom = m + x + ins + del;
    if denom == 0 {
        None
    } else {
        Some(m as f32 / denom as f32)
    }
}

#[derive(Copy, Clone, Debug)]
enum CapPolicy {
    Clamp,
    Drop,
    Penalize,
}

impl CapPolicy {
    fn parse(s: &str) -> anyhow::Result<Self> {
        match s {
            "clamp" => Ok(Self::Clamp),
            "drop" => Ok(Self::Drop),
            "penalize" => Ok(Self::Penalize),
            _ => anyhow::bail!("--cap-policy must be one of: clamp, drop, penalize"),
        }
    }
}

#[derive(Default, Debug, Clone)]
struct Stats {
    // pipeline volumes
    windows: u64,
    // FIXME: not currently used, are these useful?
    // hits: u64,
    // refined_ok: u64,
    // ident_ok: u64,

    // after collecting + sorting
    cands_total: u64,

    // main-thread dedup stages
    stage_a_kept: u64,
    stage_a_dups: u64,
    stage_b_kept: u64,
    low_ident: u64,
    stage_b_dups: u64,

    // cap diagnostics
    clamped: u64,  // clamped due to cap policy
    near_cap: u64, // emitted intervals close to max_interval_bp (optional)
    emitted: u64,

    // diagnostics so we can see where caps come from
    arm_missing: u64,
    precap_too_big: u64,
}

#[derive(Debug, Clone, Serialize)]
struct ReportRow {
    contig: String,
    start: usize,
    end: usize,
    break_pos: usize,
    identity_est: f32,
    tir_ident: Option<f32>,
    matches: usize,
    span: usize,
    la0: usize,
    la1: usize,
    ra0: usize,
    ra1: usize,
    contig_len: usize,
    win_start: usize,
    win_end: usize,
    kept_pts: usize,
    bin: i32,
    arm_len: usize,
    spacer: isize,
    ir_class: String,
    interval_len: usize,
    left_arm_len: usize,
    right_arm_len: usize,
    arm_ratio: f32,
    dir_ident: Option<f32>,
    direct_better: bool,
    locus_id: String,
    locus_n_irs: usize,
    pass: &'static str,
}

/// Percent-encode characters GFF3 does not allow verbatim in a seqid or
/// attribute value (control characters, whitespace, `%;=&,` and other
/// non-seqid punctuation).
fn gff_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        let ok = b.is_ascii_alphanumeric() || b".:^*$@!+_?-|".contains(&b);
        if ok {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

/// Write IR calls as GFF3: one `inverted_repeat` parent per call spanning
/// both arms, with two `repeat_component` children (left arm `+`, right arm
/// `-`). Coordinates are converted from 0-based half-open to 1-based closed.
///
/// Rows for a contig are contiguous in `rows` (contigs are processed in
/// order); within a contig, IRs are numbered and features sorted by start so
/// the output is tabix-ready.
fn write_gff3<W: Write>(
    w: &mut W,
    cmdline: &str,
    contigs: &[(String, usize)],
    rows: &[ReportRow],
) -> std::io::Result<()> {
    writeln!(w, "##gff-version 3")?;
    writeln!(w, "#!irx-version {}", clap::crate_version!())?;
    writeln!(w, "#!irx-command {cmdline}")?;
    for (id, len) in contigs {
        if *len > 0 {
            writeln!(w, "##sequence-region {} 1 {len}", gff_escape(id))?;
        }
    }

    let mut n_ir: usize = 0;
    for chunk in rows.chunk_by(|a, b| a.contig == b.contig) {
        let seqid = gff_escape(&chunk[0].contig);

        // Normalised arms, then IRs in genome order (ids follow this order).
        let mut irs: Vec<(&ReportRow, (usize, usize, usize, usize))> = chunk
            .iter()
            .map(|r| (r, normalize_arms(r.la0, r.la1, r.ra0, r.ra1)))
            .collect();
        irs.sort_by_key(|&(r, (la0, la1, ra0, ra1))| {
            if la1 > la0 && ra1 > ra0 {
                (la0, ra1)
            } else {
                (r.start, r.end)
            }
        });

        // (start, kind, line): kind orders parent < left < right at equal starts.
        let mut feats: Vec<(usize, u8, String)> = Vec::with_capacity(irs.len() * 3);
        for (r, (la0, la1, ra0, ra1)) in irs {
            let has_arms = la1 > la0 && ra1 > ra0;
            let (p0, p1) = if has_arms {
                (la0, la1.max(ra1))
            } else {
                (r.start, r.end)
            };
            if p1 <= p0 {
                continue;
            }
            n_ir += 1;
            let ir_id = format!("IR{n_ir:06}");
            let score = match r.tir_ident {
                Some(t) => format!("{t:.3}"),
                None => ".".to_string(),
            };
            let dir_attr = match r.dir_ident {
                Some(t) => format!("{t:.3}"),
                None => "NA".to_string(),
            };
            let tir_attr = match r.tir_ident {
                Some(t) => format!("{t:.3}"),
                None => "NA".to_string(),
            };
            feats.push((
                p0,
                0,
                format!(
                    "{seqid}\tirx\tinverted_repeat\t{}\t{p1}\t{score}\t.\t.\t\
ID={ir_id};arm_len={};left_arm_len={};right_arm_len={};spacer={};ir_class={};\
tir_ident={tir_attr};dir_ident={dir_attr};direct_better={};identity_est={:.3};\
break_pos={};matches={};kept_pts={};locus={};locus_n_irs={};pass={}",
                    p0 + 1,
                    r.arm_len,
                    r.left_arm_len,
                    r.right_arm_len,
                    r.spacer,
                    gff_escape(&r.ir_class),
                    u8::from(r.direct_better),
                    r.identity_est,
                    r.break_pos + 1,
                    r.matches,
                    r.kept_pts,
                    gff_escape(&r.locus_id),
                    r.locus_n_irs,
                    r.pass,
                ),
            ));
            if has_arms {
                feats.push((
                    la0,
                    1,
                    format!(
                        "{seqid}\tirx\trepeat_component\t{}\t{la1}\t.\t+\t.\tID={ir_id}.L;Parent={ir_id};Name=left_arm",
                        la0 + 1
                    ),
                ));
                feats.push((
                    ra0,
                    2,
                    format!(
                        "{seqid}\tirx\trepeat_component\t{}\t{ra1}\t.\t-\t.\tID={ir_id}.R;Parent={ir_id};Name=right_arm",
                        ra0 + 1
                    ),
                ));
            }
        }

        // Stable sort keeps IR-id order for features sharing (start, kind).
        feats.sort_by_key(|f| (f.0, f.1));
        for (_, _, line) in feats {
            writeln!(w, "{line}")?;
        }
    }
    Ok(())
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
    rank: u64,
    key_a: u64, // Stage A key: (break_bin, diag_bin)

    capped: bool,
    had_arms: bool,
    precap_too_big: bool,

    // ---- classification ----
    class: &'static str,
    spacer_bp: i32, // or isize/i64

    // ---- which pass found it (false = local windows, true = large windows) ----
    large: bool,
}

fn write_html_report(
    out_path: &std::path::Path,
    input: &std::path::Path,
    rows: &[ReportRow],
) -> Result<()> {
    const IRX_HTML_TEMPLATE: &str = include_str!("irx.html");
    const IRX_DASHBOARD_JS: &str = include_str!("irx.js");

    let data_json = serde_json::to_string(rows)?;
    let species_name = input
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("irx report");

    let html = IRX_HTML_TEMPLATE
        .replace("__TITLE__", species_name)
        .replace("__DATA_JSON__", &data_json)
        .replace("__DASHBOARD_JS__", IRX_DASHBOARD_JS);

    std::fs::write(out_path, html)?;
    Ok(())
}

#[inline]
const fn build_invalid_dna_table() -> [u8; 256] {
    let mut t = [1u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 0;
    t[b'G' as usize] = 0;
    t[b'T' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 0;
    t[b'g' as usize] = 0;
    t[b't' as usize] = 0;
    t
}

const INVALID_DNA: [u8; 256] = build_invalid_dna_table();

/// Very fast yes/no check.
#[inline]
fn is_valid_dna_acgt(seq: &[u8]) -> bool {
    let mut bad = 0u8;
    for &b in seq {
        bad |= INVALID_DNA[b as usize];
    }
    bad == 0
}

#[inline]
fn is_valid_dna_byte(b: u8) -> bool {
    INVALID_DNA[b as usize] == 0
}

/// Return the index and byte value of the first invalid character.
#[inline]
fn first_invalid_dna_byte(seq: &[u8]) -> Option<(usize, u8)> {
    for (i, &b) in seq.iter().enumerate() {
        if INVALID_DNA[b as usize] != 0 {
            return Some((i, b));
        }
    }
    None
}

/// Copy of `seq` with every non-ACGT byte replaced by a pseudo-random base that
/// depends only on its position (splitmix64), so output is identical for any
/// thread count or window layout.
fn fill_invalid_bases(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .map(|(i, &b)| {
            if is_valid_dna_byte(b) {
                b
            } else {
                let mut z = (i as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                b"ACGT"[(z >> 62) as usize]
            }
        })
        .collect()
}

fn main() -> Result<()> {
    let args = build_cli();

    let input = args.get_one::<PathBuf>("input").unwrap();
    let bed = args.get_one::<PathBuf>("bed").unwrap();
    let fasta = args.get_one::<PathBuf>("fasta").cloned();
    let fasta_wrap = *args.get_one::<usize>("fasta_wrap").unwrap();

    let window_len = *args.get_one::<usize>("window_len").unwrap();
    let step = *args.get_one::<usize>("step").unwrap();
    let hits_per_window = *args.get_one::<usize>("hits_per_window").unwrap();
    let min_arm = *args.get_one::<usize>("min_arm").unwrap();
    let min_matches = *args.get_one::<usize>("min_matches").unwrap();
    let dedup_bp = *args.get_one::<usize>("dedup_bp").unwrap();
    let threads = *args.get_one::<usize>("threads").unwrap();
    let max_interval_bp = *args.get_one::<usize>("max_interval_bp").unwrap();
    let large_window_len = *args.get_one::<usize>("large_window_len").unwrap();
    let large_step = *args.get_one::<usize>("large_step").unwrap();
    let large_max_interval_bp = *args.get_one::<usize>("large_max_interval_bp").unwrap();
    // The local pass sees an IR whole only up to this footprint; the large pass
    // takes over above it.
    let large_min_footprint = window_len.saturating_sub(step.max(1));
    let diag_widen = *args.get_one::<i32>("diag_widen").unwrap();
    let arm_gap_bp = *args.get_one::<i32>("arm_gap_bp").unwrap();
    let min_tir_ident = *args.get_one::<f32>("min_tir_ident").unwrap();
    let cap_policy = CapPolicy::parse(args.get_one::<String>("cap_policy").unwrap())?;

    // IR classification
    let classify_ir = args.get_flag("classify_ir");
    let immediate_bp = *args.get_one::<usize>("immediate_bp").unwrap();
    let wide_bp = *args.get_one::<usize>("wide_bp").unwrap();

    // HTML output
    let html = args.get_one::<PathBuf>("html").cloned();

    // GFF3 output
    let gff = args.get_one::<PathBuf>("gff").cloned();

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
            min_matches,
            end_guard: 0,
            refine: RefineCfg {
                window: 500,
                arm: 1200,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.25,
                max_jump_clip: 1000,
            },
            fold_diag_tol: 120,
        },
        // Coarse gate only: mdax's bin span is the IR footprint (arm + spacer + arm),
        // which is >= one arm. The per-arm check happens after arm bounds are derived.
        fold: FoldOnlyCfg { min_arm },
        fold2: FoldSecondPassCfg {
            min_support: 0,
            min_identity: 0.0,
            min_support_ident: 0.0,
            cut_low_ident: false,
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
        // TSV writer
        let bed_is_stdout = bed.to_string_lossy() == "-";
        let mut out: Box<dyn Write> = if bed_is_stdout {
            Box::new(std::io::BufWriter::new(std::io::stdout()))
        } else {
            Box::new(std::io::BufWriter::new(std::fs::File::create(bed)?))
        };

        // Optional FASTA writer
        let fasta_is_stdout = fasta
            .as_ref()
            .map(|p| p.to_string_lossy() == "-")
            .unwrap_or(false);

        let gff_is_stdout = gff
            .as_ref()
            .map(|p| p.to_string_lossy() == "-")
            .unwrap_or(false);

        // at most one output can go to stdout
        if [bed_is_stdout, fasta_is_stdout, gff_is_stdout]
            .iter()
            .filter(|&&b| b)
            .count()
            > 1
        {
            anyhow::bail!(
                "Only one of TSV (--bed -), FASTA (--fasta -) and GFF (--gff -) can go to stdout. Send the others to files."
            );
        }

        let mut fasta_out: Option<Box<dyn Write>> = match fasta.as_ref() {
            None => None,
            Some(p) if p.to_string_lossy() == "-" => {
                Some(Box::new(std::io::BufWriter::new(std::io::stdout())))
            }
            Some(p) => Some(Box::new(std::io::BufWriter::new(std::fs::File::create(p)?))),
        };

        // TSV header
        let header = if classify_ir {
            writeln!(
                out,
                "#contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\ttir_ident\tmatches\tspan\tla0\tla1\tra0\tra1\tcontig_len\twin_start\twin_end\tkept_pts\tbin\tarm_len\tspacer\tir_class\tdir_ident\tdirect_better\tlocus_id\tlocus_n_irs\tpass"
            )
        } else {
            writeln!(
                out,
                "#contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\ttir_ident\tmatches\tspan\tla0\tla1\tra0\tra1\tcontig_len\twin_start\twin_end\tkept_pts\tbin\tdir_ident\tdirect_better\tlocus_id\tlocus_n_irs\tpass"
            )
        };

        match header {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
            Err(e) => return Err(e.into()),
        }

        let mut report_rows: Vec<ReportRow> = Vec::new();
        // Running locus number across contigs (locus ids are unique per file).
        let mut locus_counter: usize = 0;
        // Every contig scanned (for GFF ##sequence-region lines).
        let mut contigs_seen: Vec<(String, usize)> = Vec::new();

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
            if gff.is_some() {
                contigs_seen.push((id.to_string(), len));
            }

            let q = dedup_bp.max(1);
            let mm_k = cfg.shared.minimizer.k;
            let fold_diag_tol = cfg.shared.fold_diag_tol;

            // Scan the whole contig. Non-ACGT bases (assembly gaps, IUPAC codes) are
            // replaced, in a private copy used only for detection, by deterministic
            // pseudo-random bases: random sequence cannot create k=17 IR matches, so
            // gaps add no signal, but an IR whose spacer or arm contains a gap is no
            // longer split apart. (Previously windows only covered ACGT-only spans, so
            // an IR with an N-gap in its spacer was invisible.) Arm identity and the
            // FASTA output still use the original sequence.
            let filled: Vec<u8>;
            let scan: &[u8] = match first_invalid_dna_byte(&seq) {
                None => &seq,
                Some((pos, bad)) => {
                    elog!(
                        "irx",
                        "{id}: filling non-ACGT bases for detection (first '{}' (ASCII {}) at offset {})",
                        bad as char,
                        bad,
                        pos,
                    );
                    filled = fill_invalid_bases(&seq);
                    &filled
                }
            };

            if len < cfg.fold.min_arm.saturating_mul(2) {
                continue;
            }
            let wins: Vec<(usize, usize)> = windows(len, window_len, step).collect();

            // For each window:
            // - detect top-k coarse foldback bins
            // - refine split position per hit
            // - compute approximate arm interval bounds
            // - package into `Cand` records
            //
            // We do *no* dedup here (that happens serially per contig).

            let gen_cands = |wins: &[(usize, usize)],
                             fold: &FoldOnlyCfg,
                             max_interval_bp: usize,
                             large: bool|
             -> Vec<Cand> {
            wins
                .par_iter()
                .flat_map_iter(|&(wstart, wend)| {
                    let sub = &scan[wstart..wend];
                    let win_len = sub.len();

                    // Defensive check: in the no-invalid path this is always true;
                    // in the split path it should also always be true because windows
                    // are generated only from valid spans.
                    debug_assert!(is_valid_dna_acgt(sub));

                    // Thread-local scratch reuse.
                    let scratch_cell = tls.get_or(|| RefCell::new(FoldScratch::new()));
                    let mut scratch = scratch_cell.borrow_mut();

                    // Find up to K foldback candidates in this window.
                    let hits = foldback::detect_foldbacks_topk(
                        sub,
                        &cfg.shared,
                        fold,
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
                        if hit.bin_span < fold.min_arm {
                            continue;
                        }
                        // Cheap pre-refinement gate: the arm bounds derived below come
                        // from a subset of these points, so this is an upper bound on
                        // the reported arm length.
                        if arm_extent_upper_bound(&hit.pts, mm_k) < min_arm {
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
                            fb.span,
                            fb.arm_end,
                            fb.arm_start_right,
                        )
                        .ok()
                        .flatten() else {
                            continue;
                        };

                        // Global break position (contig coords).
                        let break_pos = wstart + rf.split_pos;

                        // Determine arm bounds:
                        //   1) prefer points near the bin's anti-diagonal (gap-clustered)
                        //   2) fall back to unfiltered bounds
                        //   3) fall back to break±span/2
                        //
                        // Centre the filter on the bin's own anti-diagonal, not on
                        // `rf.split_pos`: refinement compares `refine.arm` bp either side
                        // of the split, which for an IR with a spacer is spacer vs spacer,
                        // so the refined split can wander by hundreds of bp. When it lands
                        // > tol/2 from the true centre every point is rejected and the
                        // unfiltered fallback inflates arms with stray matches (48% of
                        // simulated IRs with a 10 kb spacer did this).
                        let diag_split = median_diag_split(&hit.pts).unwrap_or(rf.split_pos);
                        let bounds = arm_bounds_near_split_diag_clustered(
                            &hit.pts,
                            diag_split,
                            mm_k,
                            win_len,
                            fold_diag_tol,
                            diag_widen,
                            arm_gap_bp,
                            ARM_GAP_FACTOR,
                        );
                        let (start, end, la0g, la1g, ra0g, ra1g, kept_pts) =
                            if let Some((la0, la1, ra0, ra1, kept)) = bounds
                            {
                                let la0g = wstart + la0;
                                let la1g = wstart + la1;
                                let ra0g = wstart + ra0;
                                let ra1g = wstart + ra1;
                                let start = la0g.min(ra0g);
                                let end = la1g.max(ra1g);
                                (start, end, la0g, la1g, ra0g, ra1g, kept)
                            } else if let Some((la0, la1, ra0, ra1)) =
                                // TODO: fall back to arm_bounds_near_split_diag?
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

                        // Hard cap interval size to avoid window-scale inflation.
                        let mut start = start;
                        let mut end = end;

                        // If we have real arm bounds, override the interval with a conservative one.
                        // This prevents huge “union inflation” when one arm is poorly bounded.
                        if la1g > la0g && ra1g > ra0g {
                            if let Some((s2, e2)) =
                                interval_from_arms_conservative(la0g, la1g, ra0g, ra1g, len)
                            {
                                start = s2;
                                end = e2;
                            }
                        }

                        // Per-arm length gate, applied before dedup so a short-armed
                        // candidate can never suppress an overlapping long-armed one.
                        // Uses the same arm_len definition as the TSV column.
                        let (n_la0, n_la1, n_ra0, n_ra1) = normalize_arms(la0g, la1g, ra0g, ra1g);
                        let cand_arm_len =
                            (n_la1.saturating_sub(n_la0)).min(n_ra1.saturating_sub(n_ra0));
                        if cand_arm_len < min_arm {
                            continue;
                        }

                        let raw_len = end.saturating_sub(start);
                        let had_arms = la1g > la0g && ra1g > ra0g;
                        let precap_too_big = raw_len > max_interval_bp && raw_len > 0;

                        let mut capped_this = false;
                        if raw_len > max_interval_bp && raw_len > 0 {
                            match cap_policy {
                                CapPolicy::Drop => {
                                    // skip this hit entirely
                                    continue;
                                }
                                CapPolicy::Clamp => {
                                    // Clamp, but do NOT mark as "capped" (no ranking penalty)
                                    let half = max_interval_bp / 2;
                                    start = break_pos.saturating_sub(half);
                                    end = (start + max_interval_bp).min(len);
                                    if end - start < max_interval_bp {
                                        start = end.saturating_sub(max_interval_bp);
                                    }
                                }
                                CapPolicy::Penalize => {
                                    // Clamp AND mark as capped (absolute ranking penalty)
                                    capped_this = true;

                                    let half = max_interval_bp / 2;
                                    start = break_pos.saturating_sub(half);
                                    end = (start + max_interval_bp).min(len);
                                    if end - start < max_interval_bp {
                                        start = end.saturating_sub(max_interval_bp);
                                    }
                                }
                            }
                        }

                        // Score is 0..1000 scale of identity estimate (BED-ish).
                        let score = (rf.identity_est.clamp(0.0, 1.0) * 1000.0).round() as i32;

                        // ---- Step 2: absolute "penalize" ranking ----
                        // Any non-capped hit MUST outrank any capped hit.
                        // We'll encode that as a top "group bit": 1 = non-capped, 0 = capped.
                        let score_u: u64 = score.clamp(0, 1000) as u64;
                        let bin_rank_u: u64 = (hit.bin_rank as u64).min(0xFFFF_FFFF);
                        let matches_u: u64 = (fb.matches as u64).min(0xFFFF);

                        // If we're in Penalize mode, capped_this means "belongs to the losing group".
                        // Otherwise (Clamp), capped_this is false and everyone stays in the winning group.
                        let group_bit: u64 = if capped_this { 0 } else { 1 };

                        // Layout (high -> low):
                        // [ group:1 ][ score:11 ][ bin_rank:32 ][ matches:16 ][ spare:4 ]
                        let rank: u64 = (group_bit << 63)
                            | (score_u << 52)
                            | (bin_rank_u << 20)
                            | (matches_u << 4);

                        // Stage A key: (break_bin, diag_bin) packed into u64.
                        // Including `diag_bin` makes this less aggressive (keeps distinct anti-diagonals).
                        let break_bin: u32 = (break_pos / q) as u32;
                        let diag_u32: u32 = bin as u32; // OK even if negative -> wraps; only used as a hash-ish key
                        let key_a: u64 = pack_u32_pair(break_bin, diag_u32);


                        // Use the same classifier (and the same user-configurable
                        // --immediate-bp/--wide-bp thresholds) as the TSV's `ir_class`
                        // column, so the FASTA header's `class=` field can't disagree
                        // with it for the same candidate.
                        let (class, spacer_bp) = if la1g > la0g && ra1g > ra0g {
                            let (_n_la0, n_la1, n_ra0, _n_ra1) =
                                normalize_arms(la0g, la1g, ra0g, ra1g);
                            let spacer = n_ra0 as isize - n_la1 as isize;
                            let cls = classify_ir_str(spacer, immediate_bp, wide_bp);
                            (cls, spacer.clamp(i32::MIN as isize, i32::MAX as isize) as i32)
                        } else {
                            ("unknown", 0)
                        };

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
                            capped: capped_this,
                            had_arms,
                            precap_too_big,
                            class,
                            spacer_bp,
                            large,
                        });
                    }

                    out_local.into_iter()
                })
                .collect()
            };

            // Pass 1 (local): complete for footprints <= window_len - step.
            let mut cands: Vec<Cand> = gen_cands(&wins, &cfg.fold, max_interval_bp, false);
            let mut n_wins = wins.len();

            // Pass 2 (large): longer windows for footprints the local pass cannot see
            // whole. The coarse bin-span gate (which measures footprint) is raised to the
            // local completeness limit so small IRs do not take the large windows' top-K
            // slots; only candidates whose arms span more than that are kept, and both
            // passes then go through the same dedup.
            if large_window_len > 0 && len > large_min_footprint {
                let large_fold = FoldOnlyCfg {
                    min_arm: large_min_footprint.max(cfg.fold.min_arm),
                };
                let large_wins: Vec<(usize, usize)> =
                    windows(len, large_window_len, large_step).collect();
                n_wins += large_wins.len();
                let large_cands =
                    gen_cands(&large_wins, &large_fold, large_max_interval_bp, true);
                cands.extend(large_cands.into_iter().filter(|c| {
                    c.had_arms && {
                        let (la0, _, _, ra1) = normalize_arms(c.la0, c.la1, c.ra0, c.ra1);
                        ra1.saturating_sub(la0) > large_min_footprint
                    }
                }));
            }

            let mut st = Stats::default();
            st.windows = n_wins as u64;
            st.cands_total = cands.len() as u64;

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

            // Dedup + identity filter. Semantics: walk candidates best-first and accept
            // one unless (A) it repeats an accepted call's (break_bin, diag_bin) key,
            // (B) both its arms overlap an accepted call's arms (the same IR from another
            // window -- whole-interval overlap is not enough, that also dropped nested and
            // merely overlapping IRs), or its arm identity is below --min-tir-ident (so a
            // non-homologous call can never suppress a genuine one).
            //
            // Arm identities (block alignment of whole arms) were ~70% of runtime when
            // computed one at a time on this (serial) thread, so: run the cheap serial
            // dedup skipping candidates already known to be low-identity, compute
            // identities for the newly accepted ones in parallel, mark failures, repeat
            // until no accepted candidate fails (usually 1-2 rounds). The fixed point is
            // the sequential result, except that a low-identity candidate no longer
            // reserves its Stage-A key.
            let cap_thresh = (max_interval_bp as f64 * 0.95).round() as usize; // “near cap” threshold
            let ov = |a0: usize, a1: usize, b0: usize, b1: usize| a0 < b1 && b0 < a1;
            let cand_arms: Vec<(usize, usize, usize, usize)> = cands
                .iter()
                .map(|c| {
                    if c.had_arms {
                        normalize_arms(c.la0, c.la1, c.ra0, c.ra1)
                    } else {
                        (c.start, c.end, c.start, c.end)
                    }
                })
                .collect();
            let mut ident: Vec<Option<(f32, f32)>> = vec![None; cands.len()];
            let mut low = vec![false; cands.len()];
            let (accepted_idx, accepted_arms) = loop {
                let mut seen_break_bins: HashSet<u64> = HashSet::new();
                let mut acc_arms: Vec<(usize, usize, usize, usize)> = Vec::new();
                let mut arm_index: HashMap<usize, Vec<usize>> = HashMap::new();
                let mut acc: Vec<usize> = Vec::new();
                st.stage_a_dups = 0;
                st.stage_a_kept = 0;
                st.stage_b_dups = 0;
                st.stage_b_kept = 0;
                for (i, c) in cands.iter().enumerate() {
                    if low[i] {
                        continue;
                    }
                    if !seen_break_bins.insert(c.key_a) {
                        st.stage_a_dups += 1;
                        continue;
                    }
                    st.stage_a_kept += 1;
                    let arms = cand_arms[i];
                    let duplicate =
                        (arms.0 / DEDUP_BIN..=arms.1.saturating_sub(1) / DEDUP_BIN).any(|bin| {
                            arm_index.get(&bin).is_some_and(|ids| {
                                ids.iter().any(|&j| {
                                    let a = acc_arms[j];
                                    ov(arms.0, arms.1, a.0, a.1) && ov(arms.2, arms.3, a.2, a.3)
                                })
                            })
                        });
                    if duplicate {
                        st.stage_b_dups += 1;
                        continue;
                    }
                    st.stage_b_kept += 1;
                    for bin in arms.0 / DEDUP_BIN..=arms.1.saturating_sub(1) / DEDUP_BIN {
                        arm_index.entry(bin).or_default().push(acc_arms.len());
                    }
                    acc_arms.push(arms);
                    acc.push(i);
                }

                let todo: Vec<usize> = acc.iter().copied().filter(|&i| ident[i].is_none()).collect();
                let results: Vec<(usize, f32, f32)> = todo
                    .par_iter()
                    .map(|&i| {
                        let c = &cands[i];
                        let (la0n, la1n, ra0n, ra1n) = normalize_arms(c.la0, c.la1, c.ra0, c.ra1);
                        if !(la1n > la0n && ra1n > ra0n && la1n <= seq.len() && ra1n <= seq.len()) {
                            return (i, f32::NAN, f32::NAN);
                        }
                        let (left, right) = (&seq[la0n..la1n], &seq[ra0n..ra1n]);
                        let mut rc_buf: Vec<u8> = Vec::with_capacity(right.len());
                        revcomp_into(&mut rc_buf, right);
                        let tir = arm_identity(left, &rc_buf);
                        // Direct identity only for calls that will be kept.
                        let dir = if min_tir_ident > 0.0 && !(tir >= min_tir_ident) {
                            f32::NAN
                        } else {
                            arm_identity(left, right)
                        };
                        (i, tir, dir)
                    })
                    .collect();
                let mut any_low = false;
                for (i, tir, dir) in results {
                    ident[i] = Some((tir, dir));
                    if min_tir_ident > 0.0 && !(tir >= min_tir_ident) {
                        low[i] = true;
                        any_low = true;
                    }
                }
                if !any_low {
                    break (acc, acc_arms);
                }
            };
            st.low_ident = low.iter().filter(|&&l| l).count() as u64;

            // Accepted calls for this contig, buffered so loci can be assigned first.
            let mut accepted: Vec<(Cand, ReportRow)> = Vec::with_capacity(accepted_idx.len());
            for &i in &accepted_idx {
                let c = cands[i].clone();
                let (tir_ident, dir_ident) = ident[i].unwrap_or((f32::NAN, f32::NAN));
                let (la0n, la1n, ra0n, ra1n) = normalize_arms(c.la0, c.la1, c.ra0, c.ra1);
                let arm_len = (la1n.saturating_sub(la0n)).min(ra1n.saturating_sub(ra0n));
                let spacer: isize = (ra0n as isize) - (la1n as isize);
                let ir_class = classify_ir_str(spacer, immediate_bp, wide_bp);

                // Cap-hitter diagnostic (use the *actual* emitted length)
                st.emitted += 1;

                // actually clamped (only happens in clamp/penalize)
                if c.capped {
                    st.clamped += 1;
                }

                // “near cap” by emitted span (diagnostic)
                let span = c.end.saturating_sub(c.start);
                if span >= cap_thresh {
                    st.near_cap += 1;
                }

                // where caps come from
                if !c.had_arms {
                    st.arm_missing += 1;
                }
                if c.precap_too_big {
                    st.precap_too_big += 1;
                }

                let left_arm_len = la1n.saturating_sub(la0n);
                let right_arm_len = ra1n.saturating_sub(ra0n);
                let arm_ratio = if left_arm_len > 0 && right_arm_len > 0 {
                    let a = left_arm_len.max(right_arm_len) as f32;
                    let b = left_arm_len.min(right_arm_len) as f32;
                    a / b
                } else {
                    f32::NAN
                };

                let pass = if c.large { "large" } else { "local" };
                accepted.push((
                    c,
                    ReportRow {
                        contig: id.to_string(),
                        start: 0,
                        end: 0,
                        break_pos: 0,
                        identity_est: 0.0,
                        tir_ident: tir_ident.is_finite().then_some(tir_ident),
                        matches: 0,
                        span: 0,
                        la0: 0,
                        la1: 0,
                        ra0: 0,
                        ra1: 0,
                        contig_len: len,
                        win_start: 0,
                        win_end: 0,
                        kept_pts: 0,
                        bin: 0,
                        arm_len,
                        spacer,
                        ir_class: ir_class.to_string(),
                        interval_len: 0,
                        left_arm_len,
                        right_arm_len,
                        arm_ratio,
                        dir_ident: dir_ident.is_finite().then_some(dir_ident),
                        direct_better: dir_ident.is_finite()
                            && tir_ident.is_finite()
                            && dir_ident > tir_ident,
                        locus_id: String::new(),
                        locus_n_irs: 0,
                        pass,
                    },
                ));
            }

            // Loci: IRs whose arms share sequence (overlap, transitively). Numbered in
            // genome order of each locus' leftmost arm, continuing across contigs.
            let (locus_of, locus_sizes) = assign_loci(&accepted_arms);
            let mut order: Vec<usize> = (0..locus_sizes.len()).collect();
            let mut first_pos = vec![usize::MAX; locus_sizes.len()];
            for (i, a) in accepted_arms.iter().enumerate() {
                first_pos[locus_of[i]] = first_pos[locus_of[i]].min(a.0);
            }
            order.sort_by_key(|&l| (first_pos[l], l));
            let mut locus_name = vec![String::new(); locus_sizes.len()];
            for l in order {
                locus_counter += 1;
                locus_name[l] = format!("L{locus_counter:07}");
            }

            for (i, (c, mut row)) in accepted.into_iter().enumerate() {
                row.start = c.start;
                row.end = c.end;
                row.break_pos = c.break_pos;
                row.identity_est = c.ident;
                row.matches = c.matches;
                row.span = c.span;
                row.la0 = c.la0;
                row.la1 = c.la1;
                row.ra0 = c.ra0;
                row.ra1 = c.ra1;
                row.win_start = c.win_start;
                row.win_end = c.win_end;
                row.kept_pts = c.kept_pts;
                row.bin = c.bin;
                row.interval_len = c.end.saturating_sub(c.start);
                row.locus_id = locus_name[locus_of[i]].clone();
                row.locus_n_irs = locus_sizes[locus_of[i]];

                let tir_ident = row.tir_ident.unwrap_or(f32::NAN);
                let dir_ident = row.dir_ident.unwrap_or(f32::NAN);
                let extra = format!(
                    "\t{:.3}\t{}\t{}\t{}\t{}",
                    dir_ident,
                    u8::from(row.direct_better),
                    row.locus_id,
                    row.locus_n_irs,
                    row.pass
                );

                let r = if classify_ir {
                    writeln!(
                        out,
                        "{id}\t{}\t{}\tIR\t{}\t.\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{extra}",
                        c.start,
                        c.end,
                        c.score,
                        c.break_pos,
                        c.ident,
                        tir_ident,
                        c.matches,
                        c.span,
                        c.la0,
                        c.la1,
                        c.ra0,
                        c.ra1,
                        len,
                        c.win_start,
                        c.win_end,
                        c.kept_pts,
                        c.bin,
                        row.arm_len,
                        row.spacer,
                        row.ir_class
                    )
                } else {
                    writeln!(
                        out,
                        "{id}\t{}\t{}\tIR\t{}\t.\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{extra}",
                        c.start,
                        c.end,
                        c.score,
                        c.break_pos,
                        c.ident,
                        tir_ident,
                        c.matches,
                        c.span,
                        c.la0,
                        c.la1,
                        c.ra0,
                        c.ra1,
                        len,
                        c.win_start,
                        c.win_end,
                        c.kept_pts,
                        c.bin
                    )
                };
                match r {
                    Ok(_) => {}
                    Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
                    Err(e) => return Err(e.into()),
                }

                if let Some(fw) = fasta_out.as_mut() {
                    if c.start < c.end && c.end <= seq.len() {
                        let ident_milli = (c.ident.clamp(0.0, 1.0) * 1000.0).round() as u16;

                        match write_fasta_record_bytes(
                            fw,
                            id,
                            c.start,
                            c.end,
                            c.break_pos,
                            ident_milli,
                            tir_ident,
                            c.matches,
                            c.span,
                            c.bin,
                            c.class,
                            c.spacer_bp,
                            &seq[c.start..c.end],
                            fasta_wrap,
                        ) {
                            Ok(_) => {}
                            Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
                            Err(e) => return Err(e.into()),
                        }
                    }
                }

                report_rows.push(row);
            }

            if st.emitted > 0 {
                let clamped_pct = 100.0 * (st.clamped as f64) / (st.emitted as f64);
                let near_cap_pct = 100.0 * (st.near_cap as f64) / (st.emitted as f64);

                elog!(
                    "irx",
                    "{id}: wins={} cands={} emitted={} clamped={} ({:.1}%) near_cap={} ({:.1}%) arm_missing={} precap_too_big={} stageA_dups={} stageB_dups={} low_ident={}",
                    st.windows,
                    st.cands_total,
                    st.emitted,
                    st.clamped,
                    clamped_pct,
                    st.near_cap,
                    near_cap_pct,
                    st.arm_missing,
                    st.precap_too_big,
                    st.stage_a_dups,
                    st.stage_b_dups,
                    st.low_ident,
                );
            }
        }

        if let Some(html_path) = html.as_ref() {
            write_html_report(html_path, input, &report_rows)?;
        }

        if let Some(gff_path) = gff.as_ref() {
            let mut gw: Box<dyn Write> = if gff_is_stdout {
                Box::new(std::io::BufWriter::new(std::io::stdout()))
            } else {
                Box::new(std::io::BufWriter::new(std::fs::File::create(gff_path)?))
            };
            let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
            match write_gff3(&mut gw, &cmdline, &contigs_seen, &report_rows)
                .and_then(|_| gw.flush())
            {
                Ok(_) => {}
                Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
                Err(e) => return Err(e.into()),
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Mirrored matchpoints for an IR centred at `c`: one point every `step` bp
    /// along an arm whose inner edge is `inner` bp from the centre, `len` bp long.
    fn arm_points(c: i32, inner: i32, len: i32, step: i32) -> Vec<(i32, i32)> {
        let mut v = Vec::new();
        let mut u = inner;
        while u < inner + len {
            let (l, r) = (c - u, c + u - 17);
            v.push((l, r));
            v.push((r, l));
            u += step;
        }
        v
    }

    #[test]
    fn loci_join_irs_sharing_arm_sequence_transitively() {
        // IR0 and IR1 share a left arm; IR2 shares a right arm with IR1; IR3 is alone.
        let arms = [
            (1_000, 3_000, 10_000, 12_000),
            (1_500, 3_500, 20_000, 22_000),
            (30_000, 32_000, 21_000, 23_000),
            (50_000, 52_000, 60_000, 62_000),
        ];
        let (locus_of, sizes) = assign_loci(&arms);
        assert_eq!(locus_of[0], locus_of[1]);
        assert_eq!(locus_of[1], locus_of[2]);
        assert_ne!(locus_of[0], locus_of[3]);
        let mut s = sizes.clone();
        s.sort();
        assert_eq!(s, vec![1, 3]);
    }

    #[test]
    fn cluster_keeps_full_clean_arm() {
        // 20 kb arms, 5 kb spacer (inner edge 2.5 kb from centre), a point every 11 bp.
        let c = 100_000;
        let pts = arm_points(c, 2_500, 20_000, 11);
        let (la0, la1, ra0, ra1, kept) =
            arm_bounds_near_split_diag_clustered(&pts, c as usize, 17, 250_000, 120, 4, 1000, 10.0)
                .unwrap();
        assert_eq!(kept, pts.len());
        // Full extent, not trimmed: within one point spacing of 20 kb.
        assert!((la1 - la0) as i32 >= 20_000 - 11 && (ra1 - ra0) as i32 >= 20_000 - 11);
    }

    #[test]
    fn cluster_excludes_stray_matches_beyond_a_gap() {
        // Dense 5 kb arm plus a few stray matches 30 kb further out on the same anti-diagonal.
        let c = 100_000;
        let mut pts = arm_points(c, 1_000, 5_000, 11);
        pts.extend(arm_points(c, 36_000, 200, 50));
        let (la0, la1, _, _, _) =
            arm_bounds_near_split_diag_clustered(&pts, c as usize, 17, 250_000, 120, 4, 1000, 10.0)
                .unwrap();
        let arm = (la1 - la0) as i32;
        assert!(arm < 6_000, "stray matches inflated the arm to {arm} bp");
        assert!(arm >= 5_000 - 11, "core arm truncated to {arm} bp");
    }

    #[test]
    fn cluster_prefers_longer_sparse_arm_over_dense_short_copy() {
        // Genuine 2 kb arm at ~90% identity (a point every ~70 bp) vs a perfect
        // 300 bp stray copy (a point every 11 bp, more points in total) 20 kb out.
        let c = 100_000;
        let mut pts = arm_points(c, 5_000, 2_000, 70);
        pts.extend(arm_points(c, 27_000, 300, 5));
        let (la0, la1, _, _, _) =
            arm_bounds_near_split_diag_clustered(&pts, c as usize, 17, 250_000, 120, 4, 1000, 10.0)
                .unwrap();
        assert!((la1 - la0) as i32 >= 2_000 - 70, "picked the stray copy: {} bp", la1 - la0);
    }

    #[test]
    fn cluster_tolerates_sparse_divergent_arm() {
        // ~85%-identity arm: irregular gaps of 50-400 bp must not split it.
        let c = 100_000;
        let mut pts = Vec::new();
        let mut u = 2_000;
        let mut i = 0;
        while u < 2_000 + 8_000 {
            let (l, r) = (c - u, c + u - 17);
            pts.push((l, r));
            pts.push((r, l));
            u += [50, 120, 400, 90, 200, 60][i % 6];
            i += 1;
        }
        let (la0, la1, _, _, _) =
            arm_bounds_near_split_diag_clustered(&pts, c as usize, 17, 250_000, 120, 4, 1000, 10.0)
                .unwrap();
        assert!((la1 - la0) as i32 >= 7_500, "divergent arm split: {} bp", la1 - la0);
    }
}
