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
  - emit BED-like TSV rows\n\n\
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
  kept_pts       Number of matchpoints retained near refined anti-diagonal (0 if unfiltered fallback)\n\
  bin            Anti-diagonal bin id (debug / dedup key component)\n\n\
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
Higher values may find multiple IRs within the same window, at increased cost.",
                )
                .long("hits-per-window")
                .default_value("4")
                .value_parser(value_parser!(usize)),
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
            Arg::new("trim_lo")
                .help("Lower quantile for trimmed arm bounds (0..1).")
                .long_help(
                    "After filtering matchpoints near the split, we derive arm bounds using\n\
trimmed quantiles rather than min/max to avoid single-point inflation.\n\n\
Example: --trim-lo 0.10 keeps the 10th percentile as the lower bound.",
                )
                .long("trim-lo")
                .default_value("0.10")
                .value_parser(value_parser!(f32)),
        )
        .arg(
            Arg::new("trim_hi")
                .help("Upper quantile for trimmed arm bounds (0..1).")
                .long_help(
                    "Upper quantile used for trimmed arm bounds.\n\
Example: --trim-hi 0.90 uses the 90th percentile as the upper bound.",
                )
                .long("trim-hi")
                .default_value("0.90")
                .value_parser(value_parser!(f32)),
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

#[inline]
fn ir_class_from_arms(
    la0: usize,
    la1: usize,
    ra0: usize,
    ra1: usize,
) -> Option<(&'static str, isize)> {
    // require sane, ordered arms
    if la0 >= la1 || ra0 >= ra1 {
        return None;
    }
    // ensure left is left
    let (_la0, la1, ra0, _ra1) = normalize_arms(la0, la1, ra0, ra1);

    // spacer between arms (can be negative if overlapping)
    let spacer = ra0 as isize - la1 as isize;

    let class = if spacer <= 0 {
        "immediate"
    } else if spacer <= 10_000 {
        "spaced"
    } else {
        "wide"
    };

    Some((class, spacer))
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

/// Robust bounds near split: keep split-consistent points, then take trimmed quantiles
/// on left/right arms to avoid single-point inflation.
///
/// q_lo/q_hi are fractions in [0,1], e.g. 0.05 and 0.95.
///
/// Returns (la0, la1, ra0, ra1, kept_pts) in window-local coords.
fn arm_bounds_near_split_diag_trimmed(
    best_pts: &[(i32, i32)],
    split: usize,
    k: usize,
    win_len: usize,
    fold_diag_tol: i32,
    widen: i32,
    q_lo: f32,
    q_hi: f32,
) -> Option<(usize, usize, usize, usize, usize)> {
    if best_pts.is_empty() {
        return None;
    }

    let split2 = 2 * (split as i32);
    let tol = fold_diag_tol.max(1) * widen.max(1);

    let mut lefts: Vec<i32> = Vec::new();
    let mut rights: Vec<i32> = Vec::new();

    for &(p1, p2) in best_pts {
        let d = p1 + p2;
        if (d - split2).abs() > tol {
            continue;
        }
        lefts.push(p1.min(p2));
        rights.push(p1.max(p2));
    }

    let kept = lefts.len();
    if kept == 0 {
        return None;
    }

    lefts.sort_unstable();
    rights.sort_unstable();

    // Quantile indices (clamped)
    let n = kept;
    let lo = ((q_lo.clamp(0.0, 1.0) * (n.saturating_sub(1) as f32)).round() as usize).min(n - 1);
    let hi = ((q_hi.clamp(0.0, 1.0) * (n.saturating_sub(1) as f32)).round() as usize).min(n - 1);

    let l0 = lefts[lo];
    let l1 = lefts[hi];
    let r0 = rights[lo];
    let r1 = rights[hi];

    let la0 = l0.max(0) as usize;
    let la1 = (l1.max(0) as usize).saturating_add(k).min(win_len);
    let ra0 = r0.max(0) as usize;
    let ra1 = (r1.max(0) as usize).saturating_add(k).min(win_len);

    if la0 < la1 && ra0 < ra1 {
        Some((la0, la1, ra0, ra1, kept))
    } else {
        None
    }
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
    qs: usize,  // Stage B quantised interval start
    qe: usize,  // Stage B quantised interval end

    capped: bool,
    had_arms: bool,
    precap_too_big: bool,

    // ---- classification ----
    class: &'static str,
    spacer_bp: i32, // or isize/i64
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

/// Return maximal half-open ACGT-only spans in `seq`.
///
/// If the whole contig is valid, this returns exactly one span: `(0, seq.len())`.
fn valid_acgt_spans(seq: &[u8]) -> Vec<(usize, usize)> {
    let mut spans = Vec::new();
    let mut i = 0usize;
    let n = seq.len();

    while i < n {
        // skip invalid run
        while i < n && !is_valid_dna_byte(seq[i]) {
            i += 1;
        }
        if i >= n {
            break;
        }

        let start = i;

        // extend valid run
        while i < n && is_valid_dna_byte(seq[i]) {
            i += 1;
        }

        let end = i;
        if start < end {
            spans.push((start, end));
        }
    }

    spans
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
    let dedup_bp = *args.get_one::<usize>("dedup_bp").unwrap();
    let threads = *args.get_one::<usize>("threads").unwrap();
    let max_interval_bp = *args.get_one::<usize>("max_interval_bp").unwrap();
    let diag_widen = *args.get_one::<i32>("diag_widen").unwrap();
    let trim_lo = *args.get_one::<f32>("trim_lo").unwrap();
    let trim_hi = *args.get_one::<f32>("trim_hi").unwrap();
    let cap_policy = CapPolicy::parse(args.get_one::<String>("cap_policy").unwrap())?;

    // IR classification
    let classify_ir = args.get_flag("classify_ir");
    let immediate_bp = *args.get_one::<usize>("immediate_bp").unwrap();
    let wide_bp = *args.get_one::<usize>("wide_bp").unwrap();

    // HTML output
    let html = args.get_one::<PathBuf>("html").cloned();

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
                window: 500,
                arm: 1200,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.25,
                max_jump_clip: 1000,
            },
            fold_diag_tol: 120,
        },
        fold: FoldOnlyCfg { min_arm: 2000 },
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

        // can't have both outputs going to stdout
        if bed_is_stdout && fasta_is_stdout {
            anyhow::bail!(
                "Cannot write both TSV (--bed -) and FASTA (--fasta -) to stdout. Send one of them to a file."
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
                "#contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\ttir_ident\tmatches\tspan\tla0\tla1\tra0\tra1\tcontig_len\twin_start\twin_end\tkept_pts\tbin\tarm_len\tspacer\tir_class"
            )
        } else {
            writeln!(
                out,
                "#contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\ttir_ident\tmatches\tspan\tla0\tla1\tra0\tra1\tcontig_len\twin_start\twin_end\tkept_pts\tbin"
            )
        };

        match header {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
            Err(e) => return Err(e.into()),
        }

        let mut report_rows: Vec<ReportRow> = Vec::new();

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

            let q = dedup_bp.max(1);
            let mm_k = cfg.shared.minimizer.k;
            let fold_diag_tol = cfg.shared.fold_diag_tol;

            // Split contig into maximal valid ACGT-only spans.
            // If the contig has no invalid bases, this is just [(0, len)].
            let spans = valid_acgt_spans(&seq);
            if spans.len() > 1
                || spans
                    .first()
                    .map(|&(s, e)| s != 0 || e != len)
                    .unwrap_or(false)
            {
                if let Some((pos, bad)) = first_invalid_dna_byte(&seq) {
                    elog!(
                        "irx",
                        "{id}: splitting contig on invalid base '{}' (ASCII {}) first seen at offset {} into {} valid spans",
                        bad as char,
                        bad,
                        pos,
                        spans.len(),
                    );
                }
            }

            if spans.is_empty() {
                elog!("irx", "{id}: no valid ACGT spans; skipping contig");
                continue;
            }

            // Since IRs currently require min_arm on each side, tiny valid fragments are not useful.
            // This also avoids wasting time on tiny islands between Ns.
            let min_span_len = cfg.fold.min_arm.saturating_mul(2);

            // Build contig-global windows from valid spans.
            let mut wins: Vec<(usize, usize)> = Vec::new();
            let mut had_invalid = false;

            for &(span_start, span_end) in &spans {
                if span_start > 0 || span_end < len {
                    had_invalid = true;
                }

                let span_len = span_end - span_start;
                if span_len < min_span_len {
                    continue;
                }

                for (rel_start, rel_end) in windows(span_len, window_len, step) {
                    wins.push((span_start + rel_start, span_start + rel_end));
                }
            }

            if wins.is_empty() {
                if had_invalid {
                    elog!(
                        "irx",
                        "{id}: contig contains invalid bases; no valid spans >= {} bp after splitting",
                        min_span_len
                    );
                }
                continue;
            }

            // For each window:
            // - detect top-k coarse foldback bins
            // - refine split position per hit
            // - compute approximate arm interval bounds
            // - package into `Cand` records
            //
            // We do *no* dedup here (that happens serially per contig).

            let mut cands: Vec<Cand> = wins
                .par_iter()
                .flat_map_iter(|&(wstart, wend)| {
                    let sub = &seq[wstart..wend];
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
                            fb.span,
                        )
                        .ok()
                        .flatten() else {
                            continue;
                        };

                        // Global break position (contig coords).
                        let break_pos = wstart + rf.split_pos;

                        // Determine arm bounds:
                        //   1) prefer points near refined split (anti-diagonal proximity)
                        //   2) fall back to unfiltered bounds
                        //   3) fall back to break±span/2
                        let (start, end, la0g, la1g, ra0g, ra1g, kept_pts) =
                            if let Some((la0, la1, ra0, ra1, kept)) =
                                arm_bounds_near_split_diag_trimmed(
                                    &hit.pts,
                                    rf.split_pos,
                                    mm_k,
                                    win_len,
                                    fold_diag_tol,
                                    diag_widen,
                                    trim_lo,
                                    trim_hi,
                                )
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

                        // Quantised interval for Stage B overlap dedup.
                        let qs: usize = (start / q) * q;
                        let qe: usize = (end / q) * q;

                        let (class, spacer_bp) = if la1g > la0g && ra1g > ra0g {
                            if let Some((cls, sp)) = ir_class_from_arms(la0g, la1g, ra0g, ra1g) {
                                (cls, sp.clamp(i32::MIN as isize, i32::MAX as isize) as i32)
                            } else {
                                ("unknown", 0)
                            }
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
                            qs,
                            qe,
                            capped: capped_this,
                            had_arms,
                            precap_too_big,
                            class,
                            spacer_bp,
                        });
                    }

                    out_local.into_iter()
                })
                .collect();

            let mut st = Stats::default();
            st.windows = wins.len() as u64;
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

            // Stage A: cheap dedup on break/bin.
            let mut seen_break_bins: HashSet<u64> = HashSet::new();
            // Stage B: greedy overlap dedup on quantised intervals.
            let mut accepted_intervals: Vec<(usize, usize)> = Vec::new();

            let cap_thresh = (max_interval_bp as f64 * 0.95).round() as usize; // “near cap” threshold

            let mut rc_buf: Vec<u8> = Vec::new();

            for c in cands {
                // A: check for duplicate break/bin key. This is very fast and removes exact duplicates and many near-duplicates.
                if !seen_break_bins.insert(c.key_a) {
                    st.stage_a_dups += 1;
                    continue;
                }
                st.stage_a_kept += 1;

                // B: check for overlap with accepted intervals. This is O(N^2) in the worst case but should be fine for small N.
                let mut overlap = false;
                for &(s, e) in &accepted_intervals {
                    if !(c.qe <= s || c.qs >= e) {
                        overlap = true;
                        break;
                    }
                }
                if overlap {
                    st.stage_b_dups += 1;
                    continue;
                }
                st.stage_b_kept += 1;

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

                accepted_intervals.push((c.qs, c.qe));

                let (la0n, la1n, ra0n, ra1n) = normalize_arms(c.la0, c.la1, c.ra0, c.ra1);
                let arm_len = (la1n.saturating_sub(la0n)).min(ra1n.saturating_sub(ra0n));
                let spacer: isize = (ra0n as isize) - (la1n as isize);
                let ir_class = classify_ir_str(spacer, immediate_bp, wide_bp);
                let tir_ident =
                    if la1n > la0n && ra1n > ra0n && la1n <= seq.len() && ra1n <= seq.len() {
                        let left = &seq[la0n..la1n];
                        let right = &seq[ra0n..ra1n];

                        // Reuse a preallocated reverse-complement buffer.
                        revcomp_into(&mut rc_buf, right);

                        // Fast path
                        let proxy =
                            arm_identity_proxy_blockalign(left, &rc_buf).unwrap_or(f32::NAN);

                        // Only do expensive traceback for promising candidates
                        if proxy.is_finite() && proxy >= 0.65 {
                            arm_identity_trace_blockalign(left, &rc_buf).unwrap_or(proxy)
                        } else {
                            proxy
                        }
                    } else {
                        f32::NAN
                    };

                let left_arm_len = la1n.saturating_sub(la0n);
                let right_arm_len = ra1n.saturating_sub(ra0n);
                let arm_ratio = if left_arm_len > 0 && right_arm_len > 0 {
                    let a = left_arm_len.max(right_arm_len) as f32;
                    let b = left_arm_len.min(right_arm_len) as f32;
                    a / b
                } else {
                    f32::NAN
                };

                report_rows.push(ReportRow {
                    contig: id.to_string(),
                    start: c.start,
                    end: c.end,
                    break_pos: c.break_pos,
                    identity_est: c.ident,
                    tir_ident: if tir_ident.is_finite() {
                        Some(tir_ident)
                    } else {
                        None
                    },
                    matches: c.matches,
                    span: c.span,
                    la0: c.la0,
                    la1: c.la1,
                    ra0: c.ra0,
                    ra1: c.ra1,
                    contig_len: len,
                    win_start: c.win_start,
                    win_end: c.win_end,
                    kept_pts: c.kept_pts,
                    bin: c.bin,
                    arm_len,
                    spacer,
                    ir_class: ir_class.to_string(),
                    interval_len: c.end.saturating_sub(c.start),
                    left_arm_len,
                    right_arm_len,
                    arm_ratio,
                });

                let row = if classify_ir {
                    writeln!(
                        out,
                        "{id}\t{}\t{}\tIR\t{}\t.\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
                        arm_len,
                        spacer,
                        ir_class
                    )
                } else {
                    writeln!(
                        out,
                        "{id}\t{}\t{}\tIR\t{}\t.\t{}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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

                // emit TSV
                match row {
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
            }

            if st.emitted > 0 {
                let clamped_pct = 100.0 * (st.clamped as f64) / (st.emitted as f64);
                let near_cap_pct = 100.0 * (st.near_cap as f64) / (st.emitted as f64);

                elog!(
                    "irx",
                    "{id}: wins={} cands={} emitted={} clamped={} ({:.1}%) near_cap={} ({:.1}%) arm_missing={} precap_too_big={} stageA_dups={} stageB_dups={}",
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
                );
            }
        }

        if let Some(html_path) = html.as_ref() {
            write_html_report(html_path, input, &report_rows)?;
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
