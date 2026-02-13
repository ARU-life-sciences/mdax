use anyhow::Result;
use clap::{Arg, Command, value_parser};
use std::collections::HashSet;
use std::{io::Write, path::PathBuf};

use mdax::{
    cfg::{FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg, SharedCfg, SigCfg},
    foldback, io,
    scratch::FoldScratch,
    utils::RefineMode,
};

fn build_cli() -> clap::ArgMatches {
    Command::new("irx")
        .version(clap::crate_version!())
        .about("Detect candidate inverted repeats in assemblies using mdax foldback logic")
        .arg(
            Arg::new("input")
                .help("Input assembly FASTA(.gz)")
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("bed")
                .help("Output BED/TSV file ('-' for stdout)")
                .long("bed")
                .short('o')
                .default_value("-")
                .value_parser(value_parser!(PathBuf)),
        )
        .arg(
            Arg::new("window_len")
                .help("Window length for synthetic reads (0 = whole contig)")
                .long("window-len")
                .default_value("250000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("step")
                .help("Step size between windows")
                .long("step")
                .default_value("50000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("hits_per_window")
                .help("How many foldback candidates to emit per window (top-k bins)")
                .long("hits-per-window")
                .default_value("4")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_ident")
                .help("Minimum identity of the left/right arms of the putative inverted repeat.")
                .long("min-ident")
                .default_value("0.2")
                .value_parser(value_parser!(f32))
        )
        .arg(
            Arg::new("dedup_bp")
                .help("Quantise starts/ends for faster dedup (bp). Used by the interval-overlap dedupper.")
                .long("dedup-bp")
                .default_value("200")
                .value_parser(value_parser!(usize)),
        )
        .get_matches()
}

/// Generate (start,end) windows over a contig. If window_len==0, yields the whole contig once.
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
/// For a foldback, matchpoints satisfy roughly `p1 + p2 ≈ 2*split`.
///
/// We keep points within `tol = fold_diag_tol * widen` of `2*split` and then
/// take min/max on left/right arms.
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

#[inline]
fn pack_u32_pair(a: u32, b: u32) -> u64 {
    ((a as u64) << 32) | (b as u64)
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

    // Minimal cfg: reuse mdax detector + refinement.
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

    let mut out: Box<dyn Write> = if bed.to_string_lossy() == "-" {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    } else {
        Box::new(std::io::BufWriter::new(std::fs::File::create(bed)?))
    };

    // TSV header (BED6+extras); keep as comment so BED tools can ignore.
    match writeln!(
        out,
        "#contig\tstart\tend\tname\tscore\tstrand\tbreak_pos\tidentity_est\tmatches\tspan\tla0\tla1\tra0\tra1\twin_start\twin_end\tkept_pts\tbin"
    ) {
        Ok(_) => {}
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
        Err(e) => return Err(e.into()),
    }

    let mut r = io::open_fasta_reader(input)?;
    let mut fold_scratch = FoldScratch::new();

    // Per-contig dedup: keep emitted intervals (quantised) and skip overlaps.
    while let Some(rec) = r.next() {
        let rec = rec?;
        let id = first_token_utf8(rec.id());
        let seq = rec.seq();
        let len = seq.len();

        // cheap dedup on breakpoint bins (very fast).
        // Use u64 key to avoid tuple hashing overhead.
        let mut seen_break_bins: HashSet<u64> = HashSet::new();

        // Stage B: overlap dedup on quantised intervals (more expensive).
        let mut seen_intervals_q: Vec<(usize, usize)> = Vec::new();

        for (wstart, wend) in windows(len, window_len, step) {
            let sub = &seq[wstart..wend];
            let win_len = sub.len();
            let k = cfg.shared.minimizer.k;

            // MULTI-HIT: top-k bins per window
            let hits = foldback::detect_foldbacks_topk(
                sub,
                &cfg.shared,
                &cfg.fold,
                &mut fold_scratch,
                hits_per_window,
            );

            if hits.is_empty() {
                continue;
            }

            for hit in hits {
                let fb = hit.fb;
                let bin = hit.bin;

                let Some(rf) = foldback::refine_breakpoint(
                    sub,
                    fb.split_pos,
                    &cfg.shared,
                    &mut fold_scratch.refine,
                )?
                else {
                    continue;
                };

                let break_pos = wstart + rf.split_pos;

                // 1) prefer split-consistent bounds from hit.pts
                // 2) fallback to unfiltered bounds from hit.pts
                // 3) last resort: break±span/2
                let (start, end, la0g, la1g, ra0g, ra1g, kept_pts) =
                    if let Some((la0, la1, ra0, ra1, kept)) = arm_bounds_near_split_diag(
                        &hit.pts,
                        rf.split_pos,
                        k,
                        win_len,
                        cfg.shared.fold_diag_tol,
                        20, // widen factor; 20*120=2400bp tolerance
                    ) {
                        let la0g = wstart + la0;
                        let la1g = wstart + la1;
                        let ra0g = wstart + ra0;
                        let ra1g = wstart + ra1;
                        let start = la0g.min(ra0g);
                        let end = la1g.max(ra1g);
                        (start, end, la0g, la1g, ra0g, ra1g, kept)
                    } else if let Some((la0, la1, ra0, ra1)) =
                        arm_bounds_from_best_pts(&hit.pts, k, win_len)
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

                // ------------------------------
                // Two-stage dedup
                // ------------------------------
                let q = dedup_bp.max(1);

                // TODO: check, is this where to reserve?
                seen_break_bins.reserve((len / q.max(1)).min(1_000_000));
                seen_intervals_q.reserve(1024);

                // Stage A: dedup by breakpoint bin (cheap).
                let break_bin = (break_pos / q) as u32;
                let diag_bin = bin as i32; // from hit.bin; may be negative in principle
                let diag_bin_u32 = diag_bin as u32; // stable bit-pattern

                // Key: (break_bin, diag_bin) packed into u64
                // If you want *more* aggressive merging, drop diag_bin from the key and only use break_bin.
                let key_a = pack_u32_pair(break_bin, diag_bin_u32);

                if !seen_break_bins.insert(key_a) {
                    continue; // already emitted something for this (break_bin, diag_bin)
                }

                // Stage B: overlap dedup on quantised intervals (more expensive).
                let qs = (start / q) * q;
                let qe = (end / q) * q;

                let mut duplicate = false;
                for &(s, e) in &seen_intervals_q {
                    if !(qe <= s || qs >= e) {
                        duplicate = true;
                        break;
                    }
                }
                if duplicate {
                    continue;
                }
                seen_intervals_q.push((qs, qe));

                // skip if identity is too low
                if rf.identity_est < min_ident {
                    continue;
                }
                // TODO: why do we need this?
                let score = (rf.identity_est.clamp(0.0, 1.0) * 1000.0).round() as i32;

                match writeln!(
                    out,
                    "{id}\t{start}\t{end}\tIR\t{score}\t.\t{break_pos}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    rf.identity_est,
                    fb.matches,
                    fb.span,
                    la0g,
                    la1g,
                    ra0g,
                    ra1g,
                    wstart,
                    wend,
                    kept_pts,
                    bin,
                ) {
                    Ok(_) => {}
                    Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => return Ok(()),
                    Err(e) => return Err(e.into()),
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arm_bounds_near_split_diag_keeps_consistent_points() {
        // p1+p2 == 2000 for split=1000
        let best_pts = vec![(900, 1100), (850, 1150), (800, 1200)];
        let out = arm_bounds_near_split_diag(&best_pts, 1000, 17, 5000, 120, 20);
        assert!(out.is_some());
        let (_la0, _la1, _ra0, _ra1, kept) = out.unwrap();
        assert_eq!(kept, 3);
    }

    #[test]
    fn test_arm_bounds_near_split_diag_rejects_off_diag_points() {
        // One good, one bad (sum far from 2*split)
        let best_pts = vec![(900, 1100), (0, 4000)];
        let out = arm_bounds_near_split_diag(&best_pts, 1000, 17, 5000, 120, 1);
        assert!(out.is_some());
        let (_la0, _la1, _ra0, _ra1, kept) = out.unwrap();
        assert_eq!(kept, 1);
    }

    #[test]
    fn test_windows_generation() {
        let wins: Vec<_> = windows(1000, 200, 100).collect();
        assert!(wins.len() > 1);
        assert_eq!(wins[0], (0, 200));
    }

    #[test]
    fn test_interval_overlap_dedup() {
        let first = (100, 200);
        let second = (150, 250); // overlaps
        let third = (300, 400); // no overlap
        assert!(!(second.1 <= first.0 || second.0 >= first.1));
        assert!(third.1 <= first.0 || third.0 >= first.1);
    }
}
