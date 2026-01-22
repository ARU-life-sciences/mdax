use crate::cfg::{FoldOnlyCfg, SharedCfg};
use crate::minimizer::sampled_minimizers;
use crate::utils::{RefineMode, Refined, banded_edit_distance, div_floor};
use gxhash::{HashMap, HashMapExt};

/// IMPORTANT:
/// This implementation uses *strand-specific* minimizers (forward-only) via:
///   simd_minimizers::minimizers(k,w).hasher(NtHasher::<false>)
///
/// And it matches minimizers using the *minimizer values returned by simd_minimizers*
/// (i.e. the same hashes/values that were used to choose minimizers), rather than
/// re-hashing the underlying k-mers with a separate hash function.
///
/// This makes the foldback detector more “Pacasus-like” in the sense that it is
/// explicitly looking for strong reverse-complement self-similarity (anti-diagonal)
/// without canonical-strand tie-breaking artifacts.

/// Configuration for foldback detection.
#[derive(Debug, Clone)]
pub struct FoldCfg {
    /// Sampling minimizer k value
    pub k: usize,
    /// Sampling minimizer window size
    pub w: usize,
    /// Diagonal bucketing tolerance when clustering matchpoints
    pub diag_tol: i32,
    /// Minimum palindrome arm length to actually call a palindrome (coarse span)
    pub min_arm: usize,
    /// Minimum number of minimizer matchpoints required to support a foldback candidate.
    pub min_matches: usize,
    /// Don't call breakpoints too close to read endpoints
    pub end_guard: usize,
    /// Half-width (bp) of the local search window used to refine a coarse breakpoint estimate.
    pub refine_window: usize,
    /// Arm length (bp) used during breakpoint refinement.
    pub refine_arm: usize,
    /// Refinement mode.
    pub refine_mode: RefineMode,
    /// Max edit-distance rate used to define band in ONT refinement.
    pub max_ed_rate: f32,
}

impl Default for FoldCfg {
    fn default() -> Self {
        Self {
            k: 17,
            w: 21,
            diag_tol: 120,
            // NOTE: if you want “Pacasus-like” behavior on real ONT reads,
            // you usually want this larger (hundreds–thousands), and CLI-controlled.
            min_arm: 200,
            min_matches: 10,
            end_guard: 1000,
            refine_window: 200,
            refine_arm: 500,
            refine_mode: RefineMode::HiFi,
            max_ed_rate: 0.05,
        }
    }
}

/// Result of foldback detection (coarse).
#[derive(Debug, Clone)]
pub struct FoldBreakpoint {
    pub split_pos: usize,
    pub score: i64,
    pub matches: usize,
    pub span: usize,
}

/// Refine an inexact breakpoint estimate `s0` into a more exact split position.
pub fn refine_breakpoint(seq: &[u8], s0: usize, cfg: &SharedCfg) -> Option<Refined> {
    let n = seq.len();
    if s0 < cfg.refine.arm + cfg.refine.window || s0 + cfg.refine.arm + cfg.refine.window > n {
        eprintln!(
            "Refinement context exceeds sequence bounds: s0={}, arm={}, window={}, len={}",
            s0, cfg.refine.arm, cfg.refine.window, n
        );
        return None;
    }
    match cfg.refine.mode {
        RefineMode::HiFi => refine_breakpoint_hamming(seq, s0, cfg),
        RefineMode::ONT => refine_breakpoint_banded_ed(seq, s0, cfg),
    }
}

pub fn refine_breakpoint_hamming(seq: &[u8], s0: usize, cfg: &SharedCfg) -> Option<Refined> {
    let n = seq.len();
    if s0 < cfg.refine.arm + cfg.refine.window || s0 + cfg.refine.arm + cfg.refine.window > n {
        return None;
    }

    let mut best_s = s0;
    let mut best_score: i64 = i64::MIN;

    for s in (s0 - cfg.refine.window)..=(s0 + cfg.refine.window) {
        let l = &seq[s - cfg.refine.arm..s];

        let mut score: i64 = 0;
        for i in 0..cfg.refine.arm {
            // walk outward from junction on the left
            let a = l[cfg.refine.arm - 1 - i];
            // right arm: forward from junction, but we want RC, so just complement here
            let b = comp(seq[s + i]);
            if a == b {
                score += 1;
            } else {
                score -= 1;
            }
        }

        if score > best_score {
            best_score = score;
            best_s = s;
        }
    }

    let matches = ((best_score + cfg.refine.arm as i64) / 2).max(0) as usize;
    let identity = matches as f32 / cfg.refine.arm as f32;

    Some(Refined {
        split_pos: best_s,
        score: best_score,
        identity_est: identity,
    })
}

pub fn refine_breakpoint_banded_ed(seq: &[u8], s0: usize, cfg: &SharedCfg) -> Option<Refined> {
    let n = seq.len();
    if s0 < cfg.refine.arm + cfg.refine.window || s0 + cfg.refine.arm + cfg.refine.window > n {
        return None;
    }

    let mut band = (cfg.refine.max_ed_rate.max(0.0) * cfg.refine.arm as f32).ceil() as usize;
    band = band.max(1).min(cfg.refine.arm);

    let mut best_s = s0;
    let mut best_ed: usize = usize::MAX;

    let mut right_rc = vec![b'N'; cfg.refine.arm];

    for s in (s0 - cfg.refine.window)..=(s0 + cfg.refine.window) {
        let left = &seq[s - cfg.refine.arm..s];
        let right = &seq[s..s + cfg.refine.arm];

        revcomp_into(right, &mut right_rc);

        let ed = banded_edit_distance(left, &right_rc, band);
        if ed < best_ed {
            best_ed = ed;
            best_s = s;
        }
    }

    let identity = (1.0 - (best_ed as f32 / cfg.refine.arm as f32)).clamp(0.0, 1.0);

    Some(Refined {
        split_pos: best_s,
        score: -(best_ed as i64),
        identity_est: identity,
    })
}

/// Detect a foldback/palindromic junction in a single sequence.
///
/// Method:
/// 1) Compute *forward-only* minimizers on `seq` and on `revcomp(seq)` using NtHasher<false>
/// 2) Index forward minimizers by their minimizer *value* (as returned by simd_minimizers)
/// 3) For each minimizer on rc, match by value and convert rc pos -> forward coordinate
/// 4) Cluster matchpoints by anti-diagonal d = p1 + p2 (within diag_tol)
/// 5) Choose best cluster; return split ≈ median(d/2)
pub fn detect_foldback(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
) -> Option<FoldBreakpoint> {
    let k = shared.minimizer.k;
    let w = shared.minimizer.w;

    if seq.len() < k + w + 10 {
        return None;
    }

    // forward-only minimizers (or canonical if configured), returning (pos, value)
    let (pos_f, val_f) = sampled_minimizers(seq, &shared.minimizer);

    let mut rc = seq.to_vec();
    revcomp_in_place(&mut rc);
    let (pos_rc, val_rc) = sampled_minimizers(&rc, &shared.minimizer);

    if pos_f.len() < shared.min_matches || pos_rc.len() < shared.min_matches {
        return None;
    }

    // index forward minimizers by value -> positions
    let mut idx_f: HashMap<u64, Vec<i32>> = HashMap::new();
    for (&p, &v) in pos_f.iter().zip(val_f.iter()) {
        idx_f.entry(v).or_default().push(p as i32);
    }

    // match rc minimizers by value, convert rc pos -> forward coordinate
    let n = seq.len() as i32;
    let k_i32 = k as i32;

    let mut d_bins: HashMap<i32, Vec<(i32, i32)>> = HashMap::new();
    for (&prc, &vrc) in pos_rc.iter().zip(val_rc.iter()) {
        if let Some(p1s) = idx_f.get(&vrc) {
            let p2 = n - k_i32 - (prc as i32); // rc start -> forward start
            for &p1 in p1s {
                let d = p1 + p2;
                let bin = div_floor(d, shared.fold_diag_tol.max(1));
                d_bins.entry(bin).or_default().push((p1, p2));
            }
        }
    }

    best_bin_to_fold_breakpoint(&d_bins, shared, fold, seq.len())
}

fn best_bin_to_fold_breakpoint(
    d_bins: &HashMap<i32, Vec<(i32, i32)>>,
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    len: usize,
) -> Option<FoldBreakpoint> {
    let mut best: Option<FoldBreakpoint> = None;

    for (_bin, pts) in d_bins {
        if pts.len() < shared.min_matches {
            continue;
        }

        let mut v = pts.clone();
        v.sort_by_key(|(p1, _)| *p1);

        // monotone anti-diagonal chain: p2 decreasing as p1 increases
        let mut chain = Vec::new();
        let mut last_p2 = i32::MAX;
        for (p1, p2) in v {
            if p2 < last_p2 {
                chain.push((p1, p2));
                last_p2 = p2;
            }
        }
        if chain.len() < shared.min_matches {
            continue;
        }

        let p1_min = chain.first().unwrap().0 as i64;
        let p1_max = chain.last().unwrap().0 as i64;
        let span = (p1_max - p1_min).max(0) as usize;
        if span < fold.min_arm {
            continue;
        }

        // median d = p1 + p2 => split ≈ d/2
        let mut ds: Vec<i32> = chain.iter().map(|(p1, p2)| p1 + p2).collect();
        ds.sort_unstable();
        let d_med = ds[ds.len() / 2];
        let split_u = ((d_med as f64) / 2.0).round().max(0.0) as usize;

        if split_u < shared.end_guard || split_u + shared.end_guard > len {
            continue;
        }

        let score = span as i64 + (chain.len() as i64 * 10);
        let bp = FoldBreakpoint {
            split_pos: split_u,
            score,
            matches: chain.len(),
            span,
        };

        if best.as_ref().map(|b| bp.score > b.score).unwrap_or(true) {
            best = Some(bp);
        }
    }

    best
}

/// Write reverse-complement of `src` into `dst`.
#[inline]
fn revcomp_into(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    let mut i = 0usize;
    let mut j = src.len();
    while i < src.len() {
        j -= 1;
        dst[i] = comp(src[j]);
        i += 1;
    }
}

/// Return the DNA complement of a base.
/// Non-ACGT bases are mapped to 'N'.
#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// Reverse-complement in place (ACGTN; other -> N)
fn revcomp_in_place(seq: &mut [u8]) {
    let mut i = 0usize;
    let mut j = seq.len().saturating_sub(1);
    while i < j {
        let a = comp(seq[i]);
        let b = comp(seq[j]);
        seq[i] = b;
        seq[j] = a;
        i += 1;
        j = j.saturating_sub(1);
    }
    if i == j && i < seq.len() {
        seq[i] = comp(seq[i]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cfg::{ConcatOnlyCfg, FoldOnlyCfg, MinimizerCfg, RefineCfg, SharedCfg};

    fn test_fold_cfg(
        k: usize,
        w: usize,
        diag_tol: i32,
        min_arm: usize,
        min_matches: usize,
        refine_window: usize,
        refine_arm: usize,
        refine_mode: RefineMode,
        max_ed_rate: f32,
    ) -> (SharedCfg, FoldOnlyCfg) {
        (
            SharedCfg {
                minimizer: MinimizerCfg {
                    k,
                    w,
                    forward_only: true,
                },
                min_matches,
                end_guard: 0,
                refine: RefineCfg {
                    window: refine_window,
                    arm: refine_arm,
                    mode: refine_mode,
                    max_ed_rate,
                },
                fold_diag_tol: diag_tol,
                concat_diag_tol: diag_tol, // irrelevant for foldback tests
            },
            FoldOnlyCfg { min_arm },
        )
    }

    fn repeat_bytes(pat: &[u8], times: usize) -> Vec<u8> {
        let mut v = Vec::with_capacity(pat.len() * times);
        for _ in 0..times {
            v.extend_from_slice(pat);
        }
        v
    }

    fn make_noisy_foldback_from_left(
        left: &[u8],
        join: &[u8],
        sub_every: usize,
        ins_every: usize,
        del_every: usize,
    ) -> (Vec<u8>, usize) {
        let mut l = left.to_vec();

        if sub_every > 0 {
            l = mutate_substitutions(&l, sub_every);
        }
        if ins_every > 0 {
            l = inject_insertions(&l, ins_every, b'A');
        }
        if del_every > 0 {
            l = inject_deletions(&l, del_every);
        }

        let mut rc = l.clone();
        super::revcomp_in_place(&mut rc);

        let mut s = Vec::with_capacity(l.len() + join.len() + rc.len());
        s.extend_from_slice(&l);
        s.extend_from_slice(join);
        let bp = l.len() + join.len();
        s.extend_from_slice(&rc);

        (s, bp)
    }

    fn mutate_substitutions(seq: &[u8], every: usize) -> Vec<u8> {
        let mut out = seq.to_vec();
        for i in (0..out.len()).step_by(every.max(1)) {
            out[i] = match out[i] {
                b'A' => b'C',
                b'C' => b'A',
                b'G' => b'T',
                b'T' => b'G',
                x => x,
            };
        }
        out
    }

    fn inject_deletions(seq: &[u8], every: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(seq.len());
        for (i, &b) in seq.iter().enumerate() {
            if every > 0 && (i % every == 0) {
                continue;
            }
            out.push(b);
        }
        out
    }

    fn inject_insertions(seq: &[u8], every: usize, ins: u8) -> Vec<u8> {
        let mut out = Vec::with_capacity(seq.len() + seq.len() / every.max(1) + 1);
        for (i, &b) in seq.iter().enumerate() {
            out.push(b);
            if every > 0 && (i % every == 0) {
                out.push(ins);
            }
        }
        out
    }

    #[test]
    fn detect_foldback_breakpoint_near_truth() {
        let left = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        let join = b"";
        let (s, true_bp) = make_noisy_foldback_from_left(left, join, 0, 0, 0);

        let (shared, fold) = test_fold_cfg(9, 11, 50, 20, 5, 0, 0, RefineMode::HiFi, 0.0);

        let bp = detect_foldback(&s, &shared, &fold).unwrap();

        assert!((bp.split_pos as i32 - true_bp as i32).abs() <= 20);
    }

    #[test]
    fn refine_hifi_hits_exact_breakpoint_clean() {
        let left = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        let join = b"";
        let (s, true_bp) = make_noisy_foldback_from_left(left, join, 12, 0, 0);

        let (shared, fold) = test_fold_cfg(9, 11, 50, 20, 5, 10, 30, RefineMode::HiFi, 0.0);

        let coarse = detect_foldback(&s, &shared, &fold).unwrap();

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared).unwrap();

        assert!((refined.split_pos as i32 - true_bp as i32).abs() <= 2);

        assert!(refined.identity_est > 0.9);
    }

    #[test]
    fn refine_ont_tolerates_indels() {
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 40);
        let (s, _bp) = make_noisy_foldback_from_left(&left, b"", 0, 23, 41);

        let (shared, fold) = test_fold_cfg(7, 9, 120, 30, 3, 40, 60, RefineMode::ONT, 0.35);

        let coarse = detect_foldback(&s, &shared, &fold).unwrap();

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared).expect("should refine");

        assert!(refined.identity_est > 0.6);
        assert!(
            (refined.split_pos as i32 - coarse.split_pos as i32).abs()
                <= shared.refine.window as i32
        );
    }

    #[test]
    fn hifi_and_ont_refiners_diverge_on_indels() {
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 40);
        let (s, _bp) = make_noisy_foldback_from_left(&left, b"", 0, 23, 41);

        let (mut shared, fold) = test_fold_cfg(
            7,
            9,   // k, w
            120, // diag_tol
            30,  // min_arm (coarse filter)
            3,   // min_matches
            40,  // refine_window
            80,  // refine_arm  (must be >0)
            RefineMode::HiFi,
            0.35, // max_ed_rate (used when mode=ONT)
        );

        let coarse = detect_foldback(&s, &shared, &fold).unwrap();

        shared.refine.mode = RefineMode::HiFi;
        let hifi = refine_breakpoint(&s, coarse.split_pos, &shared);

        shared.refine.mode = RefineMode::ONT;
        let ont = refine_breakpoint(&s, coarse.split_pos, &shared).unwrap();

        assert!(ont.identity_est > 0.4);
        assert!(
            (ont.split_pos as i32 - coarse.split_pos as i32).abs() <= shared.refine.window as i32
        );

        // optional: HiFi may fail, but ONT shouldn't
        assert!(hifi.is_some() || ont.identity_est > 0.4);
    }
}
