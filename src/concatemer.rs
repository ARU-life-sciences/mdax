use gxhash::{HashMap, HashMapExt};

use crate::cfg::{ConcatOnlyCfg, SharedCfg};
use crate::minimizer::sampled_minimizers_into;
use crate::utils::{RefineMode, Refined, banded_edit_distance, div_floor};

// NOTE: this code is very similar to `foldback.rs`
// but kept distinct in case we want to delete later

#[derive(Debug, Clone)]
pub struct ConcatBreakpoint {
    pub split_pos: usize, // coarse junction estimate
    pub delta: usize,     // repeat offset (median p2 - p1)
    pub score: i64,
    pub matches: usize,
    pub span: usize,

    // diagnostics
    pub cross_frac: f32,
    pub p2_span: usize,
}

pub fn refine_concatemer_breakpoint(
    seq: &[u8],
    s0: usize,
    delta: usize,
    cfg: &SharedCfg,
) -> Option<Refined> {
    if delta == 0 {
        return None;
    }
    match cfg.refine.mode {
        RefineMode::HiFi => refine_concat_hamming(seq, s0, delta, cfg),
        RefineMode::ONT => refine_concat_banded_ed(seq, s0, delta, cfg),
    }
}

/// Refine breakpoint using indel-tolerant banded edit distance.
///
/// For each candidate split position `s` in `[s0-window, s0+window]`:
/// - left arm  = seq[s-arm .. s]
/// - right arm = seq[s .. s+arm]
/// We compute the banded Levenshtein distance between the two arms and choose
/// the split with the smallest distance.
///
/// `max_ed_rate` controls the bandwidth: band = ceil(max_ed_rate * arm).
/// Typical ONT values: 0.15–0.25 depending on error rate.
///
/// Returns `None` if the refinement window/arms exceed sequence bounds.
pub fn refine_concat_banded_ed(
    seq: &[u8],
    s0: usize,
    delta: usize,
    cfg: &SharedCfg,
) -> Option<Refined> {
    let n = seq.len();
    if n == 0 || cfg.refine.arm == 0 {
        return None;
    }

    let mut band = ((cfg.refine.max_ed_rate.max(0.0)) * cfg.refine.arm as f32).ceil() as usize;
    band = band.max(1).min(cfg.refine.arm);

    let lo = s0.saturating_sub(cfg.refine.window);
    let hi = (s0 + cfg.refine.window).min(n);

    let mut best_s: usize = 0;
    let mut best_ed: usize = usize::MAX;
    let mut any = false;

    for s in lo..=hi {
        if s < delta {
            continue;
        }
        let a0 = s;
        let b0 = s - delta;

        if a0 + cfg.refine.arm > n {
            continue;
        }
        if b0 + cfg.refine.arm > n {
            continue;
        }

        let a = &seq[a0..a0 + cfg.refine.arm];
        let b = &seq[b0..b0 + cfg.refine.arm];

        let ed = banded_edit_distance(a, b, band);
        if ed < best_ed {
            best_ed = ed;
            best_s = s;
        }
        any = true;
    }

    if !any {
        return None;
    }

    let identity = (1.0 - (best_ed as f32 / cfg.refine.arm as f32)).clamp(0.0, 1.0);
    Some(Refined {
        split_pos: best_s,
        score: -(best_ed as i64),
        identity_est: identity,
    })
}

fn refine_concat_hamming(seq: &[u8], s0: usize, delta: usize, cfg: &SharedCfg) -> Option<Refined> {
    let n = seq.len();
    if n == 0 || cfg.refine.arm == 0 {
        return None;
    }

    let lo = s0.saturating_sub(cfg.refine.window);
    let hi = (s0 + cfg.refine.window).min(n); // note: s itself can be up to n-1; bounds checked inside

    let mut best_s: usize = 0;
    let mut best_score: i64 = i64::MIN;
    let mut any = false;

    for s in lo..=hi {
        // Need:
        //   a = seq[s .. s+arm]
        //   b = seq[s-delta .. s-delta+arm]
        // so require s>=delta and both intervals within bounds.
        if s < delta {
            continue;
        }
        let a0 = s;
        let b0 = s - delta;

        if a0 + cfg.refine.arm > n {
            continue;
        }
        if b0 + cfg.refine.arm > n {
            continue;
        }

        let a = &seq[a0..a0 + cfg.refine.arm];
        let b = &seq[b0..b0 + cfg.refine.arm];

        let mut score: i64 = 0;
        for i in 0..cfg.refine.arm {
            if a[i] == b[i] {
                score += 1;
            } else {
                score -= 1;
            }
        }

        if score > best_score {
            best_score = score;
            best_s = s;
        }
        any = true;
    }

    if !any {
        return None;
    }

    let matches = ((best_score + cfg.refine.arm as i64) / 2).max(0) as usize;
    let identity = matches as f32 / cfg.refine.arm as f32;

    Some(Refined {
        split_pos: best_s,
        score: best_score,
        identity_est: identity,
    })
}

pub fn detect_concatemer(
    seq: &[u8],
    shared: &SharedCfg,
    con: &ConcatOnlyCfg,
) -> Option<ConcatBreakpoint> {
    let k = shared.minimizer.k;
    let w = shared.minimizer.w;

    if seq.len() < k + w + 10 {
        return None;
    }

    let (pos_f, val_f) = sampled_minimizers_into(seq, &shared.minimizer, &mut vec![], &mut vec![]);

    if pos_f.len() < shared.min_matches {
        return None;
    }

    // index by minimizer value -> positions (sorted)
    let mut idx: HashMap<u64, Vec<i32>> = HashMap::new();
    for (&p, &v) in pos_f.iter().zip(val_f.iter()) {
        idx.entry(v).or_default().push(p as i32);
    }
    for v in idx.values_mut() {
        v.sort_unstable();
    }

    // delta bins: delta = p2 - p1
    let mut bins: HashMap<i32, Vec<(i32, i32)>> = HashMap::new();

    for (&p2_u, &v2) in pos_f.iter().zip(val_f.iter()) {
        let p2 = p2_u as i32;

        if let Some(p1s) = idx.get(&v2) {
            for &p1 in p1s {
                if p1 >= p2 {
                    break;
                }
                let delta = p2 - p1;
                if delta < con.min_delta as i32 {
                    continue;
                }
                let bin = div_floor(delta, shared.concat_diag_tol.max(1));
                bins.entry(bin).or_default().push((p1, p2));
            }
        }
    }

    best_bin_to_concat_breakpoint(&bins, shared, con, seq.len())
}

fn best_bin_to_concat_breakpoint(
    bins: &HashMap<i32, Vec<(i32, i32)>>,
    shared: &SharedCfg,
    con: &ConcatOnlyCfg,
    len: usize,
) -> Option<ConcatBreakpoint> {
    let mut best: Option<ConcatBreakpoint> = None;

    for (_bin, pts) in bins {
        if pts.len() < shared.min_matches {
            continue;
        }

        let mut v = pts.clone();
        v.sort_by_key(|(p1, _)| *p1);

        // monotone chain: p2 increasing
        let mut chain = Vec::new();
        let mut last_p2 = i32::MIN;
        for (p1, p2) in v {
            if p2 > last_p2 {
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
        if span < con.min_span {
            continue;
        }

        // delta median
        let mut ds: Vec<i32> = chain.iter().map(|(p1, p2)| p2 - p1).collect();
        ds.sort_unstable();
        let delta_med = ds[ds.len() / 2].max(0) as usize;

        let split_u = delta_med;
        if split_u < shared.end_guard || split_u + shared.end_guard > len {
            continue;
        }

        // crossing + p2 span checks
        let mut cross = 0usize;
        let mut p2_min = i32::MAX;
        let mut p2_max = i32::MIN;

        for &(p1, p2) in &chain {
            if (p1 as usize) < split_u && (p2 as usize) >= split_u {
                cross += 1;
            }
            p2_min = p2_min.min(p2);
            p2_max = p2_max.max(p2);
        }

        let cross_frac = cross as f32 / chain.len().max(1) as f32;
        if cross_frac < con.cross_frac {
            continue;
        }

        let p2_span = (p2_max as i64 - p2_min as i64).max(0) as usize;
        if p2_span < con.min_span {
            continue;
        }

        let score = span as i64 + (chain.len() as i64 * 10);
        let bp = ConcatBreakpoint {
            split_pos: split_u,
            delta: delta_med,
            score,
            matches: chain.len(),
            span,
            cross_frac,
            p2_span,
        };

        if best.as_ref().map(|b| bp.score > b.score).unwrap_or(true) {
            best = Some(bp);
        }
    }

    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cfg::{ConcatOnlyCfg, MinimizerCfg, RefineCfg, SharedCfg};

    fn test_concat_cfg(
        k: usize,
        w: usize,
        diag_tol: i32,
        min_span: usize,
        min_matches: usize,
        min_delta: usize,
        refine_window: usize,
        refine_arm: usize,
        refine_mode: RefineMode,
        max_ed_rate: f32,
    ) -> (SharedCfg, ConcatOnlyCfg) {
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
                concat_diag_tol: diag_tol,
            },
            ConcatOnlyCfg {
                min_span,
                min_delta,
                cross_frac: 0.8,
            },
        )
    }

    fn pseudo_dna(len: usize, seed: u64) -> Vec<u8> {
        // Simple LCG -> 2-bit base. Deterministic but not periodic like mod 4.
        let mut x = seed;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let b = match (x >> 62) & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            out.push(b);
        }
        out
    }

    fn repeat_bytes(pat: &[u8], times: usize) -> Vec<u8> {
        let mut v = Vec::with_capacity(pat.len() * times);
        for _ in 0..times {
            v.extend_from_slice(pat);
        }
        v
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

    /// Construct a concatemer: left + join + left(with optional noise)
    /// True split is at left.len() + join.len()
    fn make_concatemer(left: &[u8], join: &[u8], sub_every: usize) -> (Vec<u8>, usize) {
        let mut l2 = left.to_vec();
        if sub_every > 0 {
            l2 = mutate_substitutions(&l2, sub_every);
        }
        let mut s = Vec::with_capacity(left.len() * 2 + join.len());
        s.extend_from_slice(left);
        s.extend_from_slice(join);
        let bp = left.len() + join.len();
        s.extend_from_slice(&l2);
        (s, bp)
    }

    #[test]
    fn detects_simple_concatemer() {
        let left = pseudo_dna(6400, 1);
        let (s, true_bp) = make_concatemer(&left, b"", 0);

        let (shared, con) =
            test_concat_cfg(17, 21, 200, 2000, 20, 2000, 0, 0, RefineMode::HiFi, 0.25);

        let bp = detect_concatemer(&s, &shared, &con).unwrap();
        assert!((bp.split_pos as i32 - true_bp as i32).abs() <= 200);
    }

    #[test]
    fn respects_min_delta() {
        // local repeat only: delta ~ 256 < min_delta
        let left = repeat_bytes(b"ACGTTGCA", 100); // 800bp
        let (s, _true_bp) = make_concatemer(&left, b"", 0);

        let (shared, con) =
            test_concat_cfg(17, 21, 200, 2000, 20, 2000, 0, 0, RefineMode::HiFi, 0.25);

        let bp = detect_concatemer(&s, &shared, &con);

        assert!(bp.is_none(), "should ignore small deltas");
    }

    #[test]
    fn end_guard_blocks_near_ends() {
        let left = repeat_bytes(b"ACGTTGCA", 300); // 2400
        let (s, true_bp) = make_concatemer(&left, b"", 0);

        let (mut shared, con) =
            test_concat_cfg(17, 21, 200, 2000, 20, 2000, 0, 0, RefineMode::HiFi, 0.25);

        shared.end_guard = true_bp + 1;

        assert!(detect_concatemer(&s, &shared, &con).is_none());
    }

    #[test]
    fn refine_hifi_is_high_identity_when_seeded_near_truth() {
        let left = pseudo_dna(6400, 2);
        let (s, true_bp) = make_concatemer(&left, b"", 19);

        let (shared, con) = test_concat_cfg(
            15,
            19,
            200,
            2000,
            15,
            2000,
            200,
            300,
            RefineMode::HiFi,
            0.25,
        );

        let bp = detect_concatemer(&s, &shared, &con).unwrap();

        let refined = refine_concatemer_breakpoint(&s, true_bp, bp.delta, &shared).unwrap();
        assert!((refined.split_pos as i32 - true_bp as i32).abs() <= shared.refine.window as i32);
        assert!(refined.identity_est > 0.85);
    }

    #[test]
    fn refine_hifi_improves_over_coarse_seed() {
        let left = pseudo_dna(6400, 2);
        let (s, true_bp) = make_concatemer(&left, b"", 19);

        let (shared, con) = test_concat_cfg(
            15,
            19,
            200,
            2000,
            15,
            2000,
            200,
            300,
            RefineMode::HiFi,
            0.25,
        );

        let coarse = detect_concatemer(&s, &shared, &con).unwrap();
        let refined =
            refine_concatemer_breakpoint(&s, coarse.split_pos, coarse.delta, &shared).unwrap();

        assert!((refined.split_pos as i32 - true_bp as i32).abs() <= shared.refine.window as i32);
        assert!(refined.identity_est > 0.5);
    }

    #[test]
    fn refine_ont_tolerates_indels() {
        let left = repeat_bytes(b"ACGTTGCA", 800); // 6400
        let (mut s, true_bp) = make_concatemer(&left, b"", 0);

        // make the second copy indel-y
        let second = &s[true_bp..].to_vec();
        let mut noisy = inject_insertions(second, 37, b'A');
        noisy = inject_deletions(&noisy, 53);

        s.truncate(true_bp);
        s.extend_from_slice(&noisy);

        let (shared, con) =
            test_concat_cfg(13, 17, 250, 2000, 12, 2000, 200, 200, RefineMode::ONT, 0.35);

        let coarse = detect_concatemer(&s, &shared, &con).unwrap();

        let refined = refine_concatemer_breakpoint(&s, coarse.split_pos, coarse.delta, &shared)
            .expect("refine");
        assert!(
            (refined.split_pos as i32 - coarse.split_pos as i32).abs()
                <= shared.refine.window as i32
        );
        assert!(refined.identity_est > 0.4);
    }

    #[test]
    fn random_sequence_is_usually_none() {
        let s = pseudo_dna(20000, 999);

        let (shared, con) =
            test_concat_cfg(17, 21, 200, 2000, 20, 2000, 0, 0, RefineMode::HiFi, 0.25);

        let bp = detect_concatemer(&s, &shared, &con);

        assert!(bp.is_none());
    }
}
