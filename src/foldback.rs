use crate::cfg::{FoldOnlyCfg, MdaxCfg, SharedCfg};
use crate::fingerprint::{
    SupportStats, foldback_signature, foldback_signature_from_local_matches,
    foldback_signature_from_matches, is_real_foldback,
};
use crate::minimizer::sampled_minimizers_into;
use crate::scratch::{FoldScratch, RefineScratch, SigScratch};
use crate::utils::{
    RefineMode, Refined, banded_edit_distance_scratch, comp, div_floor, revcomp_in_place,
    revcomp_into,
};
use gxhash::HashMap;

#[derive(Copy, Clone, Debug)]
pub struct BinStat {
    count: usize,
    min_p1: i32,
    max_p1: i32,
}
impl BinStat {
    fn new(p1: i32) -> Self {
        Self {
            count: 1,
            min_p1: p1,
            max_p1: p1,
        }
    }
    fn update(&mut self, p1: i32) {
        self.count += 1;
        self.min_p1 = self.min_p1.min(p1);
        self.max_p1 = self.max_p1.max(p1);
    }
    fn span(&self) -> usize {
        (self.max_p1 as i64 - self.min_p1 as i64).max(0) as usize
    }
    fn score_proxy(&self) -> i64 {
        self.span() as i64 + (self.count as i64 * 10)
    }
}

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

/// Result of foldback detection (coarse).
#[derive(Debug, Clone)]
pub struct FoldBreakpoint {
    pub split_pos: usize,
    pub score: i64,
    pub matches: usize,
    pub span: usize,
}

/// Refine an inexact breakpoint estimate `s0` into a more exact split position.
pub fn refine_breakpoint(
    seq: &[u8],
    s0: usize,
    cfg: &SharedCfg,
    refine_scratch: &mut RefineScratch,
) -> Option<Refined> {
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
        RefineMode::ONT => refine_breakpoint_banded_ed(seq, s0, cfg, refine_scratch),
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

// for ONT data
pub fn refine_breakpoint_banded_ed(
    seq: &[u8],
    s0: usize,
    cfg: &SharedCfg,
    scratch: &mut crate::scratch::RefineScratch,
) -> Option<Refined> {
    // ... bounds checks ...

    let mut band = (cfg.refine.max_ed_rate.max(0.0) * cfg.refine.arm as f32).ceil() as usize;
    band = band.max(1).min(cfg.refine.arm);

    let mut best_s = s0;
    let mut best_ed: usize = usize::MAX;

    scratch.right_rc.resize(cfg.refine.arm, b'N');
    let mut found = false;

    for s in (s0 - cfg.refine.window)..=(s0 + cfg.refine.window) {
        let left = &seq[s - cfg.refine.arm..s];
        let right = &seq[s..s + cfg.refine.arm];

        revcomp_into(right, &mut scratch.right_rc);

        let ed = banded_edit_distance_scratch(
            left,
            &scratch.right_rc,
            band,
            &mut scratch.prev,
            &mut scratch.curr,
        );

        // If ed == band+1, true distance is > band => treat as failure
        if ed > band {
            continue;
        }
        found = true;

        if ed < best_ed {
            best_ed = ed;
            best_s = s;
        }
    }

    if !found {
        return None; // or identity_est=0.0
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
/// 2) Index forward minimizers by their minimizer *value*
/// 3) For each minimizer on rc, match by value and convert rc pos -> forward coordinate
/// 4) Cluster matchpoints by anti-diagonal d = p1 + p2 (within diag_tol)
/// 5) Choose best cluster; return split ≈ median(d/2)
pub fn detect_foldback(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    scratch: &mut FoldScratch,
) -> Option<FoldBreakpoint> {
    let k = shared.minimizer.k;
    let w = shared.minimizer.w;

    if seq.len() < k + w + 10 {
        return None;
    }

    // forward-only minimizers (or canonical if configured), returning (pos, value)
    sampled_minimizers_into(
        seq,
        &shared.minimizer,
        &mut scratch.pos_f,
        &mut scratch.val_f,
    );

    scratch.rc.clear();
    scratch.rc.extend_from_slice(seq);
    revcomp_in_place(&mut scratch.rc);

    sampled_minimizers_into(
        &scratch.rc,
        &shared.minimizer,
        &mut scratch.pos_rc,
        &mut scratch.val_rc,
    );

    let pos_f = &scratch.pos_f;
    let val_f = &scratch.val_f;
    let pos_rc = &scratch.pos_rc;
    let val_rc = &scratch.val_rc;

    if pos_f.len() < shared.min_matches || pos_rc.len() < shared.min_matches {
        return None;
    }

    // ---- 1) index forward minimizers by value -> positions (with cap) ----
    // Avoid pathological buckets from low-complexity sequence.
    const MAX_BUCKET: usize = 64;

    scratch.idx_f.clear();
    scratch.repetitive.clear();

    for (&p, &v) in pos_f.iter().zip(val_f.iter()) {
        if scratch.repetitive.contains_key(&v) {
            continue;
        }
        let e = scratch.idx_f.entry(v).or_default();
        if e.len() < MAX_BUCKET {
            e.push(p as i32);
        } else {
            // mark as repetitive; drop positions to stop quadratic blowups
            scratch.idx_f.remove(&v);
            scratch.repetitive.insert(v, ());
        }
    }

    // ---- 2) first phase: lightweight per-bin stats, choose best bin ----
    let n = seq.len() as i32;
    let k_i32 = k as i32;

    scratch.stats.clear();

    for (&prc, &vrc) in pos_rc.iter().zip(val_rc.iter()) {
        if scratch.repetitive.contains_key(&vrc) {
            continue;
        }
        let Some(p1s) = scratch.idx_f.get(&vrc) else {
            continue;
        };

        let p2 = n - k_i32 - (prc as i32);

        for &p1 in p1s {
            let d = p1 + p2;
            let bin = div_floor(d, shared.fold_diag_tol.max(1));
            scratch
                .stats
                .entry(bin)
                .and_modify(|st| st.update(p1))
                .or_insert_with(|| BinStat::new(p1));
        }
    }

    // Choose best bin using gates + proxy score
    let mut best_bin: Option<i32> = None;
    let mut best_score: i64 = i64::MIN;
    let mut best_count: usize = 0;

    for (&bin, st) in scratch.stats.iter() {
        if st.count < shared.min_matches {
            continue;
        }
        if st.span() < fold.min_arm {
            continue;
        }
        let sc = st.score_proxy();
        if sc > best_score {
            best_score = sc;
            best_bin = Some(bin);
            best_count = st.count;
        }
    }

    let best_bin = best_bin?;

    if best_count < shared.min_matches {
        return None;
    }

    // ---- 3) second phase: collect points only for best bin ----
    let cap = scratch.stats.get(&best_bin).map(|s| s.count).unwrap_or(0);
    scratch.best_pts.clear();
    scratch.best_matches.clear();
    scratch.best_pts.reserve(cap);
    scratch.best_matches.reserve(cap);

    for (&prc, &vrc) in pos_rc.iter().zip(val_rc.iter()) {
        if scratch.repetitive.contains_key(&vrc) {
            continue;
        }
        let Some(p1s) = scratch.idx_f.get(&vrc) else {
            continue;
        };

        let p2 = n - k_i32 - (prc as i32);

        for &p1 in p1s {
            let d = p1 + p2;
            let bin = div_floor(d, shared.fold_diag_tol.max(1));
            if bin == best_bin {
                scratch.best_pts.push((p1, p2));
                // TODO: is this correct?
                scratch.best_matches.push((p2 as usize, vrc));
            }
        }
    }

    if scratch.best_pts.len() < shared.min_matches {
        return None;
    }

    // reuse existing scoring logic
    fold_breakpoint_from_pts(&mut scratch.best_pts, shared, fold, seq.len())
}

fn fold_breakpoint_from_pts(
    pts: &mut Vec<(i32, i32)>,
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    len: usize,
) -> Option<FoldBreakpoint> {
    if pts.len() < shared.min_matches {
        return None;
    }

    pts.sort_by_key(|(p1, _)| *p1);

    // monotone anti-diagonal chain: p2 decreasing as p1 increases
    let mut chain: Vec<(i32, i32)> = Vec::with_capacity(pts.len());
    let mut last_p2 = i32::MAX;
    for &(p1, p2) in pts.iter() {
        if p2 < last_p2 {
            chain.push((p1, p2));
            last_p2 = p2;
        }
    }
    if chain.len() < shared.min_matches {
        return None;
    }

    let p1_min = chain.first().unwrap().0 as i64;
    let p1_max = chain.last().unwrap().0 as i64;
    let span = (p1_max - p1_min).max(0) as usize;
    if span < fold.min_arm {
        return None;
    }

    // median d = p1 + p2 => split ≈ d/2
    let mut ds: Vec<i32> = chain.iter().map(|(p1, p2)| p1 + p2).collect();
    ds.sort_unstable();
    let d_med = ds[ds.len() / 2];
    let split_u = ((d_med as f64) / 2.0).round().max(0.0) as usize;

    if split_u < shared.end_guard || split_u + shared.end_guard > len {
        return None;
    }

    let score = span as i64 + (chain.len() as i64 * 10);
    Some(FoldBreakpoint {
        split_pos: split_u,
        score,
        matches: chain.len(),
        span,
    })
}

pub fn recursive_foldback_cut<'a>(
    mut seq: &'a [u8],
    cfg: &MdaxCfg,
    support: &HashMap<u64, SupportStats>,
    max_depth: usize,
    scratch: &mut FoldScratch,
    sig_scratch: &mut SigScratch,
) -> anyhow::Result<&'a [u8]> {
    for _ in 0..max_depth {
        let Some(fb) = detect_foldback(seq, &cfg.shared, &cfg.fold, scratch) else {
            break;
        };
        let Some(rf) = refine_breakpoint(seq, fb.split_pos, &cfg.shared, &mut scratch.refine)
        else {
            break;
        };
        if rf.identity_est < cfg.fold2.min_identity {
            break;
        }

        // Prefer evidence-based signature (more robust to refinement jitter).
        let sig = foldback_signature_from_local_matches(
            &mut scratch.best_matches,
            rf.split_pos,
            cfg.sig.flank_bp,
            cfg.sig.take,
            cfg.sig.value_shift,
        )
        .or_else(|| {
            // fallback: flank signature, but quantize split
            let q = 150usize;
            let split_q = (rf.split_pos / q) * q;

            foldback_signature(
                seq,
                split_q,
                &cfg.shared,
                cfg.sig.flank_bp,
                cfg.sig.take,
                sig_scratch,
                cfg.sig.value_shift,
            )
        });

        let Some(sig) = sig else {
            break;
        };

        if is_real_foldback(sig, support, cfg) {
            // don't cut
            break;
        }

        // artefact -> cut left

        let split = rf.split_pos.min(seq.len());
        seq = &seq[..split];
    }
    Ok(seq)
}

pub fn recursive_foldback_cut_from_first<'a>(
    mut seq: &'a [u8],
    first_fb: FoldBreakpoint,
    first_rf: Refined,
    cfg: &MdaxCfg,
    support: &HashMap<u64, SupportStats>,
    max_depth: usize,
    scratch: &mut FoldScratch,
    sig_scratch: &mut SigScratch,
) -> anyhow::Result<&'a [u8]> {
    // We already have the first hit/refine. After the first cut, detect/refine normally.
    let mut fb_opt = Some(first_fb);
    let mut rf_opt = Some(first_rf);

    for _ in 0..max_depth {
        // --- 1) detect foldback (skip for first iteration) ---
        let fb = match fb_opt.take() {
            Some(fb) => fb,
            None => match detect_foldback(seq, &cfg.shared, &cfg.fold, scratch) {
                Some(fb) => fb,
                None => break,
            },
        };

        // --- 2) refine breakpoint (skip for first iteration) ---
        let rf = match rf_opt.take() {
            Some(rf) => rf,
            None => match refine_breakpoint(seq, fb.split_pos, &cfg.shared, &mut scratch.refine) {
                Some(rf) => rf,
                None => break,
            },
        };

        if rf.identity_est < cfg.fold2.min_identity {
            break;
        }

        // --- 3) signature (prefer local matches computed by detect_foldback) ---
        // Important: on the first iteration, scratch.best_matches is already populated
        // because pass2 called detect_foldback before calling this function.
        let sig = foldback_signature_from_local_matches(
            &mut scratch.best_matches,
            rf.split_pos,
            cfg.sig.flank_bp,
            cfg.sig.take,
            cfg.sig.value_shift,
        )
        .or_else(|| {
            // fallback: flank signature, quantize split
            let q = 150usize;
            let split_q = (rf.split_pos / q) * q;

            foldback_signature(
                seq,
                split_q,
                &cfg.shared,
                cfg.sig.flank_bp,
                cfg.sig.take,
                sig_scratch,
                cfg.sig.value_shift,
            )
        });

        let Some(sig) = sig else {
            break;
        };

        if is_real_foldback(sig, support, cfg) {
            // real -> don't cut
            break;
        }

        // artefact -> cut left
        let split = rf.split_pos.min(seq.len());
        seq = &seq[..split];

        // next iteration will (re)detect/refine on the shortened seq
    }

    Ok(seq)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cfg::{FoldOnlyCfg, MinimizerCfg, RefineCfg, SharedCfg},
        scratch,
    };

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
        let mut scratch = scratch::FoldScratch::new();

        let bp = detect_foldback(&s, &shared, &fold, &mut scratch).unwrap();

        assert!((bp.split_pos as i32 - true_bp as i32).abs() <= 20);
    }

    #[test]
    fn refine_hifi_hits_exact_breakpoint_clean() {
        let left = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        let join = b"";
        let (s, true_bp) = make_noisy_foldback_from_left(left, join, 12, 0, 0);

        let (shared, fold) = test_fold_cfg(9, 11, 50, 20, 5, 10, 30, RefineMode::HiFi, 0.0);
        let mut scratch = scratch::FoldScratch::new();

        let coarse = detect_foldback(&s, &shared, &fold, &mut scratch).unwrap();
        let refine_scratch = &mut RefineScratch::default();

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch).unwrap();

        assert!((refined.split_pos as i32 - true_bp as i32).abs() <= 2);

        assert!(refined.identity_est > 0.9);
    }

    #[test]
    fn refine_ont_tolerates_indels() {
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 40);
        let (s, _bp) = make_noisy_foldback_from_left(&left, b"", 0, 23, 41);

        let (shared, fold) = test_fold_cfg(7, 9, 120, 30, 3, 40, 60, RefineMode::ONT, 0.35);

        let mut scratch = scratch::FoldScratch::new();
        let coarse = detect_foldback(&s, &shared, &fold, &mut scratch).unwrap();
        let refine_scratch = &mut RefineScratch::default();

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch)
            .expect("should refine");

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
        let mut scratch = scratch::FoldScratch::new();

        let coarse = detect_foldback(&s, &shared, &fold, &mut scratch).unwrap();
        let refine_scratch = &mut RefineScratch::default();

        shared.refine.mode = RefineMode::HiFi;
        let hifi = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch);

        shared.refine.mode = RefineMode::ONT;
        let ont = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch).unwrap();

        assert!(ont.identity_est > 0.4);
        assert!(
            (ont.split_pos as i32 - coarse.split_pos as i32).abs() <= shared.refine.window as i32
        );

        // optional: HiFi may fail, but ONT shouldn't
        assert!(hifi.is_some() || ont.identity_est > 0.4);
    }
}
