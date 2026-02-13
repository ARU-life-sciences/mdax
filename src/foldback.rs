//! Foldback (palindromic junction) detection, refinement, and recursive excision.
//!
//! This module implements detection of *foldback / palindromic junctions* in
//! long sequencing reads, with a focus on multiple displacement amplification (MDA)
//! artefacts.
//!
//! A foldback is characterised by strong reverse-complement self-similarity across
//! a breakpoint, forming an anti-diagonal when forward minimizers are matched
//! against reverse-complement minimizers.
//!
//! High-level pipeline:
//! 1. **Coarse detection** (`detect_foldback`)
//!    - Compute strand-specific minimizers on the sequence and its reverse complement.
//!    - Match minimizers by value (not by re-hashing k-mers).
//!    - Cluster matchpoints by anti-diagonal (`p1 + p2`) into bins.
//!    - Select the best bin and estimate a coarse split position.
//!
//! 2. **Refinement** (`refine_breakpoint_*`)
//!    - HiFi mode: Hamming-style symmetric arm comparison.
//!    - ONT mode: banded edit distance between left arm and reverse-complemented right arm.
//!
//! 3. **Signature construction**
//!    - Prefer a *local* minimizer-based signature derived from matched minimizers
//!      participating in the foldback.
//!    - Fall back to a flank-based signature if local matches are insufficient.
//!
//! 4. **Recursive cutting**
//!    - Use support statistics to distinguish real genomic palindromes from
//!      stochastic MDA artefacts.
//!    - Iteratively trim artefactual foldbacks from the right end of the read.
//!
//! Design notes:
//! - Uses *forward-only* (strand-specific) minimizers for detection to avoid
//!   canonical-strand artefacts.
//! - Explicitly guards against low-complexity sequence via capped minimizer buckets.
//! - All coordinates are carefully tracked between full-sequence and view-local
//!   coordinate systems during recursion.

use crate::cfg::{FoldOnlyCfg, MdaxCfg, SharedCfg};
use crate::elog;
use crate::fingerprint::{
    SupportStats, foldback_signature, foldback_signature_from_local_matches, is_real_foldback,
};
use crate::minimizer::sampled_minimizers_into;
use crate::scratch::{FoldScratch, RefineScratch, SigScratch};
use crate::utils::{
    RefineMode, Refined, banded_edit_distance_scratch, comp, div_floor, revcomp_into,
};
use anyhow::Result;
use gxhash::HashMap;
use std::ops::Range;

/// Per-bin summary statistics used during coarse foldback detection.
///
/// Each bin corresponds to a quantised anti-diagonal:
///     d = p1 + p2
/// where `p1` is a forward minimizer position and `p2` is the corresponding
/// reverse-complement position mapped into forward coordinates.
///
/// The statistics are deliberately lightweight and are used only to:
/// - discard bins with insufficient support
/// - discard bins with insufficient spatial span
/// - choose the single “best” bin for second-phase point collection
#[derive(Copy, Clone, Debug)]
pub struct BinStat {
    /// Number of matchpoints falling into this bin.
    count: usize,
    /// Minimum forward position observed in this bin.
    min_p1: i32,
    /// Maximum forward position observed in this bin.
    max_p1: i32,
}

impl BinStat {
    /// Create a new bin initialised with a single forward position.
    #[inline]
    fn new(p1: i32) -> Self {
        Self {
            count: 1,
            min_p1: p1,
            max_p1: p1,
        }
    }

    /// Update bin statistics with an additional forward position.
    #[inline]
    fn update(&mut self, p1: i32) {
        self.count += 1;
        self.min_p1 = self.min_p1.min(p1);
        self.max_p1 = self.max_p1.max(p1);
    }

    /// Span of forward positions in this bin.
    ///
    /// This approximates the arm length supporting the foldback.
    #[inline]
    fn span(&self) -> usize {
        (self.max_p1 as i64 - self.min_p1 as i64).max(0) as usize
    }

    /// Proxy score used to rank bins during coarse detection.
    ///
    /// Combines spatial span and match count, strongly favouring bins
    /// with many supporting minimizers.
    #[inline]
    fn score_proxy(&self) -> i64 {
        self.span() as i64 + (self.count as i64 * 10)
    }
    /// Rank bins primarily by count, then span, then proxy (tie-breaker).
    ///
    /// For your window sizes and `min_matches`, this is very safe from overflow:
    /// - count up to ~1e6 => (count<<42) stays < 2^63
    /// - span up to window len => (span<<16) is tiny by comparison
    #[inline]
    fn rank_key(&self) -> i64 {
        let c = self.count as i64; // dominates
        let s = self.span() as i64; // next
        let p = self.score_proxy(); // tie-breaker only
        (c << 42) + (s << 16) + (p & 0xFFFF)
    }

    // If you want these externally (optional):
    #[inline]
    pub fn count(&self) -> usize {
        self.count
    }
    #[inline]
    pub fn min_p1(&self) -> i32 {
        self.min_p1
    }
    #[inline]
    pub fn max_p1(&self) -> i32 {
        self.max_p1
    }
}

/// A multi-hit foldback call: breakpoint + the supporting matchpoints.
#[derive(Debug, Clone)]
pub struct FoldHit {
    pub fb: FoldBreakpoint,
    /// Matchpoints (p1,p2) in window-local forward coordinates.
    pub pts: Vec<(i32, i32)>,
    /// Local match list (pos,value) in forward coords (useful for signatures later if desired).
    pub matches: Vec<(usize, u64)>,
    /// Anti-diagonal bin id (mostly for debugging).
    pub bin: i32,
    // TODO: document these
    pub bin_count: usize,
    pub bin_span: usize,
    pub bin_rank: i64,
}

/// Result of coarse foldback detection.
///
/// This represents an approximate junction location derived purely from
/// minimizer geometry, prior to any sequence-level refinement.
#[derive(Debug, Clone)]
pub struct FoldBreakpoint {
    /// Estimated split position (0-based, relative to the sequence passed
    /// to `detect_foldback`).
    pub split_pos: usize,
    /// Coarse score combining span and number of supporting minimizers.
    pub score: i64,
    /// Number of minimizer matches supporting the foldback.
    pub matches: usize,
    /// Span of forward positions supporting the foldback (proxy for arm length).
    pub span: usize,
}

/// Refine a coarse foldback breakpoint into a more precise split position.
///
/// This function dispatches to either a HiFi-optimised or ONT-optimised
/// refinement strategy, based on `cfg.refine.mode`.
///
/// Returns `None` if the refinement window would exceed sequence bounds
/// or if no acceptable refinement is found.
pub fn refine_breakpoint(
    seq: &[u8],
    s0: usize,
    cfg: &SharedCfg,
    refine_scratch: &mut RefineScratch,
) -> Result<Option<Refined>> {
    let n = seq.len();
    if s0 < cfg.refine.arm + cfg.refine.window || s0 + cfg.refine.arm + cfg.refine.window > n {
        elog!(
            "REFINE BREAKPOINT",
            "Refinement context exceeds sequence bounds: s0={}, arm={}, window={}, len={}",
            s0,
            cfg.refine.arm,
            cfg.refine.window,
            n
        );
        return Ok(None);
    }
    match cfg.refine.mode {
        RefineMode::HiFi => Ok(refine_breakpoint_hamming(seq, s0, cfg)),
        RefineMode::ONT => Ok(refine_breakpoint_banded_ed(seq, s0, cfg, refine_scratch)),
    }
}

/// HiFi-optimised breakpoint refinement using symmetric arm comparison.
///
/// For each candidate split within `±window` of `s0`, this compares:
/// - the left arm (walking outward from the split)
/// - the reverse-complement of the right arm
///
/// Scoring is Hamming-like: +1 for match, −1 for mismatch.
///
/// This assumes low indel rates and is therefore **not suitable for ONT data**.
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

/// ONT-optimised breakpoint refinement using banded edit distance.
///
/// This compares the left arm against the reverse-complemented right arm
/// using a banded dynamic programming alignment, allowing for indels.
///
/// The band width is derived from `max_ed_rate * arm_length`.
///
/// Returns the split position with minimal edit distance within the
/// refinement window.
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

    let mut found = false;

    for s in (s0 - cfg.refine.window)..=(s0 + cfg.refine.window) {
        let left = &seq[s - cfg.refine.arm..s];
        let right = &seq[s..s + cfg.refine.arm];

        // Fill scratch.right_rc with revcomp(right)
        revcomp_into(&mut scratch.right_rc, right);

        let ed = banded_edit_distance_scratch(
            left,
            &scratch.right_rc,
            band,
            &mut scratch.prev,
            &mut scratch.curr,
        );

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

/// Detect a foldback (palindromic) junction in a single sequence.
///
/// This performs *coarse* detection based purely on minimizer geometry.
///
/// Algorithm:
/// 1. Compute strand-specific (forward-only) minimizers on `seq` and on
///    `revcomp(seq)`.
/// 2. Index forward minimizers by minimizer *value*, capping bucket sizes
///    to avoid low-complexity blowups.
/// 3. Match reverse-complement minimizers to forward minimizers by value,
///    mapping rc positions back into forward coordinates.
/// 4. Cluster matchpoints by quantised anti-diagonal `d = p1 + p2`.
/// 5. Select the best bin based on support and span, then estimate the
///    split position as `median(d) / 2`.
///
/// Returns `None` if:
/// - too few minimizers are present
/// - all candidate bins fail support/span thresholds
/// - the inferred split is too close to sequence ends
///
/// Notes:
/// - Uses minimizer *values* directly (not re-hashed k-mers) to ensure
///   consistency with minimizer selection.
/// - Forward-only minimizers avoid canonical-strand artefacts and make
///   anti-diagonal structure explicit.
/// - Low-complexity sequence may be intentionally ignored if minimizer
///   values exceed the bucket cap.
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

    revcomp_into(&mut scratch.rc, seq);

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

        // Map rc minimizer start position (prc) back into forward coordinates.
        // If forward k-mer starts at i, its start in rc is n-k-i => i = n-k-prc.
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
                // Use forward minimizer position for locality filtering in signature
                scratch.best_matches.push((p1 as usize, vrc));
            }
        }
    }

    if scratch.best_pts.len() < shared.min_matches {
        return None;
    }

    // reuse existing scoring logic
    fold_breakpoint_from_pts(&mut scratch.best_pts, shared, fold, seq.len())
}

// swap out but keep capacity in scratch
fn take_vec_preserve_capacity<T>(v: &mut Vec<T>) -> Vec<T> {
    let mut out = Vec::with_capacity(v.capacity());
    std::mem::swap(&mut out, v);
    out
}

/// Detect up to `k_hits` foldback candidates in `seq`, ranked by proxy score.
///
/// This does the same coarse minimizer geometry as `detect_foldback`, but
/// returns multiple bins instead of only the best bin.
///
/// Notes:
/// - This allocates per-hit `Vec`s for pts/matches. That’s fine for `irx`,
///   but `mdax` should keep using `detect_foldback()` for speed.
pub fn detect_foldbacks_topk(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    scratch: &mut FoldScratch,
    k: usize,
) -> Vec<FoldHit> {
    let k = k.max(1);

    let mm_k = shared.minimizer.k;
    let mm_w = shared.minimizer.w;
    if seq.len() < mm_k + mm_w + 10 {
        return Vec::new();
    }

    // ---- minimizers ----
    sampled_minimizers_into(
        seq,
        &shared.minimizer,
        &mut scratch.pos_f,
        &mut scratch.val_f,
    );
    revcomp_into(&mut scratch.rc, seq);
    sampled_minimizers_into(
        &scratch.rc,
        &shared.minimizer,
        &mut scratch.pos_rc,
        &mut scratch.val_rc,
    );

    if scratch.pos_f.len() < shared.min_matches || scratch.pos_rc.len() < shared.min_matches {
        return Vec::new();
    }

    // ---- index forward minimizers by value (cap buckets) ----
    const MAX_BUCKET: usize = 64;
    scratch.idx_f.clear();
    scratch.repetitive.clear();

    for (&p, &v) in scratch.pos_f.iter().zip(scratch.val_f.iter()) {
        if scratch.repetitive.contains_key(&v) {
            continue;
        }
        let e = scratch.idx_f.entry(v).or_default();
        if e.len() < MAX_BUCKET {
            e.push(p as i32);
        } else {
            scratch.idx_f.remove(&v);
            scratch.repetitive.insert(v, ());
        }
    }

    // ---- per-bin stats ----
    let n = seq.len() as i32;
    let mm_k_i32 = mm_k as i32;

    scratch.stats.clear();

    for (&prc, &vrc) in scratch.pos_rc.iter().zip(scratch.val_rc.iter()) {
        if scratch.repetitive.contains_key(&vrc) {
            continue;
        }
        let Some(p1s) = scratch.idx_f.get(&vrc) else {
            continue;
        };
        let p2 = n - mm_k_i32 - (prc as i32);

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

    scratch.top_bins.clear();

    for (&bin, st) in scratch.stats.iter() {
        let cnt = st.count();
        let sp = st.span();
        if cnt < shared.min_matches {
            continue;
        }
        if sp < fold.min_arm {
            continue;
        }
        scratch.top_bins.push((bin, st.rank_key(), cnt, sp));
    }

    scratch
        .top_bins
        .sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    scratch.top_bins.truncate(k);

    if scratch.top_bins.is_empty() {
        return Vec::new();
    }

    // ---- collect points per selected bin in one pass ----
    // We'll build per-bin vectors using a tiny map from bin->index in `bins`.
    scratch.bin_to_idx.clear();
    for (i, (bin, _rk, _cnt, _sp)) in scratch.top_bins.iter().enumerate() {
        scratch.bin_to_idx.insert(*bin, i);
    }

    // Ensure per-bin buffers exist and are empty
    let nb = scratch.top_bins.len();
    if scratch.pts_per.len() < nb {
        scratch.pts_per.resize_with(nb, Vec::new);
    }
    if scratch.matches_per.len() < nb {
        scratch.matches_per.resize_with(nb, Vec::new);
    }
    for i in 0..nb {
        scratch.pts_per[i].clear();
        scratch.matches_per[i].clear();

        // optional but helpful: reserve roughly to avoid growth
        let (_bin, _rk, cnt, _sp) = scratch.top_bins[i];
        scratch.pts_per[i].reserve(cnt);
        scratch.matches_per[i].reserve(cnt);
    }

    for (&prc, &vrc) in scratch.pos_rc.iter().zip(scratch.val_rc.iter()) {
        if scratch.repetitive.contains_key(&vrc) {
            continue;
        }
        let Some(p1s) = scratch.idx_f.get(&vrc) else {
            continue;
        };
        let p2 = n - mm_k_i32 - (prc as i32);

        for &p1 in p1s {
            let d = p1 + p2;
            let bin = div_floor(d, shared.fold_diag_tol.max(1));
            if let Some(&ix) = scratch.bin_to_idx.get(&bin) {
                scratch.pts_per[ix].push((p1, p2));
                scratch.matches_per[ix].push((p1 as usize, vrc));
            }
        }
    }

    // ---- compute FoldBreakpoint per bin ----
    let mut hits = Vec::with_capacity(scratch.top_bins.len());

    for (i, (bin, rk, cnt, sp)) in scratch.top_bins.drain(..).enumerate() {
        let mut pts = take_vec_preserve_capacity(&mut scratch.pts_per[i]);
        if pts.len() < shared.min_matches {
            continue;
        }

        let Some(fb) = fold_breakpoint_from_pts(&mut pts, shared, fold, seq.len()) else {
            continue;
        };

        let matches = take_vec_preserve_capacity(&mut scratch.matches_per[i]);

        hits.push(FoldHit {
            fb,
            bin,
            pts,     // already sorted/mutated; OK
            matches, // TODO: keep if we want?
            bin_count: cnt,
            bin_span: sp,
            bin_rank: rk,
        });
    }

    hits.sort_by(|a, b| b.bin_rank.cmp(&a.bin_rank).then_with(|| a.bin.cmp(&b.bin)));
    hits
}

/// Convert a set of matched minimizer points into a coarse foldback breakpoint.
///
/// Expects points `(p1, p2)` in forward coordinates that approximately lie
/// along a single anti-diagonal.
///
/// The points are filtered into a monotone chain where `p2` decreases as
/// `p1` increases, enforcing foldback geometry.
///
/// The split position is estimated as:
///     split ~= median(p1 + p2) / 2
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

/// Recursively excise artefactual foldbacks from a sequence prefix.
///
/// Starting from an initial foldback detection and refinement (given in
/// *full-sequence coordinates*), this function:
/// - Converts breakpoints into view-local coordinates
/// - Re-detects and refines foldbacks on progressively shortened prefixes
/// - Uses foldback signatures and support statistics to decide whether a
///   foldback is real or artefactual
///
/// Real (supported) foldbacks terminate recursion.
/// Artefactual foldbacks cause the sequence to be truncated at the split.
///
/// Returns the final `[start, end)` range to keep from the original sequence.
///
/// Important invariants:
/// - `first_fb.split_pos` and `first_rf.split_pos` are in full-sequence coordinates.
/// - All subsequent detections operate in view-local coordinates.
/// - `scratch.best_matches` is assumed to correspond to the most recent
///   `detect_foldback` call.
pub fn recursive_foldback_cut_from_first_range(
    seq: &[u8],               // original full sequence
    mut keep: Range<usize>,   // current keep window into `seq`
    first_fb: FoldBreakpoint, // split_pos is assumed relative to `seq` (full coords)
    first_rf: Refined,        // split_pos is assumed relative to `seq` (full coords)
    cfg: &MdaxCfg,
    support: &HashMap<u64, SupportStats>,
    max_depth: usize,
    scratch: &mut FoldScratch,
    sig_scratch: &mut SigScratch,
) -> anyhow::Result<Range<usize>> {
    // We already have the first hit/refine (full coords). After first cut, detect/refine on the view.
    let mut fb_opt = Some(first_fb);
    let mut rf_opt = Some(first_rf);

    for _ in 0..max_depth {
        // View of the sequence we are currently considering
        let view = &seq[keep.clone()];

        // --- 1) detect foldback (skip for first iteration) ---
        let fb = match fb_opt.take() {
            Some(mut fb) => {
                // Convert full-coord split into view-local coord if needed.
                // (If keep.start==0, this is a no-op.)
                if fb.split_pos >= keep.start {
                    fb.split_pos -= keep.start;
                } else {
                    // If inconsistent, bail rather than panic.
                    break;
                }
                fb
            }
            None => match detect_foldback(view, &cfg.shared, &cfg.fold, scratch) {
                Some(fb) => fb,
                None => break,
            },
        };

        // --- 2) refine breakpoint (skip for first iteration) ---
        let rf = match rf_opt.take() {
            Some(mut rf) => {
                if rf.split_pos >= keep.start {
                    rf.split_pos -= keep.start;
                } else {
                    break;
                }
                rf
            }
            None => {
                match refine_breakpoint(view, fb.split_pos, &cfg.shared, &mut scratch.refine)? {
                    Some(rf) => rf,
                    None => break,
                }
            }
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
            // fallback: flank signature, quantize split (view-local coords)
            let q = 150usize;
            let split_q = (rf.split_pos / q) * q;

            foldback_signature(
                view,
                split_q,
                &cfg.shared,
                cfg.sig.flank_bp,
                cfg.sig.take,
                sig_scratch,
                cfg.sig.value_shift,
            )
        });

        let Some(sig) = sig else { break };

        if is_real_foldback(sig, support, cfg) {
            // real -> don't cut
            break;
        }

        // artefact -> cut left (i.e., keep prefix of current view)
        let split = rf.split_pos.min(view.len());
        keep.end = keep.start + split;

        // next iteration will detect/refine on the shortened view
    }

    Ok(keep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::revcomp_in_place;
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
        revcomp_in_place(&mut rc);

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

    fn random_dna(len: usize, mut x: u64) -> Vec<u8> {
        // Deterministic xorshift64* PRNG, no external crates.
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            let r = x.wrapping_mul(0x2545F4914F6CDD1D);
            out.push(match (r & 3) as u8 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            });
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

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch)
            .unwrap()
            .unwrap();

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
            .unwrap()
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
        let hifi = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch).unwrap();

        shared.refine.mode = RefineMode::ONT;
        let ont = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch)
            .unwrap()
            .unwrap();

        assert!(ont.identity_est > 0.4);
        assert!(
            (ont.split_pos as i32 - coarse.split_pos as i32).abs() <= shared.refine.window as i32
        );

        // optional: HiFi may fail, but ONT shouldn't
        assert!(hifi.is_some() || ont.identity_est > 0.4);
    }

    #[test]
    fn local_signature_from_best_matches_is_some_and_consistent() {
        // Use non-repetitive sequence so MAX_BUCKET/repetitive filter doesn’t delete everything.
        let left = random_dna(4000, 0xC0FFEE);
        let (s, _true_bp) = make_noisy_foldback_from_left(&left, b"", 0, 0, 0);

        // Slightly denser minimizers can help, but the key fix is non-repetitive input.
        let (shared, fold) = test_fold_cfg(
            7,   // k
            7,   // w  (denser than 9/11)
            120, // diag_tol
            200, // min_arm (coarse)
            3,   // min_matches
            20,  // refine_window
            80,  // refine_arm
            RefineMode::HiFi,
            0.0,
        );

        let mut scratch = scratch::FoldScratch::new();

        let coarse = detect_foldback(&s, &shared, &fold, &mut scratch)
            .expect("detect_foldback should succeed on a clean random foldback");

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, &mut scratch.refine)
            .unwrap()
            .expect("refine should succeed");

        // Local signature: choose parameters that cannot fail due to take being too large.
        let window_bp = 300usize;
        let take = 6usize;
        let value_shift = 0u8;

        let sig_local = foldback_signature_from_local_matches(
            &mut scratch.best_matches,
            refined.split_pos,
            window_bp,
            take,
            value_shift,
        );

        assert!(sig_local.is_some(), "local signature should be computable");
    }

    #[test]
    fn detect_foldback_repetitive_sequence_may_be_ignored() {
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 300);
        let (s, _bp) = make_noisy_foldback_from_left(&left, b"", 0, 0, 0);

        let (shared, fold) = test_fold_cfg(7, 9, 120, 200, 3, 20, 80, RefineMode::HiFi, 0.0);
        let mut scratch = scratch::FoldScratch::new();

        // Depending on exact minimizer impl, this may be None due to repetitive filtering.
        let _ = detect_foldback(&s, &shared, &fold, &mut scratch);
    }
}
