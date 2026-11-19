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
use crate::fingerprint::{JunctionSupport, axis_tol, junction_keys};
use crate::minimizer::sampled_minimizers_into;
use crate::scratch::{FoldScratch, RefineScratch, SigScratch};
use crate::utils::{
    RefineMode, Refined, comp, div_floor, edit_distance_myers,
    revcomp_into,
};
use anyhow::Result;
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
    /// Matchpoints come in mirrored pairs (one with `p1` in each arm), so this
    /// approximates the whole foldback footprint (arm + gap + arm), not a
    /// single arm. It is an upper bound on arm length.
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
    /// Estimated end of the left arm (exclusive) in read coordinates.
    /// Derived from the rightmost forward minimizer position + k.
    /// Zero if unavailable.
    pub arm_end: usize,
    /// Estimated start of the right arm (inclusive) in read coordinates.
    /// Derived from the leftmost right-arm minimizer position.
    /// Zero if unavailable.
    pub arm_start_right: usize,
    /// Set only for candidates found by the alignment-verified rescue path
    /// (see [`crate::cfg::RescueCfg`]). Such a candidate already carries a
    /// base-level arm alignment, so it must not go through
    /// `refine_breakpoint` (which needs far more flanking sequence than a
    /// short or near-end foldback has); use [`refine_detected`] instead.
    pub verified: Option<Verified>,
}

/// Outcome of the rescue path's arm alignment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Verified {
    /// Identity of the aligned arm window (1 - edit distance / window).
    pub identity: f32,
    /// Unpaired sequence between the two arms (bp); 0 for a clean hairpin.
    pub gap: usize,
    /// Length of the arm window that was aligned (bp).
    pub arm: usize,
}

impl Verified {
    fn as_refined(&self, split_pos: usize) -> Refined {
        Refined {
            split_pos,
            score: (self.arm as f32 * (2.0 * self.identity - 1.0)) as i64,
            identity_est: self.identity,
            gap_est: self.gap,
        }
    }
}

/// Refine a detected foldback, whichever path found it: rescue candidates
/// reuse the alignment that verified them, everything else goes through
/// [`refine_breakpoint`] exactly as before.
pub fn refine_detected(
    seq: &[u8],
    fb: &FoldBreakpoint,
    cfg: &SharedCfg,
    refine_scratch: &mut RefineScratch,
) -> Result<Option<Refined>> {
    match fb.verified {
        Some(v) => Ok(Some(v.as_refined(fb.split_pos))),
        None => refine_breakpoint(
            seq,
            fb.split_pos,
            cfg,
            refine_scratch,
            fb.span,
            fb.arm_end,
            fb.arm_start_right,
        ),
    }
}

/// Refine a coarse foldback breakpoint into a more precise split position.
///
/// `span_hint` is the coarse detector's estimate of the palindromic arm span
/// (typically `fb.span`).  It is used to limit the identity comparison window
/// so that identity reflects only the palindromic region rather than being
/// diluted by unrelated sequence beyond the arm.  Pass `usize::MAX` to disable.
///
/// `arm_end` and `arm_start_right` are the minimizer-chain arm boundary estimates
/// from the coarse detector (`fb.arm_end`, `fb.arm_start_right`).  When the gap
/// between them exceeds the δ-scan ceiling, the identity comparison is placed
/// directly at the arm boundaries rather than using the symmetric δ-shift.
/// Pass `0, 0` to disable (falls back to δ-shift behaviour).
///
/// The configured arm is clamped to the sequence space available around `s0`
/// (keeping `window` bp clear on each side), so reads whose split sits close
/// to either end are still refined with a shorter arm rather than being dropped.
/// Returns `None` only when fewer than 10 bp are available (degenerate case).
pub fn refine_breakpoint(
    seq: &[u8],
    s0: usize,
    cfg: &SharedCfg,
    refine_scratch: &mut RefineScratch,
    span_hint: usize,
    arm_end: usize,
    arm_start_right: usize,
) -> Result<Option<Refined>> {
    let n = seq.len();
    let w = cfg.refine.window;

    // Maximum arm that fits while keeping `window` clear on each side.
    let available = s0.min(n.saturating_sub(s0)).saturating_sub(w);
    let arm = cfg.refine.arm.min(available);

    if arm < 10 {
        return Ok(None);
    }

    // Use a local cfg with the (possibly clamped) arm so the individual
    // refiners don't need to repeat this logic.
    if arm == cfg.refine.arm {
        // No clamping needed — avoid the clone on the hot path.
        return Ok(match cfg.refine.mode {
            RefineMode::HiFi => refine_breakpoint_hamming(seq, s0, cfg, refine_scratch, span_hint, arm_end, arm_start_right),
            RefineMode::ONT  => refine_breakpoint_banded_ed(seq, s0, cfg, refine_scratch, span_hint),
        });
    }

    let mut clamped = cfg.clone();
    clamped.refine.arm = arm;
    Ok(match cfg.refine.mode {
        RefineMode::HiFi => refine_breakpoint_hamming(seq, s0, &clamped, refine_scratch, span_hint, arm_end, arm_start_right),
        RefineMode::ONT  => refine_breakpoint_banded_ed(seq, s0, &clamped, refine_scratch, span_hint),
    })
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
pub fn refine_breakpoint_hamming(
    seq: &[u8],
    s0: usize,
    cfg: &SharedCfg,
    scratch: &mut RefineScratch,
    span_hint: usize,
    arm_end: usize,
    arm_start_right: usize,
) -> Option<Refined> {
    let mut best_s = s0;
    let mut best_score: i64 = i64::MIN;
    let arm = cfg.refine.arm;

    // Reformulate the inner loop to use forward-only indexing so the compiler
    // can autovectorise it.
    //
    // Original: score += (l[arm-1-i] == comp(seq[s+i]))
    //         = (comp(l[arm-1-i]) == seq[s+i])   (comp is self-inverse)
    //         = (revcomp(l)[i]     == seq[s+i])
    //
    // We maintain revcomp(seq[s-arm..s]) in `scratch.left_rc`.  As s advances
    // by 1 the new entry is comp(seq[s]) prepended; the last entry falls off:
    //   left_rc[0..arm-1] = old left_rc[1..arm]
    //   left_rc[arm-1]    = comp(seq[s])
    // Implemented as a left-rotate of 1 + overwrite of the last element.
    let s_start = s0 - cfg.refine.window;
    revcomp_into(&mut scratch.left_rc, &seq[s_start - arm..s_start]);

    for s in s_start..=(s0 + cfg.refine.window) {
        // Update left_rc for this s: shift left by 1, append comp(seq[s-1]).
        // (At s_start the buffer is already correct from the init above.)
        if s > s_start {
            // Shift right by 1: left_rc[1..arm] ← left_rc[0..arm-1].
            // Prepend comp(seq[s-1]) at [0].
            // After update: left_rc[i] = comp(seq[s-1-i])  ✓
            scratch.left_rc.copy_within(0..arm - 1, 1);
            scratch.left_rc[0] = comp(seq[s - 1]);
        }

        let right = &seq[s..s + arm];
        let left_rc = &scratch.left_rc;

        let mut score: i64 = 0;
        for i in 0..arm {
            // Forward-indexed comparison; both slices advance monotonically.
            if left_rc[i] == right[i] {
                score += 1;
            } else {
                score -= 1;
            }
            // Early exit: remaining bases cannot beat best_score even if all match.
            if score + (arm - 1 - i) as i64 <= best_score {
                score = i64::MIN / 2; // mark as dead
                break;
            }
        }

        if score > best_score {
            best_score = score;
            best_s = s;
        }
    }

    // ── symmetric gap-aware identity ─────────────────────────────────────────
    // `split_pos` is the anti-diagonal midpoint: for [A, L1][gap, g][RC(A'), L2]
    // it sits at L1 + g/2, not at the true junction (L1).  Shifting both arms
    // symmetrically by δ moves the comparison window so that:
    //   left  arm ends   at split − δ  (approaching true end of A)
    //   right arm starts at split + δ  (approaching true start of RC(A'))
    // At δ = g/2 the arms are correctly aligned.  gap_est = 2*δ_best ≈ g.
    //
    // HiFi identity is estimated with Hamming (no banded ED): for CCS reads
    // the indel rate is <0.1%, so over a 500 bp window at most 1–2 indels are
    // expected.  Hamming gives identity accurate to ±0.2%, well within the
    // margin of any classification threshold.
    //
    // ── identity via arm hint (direct placement) or δ-scan (symmetric shift) ──
    //
    // When the minimizer chain provides valid arm boundaries and the gap
    // (arm_start_right - arm_end) exceeds the δ-scan ceiling, the δ-scan
    // cannot reach the true arm position and would measure identity inside the
    // gap region instead.  In that case, skip the δ-scan entirely and place
    // the comparison windows directly at the arm boundaries from the chain.
    //
    // For small/zero gaps (or when the hint is unavailable), run the δ-scan.
    const MIN_EFF_ARM: usize = 100;
    // `arm` was bound earlier for the split scan.
    let max_delta = (cfg.refine.max_jump_clip / 2).min(arm.saturating_sub(MIN_EFF_ARM));

    let gap_from_hint = arm_start_right.saturating_sub(arm_end);
    let use_hint = arm_end > 0
        && arm_start_right > arm_end
        && gap_from_hint > max_delta * 2
        && arm_end <= seq.len()
        && arm_start_right < seq.len();

    // Cap the identity comparison window.
    const ID_WINDOW_CAP: usize = 500;

    let (identity_est, final_gap_est) = if use_hint {
        // Direct arm placement: compare the innermost `id_window` bp of each arm.
        // Clamp to available sequence so neither window overruns the slice.
        let id_window = arm.min(span_hint).min(ID_WINDOW_CAP).max(1);
        let avail_left  = arm_end.min(id_window);
        let avail_right = seq.len().saturating_sub(arm_start_right).min(id_window);
        let actual = avail_left.min(avail_right).max(1);
        let id_left  = &seq[arm_end - actual .. arm_end];
        let id_right = &seq[arm_start_right .. arm_start_right + actual];
        scratch.right_rc.clear();
        revcomp_into(&mut scratch.right_rc, id_right);
        let band = ((cfg.refine.max_ed_rate * actual as f32).ceil() as usize).max(1);
        let ed = edit_distance_myers(id_left, &scratch.right_rc, band);
        let ident = if ed <= band {
            (1.0 - ed as f32 / actual as f32).clamp(0.0, 1.0)
        } else {
            id_left.iter().zip(scratch.right_rc.iter()).filter(|(a, b)| a == b).count() as f32
                / actual as f32
        };
        (ident, gap_from_hint)
    } else {
        // Recompute left_rc for best_s (the split scan left it for s_start+window).
        // left_rc[i] = comp(seq[best_s - 1 - i])  for i in 0..arm.
        // Enables forward-indexed comparison: left_rc[δ+i] == seq[best_s+δ+i].
        revcomp_into(&mut scratch.left_rc, &seq[best_s - arm .. best_s]);

        // Symmetric δ-shift: scan for best gap estimate (Hamming, forward-indexed,
        // autovectorisable). Only used to rank deltas, not for the final identity.
        let eval = |delta: usize, lrc: &[u8]| -> f32 {
            let eff_arm = (arm - delta).min(ID_WINDOW_CAP);
            let left_slice  = &lrc[delta .. delta + eff_arm];
            let right_slice = &seq[best_s + delta .. best_s + delta + eff_arm];
            let matches = left_slice.iter().zip(right_slice).filter(|(a, b)| a == b).count();
            matches as f32 / eff_arm as f32
        };

        let base_ident = eval(0, &scratch.left_rc);
        let mut best_delta = 0usize;
        let mut best_ident  = base_ident;

        if max_delta > 0 {
            let stride = (max_delta / 200).max(1);

            // Stage 1: coarse scan.
            let mut coarse_best = 0usize;
            for delta in (1..=max_delta).step_by(stride) {
                let ident = eval(delta, &scratch.left_rc);
                if ident > best_ident {
                    best_ident  = ident;
                    coarse_best = delta;
                    best_delta  = delta;
                }
            }

            // Stage 2: exhaustive refine around the coarse best.
            // Skip when stride==1: the coarse scan was already exhaustive.
            if coarse_best > 0 && stride > 1 {
                let lo = coarse_best.saturating_sub(stride).max(1);
                let hi = (coarse_best + stride).min(max_delta);
                for delta in lo..=hi {
                    let ident = eval(delta, &scratch.left_rc);
                    if ident > best_ident {
                        best_ident = ident;
                        best_delta = delta;
                    }
                }
            }
        }
        let _ = best_ident; // used only to drive the scan above

        // Final identity: Myers' bit-parallel ED on capped window (≤ID_WINDOW_CAP bp).
        // Window is always ≤512 bp, within the 8-word bitvector limit.
        let eff_arm = arm - best_delta;
        let id_window = eff_arm.min(span_hint).min(ID_WINDOW_CAP).max(1);
        let id_left  = &seq[best_s - best_delta - id_window .. best_s - best_delta];
        let id_right = &seq[best_s + best_delta .. best_s + best_delta + id_window];
        scratch.right_rc.clear();
        revcomp_into(&mut scratch.right_rc, id_right);
        let band = ((cfg.refine.max_ed_rate * id_window as f32).ceil() as usize).max(1);
        let ed = edit_distance_myers(id_left, &scratch.right_rc, band);
        let ident = if ed <= band {
            (1.0 - ed as f32 / id_window as f32).clamp(0.0, 1.0)
        } else {
            // Band exceeded: Hamming lower bound using pre-computed left_rc.
            let left_slice  = &scratch.left_rc[best_delta .. best_delta + id_window];
            left_slice.iter().zip(id_right.iter()).filter(|(a, b)| a == b).count() as f32
                / id_window as f32
        };
        (ident, best_delta * 2)
    };

    Some(Refined {
        split_pos: best_s,
        score: best_score,
        identity_est,
        gap_est: final_gap_est,
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
    span_hint: usize,
) -> Option<Refined> {
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

        let ed = edit_distance_myers(left, &scratch.right_rc, band);

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

    // ── symmetric gap-aware identity (same model as Hamming refiner) ─────────
    // See Hamming refiner comment for full rationale.  Symmetric δ-shift:
    //   left  = seq[split − arm .. split − δ]  (length arm − δ)
    //   right = seq[split + δ  .. split + arm]  (length arm − δ)
    // Bounds identical to δ=0 case, already validated by outer check.
    const MIN_EFF_ARM: usize = 100;
    let arm = cfg.refine.arm;
    let max_delta = (cfg.refine.max_jump_clip / 2).min(arm.saturating_sub(MIN_EFF_ARM));

    // Pre-compute revcomp of the full right arm once.
    // revcomp(seq[best_s+delta..best_s+arm]) is the prefix [0..arm-delta] of this.
    revcomp_into(&mut scratch.right_rc, &seq[best_s..best_s + arm]);

    let eval_ed = |delta: usize| -> f32 {
        let eff_arm = arm - delta;
        let left     = &seq[best_s - arm .. best_s - delta];
        let right_rc = &scratch.right_rc[..eff_arm];
        let ed = edit_distance_myers(left, right_rc, band);
        if ed > band { return 0.0; }
        (1.0 - (ed as f32 / eff_arm as f32)).clamp(0.0, 1.0)
    };

    // δ=0 identity already known from best_ed.
    let mut best_delta = 0usize;
    let mut best_ident  = (1.0 - (best_ed as f32 / arm as f32)).clamp(0.0, 1.0);

    if max_delta > 0 {
        let stride = (max_delta / 200).max(1);

        // Stage 1: coarse scan.
        let mut coarse_best = 0usize;
        for delta in (1..=max_delta).step_by(stride) {
            let ident = eval_ed(delta);
            if ident > best_ident {
                best_ident  = ident;
                coarse_best = delta;
                best_delta  = delta;
            }
        }

        // Stage 2: exhaustive refine around the coarse best.
        // Skip when stride==1: the coarse scan was already exhaustive.
        if coarse_best > 0 && stride > 1 {
            let lo = coarse_best.saturating_sub(stride).max(1);
            let hi = (coarse_best + stride).min(max_delta);
            for delta in lo..=hi {
                let ident = eval_ed(delta);
                if ident > best_ident {
                    best_ident = ident;
                    best_delta = delta;
                }
            }
        }
    }

    // ── span-limited final identity ───────────────────────────────────────────
    let eff_arm = arm - best_delta;
    let id_window = eff_arm.min(span_hint).max(1);
    let left_short     = &seq[best_s - best_delta - id_window .. best_s - best_delta];
    // revcomp(seq[best_s+best_delta..best_s+best_delta+id_window]) is
    // scratch.right_rc[arm-best_delta-id_window..arm-best_delta], already computed.
    let right_rc_short = &scratch.right_rc[arm - best_delta - id_window .. arm - best_delta];
    let final_ed = edit_distance_myers(left_short, right_rc_short, band);
    let identity_est = if final_ed <= band {
        (1.0 - final_ed as f32 / id_window as f32).clamp(0.0, 1.0)
    } else {
        // Band exceeded on short window; fall back to Hamming.
        let mut m = 0i64;
        for i in 0..id_window {
            if seq[best_s - best_delta - 1 - i] == comp(seq[best_s + best_delta + i]) { m += 1; }
        }
        m as f32 / id_window as f32
    };

    Some(Refined {
        split_pos: best_s,
        score: -(best_ed as i64),
        identity_est,
        gap_est: best_delta * 2,
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
    // Anything the evidence gates accept behaves exactly as it always has;
    // the rescue path only ever sees reads the gates turned away.
    if let Some(fb) = detect_foldback_gated(seq, shared, fold, scratch) {
        return Some(fb);
    }
    if fold.rescue.enabled {
        return rescue_foldback(seq, shared, fold, scratch);
    }
    None
}

/// Rescue path: propose candidates from very little seed evidence, then
/// confirm each by aligning one arm against the reverse complement of the
/// other. See [`crate::cfg::RescueCfg`] for the rationale.
///
/// Relies on the matching state `detect_foldback_gated` left in `scratch`
/// (`stats`, `idx_f`, `repetitive`, `pos_rc`, `val_rc`) for this same `seq`.
fn rescue_foldback(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    scratch: &mut FoldScratch,
) -> Option<FoldBreakpoint> {
    if let Some(fb) = rescue_seeded(seq, shared, fold, scratch) {
        return Some(fb);
    }
    if fold.rescue.terminal {
        return rescue_terminal(seq, shared, fold, scratch);
    }
    None
}

/// Longest read-end arm the end-anchored rescue paths accept.
///
/// The evidence gates refuse any junction within `end_guard` of a read end,
/// however much evidence it has, so the rescue paths have to reach at least
/// that far in or foldbacks with arms between the two limits fall through
/// both (on ONT, with `end_guard` 1000, that gap held half of all foldbacks
/// with 500-1000 bp arms).
#[inline]
fn max_tail(shared: &SharedCfg, rc: &crate::cfg::RescueCfg) -> usize {
    let base = rc.max_window.max(shared.end_guard);
    base + base / 10
}

/// Length of a read-end arm to actually align: all of it bar the outermost
/// `end_slack` bases, unless that would leave less than `min_arm`.
#[inline]
fn trim_end_slack(tail: usize, rc: &crate::cfg::RescueCfg) -> usize {
    if tail >= rc.min_arm + rc.end_slack {
        tail - rc.end_slack
    } else {
        tail
    }
}

/// 2-bit code of a base, or `None` for anything that is not ACGT.
#[inline]
fn base2(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// A terminal candidate in window-local coordinates.
struct TermHit {
    arm_end: usize,
    arm_start: usize,
    tail: usize,
    identity: f32,
    seeds: usize,
    /// `p + q + ks` of the anchoring seed pair: twice the fold axis.
    axis2: i32,
}

/// Dense-seed scan of one read end. `win` is the first or last
/// `term_window` bp of the read; `at_right_end` says which.
fn terminal_candidate(
    win: &[u8],
    at_right_end: bool,
    shared: &SharedCfg,
    rc: &crate::cfg::RescueCfg,
    scratch: &mut FoldScratch,
) -> Option<TermHit> {
    let ks = rc.term_k.clamp(6, 15);
    let len = win.len();
    if len < 2 * rc.min_arm + ks {
        return None;
    }
    let mask: u32 = (1u32 << (2 * ks)) - 1;
    let shift = 2 * (ks as u32 - 1);

    // Forward k-mers, sorted by code.
    scratch.term_fwd.clear();
    let (mut code, mut valid) = (0u32, 0usize);
    for (i, &b) in win.iter().enumerate() {
        match base2(b) {
            Some(x) => {
                code = ((code << 2) | x) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= ks {
            scratch.term_fwd.push((code, (i + 1 - ks) as u32));
        }
    }
    scratch.term_fwd.sort_unstable();

    // For every position q, look up the reverse complement of its k-mer:
    // a hit at p means win[q..q+ks] is the inverted copy of win[p..p+ks].
    const MAX_OCC: usize = 4;
    scratch.term_pairs.clear();
    let (mut rcode, mut valid) = (0u32, 0usize);
    for (i, &b) in win.iter().enumerate() {
        match base2(b) {
            Some(x) => {
                rcode = (rcode >> 2) | ((3 - x) << shift);
                valid += 1;
            }
            None => valid = 0,
        }
        if valid < ks {
            continue;
        }
        let q = (i + 1 - ks) as i32;
        let lo = scratch.term_fwd.partition_point(|&(c, _)| c < rcode);
        let mut hi = lo;
        while hi < scratch.term_fwd.len() && scratch.term_fwd[hi].0 == rcode {
            hi += 1;
        }
        // k-mers frequent within the window are low-complexity / tandem
        // sequence, which pairs with itself on every anti-diagonal.
        if hi == lo || hi - lo > MAX_OCC {
            continue;
        }
        for j in lo..hi {
            let p = scratch.term_fwd[j].1 as i32;
            if p + ks as i32 <= q {
                scratch.term_pairs.push((p + q, p, q));
            }
        }
    }
    if scratch.term_pairs.len() < rc.term_min_seeds {
        return None;
    }
    scratch.term_pairs.sort_unstable();

    // Densest run of anti-diagonals within the drift tolerance.
    let d_w: i32 = match shared.refine.mode {
        RefineMode::HiFi => 8,
        RefineMode::ONT => 40,
    };
    let pairs = &scratch.term_pairs;
    let (mut best_i, mut best_j) = (0usize, 0usize);
    let mut j = 0usize;
    for i in 0..pairs.len() {
        if j < i {
            j = i;
        }
        while j < pairs.len() && pairs[j].0 - pairs[i].0 <= d_w {
            j += 1;
        }
        if j - i > best_j - best_i {
            best_i = i;
            best_j = j;
        }
    }
    let run = &pairs[best_i..best_j];
    if run.len() < rc.term_min_seeds {
        return None;
    }
    // Seeds must come from distinct places on the arm, not one spot hit
    // repeatedly.
    let (mut p_lo, mut p_hi) = (i32::MAX, i32::MIN);
    for &(_, p, _) in run {
        p_lo = p_lo.min(p);
        p_hi = p_hi.max(p);
    }
    if ((p_hi - p_lo) as usize) < rc.term_min_seeds * 4 {
        return None;
    }
    // Innermost seed anchors the alignment.
    let &(_, ap, aq) = run.iter().max_by_key(|&&(_, p, _)| p)?;
    let arm_end = ap as usize + ks;
    let arm_start = aq as usize;

    // The inverted copy must run to the read end, and be short enough that
    // the evidence gates could not have been expected to see it.
    let tail = if at_right_end { len - arm_start } else { arm_end };
    if tail < rc.min_arm || tail > max_tail(shared, rc) {
        return None;
    }
    let other = if at_right_end { arm_end } else { len - arm_start };
    if other < tail {
        return None;
    }
    // Align from the junction outwards, stopping `end_slack` short of the
    // read end.
    let aln = trim_end_slack(tail, rc);
    let left = &win[arm_end - aln..arm_end];
    let right = &win[arm_start..arm_start + aln];
    revcomp_into(&mut scratch.refine.right_rc, right);
    let band = ((1.0 - rc.min_identity).max(0.0) * aln as f32).ceil() as usize + 1;
    let ed = edit_distance_myers(left, &scratch.refine.right_rc, band);
    if ed > band {
        return None;
    }
    let identity = 1.0 - ed as f32 / aln as f32;
    if identity < rc.min_identity {
        return None;
    }
    Some(TermHit {
        arm_end,
        arm_start,
        tail,
        identity,
        seeds: run.len(),
        axis2: ap + aq + ks as i32,
    })
}

/// Terminal rescue: dense-seed scan of both read ends. See
/// [`crate::cfg::RescueCfg::terminal`].
fn rescue_terminal(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    scratch: &mut FoldScratch,
) -> Option<FoldBreakpoint> {
    let rc = &fold.rescue;
    let n = seq.len();
    // Both arms of the longest tail we accept, plus room for a spacer.
    let t = rc.term_window.max(2 * max_tail(shared, rc) + 200).min(n);
    let mut best: Option<(f32, usize, TermHit)> = None;
    // (window start, is this the right-hand end?)
    let ends: [(usize, bool); 2] = [(n - t, true), (0, false)];
    for &(start, at_right_end) in ends.iter() {
        let win = &seq[start..start + t];
        if let Some(hit) = terminal_candidate(win, at_right_end, shared, rc, scratch) {
            let score = hit.tail as f32 * hit.identity;
            if best.as_ref().map_or(true, |(bs, _, _)| score > *bs) {
                best = Some((score, start, hit));
            }
        }
    }
    let (_, start, hit) = best?;
    let arm_end = start + hit.arm_end;
    let arm_start = start + hit.arm_start;
    let gap = arm_start - arm_end;

    // Hand the support gate whatever minimizer matchpoints lie on this
    // fold's anti-diagonal (there may be too few for a key; then the
    // junction is simply unsupported).
    let k = shared.minimizer.k as i32;
    let target = 2 * start as i32 + hit.axis2 - k;
    let d_tol: i32 = match shared.refine.mode {
        RefineMode::HiFi => 8,
        RefineMode::ONT => 40,
    };
    collect_diag_points(n, shared, scratch, target, d_tol);

    Some(FoldBreakpoint {
        split_pos: (arm_end + arm_start) / 2,
        score: (2 * hit.tail + gap) as i64 + hit.seeds as i64 * 10,
        matches: hit.seeds,
        span: 2 * hit.tail + gap,
        arm_end,
        arm_start_right: arm_start,
        verified: Some(Verified {
            identity: hit.identity,
            gap,
            arm: hit.tail,
        }),
    })
}

/// Like `collect_bin_points`, but selects matchpoints by distance from a
/// target anti-diagonal rather than by bin. Leaves the buffers empty if the
/// read never reached the matching stage.
fn collect_diag_points(
    seq_len: usize,
    shared: &SharedCfg,
    scratch: &mut FoldScratch,
    target: i32,
    d_tol: i32,
) {
    scratch.best_pts.clear();
    scratch.best_matches.clear();
    scratch.best_trip.clear();
    if scratch.stats.is_empty() {
        return;
    }
    let n = seq_len as i32;
    let k_i32 = shared.minimizer.k as i32;
    for (&prc, &vrc) in scratch.pos_rc.iter().zip(scratch.val_rc.iter()) {
        if scratch.repetitive.contains_key(&vrc) {
            continue;
        }
        let Some(p1s) = scratch.idx_f.get(&vrc) else {
            continue;
        };
        let p2 = n - k_i32 - (prc as i32);
        for &p1 in p1s {
            if (p1 + p2 - target).abs() <= d_tol {
                scratch.best_pts.push((p1, p2));
                scratch.best_matches.push((p1 as usize, vrc));
                scratch.best_trip.push((p1, p2, vrc));
            }
        }
    }
}

/// Minimizer-seeded rescue: candidates are anti-diagonal bins the gates
/// rejected.
fn rescue_seeded(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    scratch: &mut FoldScratch,
) -> Option<FoldBreakpoint> {
    let rc = &fold.rescue;
    if scratch.stats.is_empty() {
        return None;
    }
    let k = shared.minimizer.k as i32;
    let n = seq.len() as i32;
    let tol = shared.fold_diag_tol.max(1);

    // Candidate anti-diagonals: every matched pair contributes two mirrored
    // points, so `min_pairs` pairs means 2 * min_pairs points in the bin.
    scratch.top_bins.clear();
    for (&bin, st) in scratch.stats.iter() {
        if st.count() >= 2 * rc.min_pairs {
            scratch.top_bins.push((bin, st.count() as i64, st.count(), st.span()));
        }
    }
    if scratch.top_bins.is_empty() {
        return None;
    }
    scratch
        .top_bins
        .sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    scratch.top_bins.truncate(rc.max_candidates.max(1));

    // How far a point's anti-diagonal may sit from the candidate's median
    // and still count as the same fold: near-exact for accurate reads,
    // looser where indels make it drift along the arm.
    let d_tol: i32 = match shared.refine.mode {
        RefineMode::HiFi => 8,
        RefineMode::ONT => 40,
    };

    let mut best: Option<(f32, i32, FoldBreakpoint)> = None;
    let n_cand = scratch.top_bins.len();
    for ci in 0..n_cand {
        let bin = scratch.top_bins[ci].0;
        collect_bin_points(seq.len(), shared, scratch, bin, tol);

        // Left-arm copy first, right-arm copy second.
        scratch.jk_tmp.clear();
        for &(p1, p2, _) in scratch.best_trip.iter() {
            if p1 < p2 {
                scratch.jk_tmp.push((p1 + p2, 0));
            }
        }
        if scratch.jk_tmp.len() < rc.min_pairs {
            continue;
        }
        scratch.jk_tmp.sort_unstable();
        let d_med = scratch.jk_tmp[scratch.jk_tmp.len() / 2].0;

        // Points consistent with one fold axis; innermost one anchors the
        // alignment, outermost one bounds the seeded arm.
        let mut inner: Option<(i32, i32)> = None;
        let mut p1_min = i32::MAX;
        let mut p2_max = i32::MIN;
        scratch.jk_tmp.clear();
        for &(p1, p2, v) in scratch.best_trip.iter() {
            if p1 < p2 && (p1 + p2 - d_med).abs() <= d_tol {
                scratch.jk_tmp.push((0, v));
                // The anchor must be a pair of two separate copies. In a
                // clean hairpin a k-mer straddling the axis also matches
                // its own (overlapping) reverse complement; that is real
                // evidence of the fold but cannot anchor an arm alignment.
                if p1 + k <= p2 && inner.map_or(true, |(ip1, _)| p1 > ip1) {
                    inner = Some((p1, p2));
                }
                p1_min = p1_min.min(p1);
                p2_max = p2_max.max(p2);
            }
        }
        let n_pairs = scratch.jk_tmp.len();
        // Distinct minimizer values, not just distinct positions: a short
        // tandem repeat matches itself on many anti-diagonals with one or
        // two values over and over, and is not evidence of a fold.
        scratch.jk_tmp.sort_unstable();
        scratch.jk_tmp.dedup();
        let n_distinct = scratch.jk_tmp.len();
        if n_distinct < rc.min_pairs || n_distinct * 2 < n_pairs {
            continue;
        }
        let Some((ip1, ip2)) = inner else { continue };

        let arm_end = ip1 + k; // exclusive end of the left arm
        let arm_start = ip2; // inclusive start of the right arm
        if arm_start < arm_end {
            // the two copies overlap across the axis: a k-mer that is its
            // own reverse complement, not a pair of arms
            continue;
        }
        let seeded = arm_end - p1_min;
        let avail = arm_end.min(n - arm_start);
        let win = if rc.end_only {
            // The shorter arm must be cut off by a read end, and all of it
            // has to align: `avail` is exactly that arm's length.
            if avail > max_tail(shared, rc) as i32 {
                continue;
            }
            trim_end_slack(avail as usize, rc) as i32
        } else {
            (rc.max_window as i32)
                .min(avail)
                .min(seeded.max(rc.min_arm as i32))
        };
        if win < rc.min_arm as i32 {
            continue;
        }
        let win_u = win as usize;
        let left = &seq[(arm_end - win) as usize..arm_end as usize];
        let right = &seq[arm_start as usize..(arm_start + win) as usize];
        revcomp_into(&mut scratch.refine.right_rc, right);
        let band = ((1.0 - rc.min_identity).max(0.0) * win_u as f32).ceil() as usize + 1;
        let ed = edit_distance_myers(left, &scratch.refine.right_rc, band);
        if ed > band {
            continue;
        }
        let identity = 1.0 - ed as f32 / win_u as f32;
        if identity < rc.min_identity {
            continue;
        }

        let score = win_u as f32 * identity;
        if best.as_ref().map_or(true, |(bs, _, _)| score > *bs) {
            let span = (p2_max + k - p1_min).max(0) as usize;
            best = Some((
                score,
                bin,
                FoldBreakpoint {
                    split_pos: ((arm_end + arm_start) / 2) as usize,
                    score: span as i64 + (n_pairs as i64 * 20),
                    matches: n_pairs * 2,
                    span,
                    arm_end: arm_end as usize,
                    arm_start_right: arm_start as usize,
                    verified: Some(Verified {
                        identity,
                        gap: (arm_start - arm_end) as usize,
                        arm: win_u,
                    }),
                },
            ));
        }
    }

    let (_, bin, fb) = best?;
    // Leave scratch describing the accepted candidate, as the gated path
    // does, so junction keys / signatures are computed from the right bin.
    collect_bin_points(seq.len(), shared, scratch, bin, tol);
    Some(fb)
}

/// Fill `best_pts` / `best_matches` / `best_trip` with the matchpoints of
/// one anti-diagonal bin, from the matching state left by
/// `detect_foldback_gated`.
fn collect_bin_points(
    seq_len: usize,
    shared: &SharedCfg,
    scratch: &mut FoldScratch,
    want_bin: i32,
    tol: i32,
) {
    let n = seq_len as i32;
    let k_i32 = shared.minimizer.k as i32;
    scratch.best_pts.clear();
    scratch.best_matches.clear();
    scratch.best_trip.clear();
    for (&prc, &vrc) in scratch.pos_rc.iter().zip(scratch.val_rc.iter()) {
        if scratch.repetitive.contains_key(&vrc) {
            continue;
        }
        let Some(p1s) = scratch.idx_f.get(&vrc) else {
            continue;
        };
        let p2 = n - k_i32 - (prc as i32);
        for &p1 in p1s {
            if div_floor(p1 + p2, tol) == want_bin {
                scratch.best_pts.push((p1, p2));
                scratch.best_matches.push((p1 as usize, vrc));
                scratch.best_trip.push((p1, p2, vrc));
            }
        }
    }
}

/// The evidence-gated coarse detector (the historical `detect_foldback`).
fn detect_foldback_gated(
    seq: &[u8],
    shared: &SharedCfg,
    fold: &FoldOnlyCfg,
    scratch: &mut FoldScratch,
) -> Option<FoldBreakpoint> {
    // Cleared up front (not just before use) so that an early return below
    // can never leave a previous read's bins for the rescue path to find.
    scratch.stats.clear();

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

    // Early exit before the more expensive RC minimizer computation.
    if scratch.pos_f.len() < shared.min_matches {
        return None;
    }

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

    if pos_rc.len() < shared.min_matches {
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
        // Deterministic tie-break: `scratch.stats` is a HashMap, whose
        // iteration order varies run-to-run (random per-process hash seed).
        // Without a tie-break, an exact score tie between two bins would
        // non-deterministically pick whichever happened to be visited first,
        // making detection results (and anything downstream: split position,
        // signature, cut decision) irreproducible across runs on identical
        // input. Prefer the smaller bin id on ties -- arbitrary but stable.
        let better = sc > best_score || (sc == best_score && Some(bin) < best_bin);
        if better {
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
    scratch.best_trip.clear();
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
                scratch.best_trip.push((p1, p2, vrc));
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

    let k_i32 = shared.minimizer.k as i32;
    let p1_min = chain.first().unwrap().0 as i64;
    let p1_max = chain.last().unwrap().0 as i64;
    let span = (p1_max - p1_min).max(0) as usize;
    if span < fold.min_arm {
        return None;
    }

    // Weighted mean of d = p1+p2, using inverse junction-distance as weight.
    //
    // For a foldback at split s, every matchpoint satisfies d ≈ 2s.  But
    // near-junction points (p1 ≈ p2 ≈ s) accumulate fewer indel errors along
    // the arm and give a more accurate estimate than far-from-junction points.
    // Weighting by 1/max(|p1-p2|, k) emphasises these near-junction matches
    // without discarding far matches entirely.
    let (weight_sum, weighted_d_sum) = chain.iter().fold(
        (0.0f64, 0.0f64),
        |(ws, wds), &(p1, p2)| {
            let diff = (p1 - p2).abs().max(k_i32) as f64;
            let w = 1.0 / diff;
            (ws + w, wds + (p1 + p2) as f64 * w)
        },
    );
    let split_u = if weight_sum > 0.0 {
        (weighted_d_sum / weight_sum / 2.0).round().max(0.0) as usize
    } else {
        0
    };

    if split_u < shared.end_guard || split_u + shared.end_guard > len {
        return None;
    }

    // Arm boundary estimates from junction-proximal (type-1) chain entries.
    //
    // The chain can contain two types of matchpoints on the same anti-diagonal:
    //   Type-1: p1 in left arm (< split), p2 in right arm (> split) — the real matches
    //   Type-2: p1 in right arm (> split), p2 in left arm (< split) — cross-matches
    //
    // For a read with a large gap, type-2 cross-matches will have p1 near the right
    // arm END and p2 near the left arm START, making chain.last() a type-2 entry.
    // Using chain.last() would give arm_end ≈ read_end and arm_start_right ≈ read_start,
    // collapsing to gap_from_hint = 0.
    //
    // Instead, filter to type-1 entries (p1 < split, p2 > split) and take the
    // maximum p1 (closest to junction from left) and minimum p2 (closest from right).
    let split_i32 = split_u as i32;
    let mut p1_max_t1 = -1i32;
    let mut p2_min_t1 = i32::MAX;
    for &(p1, p2) in &chain {
        if p1 < split_i32 && p2 > split_i32 {
            if p1 > p1_max_t1 { p1_max_t1 = p1; }
            if p2 < p2_min_t1 { p2_min_t1 = p2; }
        }
    }
    // If no type-1 entries exist (very short arm or degenerate chain), fall back
    // to type-2 chain.last() interpretation (arm_end=0 sentinel disables hint).
    let arm_end = if p1_max_t1 >= 0 {
        (p1_max_t1 as i64 + k_i32 as i64).max(0) as usize
    } else {
        0
    };
    let arm_start_right = if p2_min_t1 < i32::MAX {
        p2_min_t1.max(0) as usize
    } else {
        0
    };

    let score = span as i64 + (chain.len() as i64 * 10);
    Some(FoldBreakpoint {
        split_pos: split_u,
        score,
        matches: chain.len(),
        span,
        arm_end,
        arm_start_right,
        verified: None,
    })
}

/// Where a cut may end on each side of a detected foldback.
///
/// Returns `(left_end, right_start)` in the same coordinates as the inputs:
/// keeping the left side means keeping `..left_end`, keeping the right side
/// means keeping `right_start..`.
///
/// The policy is chosen once per run by `set_cut_spacer_policy`; see
/// `CutSpacer`.
pub fn cut_boundaries(
    seq: &[u8],
    trips: &[(i32, i32, u64)],
    k: usize,
    split: usize,
    gap_est: usize,
    chain_arm_end: usize,
    chain_arm_start: usize,
) -> (usize, usize) {
    let len = seq.len();
    let half = gap_est / 2;
    let lo = split.saturating_sub(half);
    let hi = (split + half).min(len);
    let policy = cut_spacer_policy();
    if policy == CutSpacer::Resolve {
        if let Some((m, side)) = crate::spacer::resolve(seq, trips, k) {
            if (m.a + m.b).abs_diff(2 * split) <= 200 {
                return match side {
                    crate::spacer::SpacerSide::Left => (m.b, m.b),
                    crate::spacer::SpacerSide::Right => (m.a, m.a),
                    crate::spacer::SpacerSide::Unknown => (m.a, m.b),
                };
            }
        }
    }
    match policy {
        CutSpacer::Keep => (hi, lo),
        CutSpacer::Drop => (lo, hi),
        CutSpacer::Chain | CutSpacer::Resolve => {
            if let Some((a, b)) = chain_boundaries(split, chain_arm_end, chain_arm_start, len) {
                (a, b)
            } else {
                (lo, hi)
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CutSpacer {
    /// Legacy: keep the un-mirrored stretch between the arms on whichever
    /// side is kept.
    Keep,
    /// Drop the un-mirrored stretch; boundaries from split +- gap_est/2.
    Drop,
    /// Drop it; boundaries from the innermost matched minimizers of the
    /// coarse chain, falling back to `Drop` when those look inconsistent.
    Chain,
    /// Use exact mirror boundaries, and keep the spacer with whichever side
    /// the rest of the library says it belongs to (see `crate::spacer`);
    /// drop it when there is no evidence. Falls back to `Chain`.
    Resolve,
}

static CUT_SPACER: std::sync::OnceLock<CutSpacer> = std::sync::OnceLock::new();

/// Developer override for benchmarking: `MDAX_CUT_SPACER=keep|drop|chain|resolve`.
fn cut_spacer_env() -> Option<CutSpacer> {
    match std::env::var("MDAX_CUT_SPACER").as_deref() {
        Ok("keep") => Some(CutSpacer::Keep),
        Ok("drop") => Some(CutSpacer::Drop),
        Ok("chain") => Some(CutSpacer::Chain),
        Ok("resolve") => Some(CutSpacer::Resolve),
        _ => None,
    }
}

/// Choose the spacer policy for this run (first call wins). `Resolve` needs
/// exact mirror boundaries, which HiFi reads give and noisy reads do not:
/// on ONT its boundary k-mers were mis-assigned often enough to keep 8.7% of
/// cut reads' bases as chimera (vs 0.5% for `Chain`), so noisy reads use
/// `Chain`. `MDAX_CUT_SPACER` overrides.
pub fn set_cut_spacer_policy(mode: RefineMode) -> CutSpacer {
    *CUT_SPACER.get_or_init(|| {
        cut_spacer_env().unwrap_or(match mode {
            RefineMode::HiFi => CutSpacer::Resolve,
            RefineMode::ONT => CutSpacer::Chain,
        })
    })
}

/// The policy in force: whatever `set_cut_spacer_policy` chose, else the
/// env override, else `Chain` (no library evidence is available in that
/// case, e.g. in unit tests calling the cutter directly).
pub fn cut_spacer_policy() -> CutSpacer {
    *CUT_SPACER.get_or_init(|| cut_spacer_env().unwrap_or(CutSpacer::Chain))
}

/// Arm boundaries from the coarse chain's innermost type-1 matchpoints, if
/// they bracket the refined split plausibly.
fn chain_boundaries(split: usize, arm_end: usize, arm_start: usize, len: usize) -> Option<(usize, usize)> {
    const AXIS_TOL: usize = 100;
    if arm_end == 0 || arm_start == 0 || arm_start > len || arm_end > len {
        return None;
    }
    // arm_end includes the last matched k-mer, so it can sit a few bases past
    // arm_start on a clean hairpin; that is still a valid boundary pair.
    if arm_end > split + AXIS_TOL || arm_start + AXIS_TOL < split {
        return None;
    }
    let axis2 = arm_end + arm_start;
    if axis2.abs_diff(2 * split) > 2 * AXIS_TOL {
        return None;
    }
    if arm_end > arm_start {
        // overlapping innermost k-mers: an undisplaced hairpin, cut at the axis
        return Some((split, split));
    }
    Some((arm_end, arm_start))
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
/// - `scratch.best_trip` is assumed to correspond to the most recent
///   `detect_foldback` call.
pub fn recursive_foldback_cut_from_first_range(
    seq: &[u8],               // original full sequence
    mut keep: Range<usize>,   // current keep window into `seq`
    first_fb: FoldBreakpoint, // split_pos is assumed relative to `seq` (full coords)
    first_rf: Refined,        // split_pos is assumed relative to `seq` (full coords)
    cfg: &MdaxCfg,
    support: &JunctionSupport,
    max_depth: usize,
    scratch: &mut FoldScratch,
    _sig_scratch: &mut SigScratch,
) -> anyhow::Result<Range<usize>> {
    // We already have the first hit/refine (full coords). After first cut, detect/refine on the view.
    let mut fb_opt = Some(first_fb);
    let mut rf_opt = Some(first_rf);

    // Bound passed to the recursive detector on *both* sides. It tracks the
    // true arm boundary (excluding any gap) on whichever side was last cut,
    // since the gap region can look reverse-complement similar to the arm
    // and would otherwise confuse detection on the next pass. On the first
    // iteration this equals `keep` exactly.
    let mut detect: Range<usize> = keep.clone();

    for _ in 0..max_depth {
        // View of the sequence we are currently considering.
        let view = &seq[detect.start..detect.end];

        // --- 1) detect foldback (skip for first iteration) ---
        let fb = match fb_opt.take() {
            Some(mut fb) => {
                // Convert full-coord positions into view-local coords.
                // (If detect.start==0, these are no-ops.)
                if fb.split_pos >= detect.start {
                    fb.split_pos -= detect.start;
                } else {
                    // If inconsistent, bail rather than panic.
                    break;
                }
                fb.arm_end = fb.arm_end.saturating_sub(detect.start);
                fb.arm_start_right = fb.arm_start_right.saturating_sub(detect.start);
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
                if rf.split_pos >= detect.start {
                    rf.split_pos -= detect.start;
                } else {
                    break;
                }
                rf
            }
            None => {
                // `refine_detected`, not `refine_breakpoint`: a foldback
                // found by the rescue path on the shortened view is short
                // and sits at the view's end, which `refine_breakpoint`
                // cannot refine (it needs `arm + window` bp on both sides).
                match refine_detected(view, &fb, &cfg.shared, &mut scratch.refine)? {
                    Some(rf) => rf,
                    None => break,
                }
            }
        };

        if rf.identity_est < cfg.fold2.min_identity {
            break;
        }

        // --- 3) junction-specific support ---
        // `scratch.best_trip` always corresponds to the most recent
        // `detect_foldback` call (pass2's own call on the first iteration),
        // and its keys are coordinate-free (minimizer value + arm
        // separation), so they are valid in any view of the read.
        junction_keys(
            &scratch.best_trip,
            axis_tol(cfg),
            &mut scratch.jk_tmp,
            &mut scratch.jk_out,
        );
        if support.support_n(&scratch.jk_out) >= cfg.fold2.min_support {
            // real -> don't cut
            break;
        }

        // artefact -> cut. Keep whichever side of the split is longer, rather
        // than always keeping the left arm.
        //
        // Real template-switch folds are usually displaced: the two arms are
        // offset on the genome, leaving an un-mirrored stretch between them
        // that belongs to one strand or the other. A single read cannot tell
        // which, and keeping it when it belongs to the dropped side leaves a
        // strand-switched stub (measured median ~800 bp) on the kept read, so
        // drop it: cut at the kept arm's own boundary. Boundaries come from
        // the innermost matched minimizers of the coarse chain when those are
        // consistent with the refined split (gap_est alone is off by >100 bp
        // in about a quarter of reads). All in `view`-local coordinates,
        // converted to full-sequence ones to compare against `keep`.
        let (left_end, right_start) = cut_boundaries(
            view,
            &scratch.best_trip,
            cfg.shared.minimizer.k,
            rf.split_pos,
            rf.gap_est,
            fb.arm_end,
            fb.arm_start_right,
        );
        let abs_arm_end = detect.start + left_end.min(right_start);
        let abs_left_end = detect.start + left_end;
        let abs_right_start = detect.start + right_start;

        // A foldback near the *start* of a long read means "keep the left
        // arm" throws away almost the entire read, even though everything
        // from the right arm onward (which may include a perfectly good
        // trailing tail beyond the palindrome, not just the artefactual
        // reverse-complement duplicate itself) is ordinary, assemblable
        // sequence -- the RC(A) arm is exactly as valid DNA content as A is,
        // just mirrored, so keeping *a* single non-self-overlapping copy is
        // what matters, not specifically the left one. Compare both
        // candidate lengths against the current `keep` bounds and keep the
        // longer side, then keep recursing into *that* side on the next
        // iteration (nested/compound artefacts can occur on either side).
        let left_len = abs_left_end.saturating_sub(keep.start);
        let right_len = keep.end.saturating_sub(abs_right_start);

        if right_len > left_len {
            keep.start = abs_right_start.max(keep.start).min(keep.end);
            detect = keep.start..keep.end;
        } else {
            keep.end = abs_left_end.min(keep.end).max(keep.start);
            // Re-detect only up to the inner edge of the mirrored region:
            // anything beyond it can look reverse-complement similar to the
            // arm and confuse the next pass.
            detect = detect.start..abs_arm_end.min(keep.end).max(detect.start);
        }

        // next iteration will detect/refine on the shortened view
    }

    Ok(keep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint::foldback_signature_from_local_matches;
    use crate::utils::revcomp_in_place;
    use crate::{
        cfg::{FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg, SharedCfg, SigCfg},
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
                    max_jump_clip: 1000,
                },
                fold_diag_tol: diag_tol,
            },
            FoldOnlyCfg { min_arm, ..Default::default() },
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

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch, coarse.span, coarse.arm_end, coarse.arm_start_right)
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

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch, coarse.span, coarse.arm_end, coarse.arm_start_right)
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
        let hifi = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch, coarse.span, coarse.arm_end, coarse.arm_start_right).unwrap();

        shared.refine.mode = RefineMode::ONT;
        let ont = refine_breakpoint(&s, coarse.split_pos, &shared, refine_scratch, coarse.span, coarse.arm_end, coarse.arm_start_right)
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

        let refined = refine_breakpoint(&s, coarse.split_pos, &shared, &mut scratch.refine, coarse.span, coarse.arm_end, coarse.arm_start_right)
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

    // ── identity_est robustness tests ─────────────────────────────────────────
    // These tests exist specifically because identity_est was systematically
    // wrong in earlier versions (comparisons extending far beyond the palindromic
    // arm).  They must remain: if identity_est regresses, callers will silently
    // filter out genuine artefacts.

    /// Perfect clean foldback [A][RC(A)] — identity must be ≥ 0.99.
    #[test]
    fn identity_est_is_high_for_clean_palindrome() {
        let arm = random_dna(3000, 0xCAFE_0001);
        let mut rc = arm.clone();
        revcomp_in_place(&mut rc);
        let seq: Vec<u8> = [arm.as_slice(), rc.as_slice()].concat();

        let (shared, fold) = test_fold_cfg(17, 21, 120, 200, 5, 200, 1200, RefineMode::HiFi, 0.25);
        let mut scratch = scratch::FoldScratch::new();

        let fb = detect_foldback(&seq, &shared, &fold, &mut scratch)
            .expect("clean palindrome should be detected");
        let rf = refine_breakpoint(&seq, fb.split_pos, &shared, &mut scratch.refine, fb.span, fb.arm_end, fb.arm_start_right)
            .unwrap()
            .expect("clean palindrome should refine");

        assert!(
            rf.identity_est >= 0.99,
            "clean palindrome: expected identity ≥ 0.99, got {:.4}  \
             (fb.span={}, split_pos={})",
            rf.identity_est, fb.span, rf.split_pos
        );
    }

    /// Foldback with a 300 bp template-switch gap: [A][gap][RC(A)].
    /// Gap correction must find δ ≈ 150 and identity must still be ≥ 0.90.
    #[test]
    fn identity_est_is_high_for_palindrome_with_gap() {
        let arm = random_dna(2000, 0xCAFE_0002);
        let gap = random_dna(300, 0xDEAD_BEEF); // unrelated gap sequence
        let mut rc = arm.clone();
        revcomp_in_place(&mut rc);
        let seq: Vec<u8> = [arm.as_slice(), gap.as_slice(), rc.as_slice()].concat();

        let (shared, fold) = test_fold_cfg(17, 21, 120, 200, 5, 200, 1200, RefineMode::HiFi, 0.25);
        let mut scratch = scratch::FoldScratch::new();

        let fb = detect_foldback(&seq, &shared, &fold, &mut scratch)
            .expect("gapped palindrome should be detected");
        let rf = refine_breakpoint(&seq, fb.split_pos, &shared, &mut scratch.refine, fb.span, fb.arm_end, fb.arm_start_right)
            .unwrap()
            .expect("gapped palindrome should refine");

        assert!(
            rf.identity_est >= 0.90,
            "gapped palindrome: expected identity ≥ 0.90, got {:.4}  \
             (fb.span={}, split_pos={}, read_len={})",
            rf.identity_est, fb.span, rf.split_pos, seq.len()
        );
    }

    /// Noisy foldback with ~5% substitutions: [A][noise_RC(A)].
    /// Simulates ONT-level error.  Identity must match the noise rate (≥ 0.90).
    #[test]
    fn identity_est_tracks_arm_error_rate() {
        let arm = random_dna(3000, 0xCAFE_0003);
        // sub_every=20 → ~5% substitution rate in the RC arm
        let (seq, _bp) = make_noisy_foldback_from_left(&arm, b"", 20, 0, 0);

        let (shared, fold) = test_fold_cfg(17, 21, 120, 200, 5, 200, 1200, RefineMode::HiFi, 0.25);
        let mut scratch = scratch::FoldScratch::new();

        let fb = detect_foldback(&seq, &shared, &fold, &mut scratch)
            .expect("noisy palindrome should be detected");
        let rf = refine_breakpoint(&seq, fb.split_pos, &shared, &mut scratch.refine, fb.span, fb.arm_end, fb.arm_start_right)
            .unwrap()
            .expect("noisy palindrome should refine");

        // With ~5% error rate, expect identity ≥ 0.90 (leaving margin for
        // the coarse span estimate covering slightly more than the palindromic arm).
        assert!(
            rf.identity_est >= 0.90,
            "noisy palindrome (~5% sub): expected identity ≥ 0.90, got {:.4}  \
             (fb.span={}, split_pos={})",
            rf.identity_est, fb.span, rf.split_pos
        );
    }

    // ── keep-longer-side cut policy ─────────────────────────────────────────

    fn cfg_for_cut_policy(sig_flank_bp: usize, sig_take: usize) -> MdaxCfg {
        MdaxCfg {
            shared: SharedCfg {
                minimizer: MinimizerCfg { k: 17, w: 21, forward_only: true },
                min_matches: 10,
                end_guard: 100,
                refine: RefineCfg {
                    window: 100,
                    arm: 900,
                    mode: RefineMode::HiFi,
                    max_ed_rate: 0.25,
                    max_jump_clip: 1000,
                },
                fold_diag_tol: 120,
            },
            fold: FoldOnlyCfg { min_arm: 200, ..Default::default() },
            fold2: FoldSecondPassCfg {
                min_support: 3,
                min_identity: 0.6,
                min_support_ident: 0.0,
                cut_low_ident: false,
            },
            sig: SigCfg { flank_bp: sig_flank_bp, take: sig_take, value_shift: 0 },
        }
    }

    /// An artefact foldback sitting near the *start* of a long read must keep
    /// the much longer right-hand remainder, not the tiny left arm -- this is
    /// the core behaviour change from "always keep left".
    #[test]
    fn artefact_near_read_start_keeps_longer_right_side() {
        let lead = random_dna(500, 0xA11CE);
        let arm = random_dna(1000, 0xB0B01);
        let mut rc_arm = arm.clone();
        revcomp_in_place(&mut rc_arm);
        let trail = random_dna(8000, 0xC0DE1);

        let mut seq = Vec::new();
        seq.extend_from_slice(&lead);
        seq.extend_from_slice(&arm);
        seq.extend_from_slice(&rc_arm);
        seq.extend_from_slice(&trail);

        let cfg = cfg_for_cut_policy(300, 8);
        let mut scratch = scratch::FoldScratch::new();

        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch)
            .expect("should detect foldback near read start");
        let mut refine_scratch = RefineScratch::default();
        let rf = refine_breakpoint(
            &seq, fb.split_pos, &cfg.shared, &mut refine_scratch, fb.span, fb.arm_end, fb.arm_start_right,
        )
        .unwrap()
        .expect("should refine");

        // Empty support map: nothing is "real", so this is treated as an artefact.
        let support = JunctionSupport::empty();
        let mut sig_scratch = SigScratch::default();

        let keep = recursive_foldback_cut_from_first_range(
            &seq, 0..seq.len(), fb, rf, &cfg, &support, 5, &mut scratch, &mut sig_scratch,
        )
        .unwrap();

        let kept_len = keep.end - keep.start;
        assert!(
            kept_len > seq.len() / 2,
            "expected the longer right-hand side (~9000bp) to be kept over \
             the short left arm (~1500bp), got kept_len={kept_len} of {} \
             (keep={:?})",
            seq.len(), keep,
        );
        assert_eq!(
            keep.end, seq.len(),
            "expected the kept range to extend to the end of the read (right-keep), got {keep:?}"
        );
    }

    /// Cut a displaced fold built from one genome stretch, return the kept range.
    /// Read = G[0..fwd_end] + RC(G[rev_low..rev_high]): the forward copy runs to
    /// `fwd_end`, the polymerase then switches strand at genome `rev_high`.
    fn cut_displaced(fwd_end: usize, rev_low: usize, rev_high: usize, tail: usize) -> (Vec<u8>, Range<usize>) {
        let g = random_dna(rev_high.max(fwd_end) + 10, 0x5EED_D15B);
        let mut rc = g[rev_low..rev_high].to_vec();
        revcomp_in_place(&mut rc);
        let mut seq = g[..fwd_end].to_vec();
        seq.extend_from_slice(&rc);
        seq.extend_from_slice(&random_dna(tail, 0x7A11));

        let cfg = cfg_for_cut_policy(300, 8);
        let mut scratch = scratch::FoldScratch::new();
        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch).expect("detect");
        let mut rs = RefineScratch::default();
        let rf = refine_breakpoint(&seq, fb.split_pos, &cfg.shared, &mut rs, fb.span, fb.arm_end, fb.arm_start_right)
            .unwrap()
            .expect("refine");
        let support = JunctionSupport::empty();
        let mut sig = SigScratch::default();
        let keep = recursive_foldback_cut_from_first_range(
            &seq, 0..seq.len(), fb, rf, &cfg, &support, 5, &mut scratch, &mut sig,
        )
        .unwrap();
        (seq, keep)
    }

    /// Spacer on the dropped (reverse) side: the forward copy ends at 3000,
    /// the reverse copy starts at genome 3600, so read 3000..3600 is
    /// reverse-strand sequence with no mirror. Keeping it (the old policy)
    /// left a 600 bp strand-switched stub on the kept read.
    #[test]
    fn displaced_fold_spacer_on_dropped_side_is_not_kept() {
        let (_seq, keep) = cut_displaced(3000, 1000, 3600, 0);
        assert_eq!(keep.start, 0, "left side is longer and should be kept: {keep:?}");
        assert!(
            keep.end <= 3000 + 40,
            "kept read must end at the strand switch (3000), not carry the reverse-strand spacer: {keep:?}"
        );
        assert!(keep.end >= 3000 - 200, "should not over-trim the forward copy: {keep:?}");
    }

    /// Spacer on the kept (forward) side: the forward copy runs on to 3000
    /// but the reverse copy starts at genome 2400. Read 2400..3000 is good
    /// forward sequence, but a single read cannot tell the two cases apart,
    /// so it is dropped too; what matters is that nothing reverse-strand is
    /// kept and the loss is bounded by the spacer.
    #[test]
    fn displaced_fold_spacer_on_kept_side_costs_at_most_the_spacer() {
        let (_seq, keep) = cut_displaced(3000, 1000, 2400, 0);
        assert_eq!(keep.start, 0, "{keep:?}");
        assert!(keep.end <= 3000 + 40, "no reverse-strand sequence may be kept: {keep:?}");
        assert!(keep.end + 650 >= 3000, "loss must be bounded by the 600 bp spacer: {keep:?}");
    }

    /// Same two geometries mirrored: the long side is the right one.
    #[test]
    fn displaced_fold_keeping_right_side_drops_the_spacer() {
        // forward copy short, reverse copy long, plus a long right tail
        let (seq, keep) = cut_displaced(1500, 0, 2100, 6000);
        assert_eq!(keep.end, seq.len(), "{keep:?}");
        // strand switch is at read 1500; spacer is read 1500..2100
        assert!(keep.start + 40 >= 2100, "kept right side must start after the spacer (2100): {keep:?}");
        assert!(keep.start <= 2100 + 200, "should not over-trim the reverse copy: {keep:?}");
    }

    /// The symmetric case: a foldback near the *end* of a long read must keep
    /// the longer left-hand side, matching the pre-existing "keep left" policy.
    #[test]
    fn artefact_near_read_end_keeps_longer_left_side() {
        let lead = random_dna(8000, 0xD00D1);
        let arm = random_dna(1000, 0xE22E2);
        let mut rc_arm = arm.clone();
        revcomp_in_place(&mut rc_arm);
        let trail = random_dna(500, 0xF44F3);

        let mut seq = Vec::new();
        seq.extend_from_slice(&lead);
        seq.extend_from_slice(&arm);
        seq.extend_from_slice(&rc_arm);
        seq.extend_from_slice(&trail);

        let cfg = cfg_for_cut_policy(300, 8);
        let mut scratch = scratch::FoldScratch::new();

        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch)
            .expect("should detect foldback near read end");
        let mut refine_scratch = RefineScratch::default();
        let rf = refine_breakpoint(
            &seq, fb.split_pos, &cfg.shared, &mut refine_scratch, fb.span, fb.arm_end, fb.arm_start_right,
        )
        .unwrap()
        .expect("should refine");

        let support = JunctionSupport::empty();
        let mut sig_scratch = SigScratch::default();

        let keep = recursive_foldback_cut_from_first_range(
            &seq, 0..seq.len(), fb, rf, &cfg, &support, 5, &mut scratch, &mut sig_scratch,
        )
        .unwrap();

        let kept_len = keep.end - keep.start;
        assert!(
            kept_len > seq.len() / 2,
            "expected the longer left-hand side (~9000bp) to be kept, got kept_len={kept_len} of {} (keep={:?})",
            seq.len(), keep,
        );
        assert_eq!(
            keep.start, 0,
            "expected the kept range to start at the beginning of the read (left-keep), got {keep:?}"
        );
    }

    /// A compound (nested) artefact must still be fully corrected regardless
    /// of which side of the *outer* split happens to be longer: cutting the
    /// outer junction must not terminate recursion, so any inner nested
    /// junction within the kept side is also removed.
    #[test]
    fn nested_artefact_within_kept_longer_side_is_also_cut() {
        // outer_lead is short, so the first (outer) cut keeps the right side --
        // which itself still contains a second, inner foldback that must also
        // be cleaned up rather than left in the output untouched.
        let outer_lead = random_dna(300, 0x1234);
        let outer_arm = random_dna(2000, 0x5678);
        let mut outer_rc = outer_arm.clone();
        revcomp_in_place(&mut outer_rc);

        let inner_arm = random_dna(1500, 0x9abc);
        let mut inner_rc = inner_arm.clone();
        revcomp_in_place(&mut inner_rc);

        let mut seq = Vec::new();
        seq.extend_from_slice(&outer_lead);
        seq.extend_from_slice(&outer_arm);
        seq.extend_from_slice(&outer_rc);
        seq.extend_from_slice(&inner_arm);
        seq.extend_from_slice(&inner_rc);

        let cfg = cfg_for_cut_policy(300, 8);
        let mut scratch = scratch::FoldScratch::new();

        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch)
            .expect("should detect at least the outer foldback");
        let mut refine_scratch = RefineScratch::default();
        let rf = refine_breakpoint(
            &seq, fb.split_pos, &cfg.shared, &mut refine_scratch, fb.span, fb.arm_end, fb.arm_start_right,
        )
        .unwrap()
        .expect("should refine");

        let support = JunctionSupport::empty();
        let mut sig_scratch = SigScratch::default();

        let keep = recursive_foldback_cut_from_first_range(
            &seq, 0..seq.len(), fb, rf, &cfg, &support, 5, &mut scratch, &mut sig_scratch,
        )
        .unwrap();

        let kept_len = keep.end - keep.start;
        // Whichever side got kept first, the surviving inner foldback must
        // also have been trimmed: nothing containing a full second arm+RC(arm)
        // pair (3000bp) should remain uncut.
        assert!(
            kept_len < seq.len() - 1000,
            "expected the nested artefact to also be cut, kept_len={kept_len} of {} (keep={:?})",
            seq.len(), keep,
        );
    }

    // ---------------------------------------------------------------------
    // Alignment-verified rescue path
    // ---------------------------------------------------------------------

    /// Production-like balanced HiFi gates (k=17, w=21, 12 matches, 1000 bp
    /// footprint, 500 bp end guard), with the rescue path on or off.
    fn cfg_rescue(enabled: bool) -> MdaxCfg {
        MdaxCfg {
            shared: SharedCfg {
                minimizer: MinimizerCfg { k: 17, w: 21, forward_only: true },
                min_matches: 12,
                end_guard: 500,
                refine: RefineCfg {
                    window: 200,
                    arm: 1200,
                    mode: RefineMode::HiFi,
                    max_ed_rate: 0.25,
                    max_jump_clip: 1000,
                },
                fold_diag_tol: 120,
            },
            fold: FoldOnlyCfg {
                min_arm: 1000,
                rescue: crate::cfg::RescueCfg {
                    enabled,
                    ..crate::cfg::RescueCfg::disabled()
                },
            },
            fold2: FoldSecondPassCfg {
                min_support: 3,
                min_identity: 0.6,
                min_support_ident: 0.0,
                cut_low_ident: true,
            },
            sig: SigCfg { flank_bp: 600, take: 8, value_shift: 0 },
        }
    }

    /// `body` followed by the reverse complement of its last `arm` bases,
    /// optionally after `spacer` unrelated bases: a foldback whose junction
    /// sits `arm` bp from the read end.
    fn tail_foldback(body_len: usize, arm: usize, spacer: usize, seed: u64) -> (Vec<u8>, usize) {
        let mut seq = random_dna(body_len, seed);
        let mut tail = seq[body_len - arm..].to_vec();
        revcomp_in_place(&mut tail);
        seq.extend(random_dna(spacer, seed ^ 0x5EED));
        seq.extend(tail);
        (seq, body_len)
    }

    #[test]
    fn rescue_finds_short_tail_foldback_the_gates_reject() {
        let (seq, junction) = tail_foldback(4000, 300, 0, 0xFACE01);

        let off = cfg_rescue(false);
        let mut scratch = scratch::FoldScratch::new();
        assert!(
            detect_foldback(&seq, &off.shared, &off.fold, &mut scratch).is_none(),
            "a 300 bp arm 300 bp from the read end is below every coarse gate"
        );

        let on = cfg_rescue(true);
        let fb = detect_foldback(&seq, &on.shared, &on.fold, &mut scratch)
            .expect("rescue path should find the tail foldback");
        let v = fb.verified.expect("rescued candidates carry their alignment");
        assert!(v.identity > 0.99, "clean copy should align perfectly, got {}", v.identity);
        assert!(v.arm >= on.fold.rescue.min_arm);
        assert!(
            (fb.split_pos as i64 - junction as i64).abs() <= 40,
            "split {} should sit at the junction {}",
            fb.split_pos,
            junction
        );

        // The support gate must be able to protect a genuine short repeat,
        // which needs junction keys from the accepted candidate's bin.
        let mut tmp = Vec::new();
        let mut keys = Vec::new();
        crate::fingerprint::junction_keys(&scratch.best_trip, 32, &mut tmp, &mut keys);
        assert!(
            keys.len() >= crate::fingerprint::JUNCTION_MIN_ELEMS,
            "rescued candidate should yield junction keys, got {}",
            keys.len()
        );
    }

    #[test]
    fn rescue_handles_spacer_between_arms() {
        let (seq, junction) = tail_foldback(4000, 400, 120, 0xFACE02);
        let on = cfg_rescue(true);
        let mut scratch = scratch::FoldScratch::new();
        let fb = detect_foldback(&seq, &on.shared, &on.fold, &mut scratch)
            .expect("rescue path should find a tail foldback with a spacer");
        let v = fb.verified.unwrap();
        assert!(v.identity > 0.99);
        assert!(
            (v.gap as i64 - 120).abs() <= 80,
            "gap {} should be close to the 120 bp spacer",
            v.gap
        );
        assert!(fb.arm_end <= junction + 20 && fb.arm_start_right + 20 >= junction + 120 - 60);
    }

    #[test]
    fn rescue_cut_drops_only_the_short_tail() {
        let (seq, junction) = tail_foldback(4000, 300, 0, 0xFACE03);
        let cfg = cfg_rescue(true);
        let mut scratch = scratch::FoldScratch::new();
        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch).unwrap();
        let rf = refine_detected(&seq, &fb, &cfg.shared, &mut scratch.refine)
            .unwrap()
            .expect("rescued candidates always refine");
        let support = JunctionSupport::empty();
        let mut sig_scratch = SigScratch::default();
        let keep = recursive_foldback_cut_from_first_range(
            &seq, 0..seq.len(), fb, rf, &cfg, &support, 5, &mut scratch, &mut sig_scratch,
        )
        .unwrap();
        assert_eq!(keep.start, 0);
        assert!(
            (keep.end as i64 - junction as i64).abs() <= 40,
            "expected the read trimmed at the junction ({junction}), got {keep:?}"
        );
    }

    /// A short tail foldback that only becomes a read-end foldback once an
    /// outer artefact has been cut away must be cut on the next recursion.
    #[test]
    fn rescued_foldback_nested_inside_a_cut_read_is_recut() {
        // inner: 3700 bp body + RC of its last 300 bp (junction at 3700)
        let (inner, inner_junction) = tail_foldback(3700, 300, 0, 0xFACE08);
        // outer: inner + RC of inner[1500..] (junction at 4000, long arm)
        let mut seq = inner.clone();
        let mut tail = inner[1500..].to_vec();
        revcomp_in_place(&mut tail);
        seq.extend(tail);

        let cfg = cfg_rescue(true);
        let mut scratch = scratch::FoldScratch::new();
        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch).expect("outer fold");
        assert!(fb.verified.is_none(), "the outer fold is a gated detection");
        let rf = refine_detected(&seq, &fb, &cfg.shared, &mut scratch.refine)
            .unwrap()
            .expect("refine");
        let support = JunctionSupport::empty();
        let mut sig_scratch = SigScratch::default();
        let keep = recursive_foldback_cut_from_first_range(
            &seq, 0..seq.len(), fb, rf, &cfg, &support, 5, &mut scratch, &mut sig_scratch,
        )
        .unwrap();
        assert_eq!(keep.start, 0);
        assert!(
            (keep.end as i64 - inner_junction as i64).abs() <= 40,
            "expected both folds cut, ending at the inner junction {inner_junction}, got {keep:?}"
        );
    }

    #[test]
    fn rescue_leaves_gated_detections_untouched() {
        let (seq, _) = tail_foldback(6000, 2500, 0, 0xFACE04);
        let mut scratch = scratch::FoldScratch::new();
        let off = cfg_rescue(false);
        let a = detect_foldback(&seq, &off.shared, &off.fold, &mut scratch).expect("gated hit");
        let on = cfg_rescue(true);
        let b = detect_foldback(&seq, &on.shared, &on.fold, &mut scratch).expect("gated hit");
        assert!(b.verified.is_none(), "a candidate the gates accept must not be rescued");
        assert_eq!((a.split_pos, a.matches, a.span), (b.split_pos, b.matches, b.span));
    }

    #[test]
    fn rescue_rejects_reads_without_a_fold() {
        let on = cfg_rescue(true);
        let mut scratch = scratch::FoldScratch::new();

        let plain = random_dna(8000, 0xFACE05);
        assert!(detect_foldback(&plain, &on.shared, &on.fold, &mut scratch).is_none());

        // A dinucleotide repeat is its own reverse complement about every
        // position: plenty of anti-diagonal matches, but one or two
        // minimizer values repeated, which is not evidence of a fold.
        let mut micro = random_dna(3000, 0xFACE06);
        micro.extend(repeat_bytes(b"AT", 300));
        micro.extend(random_dna(3000, 0xFACE07));
        assert!(
            detect_foldback(&micro, &on.shared, &on.fold, &mut scratch).is_none(),
            "a microsatellite must not be rescued as a foldback"
        );

        // Stale-state guard: a read too short to be matched at all must not
        // inherit the previous read's bins.
        let (fold, _) = tail_foldback(4000, 300, 0, 0xFACE08);
        assert!(detect_foldback(&fold, &on.shared, &on.fold, &mut scratch).is_some());
        assert!(detect_foldback(&plain[..30], &on.shared, &on.fold, &mut scratch).is_none());
    }

    #[test]
    fn rescue_requires_the_arms_to_actually_align() {
        // Same seeds would be proposed either way; what differs is whether
        // the bases between them agree.
        let on = cfg_rescue(true);
        let mut scratch = scratch::FoldScratch::new();

        let body = random_dna(4000, 0xFACE09);
        let arm = &body[3600..];
        let mut light = mutate_substitutions(arm, 25); // 4% divergence
        revcomp_in_place(&mut light);
        let mut seq = body.clone();
        seq.extend(light);
        assert!(
            detect_foldback(&seq, &on.shared, &on.fold, &mut scratch).is_some(),
            "a 96%-identical 400 bp arm should verify at the 90% threshold"
        );

        // Arm shorter than the minimum verified length: seeds exist, but
        // the alignment window runs into unrelated sequence.
        let mut seq = body.clone();
        let mut stub = body[3920..].to_vec(); // 80 bp
        revcomp_in_place(&mut stub);
        seq.extend(stub);
        assert!(
            detect_foldback(&seq, &on.shared, &on.fold, &mut scratch).is_none(),
            "an 80 bp arm is below rescue.min_arm (150) and must not be called"
        );
    }

    /// Terminal dense scan only: the minimizer-seeded rescue is switched off
    /// by an unreachable `min_pairs`, so these tests exercise one path.
    fn cfg_terminal_only(mode: RefineMode) -> MdaxCfg {
        let mut cfg = cfg_rescue(true);
        cfg.shared.refine.mode = mode;
        cfg.fold.rescue.min_pairs = 1_000_000;
        cfg.fold.rescue.terminal = true;
        if mode == RefineMode::ONT {
            cfg.shared.minimizer = MinimizerCfg { k: 11, w: 11, forward_only: true };
            cfg.fold.rescue.term_k = 9;
            cfg.fold.rescue.term_min_seeds = 5;
            cfg.fold.rescue.min_identity = 0.75;
            cfg.fold.rescue.min_arm = 200;
        }
        cfg
    }

    #[test]
    fn terminal_scan_finds_noisy_tail_foldback() {
        // ~9% divergence between the arms, substitutions and indels both:
        // the regime where sampled 11-mers rarely survive in both copies.
        let body = random_dna(5000, 0x7E2A01);
        let arm = &body[4650..];
        let noisy = inject_deletions(&inject_insertions(&mutate_substitutions(arm, 23), 41, b'G'), 37);
        let mut tail = noisy;
        revcomp_in_place(&mut tail);
        let mut seq = body.clone();
        seq.extend(tail);

        let cfg = cfg_terminal_only(RefineMode::ONT);
        let mut scratch = scratch::FoldScratch::new();
        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch)
            .expect("terminal scan should find a noisy 350 bp tail foldback");
        let v = fb.verified.expect("terminal candidates are verified");
        assert!(v.identity >= 0.75 && v.identity < 0.99, "identity {}", v.identity);
        assert!(
            (fb.split_pos as i64 - 5000).abs() <= 60,
            "split {} should sit at the junction (5000)",
            fb.split_pos
        );

        // and the same read from the other strand is a *head* foldback
        let mut rev = seq.clone();
        revcomp_in_place(&mut rev);
        let fb = detect_foldback(&rev, &cfg.shared, &cfg.fold, &mut scratch)
            .expect("terminal scan should find the same foldback at the read start");
        let junction = rev.len() - 5000;
        assert!(
            (fb.split_pos as i64 - junction as i64).abs() <= 60,
            "split {} should sit at the junction ({junction})",
            fb.split_pos
        );
    }

    #[test]
    fn terminal_scan_requires_the_copy_to_reach_the_read_end() {
        let cfg = cfg_terminal_only(RefineMode::HiFi);
        let mut scratch = scratch::FoldScratch::new();

        let (mut seq, _) = tail_foldback(4000, 300, 0, 0x7E2A02);
        assert!(
            detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch).is_some(),
            "a clean 300 bp tail foldback reaches the read end and must be found"
        );

        // The same inverted copy followed by 400 bp of unrelated sequence is
        // a small hairpin that happens to sit near the end, not a strand
        // switch that ran off the fragment.
        seq.extend(random_dna(400, 0x7E2A03));
        assert!(
            detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch).is_none(),
            "an inverted copy that stops short of the read end must not be called by the terminal scan"
        );

        let plain = random_dna(6000, 0x7E2A04);
        assert!(detect_foldback(&plain, &cfg.shared, &cfg.fold, &mut scratch).is_none());

        let mut micro = random_dna(5000, 0x7E2A05);
        micro.extend(repeat_bytes(b"AT", 250));
        assert!(
            detect_foldback(&micro, &cfg.shared, &cfg.fold, &mut scratch).is_none(),
            "a read ending in a dinucleotide repeat is not a tail foldback"
        );
    }

    #[test]
    fn rescue_end_only_ignores_short_hairpins_inside_a_read() {
        // A 300 bp inverted copy followed by 2 kb of unrelated sequence: a
        // short hairpin in the middle of a read, which on real data is as
        // often a genome-encoded repeat as an artefact.
        let (mut seq, junction) = tail_foldback(4000, 300, 0, 0xE0D001);
        seq.extend(random_dna(2000, 0xE0D002));

        let mut scratch = scratch::FoldScratch::new();
        let end_only = cfg_rescue(true);
        assert!(end_only.fold.rescue.end_only, "end_only is the default");
        assert!(
            detect_foldback(&seq, &end_only.shared, &end_only.fold, &mut scratch).is_none(),
            "a mid-read short hairpin must not be rescued when end_only is set"
        );

        let mut anywhere = cfg_rescue(true);
        anywhere.fold.rescue.end_only = false;
        let fb = detect_foldback(&seq, &anywhere.shared, &anywhere.fold, &mut scratch)
            .expect("with end_only off the same hairpin is rescued");
        assert!((fb.split_pos as i64 - junction as i64).abs() <= 40);
    }
}
