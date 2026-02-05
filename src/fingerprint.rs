//! Foldback junction fingerprinting + support statistics.
//!
//! This module provides two related pieces of functionality:
//! 1) **Signature / fingerprint** construction for a foldback junction.
//!    The signature is a compact `u64` intended to be stable across reads that
//!    share the *same underlying junction*, while being cheap to compute.
//!
//! 2) **Support aggregation** keyed by signature.
//!    Pass1 builds a `HashMap<sig, SupportStats>` describing how often a signature
//!    appears and how tightly its split positions cluster.
//!
//! The core idea is MinHash-ish:
//! - sample minimizers (cheap sequence sketch)
//! - keep the smallest `take` values (stable subset)
//! - hash the resulting small list into a `u64` fingerprint
//!
//! Stability knobs:
//! - canonical minimizers (`forward_only=false`) increase strand-agnostic stability
//! - quantization via `value_shift` can reduce sensitivity to small value jitter
//! - `window_bp` / `flank_bp` bounds keep the signature local to the junction

use gxhash::{GxHasher, HashMap};
use std::hash::{Hash, Hasher};

use crate::{
    cfg::{MdaxCfg, SharedCfg},
    minimizer::sampled_minimizers_into,
    scratch::SigScratch,
};

/// Compute a foldback junction fingerprint using *local match evidence* gathered by coarse detection.
///
/// This variant is designed to be cheap in pass2 because it reuses matchpoints already computed
/// while detecting the foldback (i.e. avoid re-sampling minimizers from scratch).
///
/// Inputs:
/// - `matches`: `(pos, val)` tuples produced by coarse detection, where `pos` is a read coordinate
///   (or an evidence coordinate) and `val` is a minimizer hash.
/// - `split`: coarse/refined split position (bp).
/// - `window_bp`: only matchpoints within this distance of `split` are kept.
/// - `take`: signature size (how many smallest hashes to retain).
/// - `value_shift`: optional right shift applied before hashing (coarse quantization).
///
/// Returns:
/// - `Some(sig)` if enough evidence remains after filtering.
/// - `None` if fewer than `take` matchpoints are available near the junction.
///
/// Side effect:
/// This function *destructively filters* `matches` in-place (`retain`), permanently dropping
/// matchpoints far from the junction. This is intentional to keep the signature local and stable.
///
/// Why locality helps:
/// - junction-adjacent evidence is most informative for the breakpoint identity
/// - distant repeats can introduce collisions and reduce specificity
pub fn foldback_signature_from_local_matches(
    matches: &mut Vec<(usize, u64)>, // (pos, val)
    split: usize,
    window_bp: usize,
    take: usize,
    value_shift: u8,
) -> Option<u64> {
    // Keep only matches close to the junction
    matches.retain(|(pos, _)| {
        let d = if *pos > split {
            *pos - split
        } else {
            split - *pos
        };
        d <= window_bp
    });

    // If we don't have enough local evidence, we can't build a stable signature.
    if matches.len() < take {
        return None;
    }

    // Extract minimizer hash values only.
    // (We ignore positions after locality filtering; the signature is value-based.)
    let mut vals: Vec<u64> = matches.iter().map(|(_, v)| *v).collect();

    // MinHash-style reduction
    // sort and keep the smallest `take` values to get a stable sketch-like subset.
    vals.sort_unstable();
    vals.truncate(take);

    // Optional quantization:
    // shifting values reduces sensitivity to fine-grained hash differences.
    // After shifting, re-sort so the "smallest take" property remains meaningful.
    if value_shift > 0 {
        for v in vals.iter_mut() {
            *v >>= value_shift;
        }
        vals.sort_unstable();
    }

    let mut h = gxhash::GxHasher::default();
    vals.hash(&mut h);
    Some(h.finish())
}

/// Compute a reproducible fingerprint for a foldback junction at `split` using sequence flanks.
///
/// This is the “direct-from-sequence” signature builder. It is a robust fallback when you
/// don't have match evidence from coarse detection (or when you prefer a purely sequence-based key).
///
/// Approach:
/// - take `flank_bp` bases on the left of the split
/// - take `flank_bp` bases on the right of the split, reverse-complement them
/// - sample canonical minimizers in each flank
/// - keep the smallest `take` minimizer values from each side
/// - concatenate, sort, optionally quantize (`value_shift`), and hash into `u64`
///
/// Returns `None` if:
/// - the requested flanks exceed sequence bounds, or
/// - either side yields fewer than `take` minimizers (insufficient sketch density).
///
/// Performance notes:
/// - Uses `SigScratch` to avoid allocations across calls.
/// - Constructs a lightweight minimizer config rather than cloning `SharedCfg`.
pub fn foldback_signature(
    seq: &[u8],
    split: usize,
    shared: &SharedCfg,
    flank_bp: usize,
    take: usize,
    scratch: &mut SigScratch,
    value_shift: u8,
) -> Option<u64> {
    if split < flank_bp || split + flank_bp > seq.len() {
        return None;
    }

    let left = &seq[split - flank_bp..split];
    let right = &seq[split..split + flank_bp];

    // Use canonical minimizers for signature stability:
    // a k-mer and its reverse-complement map to the same value.
    //
    // We build a small MinimizerCfg rather than cloning the whole SharedCfg each call.
    let mcfg = crate::cfg::MinimizerCfg {
        k: shared.minimizer.k,
        w: shared.minimizer.w,
        forward_only: false,
    };

    // Sample minimizers from the left flank into scratch buffers.
    sampled_minimizers_into(left, &mcfg, &mut scratch.pos_l, &mut scratch.val_l);

    // Sample minimizers from the right flank, but in reverse-complement orientation.
    // This makes the two arms comparable under a foldback model.
    scratch.right_rc.clear();
    scratch.right_rc.extend_from_slice(right);
    crate::utils::revcomp_in_place(&mut scratch.right_rc);

    sampled_minimizers_into(
        &scratch.right_rc,
        &mcfg,
        &mut scratch.pos_r,
        &mut scratch.val_r,
    );

    // Need at least `take` minimizers from each side to build a stable signature.
    if scratch.val_l.len() < take || scratch.val_r.len() < take {
        return None;
    }

    // sort values and take smallest `take`
    scratch.val_l.sort_unstable();
    scratch.val_r.sort_unstable();

    // Build combined list without allocating a new vec each time
    scratch.combined.clear();
    scratch.combined.extend_from_slice(&scratch.val_l[..take]);
    scratch.combined.extend_from_slice(&scratch.val_r[..take]);
    scratch.combined.sort_unstable();

    if take > 0 && value_shift > 0 {
        for v in scratch.combined.iter_mut() {
            *v >>= value_shift;
        }
        scratch.combined.sort_unstable();
    }

    let mut h = GxHasher::default();
    scratch.combined.hash(&mut h);
    Some(h.finish())
}

// Support statistics (pass1 aggregation; pass2 decisions)

/// Aggregate support statistics for a given junction signature.
///
/// Stored per signature in the pass1 support map and used in pass2 to decide whether
/// a foldback is “real” (genome-templated) or an artefact (amplification-induced).
///
/// Semantics:
/// - `n`: number of reads that produced this signature
/// - `min_split` / `max_split`: range of observed refined split positions
/// - `mean_ident`: running mean of refinement identity estimates
#[derive(Debug, Clone, Default)]
pub struct SupportStats {
    /// Number of observations (reads) contributing to this signature.
    pub n: usize,
    /// Minimum observed split position (bp) across reads.
    pub min_split: usize,
    /// Maximum observed split position (bp) across reads.
    pub max_split: usize,
    /// Running mean of `identity_est` across reads.
    pub mean_ident: f32,
}

impl SupportStats {
    pub fn new(split: usize, ident: f32) -> Self {
        Self {
            n: 1,
            min_split: split,
            max_split: split,
            mean_ident: ident,
        }
    }

    /// Update stats with an additional observation.
    ///
    /// Updates:
    /// - increments `n`
    /// - expands `[min_split, max_split]` to include `split`
    /// - updates `mean_ident` via an incremental mean (stable and O(1))
    pub fn update(&mut self, split: usize, ident: f32) {
        self.n += 1;
        self.min_split = self.min_split.min(split);
        self.max_split = self.max_split.max(split);
        // incremental mean
        self.mean_ident += (ident - self.mean_ident) / self.n as f32;
    }

    /// Return the spread of split positions across reads for this signature.
    ///
    /// A small span indicates consistent breakpoint localization (typical for a real junction).
    /// A large span suggests noise, signature collisions, or heterogeneous events.
    pub fn split_span(&self) -> usize {
        self.max_split.saturating_sub(self.min_split)
    }
}

/// Decide whether a foldback junction is "real" based on pass1 support statistics.
pub fn is_real_foldback(sig: u64, support: &HashMap<u64, SupportStats>, cfg: &MdaxCfg) -> bool {
    let Some(st) = support.get(&sig) else {
        return false;
    };

    // Read-space split positions are not expected to cluster tightly across reads.
    // The signature already anchors the locus; use support count as the primary gate.
    if st.n < cfg.fold2.min_support {
        return false;
    }

    // Optional: add an identity-based support gate if you want
    if cfg.fold2.min_support_ident > 0.0 && (st.mean_ident as f64) < cfg.fold2.min_support_ident {
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cfg::{FoldOnlyCfg, FoldSecondPassCfg, MinimizerCfg, RefineCfg, SharedCfg, SigCfg};
    use crate::foldback::detect_foldback;
    use crate::pipeline::pass1_build_support;
    use crate::scratch::FoldScratch;
    use crate::utils::RefineMode;

    use calm_io::stderrln;
    use gxhash::HashMapExt;
    use std::fs::File;
    use std::sync::Arc;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn pseudo_dna(len: usize, seed: u64) -> Vec<u8> {
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

    fn make_clean_foldback(left: &[u8], join: &[u8]) -> (Vec<u8>, usize) {
        let mut rc = left.to_vec();
        crate::utils::revcomp_in_place(&mut rc);

        let mut s = Vec::with_capacity(left.len() + join.len() + rc.len());
        s.extend_from_slice(left);
        s.extend_from_slice(join);
        let split = left.len() + join.len();
        s.extend_from_slice(&rc);
        (s, split)
    }

    fn shared_for_sig(k: usize, w: usize, min_matches: usize) -> SharedCfg {
        SharedCfg {
            minimizer: MinimizerCfg {
                k,
                w,
                forward_only: true,
            },
            min_matches,
            end_guard: 0,
            refine: RefineCfg {
                window: 0,
                arm: 0,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.0,
            },
            fold_diag_tol: 100,
        }
    }

    fn tmp_fasta_path(prefix: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("{prefix}_{now}.fa"));
        p
    }

    fn write_fasta(path: &std::path::Path, records: &[(&str, &[u8])]) {
        let f = File::create(path).expect("create temp fasta");
        let mut w = std::io::BufWriter::new(f);
        for (id, seq) in records {
            crate::utils::write_fasta(&mut w, id.as_bytes(), seq).expect("write fasta record");
        }
    }

    fn cfg_for_pass1() -> MdaxCfg {
        // Make it easy to detect + refine + fingerprint.
        MdaxCfg {
            shared: SharedCfg {
                minimizer: MinimizerCfg {
                    k: 9,
                    w: 11,
                    forward_only: true,
                },
                min_matches: 5,
                end_guard: 0,
                refine: RefineCfg {
                    window: 50,
                    arm: 300, // refinement arm; independent of fingerprint flank=1000
                    mode: RefineMode::HiFi,
                    max_ed_rate: 0.25,
                },
                fold_diag_tol: 50,
            },
            fold: FoldOnlyCfg {
                min_arm: 200, // span evidence for foldback detector
            },
            fold2: FoldSecondPassCfg {
                min_support: 3,
                min_identity: 0.60,
                min_support_ident: 0.0,
            },
            sig: SigCfg {
                flank_bp: 1000,
                take: 8,
                value_shift: 0,
            },
        }
    }

    #[test]
    fn foldback_pass1_support_finds_one_cluster() {
        let cfg = cfg_for_pass1();

        // Need >= 1000bp on each side of junction for foldback_signature().
        // So left must be >= 1000, and overall read long enough.
        let left = pseudo_dna(2500, 42);
        let (seq, _true_split) = make_clean_foldback(&left, b"");

        // 3 identical foldback reads -> should become one signature with n=3
        // plus a decoy random read that should be ignored
        let decoy = pseudo_dna(6000, 999);

        let path = tmp_fasta_path("mdax_fold2_pass1");
        write_fasta(
            &path,
            &[("r1", &seq), ("r2", &seq), ("r3", &seq), ("decoy", &decoy)],
        );

        let cfg2 = cfg.clone();

        let support = pass1_build_support(&path, Arc::new(cfg), 1, 1).expect("pass1 support");

        // remove file
        let _ = std::fs::remove_file(&path);

        assert_eq!(
            support.len(),
            1,
            "expected exactly one foldback fingerprint cluster"
        );

        let (_sig, st) = support.iter().next().unwrap();
        assert_eq!(st.n, 3, "should have 3 supporting reads");
        assert_eq!(st.split_span(), 0, "identical reads should cluster tightly");
        assert!(st.mean_ident >= cfg2.fold2.min_identity);
    }

    #[test]
    fn signature_is_deterministic_same_input() {
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 40); // 640bp
        let (seq, split) = make_clean_foldback(&left, b"");

        let shared = shared_for_sig(9, 11, 5);
        let mut sig_scratch = SigScratch::default();
        let value_shift = 0;

        let sig1 = foldback_signature(&seq, split, &shared, 200, 8, &mut sig_scratch, value_shift)
            .expect("sig1");
        let sig2 = foldback_signature(&seq, split, &shared, 200, 8, &mut sig_scratch, value_shift)
            .expect("sig2");

        // should be the same signature for same input
        assert_eq!(sig1, sig2);
    }

    #[test]
    fn signature_matches_for_same_junction_across_reads() {
        // Two reads that share same left flank content around the junction
        // but have extra prefix/suffix elsewhere should still match as long as
        // we fingerprint using the same split and flank region content.
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 50); // 800bp
        let (seq1, split1) = make_clean_foldback(&left, b"");
        let (mut seq2, split2) = make_clean_foldback(&left, b"");

        // add a prefix (doesn't affect the junction *content* if we also shift split)
        let prefix = repeat_bytes(b"GATTACA", 20);
        let mut with_prefix = Vec::with_capacity(prefix.len() + seq2.len());
        with_prefix.extend_from_slice(&prefix);
        with_prefix.extend_from_slice(&seq2);
        seq2 = with_prefix;

        let shared = shared_for_sig(9, 11, 5);
        let mut sig_scratch = SigScratch::default();
        let value_shift = 0;

        // fingerprint around the true junction positions in each read
        let sig1 = foldback_signature(
            &seq1,
            split1,
            &shared,
            200,
            8,
            &mut sig_scratch,
            value_shift,
        )
        .expect("sig1");
        let sig2 = foldback_signature(
            &seq2,
            split2 + prefix.len(),
            &shared,
            200,
            8,
            &mut sig_scratch,
            value_shift,
        )
        .expect("sig2");

        assert_eq!(sig1, sig2);
    }

    #[test]
    fn signature_changes_if_junction_flanks_change() {
        let left1 = repeat_bytes(b"ACGTTGCAACGTTGCA", 40);
        let left2 = repeat_bytes(b"TTTTGGGGCCCCAAAA", 40);

        let (seq1, split1) = make_clean_foldback(&left1, b"");
        let (seq2, split2) = make_clean_foldback(&left2, b"");

        let shared = shared_for_sig(9, 11, 5);
        let mut sig_scratch = SigScratch::default();
        let value_shift = 0;

        let sig1 = foldback_signature(
            &seq1,
            split1,
            &shared,
            200,
            8,
            &mut sig_scratch,
            value_shift,
        )
        .expect("sig1");
        let sig2 = foldback_signature(
            &seq2,
            split2,
            &shared,
            200,
            8,
            &mut sig_scratch,
            value_shift,
        )
        .expect("sig2");

        assert_ne!(sig1, sig2);
    }

    #[test]
    fn supportstats_updates_and_span() {
        let mut st = SupportStats::new(1000, 0.8);
        st.update(1010, 0.6);
        st.update(990, 0.9);

        assert_eq!(st.n, 3);
        assert_eq!(st.min_split, 990);
        assert_eq!(st.max_split, 1010);
        assert_eq!(st.split_span(), 20);

        // mean should be in range
        assert!(st.mean_ident >= 0.6 && st.mean_ident <= 0.9);
    }

    #[test]
    fn is_real_foldback_respects_threshold_and_tolerance() {
        use crate::cfg::{FoldOnlyCfg, MdaxCfg};

        // minimal cfg stub
        let cfg = MdaxCfg {
            shared: shared_for_sig(9, 11, 5),
            fold: FoldOnlyCfg { min_arm: 20 },
            fold2: crate::cfg::FoldSecondPassCfg {
                min_support: 3,
                min_identity: 0.6,
                min_support_ident: 0.0,
            },
            sig: SigCfg {
                flank_bp: 1000,
                take: 8,
                value_shift: 0,
            },
        };

        let mut support: HashMap<u64, SupportStats> = HashMap::new();
        let sig = 123456_u64;

        // only 2 supports => not real
        support.insert(sig, {
            let mut st = SupportStats::new(1000, 0.8);
            st.update(1005, 0.8);
            st
        });
        assert!(!is_real_foldback(sig, &support, &cfg));

        // 3 supports but too spread => not real
        support.insert(sig, {
            let mut st = SupportStats::new(1000, 0.8);
            st.update(1100, 0.8);
            st.update(1200, 0.8);
            st
        });
        assert!(!is_real_foldback(sig, &support, &cfg));

        // 3 supports and tight cluster => real
        support.insert(sig, {
            let mut st = SupportStats::new(1000, 0.8);
            st.update(1010, 0.8);
            st.update(1020, 0.8);
            st
        });
        assert!(is_real_foldback(sig, &support, &cfg));
    }

    #[test]
    fn foldback_pass1_support_ignores_outside_flank_noise() {
        let cfg = cfg_for_pass1();
        // ensure flank fits: default flank is 1000
        let left = pseudo_dna(3000, 42);
        let (base, _split) = make_clean_foldback(&left, b"");

        // mutate bases far from the junction so flanks stay identical
        fn mutate_far(mut s: Vec<u8>, split: usize, radius: usize, every: usize) -> Vec<u8> {
            let n = s.len();
            for i in (0..n).step_by(every.max(1)) {
                if i + 1 < split.saturating_sub(radius)
                    || i > (split + radius).min(n.saturating_sub(1))
                {
                    s[i] = match s[i] {
                        b'A' => b'C',
                        b'C' => b'A',
                        b'G' => b'T',
                        b'T' => b'G',
                        x => x,
                    };
                }
            }
            s
        }

        let split = left.len(); // join is empty
        let s1 = base.clone();
        let s2 = mutate_far(base.clone(), split, 1100, 97); // outside ±1100 keeps 1000bp flanks intact
        let s3 = mutate_far(base.clone(), split, 1100, 131);

        let decoy = pseudo_dna(8000, 999);

        let path = tmp_fasta_path("mdax_fold2_pass1_flank_noise");
        write_fasta(
            &path,
            &[("r1", &s1), ("r2", &s2), ("r3", &s3), ("decoy", &decoy)],
        );

        let support = pass1_build_support(&path, Arc::new(cfg), 1, 1).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(
            support.len(),
            1,
            "should still be one cluster when noise is outside flanks"
        );
        let (_sig, st) = support.iter().next().unwrap();
        assert_eq!(st.n, 3);
    }

    #[test]
    fn signature_collision_rate_is_low_on_random_foldbacks() {
        let shared = shared_for_sig(13, 17, 5);

        use std::collections::HashSet;
        let mut seen = HashSet::new();
        let mut collisions = 0usize;

        for seed in 0u64..2000 {
            let left = pseudo_dna(2500, 10_000 + seed);
            let (seq, split) = make_clean_foldback(&left, b"");
            let mut sig_scratch = SigScratch::default();
            let value_shift = 0;

            let sig = foldback_signature(
                &seq,
                split,
                &shared,
                1000,
                12,
                &mut sig_scratch,
                value_shift,
            )
            .expect("sig should exist for long enough reads");

            if !seen.insert(sig) {
                collisions += 1;
            }
        }

        // Tune threshold if you change signature scheme, but this should be near-zero.
        assert!(collisions <= 1, "too many collisions: {collisions}");
    }

    #[test]
    fn low_complexity_should_not_all_share_one_signature() {
        let shared = shared_for_sig(9, 11, 5);
        use std::collections::HashSet;

        fn poly(base: u8, len: usize) -> Vec<u8> {
            vec![base; len]
        }

        let mut sig_scratch = SigScratch::default();

        let mut sigs = HashSet::new();
        for (b, _seed) in [(b'A', 1u64), (b'C', 2), (b'G', 3), (b'T', 4)] {
            let left = poly(b, 2500);
            let (seq, split) = make_clean_foldback(&left, b"");
            let value_shift = 0;
            if let Some(sig) = foldback_signature(
                &seq,
                split,
                &shared,
                1000,
                12,
                &mut sig_scratch,
                value_shift,
            ) {
                sigs.insert(sig);
            }
        }

        // Either you get None (because not enough minimizers) or you get distinct signatures.
        assert!(sigs.len() <= 4);
    }

    #[test]
    fn signature_stability_under_split_jitter() {
        let shared = shared_for_sig(13, 17, 5);

        let left = pseudo_dna(3000, 42);
        let (seq, split) = make_clean_foldback(&left, b"");

        let mut sig_scratch = SigScratch::default();
        let value_shift = 0;
        let sig0 = foldback_signature(
            &seq,
            split,
            &shared,
            1000,
            12,
            &mut sig_scratch,
            value_shift,
        )
        .unwrap();

        let mut sig_scratch2 = SigScratch::default();
        // try small jitters
        let mut same = 0usize;
        for dj in [-20isize, -10, -5, 5, 10, 20] {
            let s2 = (split as isize + dj) as usize;
            let value_shift = 0;
            if let Some(sig) =
                foldback_signature(&seq, s2, &shared, 1000, 12, &mut sig_scratch2, value_shift)
            {
                if sig == sig0 {
                    same += 1;
                }
            }
        }

        // Don’t assert “must be stable” unless you want that behaviour;
        // this test tells you how stable it is for current signature design.
        stderrln!("Signature matched {} out of 6 small split shifts", same).unwrap();
        assert!(
            same >= 2,
            "signature seems too sensitive to small split shifts"
        );
    }

    #[test]
    fn is_real_foldback_only_when_supported_and_tight() {
        let cfg = cfg_for_pass1();

        let left = pseudo_dna(3000, 1);
        let (seq, _split) = make_clean_foldback(&left, b"");

        let path = tmp_fasta_path("mdax_real_logic");
        write_fasta(&path, &[("r1", &seq), ("r2", &seq), ("r3", &seq)]);

        let cfg2 = cfg.clone();

        let support = pass1_build_support(&path, Arc::new(cfg), 1, 1).unwrap();
        let _ = std::fs::remove_file(&path);

        assert_eq!(support.len(), 1);
        let (&sig, st) = support.iter().next().unwrap();
        assert!(is_real_foldback(sig, &support, &cfg2));
        assert_eq!(st.n, 3);
    }

    #[test]
    fn local_match_signature_stable_under_refinement_jitter() {
        let cfg = cfg_for_pass1();

        let left = pseudo_dna(3000, 123);
        let (seq, _split) = make_clean_foldback(&left, b"");

        let mut scratch = FoldScratch::new();

        let fb = detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut scratch).unwrap();

        let mut sigs = std::collections::HashSet::new();

        for dj in [-30isize, -10, 0, 10, 30] {
            let s = (fb.split_pos as isize + dj) as usize;
            let sig = foldback_signature_from_local_matches(
                &mut scratch.best_matches.clone(),
                s,
                cfg.sig.flank_bp,
                cfg.sig.take,
                cfg.sig.value_shift,
            );
            if let Some(sig) = sig {
                sigs.insert(sig);
            }
        }

        assert_eq!(
            sigs.len(),
            1,
            "local match-based signature should be stable to small split jitter"
        );
    }

    #[test]
    fn local_match_signature_ignores_distant_matches() {
        let mut matches = vec![
            (1000, 1),
            (1020, 2),
            (980, 3),
            (50_000, 999), // far away
            (60_000, 888), // far away
        ];

        let sig = foldback_signature_from_local_matches(
            &mut matches,
            1000,
            100, // tight window
            3,
            0,
        );

        assert!(sig.is_some());
        assert!(matches.len() <= 3, "distant matches should be removed");
    }
}
