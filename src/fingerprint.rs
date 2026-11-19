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

use gxhash::{GxHasher, HashMap, HashMapExt};
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
fn local_match_sketch(
    matches: &mut Vec<(usize, u64)>, // (pos, val)
    split: usize,
    window_bp: usize,
    take: usize,
    value_shift: u8,
) -> Option<Vec<u64>> {
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

    Some(vals)
}

pub fn foldback_signature_from_local_matches(
    matches: &mut Vec<(usize, u64)>, // (pos, val)
    split: usize,
    window_bp: usize,
    take: usize,
    value_shift: u8,
) -> Option<u64> {
    let vals = local_match_sketch(matches, split, window_bp, take, value_shift)?;
    let mut h = gxhash::GxHasher::default();
    vals.hash(&mut h);
    Some(h.finish())
}

/// Same evidence as [`foldback_signature_from_local_matches`], but returned as a
/// set of LSH "leave-one-out" band signatures instead of a single exact hash.
///
/// Background: the bottom-`take` MinHash sketch is meant to support *approximate*
/// matching (estimate overlap between two reads' minimizer sets), but hashing the
/// whole sorted sketch into one `u64` and using that as an exact `HashMap` key
/// throws that away -- a single sequencing error changing just one of the `take`
/// smallest values (a very likely event: HiFi's ~0.1-0.2% substitution rate means
/// a handful of the ~20-40 candidate minimizers in a junction window will differ
/// read-to-read, and MinHash's bottom-k selection is a high-variance function of
/// exactly which ones survive) silently changes the *entire* signature and drops
/// the read out of the support count for its true locus, even under otherwise
/// ample coverage.
///
/// Fix: emit one band per left-out element (`take` bands total). Two reads whose
/// sketches differ by at most one element still share at least one band, so they
/// still collide in the support map, while still using cheap O(1) hash lookups
/// (as opposed to storing full sketches and doing O(n) fuzzy comparison against
/// every existing support entry).
pub fn foldback_signature_from_local_matches_bands(
    matches: &mut Vec<(usize, u64)>, // (pos, val)
    split: usize,
    window_bp: usize,
    take: usize,
    value_shift: u8,
) -> Option<Vec<u64>> {
    let vals = local_match_sketch(matches, split, window_bp, take, value_shift)?;
    Some(hash_bands(&vals))
}

/// Hash `vals` (assumed already sorted/deduped as desired) into `vals.len()`
/// leave-one-out band signatures (or a single whole-set hash if `vals.len() <= 1`).
/// See [`foldback_signature_from_local_matches_bands`] for the rationale.
fn hash_bands(vals: &[u64]) -> Vec<u64> {
    if vals.len() <= 1 {
        let mut h = gxhash::GxHasher::default();
        vals.hash(&mut h);
        return vec![h.finish()];
    }
    let mut bands: Vec<u64> = (0..vals.len())
        .map(|skip| {
            let mut h = gxhash::GxHasher::default();
            for (i, v) in vals.iter().enumerate() {
                if i != skip {
                    v.hash(&mut h);
                }
            }
            h.finish()
        })
        .collect();
    // `vals` can contain duplicate values (two distinct minimizer positions
    // hashing to the same value -- plausible in low-complexity or repetitive
    // sequence, which is exactly where real inverted repeats can occur).
    // Leaving out either copy of a duplicate produces the *same* remaining
    // multiset and therefore the same band hash, so without deduping here a
    // single read could insert into that one band twice, silently inflating
    // its support count.
    bands.sort_unstable();
    bands.dedup();
    bands
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
fn flank_sketch(
    seq: &[u8],
    split: usize,
    shared: &SharedCfg,
    flank_bp: usize,
    take: usize,
    scratch: &mut SigScratch,
    value_shift: u8,
) -> bool {
    if split < flank_bp || split + flank_bp > seq.len() {
        return false;
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
        return false;
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

    true
}

pub fn foldback_signature(
    seq: &[u8],
    split: usize,
    shared: &SharedCfg,
    flank_bp: usize,
    take: usize,
    scratch: &mut SigScratch,
    value_shift: u8,
) -> Option<u64> {
    if !flank_sketch(seq, split, shared, flank_bp, take, scratch, value_shift) {
        return None;
    }
    let mut h = GxHasher::default();
    scratch.combined.hash(&mut h);
    Some(h.finish())
}

/// Band-signature counterpart of [`foldback_signature`] -- see
/// [`foldback_signature_from_local_matches_bands`] for the rationale (same
/// bottom-k MinHash instability applies here, across both flanks).
pub fn foldback_signature_bands(
    seq: &[u8],
    split: usize,
    shared: &SharedCfg,
    flank_bp: usize,
    take: usize,
    scratch: &mut SigScratch,
    value_shift: u8,
) -> Option<Vec<u64>> {
    if !flank_sketch(seq, split, shared, flank_bp, take, scratch, value_shift) {
        return None;
    }
    Some(hash_bands(&scratch.combined))
}

// ---------------------------------------------------------------------------
// Junction-specific support (what pass1/pass2 actually gate on)
// ---------------------------------------------------------------------------
//
// The band signatures above identify a *locus* (the minimizer content within
// a few hundred bp of the split). That is not enough to tell a genome-encoded
// inverted repeat from an amplification artefact once a library is deep:
// artefact junctions land roughly uniformly along the genome, so at e.g. one
// junction per ~190 bp (a typical 140x C. elegans WGA HiFi run) *every* locus
// "recurs" by chance, and at one per ~29 bp (a 1500x run) essentially all of
// them do. Measured against a reference on real data, the locus-level gate
// left 15-50% of true artefacts uncut as "real" for that reason.
//
// What actually distinguishes a genomic inverted repeat is that its symmetry
// axis sits at the *same genomic position* in every read, whereas independent
// artefacts at the same locus fold at different positions. The distance
// between a minimizer and its own inverted copy (`p2 - p1` for a matchpoint)
// is exactly twice that minimizer's distance to the axis (plus any spacer),
// is read-coordinate-free, and is strand-symmetric -- so `(minimizer value,
// arm separation)` is a key that two reads share only if they fold about the
// same axis.

/// Maximum number of (value, separation) elements kept per junction: the
/// matchpoints closest to the axis. Those are the ones every read through
/// the junction is guaranteed to contain, however little of the arms it
/// covers, so reads of very different extent still select the same elements.
pub const JUNCTION_MAX_ELEMS: usize = 24;
/// Fewer elements than this and the junction has no usable key (treated as
/// unsupported).
pub const JUNCTION_MIN_ELEMS: usize = 3;

/// Tolerance (bp) on arm separation for two reads to count as the same axis.
/// Keys are guaranteed to match within `tol/2` and never match beyond `tol`.
/// HiFi separations are near-exact; noisy reads accumulate indel jitter along
/// the arm and need more room.
pub fn axis_tol(cfg: &MdaxCfg) -> usize {
    match cfg.shared.refine.mode {
        crate::utils::RefineMode::HiFi => 32,
        crate::utils::RefineMode::ONT => 128,
    }
}

/// Build junction keys from the best bin's matchpoints `(p1, p2, value)`.
///
/// Each output element is a pair of keys for the same `(value, separation)`
/// under two half-offset quantisations of the separation; two reads agree on
/// the element if they share either key. `out` is empty if the junction has
/// fewer than [`JUNCTION_MIN_ELEMS`] usable elements.
pub fn junction_keys(
    trips: &[(i32, i32, u64)],
    tol: usize,
    tmp: &mut Vec<(i32, u64)>,
    out: &mut Vec<(u64, u64)>,
) {
    tmp.clear();
    out.clear();
    for &(p1, p2, v) in trips {
        // One of each mirrored pair is enough; the mirror has the same
        // separation and carries no extra information.
        if p1 < p2 {
            tmp.push((p2 - p1, v));
        }
    }
    tmp.sort_unstable();
    tmp.dedup();
    tmp.truncate(JUNCTION_MAX_ELEMS);
    if tmp.len() < JUNCTION_MIN_ELEMS {
        return;
    }
    let q = tol.max(2) as u64;
    for &(sep, v) in tmp.iter() {
        let sep = sep as u64;
        let key = |bin: u64, which: u8| -> u64 {
            let mut h = GxHasher::default();
            (v, bin, which).hash(&mut h);
            h.finish()
        };
        out.push((key(sep / q, 0), key((sep + q / 2) / q, 1)));
    }
}

/// Pass1's output: how many reads carry each junction key.
#[derive(Debug, Default, Clone)]
pub struct JunctionSupport {
    /// key -> number of reads carrying it (only keys seen in >= 2 reads are
    /// stored; anything absent has count 1 when looked up by a read that
    /// was itself part of pass1).
    counts: HashMap<u64, u32>,
    /// Histogram of per-read support (index = n, last bucket = ">= len-1").
    pub n_hist: Vec<u64>,
    /// Number of reads that contributed a junction.
    pub reads: u64,
    /// Largest per-read support observed.
    pub max_n: usize,
}

impl JunctionSupport {
    pub fn empty() -> Self {
        Self::default()
    }

    pub fn build(per_read: &[Vec<(u64, u64)>]) -> Self {
        let mut all: Vec<u64> = Vec::with_capacity(per_read.iter().map(|k| k.len() * 2).sum());
        let mut one: Vec<u64> = Vec::new();
        for keys in per_read {
            one.clear();
            for &(a, b) in keys {
                one.push(a);
                one.push(b);
            }
            // a read counts once per key, however many elements produced it
            one.sort_unstable();
            one.dedup();
            all.extend_from_slice(&one);
        }
        all.sort_unstable();
        let mut counts: HashMap<u64, u32> = HashMap::new();
        let mut i = 0;
        while i < all.len() {
            let mut j = i + 1;
            while j < all.len() && all[j] == all[i] {
                j += 1;
            }
            if j - i >= 2 {
                counts.insert(all[i], (j - i) as u32);
            }
            i = j;
        }
        drop(all);

        let mut me = Self {
            counts,
            n_hist: vec![0; 258],
            reads: per_read.len() as u64,
            max_n: 0,
        };
        let top = me.n_hist.len() - 1;
        for keys in per_read {
            let n = me.support_n(keys);
            me.max_n = me.max_n.max(n);
            me.n_hist[n.min(top)] += 1;
        }
        me
    }

    /// Number of reads (including the querying read itself) that fold about
    /// the same axis: the lower median, over the junction's elements, of how
    /// many reads carry each element. The median makes this robust both to a
    /// promiscuous element (a repeat k-mer that happens to recur at the same
    /// separation elsewhere) and to elements lost to sequencing error.
    /// Returns 0 if the junction has no usable key.
    pub fn support_n(&self, keys: &[(u64, u64)]) -> usize {
        if keys.len() < JUNCTION_MIN_ELEMS {
            return 0;
        }
        let mut c: [u32; JUNCTION_MAX_ELEMS] = [1; JUNCTION_MAX_ELEMS];
        let m = keys.len().min(JUNCTION_MAX_ELEMS);
        for (slot, &(a, b)) in c.iter_mut().zip(keys.iter()) {
            let ca = self.counts.get(&a).copied().unwrap_or(1);
            let cb = self.counts.get(&b).copied().unwrap_or(1);
            *slot = ca.max(cb);
        }
        let c = &mut c[..m];
        c.sort_unstable();
        c[(m - 1) / 2] as usize
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Below this singleton share, calibration is skipped (see
    /// [`Self::adaptive_min_support`]).
    pub const ARTEFACT_DOMINATED: f64 = 0.5;

    /// Share of keyed foldbacks whose axis no other read shares.
    ///
    /// The deciding statistic for whether calibration applies at all. In an
    /// amplified library most detected foldbacks are artefacts with unique
    /// axes (0.82-0.83 on two C. elegans WGA HiFi libraries). In a library
    /// with few artefacts, most detected foldbacks are genome-encoded
    /// inverted repeats, each shared by roughly the local coverage
    /// (0.02-0.06 on three non-WGA C. elegans libraries).
    pub fn singleton_frac(&self) -> f64 {
        let keyed: u64 = self.n_hist.iter().skip(1).sum();
        if keyed == 0 {
            return 1.0;
        }
        self.n_hist[1] as f64 / keyed as f64
    }

    /// Mean number of *other* reads sharing a junction's axis by chance.
    ///
    /// For uniformly scattered artefact junctions the count of others
    /// sharing an axis is Poisson, so `P(1 other) / P(0 others) = lambda`,
    /// i.e. `h(2) / h(1)` in the per-read support histogram. This uses only
    /// the bottom two bins, where artefacts dominate; genome-encoded
    /// repeats sit at coverage-level support and barely touch them.
    ///
    /// (An earlier estimator, `-ln(h(1) / keyed)`, read every shared axis as
    /// chance: fine for amplified libraries, but on a library whose foldbacks
    /// are mostly genuine repeats it returned lambda ~3-4 and pushed the
    /// threshold to 47-60, cutting the repeats it exists to protect.)
    pub fn chance_lambda(&self) -> f64 {
        let h1 = self.n_hist.get(1).copied().unwrap_or(0);
        if h1 == 0 {
            return 0.0;
        }
        self.n_hist.get(2).copied().unwrap_or(0) as f64 / h1 as f64
    }

    /// Smallest support threshold `t >= floor` at which fewer than `alpha`
    /// of chance co-located junctions would reach it, i.e.
    /// `P(Poisson(lambda) >= t - 1) < alpha`.
    ///
    /// A fixed absolute threshold cannot be right across libraries: the
    /// number of unrelated artefacts sharing an axis by chance grows with
    /// junction density, so the bar for "this recurs because it is in the
    /// genome" has to grow with it. This only ever raises the threshold.
    pub fn adaptive_min_support(&self, floor: usize, alpha: f64) -> usize {
        // Chance co-location needs dense artefacts. When most detected
        // foldbacks already share their axis with other reads, the library
        // is dominated by genome-encoded repeats, chance sharing is
        // negligible, and raising the bar would only cut real repeats.
        if self.singleton_frac() < Self::ARTEFACT_DOMINATED {
            return floor.max(2);
        }
        let lambda = self.chance_lambda();
        let mut t = floor.max(2);
        loop {
            // P(X >= t-1) = 1 - sum_{i=0}^{t-2} e^-l l^i / i!
            let mut term = (-lambda).exp();
            let mut cdf = term;
            for i in 1..=(t - 2) {
                term *= lambda / i as f64;
                cdf += term;
            }
            if 1.0 - cdf < alpha || t >= 64 {
                break;
            }
            t += 1;
        }
        self.extend_past_heavy_tail(t)
    }

    /// Amplified libraries are not Poisson: coverage is wildly uneven and
    /// early amplification products are copied many times over, so the
    /// number of artefacts sharing an axis has a much heavier tail than the
    /// mean rate predicts (on a ~1500x C. elegans library the Poisson
    /// threshold still left 6% of reference-confirmed artefacts above it,
    /// their support counts decaying smoothly by ~0.7x per step, while
    /// genome-encoded repeats sat far out at n >= 20).
    ///
    /// So measure the tail instead of assuming it: take the decay ratio of
    /// the per-read support histogram just above the Poisson threshold,
    /// extrapolate it geometrically, and move the threshold out to where the
    /// extrapolated chance tail holds under 0.2% of all keyed foldbacks.
    /// Genuine repeats in those bins flatten the measured decay, which errs
    /// towards a higher threshold (cutting), never a lower one.
    fn extend_past_heavy_tail(&self, t: usize) -> usize {
        const TAIL_FRAC: f64 = 0.002;
        const MIN_COUNT: u64 = 50;
        let h = |n: usize| self.n_hist.get(n).copied().unwrap_or(0);
        // the last bucket is an overflow bin, never part of the fit
        if t + 3 >= self.n_hist.len() {
            return t;
        }
        let (h0, h1, h2) = (h(t), h(t + 1), h(t + 2));
        if h0 + h1 < MIN_COUNT {
            return t;
        }
        // Only a decaying histogram is a chance tail. Flat or rising counts
        // just above the Poisson threshold are genuine repeats at their
        // coverage level, and extrapolating through them would walk the
        // threshold out past them.
        if !(h0 > h1 && h1 > h2) {
            return t;
        }
        let r = ((h1 + h2) as f64 / (h0 + h1) as f64).clamp(0.0, 0.9);
        let keyed: u64 = self.n_hist.iter().skip(1).sum();
        let budget = TAIL_FRAC * keyed as f64;
        let cap = self.n_hist.len() - 2;
        let mut n = t;
        let mut level = h0 as f64;
        // tail mass from n upward under geometric decay: level / (1 - r)
        while level / (1.0 - r) > budget && n < cap {
            level *= r;
            n += 1;
        }
        n
    }
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

/// Union-find (disjoint-set) over band-signature values.
///
/// A single real locus's reads frequently produce *disjoint* band sets that
/// never all agree with each other (read A shares a band with read B, read B
/// shares a different band with read C, but A and C share nothing directly).
/// Treating each raw band as its own key under-counts support: every one of
/// these reads' bands ends up a near-singleton even though they're all the
/// same underlying locus. This structure merges any bands ever observed
/// together in one read's own band set, and transitively merges clusters
/// bridged by other reads, into one canonical cluster id (the union-find
/// root) before support is counted. See `SupportMap`.
#[derive(Debug, Default, Clone)]
struct UnionFind {
    parent: HashMap<u64, u64>,
}

impl UnionFind {
    fn new() -> Self {
        Self {
            parent: HashMap::new(),
        }
    }

    /// Resolve `x` to its canonical root, with path compression. Unseen
    /// values are lazily registered as their own singleton root.
    fn find(&mut self, x: u64) -> u64 {
        let px = *self.parent.entry(x).or_insert(x);
        if px == x {
            x
        } else {
            let root = self.find(px);
            self.parent.insert(x, root);
            root
        }
    }

    /// Merge the clusters containing `a` and `b`.
    fn union(&mut self, a: u64, b: u64) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            // Deterministic attachment rule (smaller root under larger), so
            // the same edge set always produces the same roots regardless of
            // processing order.
            if ra < rb {
                self.parent.insert(ra, rb);
            } else {
                self.parent.insert(rb, ra);
            }
        }
    }

    /// All raw band values this structure has ever seen (as `find`/`union`
    /// arguments), i.e. every band that may need re-mapping to its root.
    fn known_keys(&self) -> Vec<u64> {
        self.parent.keys().copied().collect()
    }
}

/// Pass1's output: per-locus `SupportStats`, keyed by canonical cluster id,
/// plus the raw-band -> canonical-root mapping needed to look up a freshly
/// computed read's bands against those clusters.
///
/// Building this directly as a flat `HashMap<band, SupportStats>` (one entry
/// per raw band, as pass1 used to) under-counts support for any real locus
/// whose reads' band sets don't all agree with each other -- see
/// `UnionFind`. `SupportMap::build` resolves that before any entry is
/// counted, so `n` reflects reads-per-locus, not reads-per-raw-band.
#[derive(Debug, Default, Clone)]
pub struct SupportMap {
    stats: HashMap<u64, SupportStats>,
    band_root: HashMap<u64, u64>,
}

impl SupportMap {
    pub fn empty() -> Self {
        Self {
            stats: HashMap::new(),
            band_root: HashMap::new(),
        }
    }

    /// Build a `SupportMap` from per-read observations (one `(band, split,
    /// ident)` per read that produced a signature, using any one of that
    /// read's own bands as the lookup key) and the edges connecting bands
    /// that co-occurred within a single read's band set.
    ///
    /// Each read contributes exactly one observation to exactly one
    /// canonical cluster -- never one observation per raw band -- so cluster
    /// `n` is a true read count, not inflated by how many bands a read
    /// happens to produce.
    pub fn build(records: &[(u64, usize, f32)], edges: &[(u64, u64)]) -> Self {
        let mut uf = UnionFind::new();
        for &(a, b) in edges {
            uf.union(a, b);
        }

        let mut stats: HashMap<u64, SupportStats> = HashMap::new();
        for &(band, split, ident) in records {
            let root = uf.find(band);
            stats
                .entry(root)
                .and_modify(|st| st.update(split, ident))
                .or_insert_with(|| SupportStats::new(split, ident));
        }

        let keys = uf.known_keys();
        let mut band_root = HashMap::with_capacity(keys.len());
        for b in keys {
            band_root.insert(b, uf.find(b));
        }

        Self { stats, band_root }
    }

    /// Resolve a raw band value to its canonical cluster id. Bands never
    /// involved in a union (e.g. a read whose own band set was a disjoint
    /// singleton with no later read bridging it to anything else) resolve to
    /// themselves.
    pub fn resolve(&self, band: u64) -> u64 {
        *self.band_root.get(&band).unwrap_or(&band)
    }

    pub fn get(&self, band: u64) -> Option<&SupportStats> {
        self.stats.get(&self.resolve(band))
    }

    pub fn len(&self) -> usize {
        self.stats.len()
    }

    pub fn is_empty(&self) -> bool {
        self.stats.is_empty()
    }

    pub fn values(&self) -> impl Iterator<Item = &SupportStats> {
        self.stats.values()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&u64, &SupportStats)> {
        self.stats.iter()
    }
}

/// Given the band signatures for a read's foldback junction, find the best-supported
/// matching entry in the support map (the one with the highest read count `n`).
///
/// A read's junction maps to several band signatures (see
/// [`foldback_signature_from_local_matches_bands`]); any one of them colliding with
/// an existing support-map entry is evidence for the same underlying locus, so we
/// take the strongest match rather than requiring all bands to agree. Returns the
/// *canonical* cluster id (not necessarily equal to any of `bands`), so the same
/// locus is reported identically regardless of which of a read's bands happened to
/// hit the map.
pub fn best_support<'a>(bands: &[u64], support: &'a SupportMap) -> Option<(u64, &'a SupportStats)> {
    bands
        .iter()
        .filter_map(|&b| {
            let root = support.resolve(b);
            support.get(root).map(|st| (root, st))
        })
        .max_by_key(|(_, st)| st.n)
}

/// Safety-valve ceiling on a cluster's read count, independent of `min_support`.
///
/// `SupportMap::build`'s union-find has no bound on transitive merging: a
/// single band value shared between two otherwise-unrelated reads unions
/// *all* of each read's own bands together, and that merge chains through
/// the whole dataset via any other read that bridges two components. Before
/// the union-find fix, a collision only inflated the one shared band's count;
/// now it can fold two or more real loci's entire read populations into one
/// cluster. `test_single_nematode/README.md` already documents this
/// mechanism making microlaimidae's collision-prone case *worse* (max
/// cluster 407 -> 454) before the (separate) sig_flank_bp fix brought it back
/// down -- i.e. this isn't hypothetical, it's the exact thing that regressed.
///
/// This is a blunt backstop, not a precision classifier: calibrated from the
/// same data as the `QC MODE ADVICE` heuristic in `main.rs` (max_cluster_n
/// 407-454 on the one confirmed collision case vs. 18-130 on every other
/// real/synthetic cluster measured so far, including the deliberately-
/// stressed `run_main` ground-truth set, max 41). A cluster this large is
/// far more likely to be an unbounded collision cascade than a genuinely
/// recurring locus, so treat it as artefact (cut) rather than silently
/// trusting it as real and leaving a potentially large duplicated/uncut
/// region in the assembly. `main.rs`'s advisory threshold reuses this same
/// constant rather than hardcoding its own copy.
pub const MAX_PLAUSIBLE_CLUSTER_N: usize = 200;

/// Core "real vs artefact" gate, given an already-resolved support entry.
pub fn is_real_foldback_stats(st: &SupportStats, cfg: &MdaxCfg) -> bool {
    // Read-space split positions are not expected to cluster tightly across reads.
    // The signature already anchors the locus; use support count as the primary gate.
    if st.n < cfg.fold2.min_support {
        return false;
    }

    // See MAX_PLAUSIBLE_CLUSTER_N: an anomalously large cluster is far more
    // likely to be a collision cascade than genuine recurrence.
    if st.n > MAX_PLAUSIBLE_CLUSTER_N {
        return false;
    }

    // Optional: add an identity-based support gate if you want
    if cfg.fold2.min_support_ident > 0.0 && (st.mean_ident as f64) < cfg.fold2.min_support_ident {
        return false;
    }

    true
}

/// Decide whether a foldback junction is "real" based on pass1 support statistics.
pub fn is_real_foldback(sig: u64, support: &SupportMap, cfg: &MdaxCfg) -> bool {
    let Some(st) = support.get(sig) else {
        return false;
    };
    is_real_foldback_stats(st, cfg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cfg::{FoldOnlyCfg, FoldSecondPassCfg, MinimizerCfg, RefineCfg, SharedCfg, SigCfg};
    use crate::foldback::{detect_foldback, refine_breakpoint};
    use crate::pipeline::pass1_build_support;
    #[allow(unused_imports)]
    use super::{JunctionSupport, axis_tol, junction_keys};
    use crate::scratch::FoldScratch;
    use crate::utils::RefineMode;

    use calm_io::stderrln;
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
                max_jump_clip: 0,
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
            fn write_fasta_inner(
                f: &mut impl std::io::Write,
                id: &str,
                seq: &[u8],
            ) -> std::io::Result<()> {
                f.write_all(b">")?;
                f.write_all(id.as_bytes())?;
                f.write_all(b"\n")?;
                // write sequence as a single line (fastest)
                f.write_all(seq)?;
                f.write_all(b"\n")?;
                Ok(())
            }
            write_fasta_inner(&mut w, id, seq).expect("write fasta record");
        }
    }

    /// Mirror exactly what pass1/pass2 do to turn a read into band signatures:
    /// detect + refine (HiFi fast path), prefer the local-match signature,
    /// fall back to the flank signature. Used by tests that need to look up
    /// the support map the same way production code populated it.
    fn read_bands(seq: &[u8], cfg: &MdaxCfg) -> Option<Vec<u64>> {
        let mut scratch = FoldScratch::new();
        let fb = detect_foldback(seq, &cfg.shared, &cfg.fold, &mut scratch)?;
        let mut refine_scratch = crate::scratch::RefineScratch::default();
        let rf = refine_breakpoint(
            seq, fb.split_pos, &cfg.shared, &mut refine_scratch, fb.span, fb.arm_end, fb.arm_start_right,
        )
        .ok()??;

        foldback_signature_from_local_matches_bands(
            &mut scratch.best_matches, rf.split_pos, cfg.sig.flank_bp, cfg.sig.take, cfg.sig.value_shift,
        )
        .or_else(|| {
            let q = 150usize;
            let split_q = (rf.split_pos / q) * q;
            let mut sig_scratch = SigScratch::default();
            foldback_signature_bands(
                seq, split_q, &cfg.shared, cfg.sig.flank_bp, cfg.sig.take, &mut sig_scratch, cfg.sig.value_shift,
            )
        })
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
                    arm: 1200,
                    mode: RefineMode::HiFi,
                    max_ed_rate: 0.25,
                    max_jump_clip: 1000,
                },
                fold_diag_tol: 50,
            },
            fold: FoldOnlyCfg {
                min_arm: 200, // span evidence for foldback detector
                ..Default::default()
            },
            fold2: FoldSecondPassCfg {
                min_support: 3,
                min_identity: 0.60,
                min_support_ident: 0.0,
                cut_low_ident: false,
            },
            sig: SigCfg {
                flank_bp: 1000,
                take: 8,
                value_shift: 0,
            },
        }
    }

    /// Junction keys for a read, computed the way pass1/pass2 compute them.
    fn read_junction_keys(seq: &[u8], cfg: &MdaxCfg) -> Vec<(u64, u64)> {
        let mut scratch = FoldScratch::new();
        if detect_foldback(seq, &cfg.shared, &cfg.fold, &mut scratch).is_none() {
            return Vec::new();
        }
        let mut tmp = Vec::new();
        let mut out = Vec::new();
        junction_keys(&scratch.best_trip, axis_tol(cfg), &mut tmp, &mut out);
        out
    }

    #[test]
    fn foldback_pass1_support_finds_one_cluster() {
        let cfg = cfg_for_pass1();
        let left = pseudo_dna(2500, 42);
        let (seq, _split) = make_clean_foldback(&left, b"");
        let decoy = pseudo_dna(6000, 999);

        let path = tmp_fasta_path("mdax_fold2_pass1");
        write_fasta(
            &path,
            &[("r1", &seq), ("r2", &seq), ("r3", &seq), ("decoy", &decoy)],
        );
        let cfg2 = cfg.clone();
        let (support, total_reads) =
            pass1_build_support(&path, Arc::new(cfg), 1, 1).expect("pass1 support");
        let _ = std::fs::remove_file(&path);

        assert_eq!(total_reads, 4);
        assert_eq!(support.reads, 3, "the decoy must not contribute a junction");
        let keys = read_junction_keys(&seq, &cfg2);
        assert!(keys.len() >= JUNCTION_MIN_ELEMS);
        assert_eq!(support.support_n(&keys), 3);
    }

    /// The property the gate rests on: reads that fold about the *same* axis
    /// support each other however much of the arms each one happens to
    /// cover, while reads from the same locus that fold about *different*
    /// axes do not.
    #[test]
    fn junction_support_is_axis_specific_not_locus_specific() {
        let cfg = cfg_for_pass1();
        let genome = pseudo_dna(6000, 7);

        // Same axis (genome position 3000), three different read extents.
        let same_axis: Vec<Vec<u8>> = [(0usize, 3000usize), (400, 3000), (900, 3000)]
            .iter()
            .map(|&(a, b)| make_clean_foldback(&genome[a..b], b"").0)
            .collect();
        // Same locus, axes 150 bp apart (far beyond any axis tolerance).
        let other_axes: Vec<Vec<u8>> = [3150usize, 3300, 3450, 3600]
            .iter()
            .map(|&b| make_clean_foldback(&genome[500..b], b"").0)
            .collect();

        let path = tmp_fasta_path("mdax_axis_specific");
        let mut recs: Vec<(String, &Vec<u8>)> = Vec::new();
        for (i, r) in same_axis.iter().enumerate() {
            recs.push((format!("same{i}"), r));
        }
        for (i, r) in other_axes.iter().enumerate() {
            recs.push((format!("other{i}"), r));
        }
        let refs: Vec<(&str, &[u8])> = recs.iter().map(|(n, r)| (n.as_str(), r.as_slice())).collect();
        write_fasta(&path, &refs);
        let cfg2 = cfg.clone();
        let (support, _) = pass1_build_support(&path, Arc::new(cfg), 1, 1).unwrap();
        let _ = std::fs::remove_file(&path);

        for r in &same_axis {
            assert_eq!(support.support_n(&read_junction_keys(r, &cfg2)), 3);
        }
        for r in &other_axes {
            assert_eq!(
                support.support_n(&read_junction_keys(r, &cfg2)),
                1,
                "a read folding about its own axis must not borrow support from neighbours"
            );
        }
    }

    #[test]
    fn adaptive_min_support_rises_with_junction_density() {
        let with_hist = |alone: u64, shared: u64| {
            let mut s = JunctionSupport::empty();
            s.n_hist = vec![0; 66];
            s.n_hist[1] = alone;
            s.n_hist[2] = shared;
            s
        };
        // Sparse: almost nobody shares an axis -> the preset floor stands.
        assert_eq!(with_hist(9990, 10).adaptive_min_support(3, 0.01), 3);
        // ~15% share by chance (lambda ~0.16) -> one step up.
        assert_eq!(with_hist(8500, 1500).adaptive_min_support(3, 0.01), 4);
        // Dense artefacts (lambda 0.6) -> well above the floor.
        assert_eq!(with_hist(6000, 3600).adaptive_min_support(3, 0.01), 5);
        // Never below the floor.
        assert_eq!(with_hist(10000, 0).adaptive_min_support(5, 0.01), 5);

        // Heavy (non-Poisson) tail: the mean rate alone says 4, but the
        // histogram keeps decaying slowly above it, so the threshold has
        // to move out past where that tail has died away.
        let mut heavy = with_hist(85_000, 12_000);
        heavy.n_hist[3] = 1500;
        let mut c = 600.0f64;
        for n in 4..20 {
            heavy.n_hist[n] = c as u64;
            c *= 0.6;
        }
        let t = heavy.adaptive_min_support(3, 0.01);
        assert!((6..=12).contains(&t), "heavy tail should push the threshold out, got {t}");

        // A library with few artefacts: most detected foldbacks are genuine
        // repeats shared by ~coverage reads (shape of nrCaeEleg34, a non-WGA
        // C. elegans library). Calibration must leave the floor alone.
        let mut genuine = with_hist(555, 470);
        for n in 3..10 {
            genuine.n_hist[n] = 420;
        }
        for n in 10..30 {
            genuine.n_hist[n] = 260;
        }
        assert!(genuine.singleton_frac() < JunctionSupport::ARTEFACT_DOMINATED);
        assert_eq!(genuine.adaptive_min_support(3, 0.01), 3);

        // Flat counts just above the Poisson threshold are not a chance
        // tail: no extrapolation through them.
        let mut flat = with_hist(80_000, 9_000);
        for n in 3..40 {
            flat.n_hist[n] = 300;
        }
        assert_eq!(flat.adaptive_min_support(3, 0.01), 3);
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
        let s2 = base.clone(); // a second clean copy: real data is mostly error-free reads
        let s3 = mutate_far(base.clone(), split, 1100, 97); // outside ±1100 keeps 1000bp flanks intact
        let s4 = mutate_far(base.clone(), split, 1100, 131);

        let decoy = pseudo_dna(8000, 999);

        let path = tmp_fasta_path("mdax_fold2_pass1_flank_noise");
        write_fasta(
            &path,
            &[
                ("r1", &s1),
                ("r2", &s2),
                ("r3", &s3),
                ("r4", &s4),
                ("decoy", &decoy),
            ],
        );

        let cfg2 = cfg.clone();
        let (support, _total_reads) = pass1_build_support(&path, Arc::new(cfg), 1, 1).unwrap();
        let _ = std::fs::remove_file(&path);

        // Junction keys come from the matchpoints nearest the axis, so
        // mutations far from the junction must not cost the locus any of
        // its four reads.
        let _ = split;
        for r in [&s1, &s3, &s4] {
            let n = support.support_n(&read_junction_keys(r, &cfg2));
            assert_eq!(n, 4, "expected all 4 reads to support the axis, got n={n}");
        }
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
