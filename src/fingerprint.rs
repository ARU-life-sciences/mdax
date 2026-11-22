//! Fold-axis junction support: the real-vs-artefact gate.
//!
//! Pass1 derives a small set of keys for every detected foldback from the
//! matchpoints of its best anti-diagonal bin ([`junction_keys`]) and counts,
//! per key, how many reads carry it ([`JunctionSupport`]). Pass2 asks how many
//! reads share a read's fold axis ([`JunctionSupport::support_n`]) and keeps
//! the foldback as genome-templated ("real") when that clears `min_support`.
//!
//! This replaced an earlier locus-level scheme (MinHash band signatures of
//! the junction flanks, merged by union-find into per-locus clusters); see
//! the git history before 2026-11-19 for that code.

use gxhash::{GxHasher, HashMap, HashMapExt};
use std::hash::{Hash, Hasher};

use crate::cfg::MdaxCfg;
// ---------------------------------------------------------------------------
// Junction-specific support (what pass1/pass2 actually gate on)
// ---------------------------------------------------------------------------
//
// The legacy locus signatures (band hashes of the minimizer content within
// a few hundred bp of the split) identified a *locus*. That is not enough to tell a genome-encoded
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cfg::{FoldOnlyCfg, FoldSecondPassCfg, MinimizerCfg, RefineCfg, SharedCfg};
    use crate::foldback::detect_foldback;
    use crate::pipeline::pass1_build_support;
    use crate::scratch::FoldScratch;
    use crate::utils::RefineMode;

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


    fn cfg_for_pass1() -> MdaxCfg {
        // Make it easy to detect a foldback and key its junction.
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
                cut_low_ident: false,
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
        for h in &mut heavy.n_hist[4..20] {
            *h = c as u64;
            c *= 0.6;
        }
        let t = heavy.adaptive_min_support(3, 0.01);
        assert!((6..=12).contains(&t), "heavy tail should push the threshold out, got {t}");

        // A library with few artefacts: most detected foldbacks are genuine
        // repeats shared by ~coverage reads (shape of nrCaeEleg34, a non-WGA
        // C. elegans library). Calibration must leave the floor alone.
        let mut genuine = with_hist(555, 470);
        genuine.n_hist[3..10].fill(420);
        genuine.n_hist[10..30].fill(260);
        assert!(genuine.singleton_frac() < JunctionSupport::ARTEFACT_DOMINATED);
        assert_eq!(genuine.adaptive_min_support(3, 0.01), 3);

        // Flat counts just above the Poisson threshold are not a chance
        // tail: no extrapolation through them.
        let mut flat = with_hist(80_000, 9_000);
        flat.n_hist[3..40].fill(300);
        assert_eq!(flat.adaptive_min_support(3, 0.01), 3);
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
}
