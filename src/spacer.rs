//! Deciding which side of a displaced foldback its spacer belongs to.
//!
//! A template-switch foldback is usually *displaced*: the read is
//! `X S RC(X)`, where `X` is the mirrored arm and `S` an un-mirrored stretch
//! between the arms. `S` is genuine sequence, but it belongs to one strand or
//! the other: either the genome reads `X S` (the polymerase ran on past the
//! end of the mirrored region before switching strand), or it reads
//! `X RC(S)` (it switched first and ran on along the other strand). A single
//! read cannot tell these apart -- both give the same `X S RC(X)` -- and
//! guessing wrong either keeps a strand-switched stub on the corrected read
//! or throws away good sequence (measured on C. elegans WGA HiFi: keeping
//! the spacer always left chimeric sequence in ~45% of cuts, dropping it
//! always discarded ~9% of all bases).
//!
//! The rest of the library can tell them apart. Exactly one of the two
//! boundaries of `S` is a real genomic junction, so the k-mer spanning it
//! occurs in roughly coverage-many other reads; the k-mer spanning the other
//! boundary exists only in this read (and its clonal copies). Collect both
//! boundary k-mers for every detected foldback, count them across the whole
//! input, and cut at whichever boundary is the strand switch.

use anyhow::Result;
use crossbeam_channel as chan;
use gxhash::{HashMap, HashMapExt, HashSet, HashSetExt};
use std::path::Path;
use std::sync::{Arc, OnceLock};
use std::thread;

use crate::{cfg::MdaxCfg, foldback, io, scratch::FoldScratch};

/// Length of the boundary k-mers (fits a 2-bit u64).
pub const BK: usize = 25;
/// Spacers shorter than this are treated as no spacer (cut at the axis).
const MIN_SPACER: usize = 30;

/// Which strand the un-mirrored stretch between the arms belongs to.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SpacerSide {
    /// Genome reads `X S`: keep `S` with the left copy.
    Left,
    /// Genome reads `X RC(S)`: keep `S` with the right copy.
    Right,
    /// No evidence either way: drop `S`.
    Unknown,
}

/// Exact mirror geometry of a detected foldback, in the coordinates of the
/// sequence the matchpoints came from: the mirrored arms end at `a` (left)
/// and start at `b` (right); `[a, b)` is the spacer (empty for a hairpin).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Mirror {
    pub a: usize,
    pub b: usize,
}

/// Recover the exact symmetry axis from the matchpoints nearest it and walk
/// the mirror inward base by base to the first sustained mismatch.
///
/// `trips` are `(p1, p2, value)` from `detect_foldback` (forward k-mer starts
/// of a k-mer and of its reverse-complement copy). For an exact mirror,
/// base `p` pairs with base `T - 1 - p` where `T = p1 + p2 + k`.
pub fn mirror(seq: &[u8], trips: &[(i32, i32, u64)], k: usize) -> Option<Mirror> {
    // Innermost type-1 matchpoints vote for T; indels far out along the arm
    // shift it, so only the closest few count.
    let mut inner: Vec<(i32, i32)> = trips
        .iter()
        .filter(|t| t.0 < t.1)
        .map(|t| (t.1 - t.0, t.0 + t.1 + k as i32))
        .collect();
    if inner.is_empty() {
        return None;
    }
    inner.sort_unstable();
    inner.truncate(12);
    let mut votes: Vec<i32> = inner.iter().map(|x| x.1).collect();
    votes.sort_unstable();
    let (mut best, mut best_n, mut i) = (votes[0], 0usize, 0usize);
    while i < votes.len() {
        let mut j = i;
        while j < votes.len() && votes[j] == votes[i] {
            j += 1;
        }
        if j - i > best_n {
            best_n = j - i;
            best = votes[i];
        }
        i = j;
    }
    if best_n < 2 || best <= 0 {
        return None;
    }
    let t = best as usize;
    if t > 2 * seq.len() {
        return None;
    }
    // Start from the end of the innermost k-mer that sits exactly on this axis.
    let start = trips
        .iter()
        .filter(|x| x.0 < x.1 && (x.0 + x.1 + k as i32) as usize == t)
        .map(|x| x.0 as usize + k)
        .max()?;
    let pairs = |j: usize| -> bool {
        let m = t - 1 - j;
        m < seq.len() && crate::utils::comp(seq[m]) == seq[j]
    };
    let mut j = start;
    // walk while j is left of its own mirror
    while 2 * j + 1 < t {
        if !pairs(j) {
            // sustained mismatch (>=3 of the next 8) ends the mirror; an
            // isolated mismatch is a sequencing error
            let mut mm = 0;
            let mut q = j;
            while q < j + 8 && 2 * q + 1 < t {
                if !pairs(q) {
                    mm += 1;
                }
                q += 1;
            }
            if mm >= 3 {
                break;
            }
        }
        j += 1;
    }
    let a = j;
    let b = t - a;
    if b < a || b > seq.len() {
        return None;
    }
    Some(Mirror { a, b })
}

#[inline]
fn enc(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Canonical 2-bit encoding of `s` (length `BK`).
fn canon(s: &[u8]) -> Option<u64> {
    let (mut f, mut r) = (0u64, 0u64);
    for &c in s {
        let x = enc(c)?;
        f = (f << 2) | x;
        r = (r >> 2) | ((3 - x) << (2 * (BK - 1)));
    }
    Some(f.min(r))
}

/// The two boundary k-mers of the spacer `[a, b)`, or `None` if the spacer
/// is too short to matter or the k-mers fall off the read.
pub fn boundary_kmers(seq: &[u8], m: Mirror) -> Option<(u64, u64)> {
    if m.b - m.a < MIN_SPACER {
        return None;
    }
    let h = BK / 2;
    if m.a < h || m.b + (BK - h) > seq.len() {
        return None;
    }
    Some((canon(&seq[m.a - h..m.a - h + BK])?, canon(&seq[m.b - h..m.b - h + BK])?))
}

/// Library-wide counts of every boundary k-mer collected.
#[derive(Debug, Default)]
pub struct SpacerEvidence {
    counts: HashMap<u64, u32>,
}

static EVIDENCE: OnceLock<SpacerEvidence> = OnceLock::new();

/// Make `ev` available to the cutters for the rest of the run.
pub fn install(ev: SpacerEvidence) {
    let _ = EVIDENCE.set(ev);
}

impl SpacerEvidence {
    fn count(&self, x: u64) -> u32 {
        self.counts.get(&x).copied().unwrap_or(0)
    }

    /// Decide which side `[a, b)` belongs to from its boundary k-mers.
    pub fn side(&self, ka: u64, kb: u64) -> SpacerSide {
        let (ca, cb) = (self.count(ka), self.count(kb));
        // The read itself contributes 1 to both. A genomic junction is seen
        // in other reads; the strand switch is (nearly) unique to this one.
        let strong = |x: u32, y: u32| x >= 3 && x >= 4 * y.max(1);
        if strong(ca, cb) {
            SpacerSide::Left
        } else if strong(cb, ca) {
            SpacerSide::Right
        } else {
            SpacerSide::Unknown
        }
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }
}

/// Resolve a foldback's spacer using installed evidence. Returns the mirror
/// geometry and the side, or `None` if the exact geometry can't be recovered.
pub fn resolve(seq: &[u8], trips: &[(i32, i32, u64)], k: usize) -> Option<(Mirror, SpacerSide)> {
    let m = mirror(seq, trips, k)?;
    let side = match (EVIDENCE.get(), boundary_kmers(seq, m)) {
        (Some(ev), Some((ka, kb))) => ev.side(ka, kb),
        _ => SpacerSide::Unknown,
    };
    Some((m, side))
}

fn stream_batches<P: AsRef<Path>>(
    input: P,
    cap: usize,
) -> (chan::Receiver<Vec<Vec<u8>>>, thread::JoinHandle<Result<()>>) {
    let (tx, rx) = chan::bounded::<Vec<Vec<u8>>>(cap);
    let path = input.as_ref().to_owned();
    let h = thread::spawn(move || -> Result<()> {
        let mut r = io::open_fasta_reader(path)?;
        let mut batch = Vec::with_capacity(256);
        while let Some(rec) = r.next() {
            batch.push(rec?.seq().to_vec());
            if batch.len() == 256 {
                tx.send(std::mem::take(&mut batch))
                    .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
            }
        }
        if !batch.is_empty() {
            tx.send(batch).map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
        }
        Ok(())
    });
    (rx, h)
}

/// Two extra passes over `input`: collect the boundary k-mers of every
/// detected displaced foldback, then count them in every read.
pub fn build_evidence<P: AsRef<Path> + Clone>(
    input: P,
    cfg: Arc<MdaxCfg>,
    threads: usize,
    chan_cap: usize,
) -> Result<SpacerEvidence> {
    // pass A: queries
    let (rx, reader) = stream_batches(input.clone(), chan_cap);
    let mut hs = Vec::new();
    for _ in 0..threads {
        let rx = rx.clone();
        let cfg = cfg.clone();
        hs.push(thread::spawn(move || -> Vec<u64> {
            let mut out = Vec::new();
            let mut sc = FoldScratch::new();
            let k = cfg.shared.minimizer.k;
            while let Ok(batch) = rx.recv() {
                for seq in batch {
                    if foldback::detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut sc).is_none() {
                        continue;
                    }
                    if let Some(m) = mirror(&seq, &sc.best_trip, k)
                        && let Some((a, b)) = boundary_kmers(&seq, m)
                    {
                        out.push(a);
                        out.push(b);
                    }
                }
            }
            out
        }));
    }
    drop(rx);
    reader.join().map_err(|_| anyhow::anyhow!("reader panicked"))??;
    let mut queries: HashSet<u64> = HashSet::new();
    for h in hs {
        queries.extend(h.join().map_err(|_| anyhow::anyhow!("worker panicked"))?);
    }
    let queries = Arc::new(queries);

    // pass B: count every canonical BK-mer of every read that is a query
    let (rx, reader) = stream_batches(input, chan_cap);
    let mut hs = Vec::new();
    for _ in 0..threads {
        let rx = rx.clone();
        let q = queries.clone();
        hs.push(thread::spawn(move || -> HashMap<u64, u32> {
            let mut c: HashMap<u64, u32> = HashMap::new();
            let mask = if BK == 32 { u64::MAX } else { (1u64 << (2 * BK)) - 1 };
            while let Ok(batch) = rx.recv() {
                for seq in batch {
                    let (mut f, mut r, mut n) = (0u64, 0u64, 0usize);
                    for &ch in &seq {
                        match enc(ch) {
                            Some(x) => {
                                f = ((f << 2) | x) & mask;
                                r = (r >> 2) | ((3 - x) << (2 * (BK - 1)));
                                n += 1;
                                if n >= BK {
                                    let key = f.min(r);
                                    if q.contains(&key) {
                                        *c.entry(key).or_insert(0) += 1;
                                    }
                                }
                            }
                            None => n = 0,
                        }
                    }
                }
            }
            c
        }));
    }
    drop(rx);
    reader.join().map_err(|_| anyhow::anyhow!("reader panicked"))??;
    let mut counts: HashMap<u64, u32> = HashMap::with_capacity(queries.len());
    for h in hs {
        for (k, v) in h.join().map_err(|_| anyhow::anyhow!("worker panicked"))? {
            *counts.entry(k).or_insert(0) += v;
        }
    }
    Ok(SpacerEvidence { counts })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x & 3) as usize]
            })
            .collect()
    }

    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| crate::utils::comp(b)).collect()
    }

    /// Trips for an exact mirror: every k-mer of `[lo, hi)` paired with its copy.
    fn trips_for(t: usize, lo: usize, hi: usize, k: usize) -> Vec<(i32, i32, u64)> {
        (lo..hi - k)
            .step_by(7)
            .map(|p1| (p1 as i32, (t - k - p1) as i32, 0))
            .collect()
    }

    #[test]
    fn mirror_finds_exact_spacer_boundaries() {
        let g = dna(5000, 9);
        // forward copy to 3000, reverse copy from genome 3600 down to 1000
        let mut seq = g[..3000].to_vec();
        seq.extend_from_slice(&rc(&g[1000..3600]));
        // mirrored arm is genome 1000..3000 at read 1000..3000, mirrored at
        // read 3600..5600; T = 3000 + 3600 = 6600
        let k = 17;
        let m = mirror(&seq, &trips_for(6600, 1000, 2900, k), k).expect("mirror");
        assert!(m.a.abs_diff(3000) <= 3, "{m:?}");
        assert!(m.b.abs_diff(3600) <= 3, "{m:?}");
    }

    #[test]
    fn genomic_boundary_kmer_decides_the_side() {
        let g = dna(6000, 3);
        let k = 17;
        // spacer belongs to the left copy: forward to 3000, reverse from 2400
        let mut left = g[..3000].to_vec();
        left.extend_from_slice(&rc(&g[1000..2400]));
        // spacer belongs to the right copy: forward to 3000, reverse from 3600
        let mut right = g[..3000].to_vec();
        right.extend_from_slice(&rc(&g[1000..3600]));

        let ml = mirror(&left, &trips_for(5400, 1000, 2300, k), k).unwrap();
        let mr = mirror(&right, &trips_for(6600, 1000, 2900, k), k).unwrap();
        let (la, lb) = boundary_kmers(&left, ml).unwrap();
        let (ra, rb) = boundary_kmers(&right, mr).unwrap();

        // "library": both chimeras plus five linear reads of the genome
        let mut counts: HashMap<u64, u32> = HashMap::new();
        let mut add = |s: &[u8]| {
            for w in s.windows(BK) {
                if let Some(x) = canon(w) {
                    *counts.entry(x).or_insert(0) += 1;
                }
            }
        };
        add(&left);
        add(&right);
        for _ in 0..5 {
            add(&g);
        }
        let ev = SpacerEvidence { counts };
        assert_eq!(ev.side(la, lb), SpacerSide::Left);
        assert_eq!(ev.side(ra, rb), SpacerSide::Right);
    }
}
