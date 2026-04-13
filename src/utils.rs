//! Utility functions and small shared types for `mdax`.
//!
//! This module focuses on hot-path helpers (FASTA writing, reverse complements,
//! and small DP utilities) that are reused across pass1/pass2 and refinement code.

use anyhow::Result;
use calm_io::stderrln;

/// Result of breakpoint refinement around a coarse split.
///
/// `Refined` captures the minimal information needed downstream:
/// - `split_pos`: refined breakpoint coordinate in read space
/// - `score`: refinement objective (higher is better)
/// - `identity_est`: quick similarity estimate (used for filtering / confidence)
///
/// This is typically produced by a local search that compares left/right arms
/// (often via banded edit distance) around candidate split positions.
#[derive(Debug, Clone)]
pub struct Refined {
    pub split_pos: usize,
    pub score: i64,
    pub identity_est: f32,
    /// Estimated full template-switch gap (bp): `2 × δ_best` where δ is the
    /// symmetric half-shift that maximises arm identity.  Zero for clean
    /// hairpin folds with no template switching.
    pub gap_est: usize,
}

/// Refinement mode selection.
///
/// This influences *how expensive / how robust* refinement is, typically by
/// changing:
/// - band size / search radius
/// - scoring parameters
/// - or assumptions about error profile
///
/// In `mdax`:
/// - `HiFi`: cheaper + stable enough for signature coherence in pass1/lookup
/// - `ONT`: more tolerant to higher indel/substitution error for real read correction (in ONT reads)
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum RefineMode {
    HiFi,
    ONT,
}

/// Integer floor-division for signed integers.
///
/// This helper returns ⌊x / d⌋ for all sign combinations, matching mathematical floor.
pub fn div_floor(x: i32, d: i32) -> i32 {
    let q = x / d;
    let r = x % d;
    if (r != 0) && ((r > 0) != (d > 0)) {
        q - 1
    } else {
        q
    }
}

/// Compute a banded Levenshtein edit distance between `a` and `b`, reusing scratch buffers.
///
/// Returns a distance in `0..=band+1`.
/// - If the true Levenshtein distance ≤ `band`, returns the exact distance.
/// - If the true distance > `band`, returns `band + 1`.
///
/// Why return `band+1` instead of infinity?
/// - Refinement often compares many candidate splits; "too far" still needs a comparable
///   sentinel value so you can pick the least-bad candidate deterministically.
///
/// Complexity:
/// - Time: `O(len(a) * (2*band + 1))` (clamped by `len(b)`)
/// - Space: `O(len(b))` for two rolling rows
///
/// `prev` and `curr` are scratch vectors to avoid reallocation in tight loops.
/// They are resized as needed and overwritten.
pub fn banded_edit_distance_scratch(
    a: &[u8],
    b: &[u8],
    band: usize,
    prev: &mut Vec<usize>,
    curr: &mut Vec<usize>,
) -> usize {
    let n = a.len();
    let m = b.len();
    let inf = band + 1;

    // Ensure scratch rows are the right length.
    prev.resize(m + 1, inf);
    curr.resize(m + 1, inf);

    // Initialize row 0: distance from empty prefix of `a` to prefix of `b`.
    // Outside the band: clamp to `inf` so it never wins.
    prev[0] = 0;
    for j in 1..=m {
        prev[j] = if j <= band { j } else { inf };
    }

    // Fill curr once before the loop.  Inside the loop we maintain the invariant
    // with a single O(1) guard write per row instead of a full O(m) reset.
    //
    // Invariant: every cell of `curr` that could be read as `curr[j-1]` at the
    // start of the inner loop (i.e. `curr[j_min - 1]`) is either:
    //   (a) `curr[0]`, which is always set explicitly before the inner loop, or
    //   (b) `curr[j_min - 1]` for j_min > 1, which is set to `inf` by the guard.
    // All other cells in [j_min..j_max] are written before they are read because
    // the loop advances left-to-right and `curr[j]` is written at step j before
    // `curr[j]` is read as `curr[j-1]` at step j+1.
    curr.fill(inf);

    for i in 1..=n {
        // Band limits: only cells j where |i - j| <= band are in scope.
        let j_min = i.saturating_sub(band).max(1);
        let j_max = (i + band).min(m);

        // Guard: the insertion path at j = j_min reads curr[j_min - 1].  When
        // j_min > 1 that cell is stale (written two rows ago); reset it to inf.
        // curr[0] is always set explicitly below.
        if j_min > 1 {
            curr[j_min - 1] = inf;
        }
        curr[0] = if i <= band { i } else { inf };

        // Hoist the row-constant character lookup.
        let ai = a[i - 1];
        let mut row_best = inf;

        for j in j_min..=j_max {
            // Values are bounded by band + n (≤ ~625 for our 500 bp window),
            // nowhere near usize overflow, so plain addition is safe.
            let cost = usize::from(ai != b[j - 1]);
            let del = prev[j] + 1;
            let ins = curr[j - 1] + 1;
            let sub = prev[j - 1] + cost;
            let v = del.min(ins).min(sub);
            curr[j] = v;
            row_best = row_best.min(v);
        }

        // Early exit: if the best cell in this row already exceeds the band,
        // the true edit distance must be > band.
        if row_best > band {
            return band + 1;
        }

        std::mem::swap(prev, curr);
    }
    prev[m].min(band + 1)
}

/// Myers' bit-parallel edit distance (Hyyrö 2001 formulation).
///
/// Computes the Levenshtein edit distance between `a` (pattern, length n) and
/// `b` (text, length m), returning the distance clamped to `band + 1`.
///
/// Runs in O(⌈n/64⌉ · m) word operations — roughly 64× fewer inner-loop
/// iterations than the scalar banded DP for typical n = m windows.
///
/// Works for any pattern length; the bitvector state is heap-allocated and
/// grows with n.
///
/// # Equations (per text character b[j])
/// ```text
/// Xv = Eq | Mv
/// X  = Xv & Pv
/// D0 = ((X + Pv) ^ Pv) | Xv    ← + is multi-word carry-propagating addition
/// Ph = Mv | ~(D0 | Pv)
/// Mh = Pv & D0
/// score += bit(n−1, Ph);  score -= bit(n−1, Mh)
/// Ph = (Ph << 1) | 1;  Mh <<= 1   ← left-shift with cross-word carry
/// Pv = Mh | ~(Xv | Ph)
/// Mv = Ph & Xv
/// ```
pub fn edit_distance_myers(a: &[u8], b: &[u8], band: usize) -> usize {
    let n = a.len();
    let m = b.len();
    if n == 0 { return m.min(band + 1); }
    if m == 0 { return n.min(band + 1); }

    let words = (n + 63) / 64;

    // Precompute PM: flat Vec indexed as pm[c * words + k] = word k of the
    // match bitvector for pattern char c.  Only the 256 ASCII byte values
    // are needed; DNA uses 5 of them.
    let mut pm = vec![0u64; 256 * words];
    for (i, &c) in a.iter().enumerate() {
        pm[c as usize * words + (i >> 6)] |= 1u64 << (i & 63);
    }

    // Mask for the last word: zero out bits above position n-1.
    let last_bits = n & 63;
    let last_mask: u64 = if last_bits == 0 { !0 } else { (1u64 << last_bits) - 1 };
    let hi_w = words - 1;       // index of the last active word
    let hi_b = (n - 1) & 63;   // bit index of pattern position n-1 within hi_w

    // DP state: Pv[k]/Mv[k] = positive/negative vertical difference bitvectors.
    // Initial state: all vertical deltas = +1 (column 0 values = 0,1,…,n).
    let mut pv = vec![!0u64; words];
    let mut mv = vec![0u64; words];
    pv[hi_w] = last_mask;

    let mut score = n;

    for (j, &bc) in b.iter().enumerate() {
        let eq = &pm[bc as usize * words .. bc as usize * words + words];

        // Single left-to-right pass per text character.
        // Carries thread state between words:
        //   `add_carry` — carry out of the multi-word addition X + Pv
        //   `ph_carry`  — bit shifted in at word boundary when left-shifting Ph
        //   `mh_carry`  — same for Mh
        let mut add_carry = 0u64;
        let mut ph_carry  = 1u64; // boundary: d[0][j] = j, so horizontal delta at row 0 = +1
        let mut mh_carry  = 0u64;

        for k in 0..words {
            let pv_k = pv[k];
            let mv_k = mv[k];

            // D0: zero-difference positions.
            let xv = eq[k] | mv_k;
            let x  = xv & pv_k;
            let (s1, c1) = x.overflowing_add(pv_k);
            let (s2, c2) = s1.overflowing_add(add_carry);
            let d0 = (s2 ^ pv_k) | xv;
            add_carry = (c1 as u64) | (c2 as u64);

            // Horizontal differences for this word.
            let ph_k = mv_k | !(d0 | pv_k);
            let mh_k = pv_k & d0;

            // Score update from the last pattern bit (only in the final word).
            if k == hi_w {
                score += ((ph_k >> hi_b) & 1) as usize;
                score  = score.saturating_sub(((mh_k >> hi_b) & 1) as usize);
            }

            // Left-shift Ph and Mh by 1 with cross-word carry, then update Pv/Mv.
            let nph = (ph_k << 1) | ph_carry;  ph_carry = ph_k >> 63;
            let nmh = (mh_k << 1) | mh_carry;  mh_carry = mh_k >> 63;
            pv[k] = nmh | !(xv | nph);
            mv[k] = nph & xv;
        }

        // Mask stray bits in the last word (produced by ~).
        pv[hi_w] &= last_mask;
        mv[hi_w] &= last_mask;

        // Early exit: score = d[n][j+1].  Even if score decreases by 1 every
        // remaining column, it cannot drop below score - (m - j - 1).
        // If that floor already exceeds band, the true distance must be > band.
        if score > band + (m - 1 - j) { return band + 1; }
    }

    score.min(band + 1)
}

/// Reverse-complement lookup table for ASCII bases.
/// Unknowns map to 'N'. Preserves case for A/C/G/T/N.
/// Includes common IUPAC ambiguity complements 
#[rustfmt::skip]
pub const RC_LUT: [u8; 256] = {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'T'; t[b'C' as usize] = b'G'; t[b'G' as usize] = b'C'; t[b'T' as usize] = b'A'; t[b'N' as usize] = b'N';
    t[b'a' as usize] = b't'; t[b'c' as usize] = b'g'; t[b'g' as usize] = b'c'; t[b't' as usize] = b'a'; t[b'n' as usize] = b'n';
    // common ambiguity codes
    t[b'R' as usize] = b'Y'; t[b'Y' as usize] = b'R';
    t[b'S' as usize] = b'S'; t[b'W' as usize] = b'W';
    t[b'K' as usize] = b'M'; t[b'M' as usize] = b'K';
    t[b'B' as usize] = b'V'; t[b'V' as usize] = b'B';
    t[b'D' as usize] = b'H'; t[b'H' as usize] = b'D';
    t[b'r' as usize] = b'y'; t[b'y' as usize] = b'r';
    t[b's' as usize] = b's'; t[b'w' as usize] = b'w';
    t[b'k' as usize] = b'm'; t[b'm' as usize] = b'k';
    t[b'b' as usize] = b'v'; t[b'v' as usize] = b'b';
    t[b'd' as usize] = b'h'; t[b'h' as usize] = b'd';
    t
};

/// Reverse-complement in place.
/// Hot path: two-pointer swap + LUT, minimal branching.
#[inline]
pub fn revcomp_in_place(seq: &mut [u8]) {
    let mut i = 0usize;
    let mut j = seq.len();

    // Safety: pointer ops are bounds-checked via i/j invariants.
    while i < j {
        j -= 1;
        if i == j {
            // middle element (odd length)
            let b = unsafe { *seq.get_unchecked(i) };
            unsafe { *seq.get_unchecked_mut(i) = RC_LUT[b as usize] };
            break;
        }

        let a = unsafe { *seq.get_unchecked(i) };
        let b = unsafe { *seq.get_unchecked(j) };
        unsafe {
            *seq.get_unchecked_mut(i) = RC_LUT[b as usize];
            *seq.get_unchecked_mut(j) = RC_LUT[a as usize];
        }

        i += 1;
    }
}

/// Return the DNA complement of a base.
///
/// Mapping:
/// - A to T
/// - C to G
/// - everything else to N
///
/// This is deliberately conservative: any ambiguity or unexpected symbol becomes `N`,
/// which keeps downstream signatures/distances well-defined.
#[inline]
pub fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// Reverse-complement `src` into preallocated `dst`.
/// `dst` will be resized to match `src` length.
/// This avoids allocations if `dst` is reused.
#[inline]
pub fn revcomp_into(dst: &mut Vec<u8>, src: &[u8]) {
    let n = src.len();
    dst.resize(n, 0);

    // Write dst[k] = RC(src[n-1-k])
    // Use raw pointers to keep bounds checks out of the inner loop.
    unsafe {
        let sp = src.as_ptr();
        let dp = dst.as_mut_ptr();
        for k in 0..n {
            let b = *sp.add(n - 1 - k);
            *dp.add(k) = RC_LUT[b as usize];
        }
    }
}

/// For logging the support clusters in `main.rs`
pub fn compact_histogram(
    buckets: &std::collections::HashMap<usize, usize>,
    width: usize,
) -> Result<()> {
    use std::collections::BTreeMap;

    // Aggregate into ranges of 10
    let mut range_buckets: BTreeMap<usize, usize> = BTreeMap::new();
    for (&key, &count) in buckets {
        let range_start = (key / 10) * 10;
        *range_buckets.entry(range_start).or_insert(0) += count;
    }

    let max_count = *range_buckets.values().max().unwrap_or(&1) as f64;
    let ratio = width as f64 / max_count; // <= 1.0 typically

    for (range_start, &count) in &range_buckets {
        let range_end = range_start + 9;
        let bar_width = ((count as f64) * ratio).round() as usize;
        let bar = "█".repeat(bar_width.min(width));
        stderrln!("{:>3}–{:>3}: {:>6} {}", range_start, range_end, count, bar)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // A tiny reference Levenshtein for small strings (no banding, no scratch).
    fn levenshtein(a: &[u8], b: &[u8]) -> usize {
        let n = a.len();
        let m = b.len();
        let mut dp = vec![vec![0usize; m + 1]; n + 1];
        for i in 0..=n {
            dp[i][0] = i;
        }
        for j in 0..=m {
            dp[0][j] = j;
        }
        for i in 1..=n {
            for j in 1..=m {
                let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
                dp[i][j] = (dp[i - 1][j] + 1)
                    .min(dp[i][j - 1] + 1)
                    .min(dp[i - 1][j - 1] + cost);
            }
        }
        dp[n][m]
    }

    #[test]
    fn div_floor_matches_mathematical_floor() {
        // Positive / positive
        assert_eq!(div_floor(7, 3), 2);
        assert_eq!(div_floor(6, 3), 2);

        // Negative numerator
        assert_eq!(div_floor(-7, 3), -3); // -2.333 -> -3
        assert_eq!(div_floor(-6, 3), -2);

        // Negative denominator
        assert_eq!(div_floor(7, -3), -3); // -2.333 -> -3
        assert_eq!(div_floor(6, -3), -2);

        // Both negative
        assert_eq!(div_floor(-7, -3), 2); // 2.333 -> 2
        assert_eq!(div_floor(-6, -3), 2);

        // Small edges around zero
        assert_eq!(div_floor(-1, 2), -1);
        assert_eq!(div_floor(1, -2), -1);
    }

    #[test]
    fn comp_is_conservative_and_case_insensitive_for_acgt() {
        assert_eq!(comp(b'A'), b'T');
        assert_eq!(comp(b'a'), b'T');
        assert_eq!(comp(b'C'), b'G');
        assert_eq!(comp(b'c'), b'G');
        assert_eq!(comp(b'G'), b'C');
        assert_eq!(comp(b'g'), b'C');
        assert_eq!(comp(b'T'), b'A');
        assert_eq!(comp(b't'), b'A');

        // Ambiguity becomes N
        assert_eq!(comp(b'R'), b'N');
        assert_eq!(comp(b'N'), b'N');
        assert_eq!(comp(b'-'), b'N');
    }

    #[test]
    fn revcomp_in_place_even_and_odd_lengths() {
        let mut s1 = b"ACGT".to_vec();
        revcomp_in_place(&mut s1);
        assert_eq!(s1, b"ACGT"); // ACGT is self-revcomp

        let mut s2 = b"AAA".to_vec();
        revcomp_in_place(&mut s2);
        assert_eq!(s2, b"TTT");

        let mut s3 = b"ACGTA".to_vec();
        revcomp_in_place(&mut s3);
        assert_eq!(s3, b"TACGT");
    }

    #[test]
    fn revcomp_into_matches_in_place() {
        let src = b"AcgTNryKMbvDh".to_vec(); // includes mixed case + ambiguity
        let mut dst = Vec::new();

        revcomp_into(&mut dst, &src);

        let mut in_place = src.clone();
        revcomp_in_place(&mut in_place);

        assert_eq!(dst, in_place);
    }

    #[test]
    fn banded_edit_distance_returns_exact_when_within_band() {
        let a = b"ACGTACGTACGT";
        let b = b"ACGTTCGTACGT"; // 1 substitution
        let true_ed = levenshtein(a, b);
        assert_eq!(true_ed, 1);

        let mut prev = Vec::new();
        let mut curr = Vec::new();
        let ed = banded_edit_distance_scratch(a, b, 2, &mut prev, &mut curr);
        assert_eq!(ed, 1);
    }

    #[test]
    fn banded_edit_distance_returns_band_plus_one_when_outside_band() {
        let a = b"AAAAAAAAAAAA";
        let b = b"TTTTTTTTTTTT"; // 12 substitutions
        let true_ed = levenshtein(a, b);
        assert_eq!(true_ed, 12);

        let mut prev = Vec::new();
        let mut curr = Vec::new();
        let band = 3;
        let ed = banded_edit_distance_scratch(a, b, band, &mut prev, &mut curr);
        assert_eq!(ed, band + 1);
    }

    #[test]
    fn banded_edit_distance_symmetric_sanity_with_large_band() {
        // With band >= max(len), it should behave like full Levenshtein.
        let a = b"GATTACA";
        let b = b"GCATGCU";

        let true_ed = levenshtein(a, b);

        let mut prev = Vec::new();
        let mut curr = Vec::new();
        let band = a.len().max(b.len());

        let ed_ab = banded_edit_distance_scratch(a, b, band, &mut prev, &mut curr);

        let mut prev2 = Vec::new();
        let mut curr2 = Vec::new();
        let ed_ba = banded_edit_distance_scratch(b, a, band, &mut prev2, &mut curr2);

        assert_eq!(ed_ab, true_ed);
        assert_eq!(ed_ba, true_ed);
    }

    // ── Myers' bit-parallel edit distance ────────────────────────────────────

    /// Exhaustively compare Myers' against the reference Levenshtein for all
    /// pairs of short DNA strings.
    #[test]
    fn myers_matches_levenshtein_exhaustive_short() {
        let alphabet = b"ACGT";
        let seqs: Vec<Vec<u8>> = (0..=8).flat_map(|len| {
            // All strings over ACGT of each length (4^len, capped at len ≤ 4 for speed)
            // For len > 4 just use a few hand-picked strings.
            if len <= 4 {
                let count = 4usize.pow(len as u32);
                (0..count).map(|mut idx| {
                    (0..len).map(|_| { let c = alphabet[idx % 4]; idx /= 4; c }).collect()
                }).collect()
            } else {
                vec![
                    b"ACGTACGT"[..len].to_vec(),
                    b"TTTTTTTT"[..len].to_vec(),
                    b"AAAAAAAA"[..len].to_vec(),
                ]
            }
        }).collect();

        for a in &seqs {
            for b in &seqs {
                let band = a.len().max(b.len());
                let expected = levenshtein(a, b);
                let got = edit_distance_myers(a, b, band);
                assert_eq!(
                    got, expected,
                    "Myers mismatch: a={:?} b={:?} expected={expected} got={got}",
                    std::str::from_utf8(a).unwrap(),
                    std::str::from_utf8(b).unwrap(),
                );
            }
        }
    }

    #[test]
    fn myers_band_exceeded_returns_band_plus_one() {
        let a = b"AAAAAAAAAAAAAAAA";
        let b = b"TTTTTTTTTTTTTTTT";
        let band = 3;
        assert_eq!(edit_distance_myers(a, b, band), band + 1);
    }

    #[test]
    fn myers_exact_match() {
        let s = b"ACGTACGTACGTACGT";
        assert_eq!(edit_distance_myers(s, s, s.len()), 0);
    }

    /// Multi-word path: n > 64. Compare Myers' against banded ED on a 200 bp window.
    #[test]
    fn myers_multiword_matches_banded_ed() {
        // Build a ~97% identity pair, 200 bp each.
        let a: Vec<u8> = (0..200usize).map(|i| b"ACGT"[i % 4]).collect();
        let mut b = a.clone();
        // Introduce 6 substitutions spread through the sequence.
        for &pos in &[10usize, 40, 70, 100, 140, 180] {
            b[pos] = if b[pos] == b'A' { b'T' } else { b'A' };
        }
        let band = 20;
        let mut prev = Vec::new();
        let mut curr = Vec::new();
        let banded = banded_edit_distance_scratch(&a, &b, band, &mut prev, &mut curr);
        let myers  = edit_distance_myers(&a, &b, band);
        assert_eq!(myers, banded, "multiword: myers={myers} banded={banded}");
    }

    /// Longer sequence (n = 500 bp).
    #[test]
    fn myers_500bp_matches_banded_ed() {
        let a: Vec<u8> = (0..500usize).map(|i| b"ACGT"[i % 4]).collect();
        let mut b = a.clone();
        for &pos in &[5usize, 50, 100, 200, 300, 400, 490] {
            b[pos] = if b[pos] == b'A' { b'T' } else { b'A' };
        }
        let band = 50;
        let mut prev = Vec::new();
        let mut curr = Vec::new();
        let banded = banded_edit_distance_scratch(&a, &b, band, &mut prev, &mut curr);
        let myers  = edit_distance_myers(&a, &b, band);
        assert_eq!(myers, banded, "500bp: myers={myers} banded={banded}");
    }

    /// n > 512: exercises the heap-allocated multi-word path beyond the old fixed limit.
    #[test]
    fn myers_large_n_matches_banded_ed() {
        let a: Vec<u8> = (0..1200usize).map(|i| b"ACGT"[i % 4]).collect();
        let mut b = a.clone();
        for &pos in &[10usize, 100, 300, 600, 900, 1100] {
            b[pos] = if b[pos] == b'A' { b'T' } else { b'A' };
        }
        let band = 50;
        let mut prev = Vec::new();
        let mut curr = Vec::new();
        let banded = banded_edit_distance_scratch(&a, &b, band, &mut prev, &mut curr);
        let myers  = edit_distance_myers(&a, &b, band);
        assert_eq!(myers, banded, "1200bp: myers={myers} banded={banded}");
    }
}
