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

    // Ensure scratch rows are the right length and default-filled.
    prev.resize(m + 1, inf);
    curr.resize(m + 1, inf);

    // Initialize row 0: distance from empty prefix of `a` to prefix of `b`.
    // Outside the band: clamp to `inf` so it never wins.
    prev[0] = 0;
    for j in 1..=m {
        prev[j] = if j <= band { j } else { inf };
    }

    for i in 1..=n {
        // reset curr row to inf
        for v in curr.iter_mut() {
            *v = inf;
        }
        curr[0] = if i <= band { i } else { inf };

        // Compute the band limits for this row:
        // only j where |i - j| <= band.
        let i_is = i as isize;
        let band_is = band as isize;
        let j_min = (1isize.max(i_is - band_is)) as usize;
        let j_max = (m as isize).min(i_is + band_is) as usize;

        // Track the best value in the row so we can early-exit if the whole row exceeds `band`.
        let mut row_best = inf;

        for j in j_min..=j_max {
            // Levenshtein cost: substitution is 0 if equal else 1.
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };

            // Standard DP transitions:
            // - deletion: drop a[i-1]
            // - insertion: insert b[j-1]
            // - substitution/match: align a[i-1] with b[j-1]
            //
            // Use saturating_add because values may be `inf` and we want to avoid overflow.

            let del = prev[j].saturating_add(1);
            let ins = curr[j - 1].saturating_add(1);
            let sub = prev[j - 1].saturating_add(cost);

            let v = del.min(ins).min(sub);
            curr[j] = v;
            row_best = row_best.min(v);
        }

        // Early exit: if the best cell in this row is already > band, true distance must be > band.
        if row_best > band {
            return band + 1;
        }

        // Roll the DP rows.
        std::mem::swap(prev, curr);
    }
    // Return final distance, clamped at band+1.
    prev[m].min(band + 1)
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
}
