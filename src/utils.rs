//! Utility functions and small shared types for `mdax`.
//!
//! This module focuses on hot-path helpers (FASTA writing, reverse complements,
//! and small DP utilities) that are reused across pass1/pass2 and refinement code.

use anyhow::Result;
use calm_io::stderrln;

/// Write one FASTA record to `f` as a single-line sequence.
///
/// This is intentionally minimal and allocation-free:
/// - writes `>` + `id` + newline
/// - writes `seq` on **one line** (fastest / simplest; avoids line wrapping)
/// - writes trailing newline
///
/// Notes:
/// - `seq` is assumed to already be ASCII bases (A/C/G/T/N) but we don't enforce it.
/// - We use `impl Write` to avoid dynamic dispatch in hot code.
pub fn write_fasta(f: &mut impl std::io::Write, id: &str, seq: &[u8]) -> std::io::Result<()> {
    f.write_all(b">")?;
    f.write_all(id.as_bytes())?;
    f.write_all(b"\n")?;
    // write sequence as a single line (fastest)
    f.write_all(seq)?;
    f.write_all(b"\n")?;
    Ok(())
}

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
/// But note that `ONT` mode is slower (~1.5x slower?)
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
/// - Time: `O(len(a) * band)` (only a diagonal band of the DP matrix is evaluated)
/// - Space: `O(len(b))` via two rolling rows (`prev`, `curr`)
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

/// Write the reverse-complement of `src` into `dst`.
///
/// This is a low-level, allocation-free primitive for situations where you want to
/// reuse an output buffer (e.g. inside inner loops).
///
/// Requirements:
/// - `src.len() == dst.len()` (enforced by debug_assert)
/// - `dst` may alias `src`? (Not safe for in-place use; use `revcomp_in_place` instead.)
///
/// Non-ACGT bases are mapped to `N` via `comp()`.
#[inline]
pub fn revcomp_into(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    let mut i = 0usize;
    let mut j = src.len();
    while i < src.len() {
        j -= 1;
        dst[i] = comp(src[j]);
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

/// Reverse-complement a sequence **in place**.
///
/// This swaps and complements from both ends toward the middle:
/// - O(n) time
/// - O(1) extra memory
///
/// Behavior:
/// - A/C/G/T (case-insensitive) are complemented as expected.
/// - Any other byte is converted to `N` (including already-ambiguous bases).
///
/// Edge cases:
/// - Empty slice: no-op.
/// - Odd length: the middle base is complemented once after the loop.
pub fn revcomp_in_place(seq: &mut [u8]) {
    let mut i = 0usize;
    let mut j = seq.len().saturating_sub(1);
    while i < j {
        // Complement first, then swap, so we don't need a temp buffer of original bases.
        let a = comp(seq[i]);
        let b = comp(seq[j]);
        seq[i] = b;
        seq[j] = a;
        i += 1;
        j = j.saturating_sub(1);
    }

    // If length is odd, complement the middle element.
    if i == j && i < seq.len() {
        seq[i] = comp(seq[i]);
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
