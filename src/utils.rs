#[derive(Debug, Clone)]
pub struct Refined {
    pub split_pos: usize,
    pub score: i64,
    pub identity_est: f32,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum RefineMode {
    HiFi,
    ONT,
}

pub fn div_floor(x: i32, d: i32) -> i32 {
    let q = x / d;
    let r = x % d;
    if (r != 0) && ((r > 0) != (d > 0)) {
        q - 1
    } else {
        q
    }
}

// Compute banded Levenshtein edit distance between `a` and `b`.
//
// Returns a distance in [0..=band+1]. If the true distance is > band, returns band+1.
// This makes refinement robust: we can still choose the best candidate split even if
// all candidates are "bad" under the current band.
//
// Complexity: O(len(a) * band)
pub fn banded_edit_distance(a: &[u8], b: &[u8], band: usize) -> usize {
    let n = a.len();
    let m = b.len();

    if n == 0 {
        return m.min(band + 1);
    }
    if m == 0 {
        return n.min(band + 1);
    }

    let maxd = band;
    let inf = maxd + 1;

    let mut prev = vec![inf; m + 1];
    let mut curr = vec![inf; m + 1];

    prev[0] = 0;
    for j in 1..=m {
        prev[j] = if j <= maxd { j } else { inf };
    }

    for i in 1..=n {
        let i_is = i as isize;
        let band_is = maxd as isize;

        let j_min = (1isize.max(i_is - band_is)) as usize;
        let j_max = (m as isize).min(i_is + band_is) as usize;

        // reset current row
        for v in curr.iter_mut() {
            *v = inf;
        }
        curr[0] = if i <= maxd { i } else { inf };

        let mut row_best = inf;

        for j in j_min..=j_max {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };

            let del = prev[j].saturating_add(1);
            let ins = curr[j - 1].saturating_add(1);
            let sub = prev[j - 1].saturating_add(cost);

            let v = del.min(ins).min(sub);
            curr[j] = v;
            row_best = row_best.min(v);
        }

        // Early exit: if best in row already exceeds band, we can stop.
        if row_best > maxd {
            return maxd + 1;
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    prev[m].min(maxd + 1)
}

// Write reverse-complement of `src` into `dst`.
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

// Return the DNA complement of a base.
// Non-ACGT bases are mapped to 'N'.
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

// Reverse-complement in place (ACGTN; other -> N)
pub fn revcomp_in_place(seq: &mut [u8]) {
    let mut i = 0usize;
    let mut j = seq.len().saturating_sub(1);
    while i < j {
        let a = comp(seq[i]);
        let b = comp(seq[j]);
        seq[i] = b;
        seq[j] = a;
        i += 1;
        j = j.saturating_sub(1);
    }
    if i == j && i < seq.len() {
        seq[i] = comp(seq[i]);
    }
}
