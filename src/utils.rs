pub fn write_fasta(f: &mut impl std::io::Write, id: &str, seq: &[u8]) -> std::io::Result<()> {
    f.write_all(b">")?;
    f.write_all(id.as_bytes())?;
    f.write_all(b"\n")?;
    // write sequence as a single line (fastest)
    f.write_all(seq)?;
    f.write_all(b"\n")?;
    Ok(())
}

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

    prev.resize(m + 1, inf);
    curr.resize(m + 1, inf);

    // init prev row
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

        let i_is = i as isize;
        let band_is = band as isize;
        let j_min = (1isize.max(i_is - band_is)) as usize;
        let j_max = (m as isize).min(i_is + band_is) as usize;

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

        if row_best > band {
            return band + 1;
        }

        std::mem::swap(prev, curr);
    }

    prev[m].min(band + 1)
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
