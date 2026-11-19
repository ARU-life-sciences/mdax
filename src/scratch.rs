//! Scratch buffers for foldback detection, refinement, and junction keying.
//!
//! These structs exist purely to reuse allocations across many reads and avoid
//! repeated heap churn in the hot path.
//!
//! Invariants / usage patterns:
//! - Buffers are reused and are not required to be empty on entry.
//! - Callers should treat contents as ephemeral: each algorithm overwrites or
//!   clears what it needs.
use gxhash::HashMapExt;

use crate::foldback::BinStat;

/// Scratch space for `detect_foldback` + subsequent operations that reuse the
/// matched minimizers.
///
/// This holds:
/// - the reverse-complemented sequence (`rc`)
/// - minimizer outputs for the forward sequence and its reverse complement
/// - hash-map indices used during minimizer matching and anti-diagonal binning
/// - the winning anti-diagonal's matchpoints, with and without their
///   minimizer values (the latter keys the junction)
///
/// Coordinate conventions:
/// - `pos_f` and `pos_rc` are minimizer *start positions* relative to the input
///   passed to `detect_foldback`.
/// - During matching, rc positions are mapped into forward coordinates via
///   `p2 = n - k - prc` (where `n = seq.len()`).
/// - `best_pts` stores `(p1, p2)` pairs in forward coordinates.
/// - `best_trip` stores `(p1, p2, value)` triples in forward coordinates.
pub struct FoldScratch {
    /// Reverse complement of the most recent sequence passed to `detect_foldback`.
    pub rc: Vec<u8>,

    // --- Minimizer buffers ---
    /// Minimizer start positions on the forward sequence.
    pub pos_f: Vec<u32>,
    /// Minimizer values on the forward sequence.
    pub val_f: Vec<u64>,
    /// Minimizer start positions on the reverse-complement sequence.
    pub pos_rc: Vec<u32>,
    /// Minimizer values on the reverse-complement sequence.
    pub val_rc: Vec<u64>,

    // --- Matching / binning state ---
    /// Index of forward minimizer value -> list of forward positions where it occurs.
    ///
    /// Used to join rc minimizers to forward minimizers by value.
    pub idx_f: gxhash::HashMap<u64, Vec<i32>>,

    /// Set of minimizer values marked as repetitive and therefore ignored.
    ///
    /// This is populated when a minimizer value exceeds a bucket cap in `idx_f`,
    /// to prevent quadratic blowups on low-complexity sequence.
    pub repetitive: gxhash::HashMap<u64, ()>,

    /// Per-bin statistics keyed by anti-diagonal bin (`(p1+p2)/diag_tol`).
    pub stats: gxhash::HashMap<i32, BinStat>,

    /// Matchpoints `(p1, p2)` for the selected best anti-diagonal bin.
    pub best_pts: Vec<(i32, i32)>,

    /// Matchpoints of the best bin with their minimizer value: `(p1, p2, value)`.
    ///
    /// Unlike `best_pts` (re-sorted by the chain builder), this is left
    /// untouched after `detect_foldback`, so junction keys can be derived
    /// from it at any later point. See `fingerprint::junction_keys`.
    pub best_trip: Vec<(i32, i32, u64)>,
    /// Scratch for `fingerprint::junction_keys`: `(arm separation, value)`.
    pub jk_tmp: Vec<(i32, u64)>,
    /// Output of `fingerprint::junction_keys` for the most recent call.
    pub jk_out: Vec<(u64, u64)>,
    /// Terminal dense scan: forward k-mers of the window as `(code, pos)`.
    pub term_fwd: Vec<(u32, u32)>,
    /// Terminal dense scan: inverted seed pairs as `(p + q, p, q)`.
    pub term_pairs: Vec<(i32, i32, i32)>,

    /// Scratch for refinement stage.
    pub refine: RefineScratch,

    // --- Multi-hit (top-k) scratch ---
    /// Candidate bins for top-k selection: (bin, rank_key, count, span)
    pub top_bins: Vec<(i32, i64, usize, usize)>,

    /// Map bin -> index in `top_bins` (for point collection).
    pub bin_to_idx: gxhash::HashMap<i32, usize>,

    /// Per-selected-bin matchpoints.
    pub pts_per: Vec<Vec<(i32, i32)>>,

    /// Per-selected-bin local matches (pos,value).
    pub matches_per: Vec<Vec<(usize, u64)>>,
}

impl FoldScratch {
    /// Construct a fresh scratch object with empty buffers.
    ///
    /// Note: buffers will grow to match the largest read processed; we reuse the
    /// same `FoldScratch` across reads to amortize allocations.
    pub fn new() -> Self {
        Self {
            rc: Vec::new(),
            pos_f: Vec::new(),
            val_f: Vec::new(),
            pos_rc: Vec::new(),
            val_rc: Vec::new(),
            idx_f: gxhash::HashMap::new(),
            repetitive: gxhash::HashMap::new(),
            stats: gxhash::HashMap::new(),
            best_pts: Vec::new(),
            best_trip: Vec::new(),
            jk_tmp: Vec::new(),
            jk_out: Vec::new(),
            term_fwd: Vec::new(),
            term_pairs: Vec::new(),
            refine: RefineScratch::default(),
            top_bins: Vec::new(),
            bin_to_idx: gxhash::HashMap::new(),
            pts_per: Vec::new(),
            matches_per: Vec::new(),
        }
    }
}

impl Default for FoldScratch {
    fn default() -> Self {
        Self::new()
    }
}

/// Scratch space for breakpoint refinement.
///
/// Holds:
/// - `right_rc`: reverse complement of the right arm
/// - `left_rc`: reverse complement of the left arm, maintained incrementally
///   across the split-position scan in the Hamming refiner so the inner
///   comparison loop uses forward indexing (compiler-autovectorisable).
#[derive(Debug, Default)]
pub struct RefineScratch {
    pub right_rc: Vec<u8>,
    pub left_rc: Vec<u8>,
}
