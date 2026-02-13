//! Scratch buffers for foldback detection, refinement, and signature construction.
//!
//! These structs exist purely to reuse allocations across many reads and avoid
//! repeated heap churn in the hot path.
//!
//! Invariants / usage patterns:
//! - Buffers are reused and are not required to be empty on entry.
//! - Callers should treat contents as ephemeral: each algorithm overwrites or
//!   clears what it needs.
//! - Some functions destructively filter buffers (e.g. `best_matches`), which is
//!   safe as long as the next call repopulates them.
use gxhash::HashMapExt;

use crate::foldback::BinStat;

/// Scratch space for `detect_foldback` + subsequent operations that reuse the
/// matched minimizers.
///
/// This holds:
/// - the reverse-complemented sequence (`rc`)
/// - minimizer outputs for the forward sequence and its reverse complement
/// - hash-map indices used during minimizer matching and anti-diagonal binning
/// - the winning anti-diagonal's matchpoints and a local match list used for
///   constructing a junction signature
///
/// Coordinate conventions:
/// - `pos_f` and `pos_rc` are minimizer *start positions* relative to the input
///   passed to `detect_foldback`.
/// - During matching, rc positions are mapped into forward coordinates via
///   `p2 = n - k - prc` (where `n = seq.len()`).
/// - `best_pts` stores `(p1, p2)` pairs in forward coordinates.
/// - `best_matches` stores `(pos, value)` pairs in forward coordinates, intended
///   for locality filtering around a refined split position.
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

    /// Local match list for signature construction: `(pos, minimizer_value)`.
    ///
    /// Important:
    /// - `pos` is in forward coordinates (relative to the same sequence view as
    ///   the refined split).
    /// - This buffer may be **destructively filtered** by signature functions
    ///   (e.g. `foldback_signature_from_local_matches`).
    pub best_matches: Vec<(usize, u64)>,

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
            best_matches: Vec::new(),
            refine: RefineScratch::default(),
            top_bins: Vec::new(),
            bin_to_idx: gxhash::HashMap::new(),
            pts_per: Vec::new(),
            matches_per: Vec::new(),
        }
    }
}

/// Scratch space for foldback signature computation.
///
/// This is used by both flank-based signatures and local-match signatures.
/// Buffers are reused across calls to avoid allocating minimizer output vectors
/// and intermediate combined-value storage.
///
/// Coordinate conventions:
/// - `pos_l/val_l` and `pos_r/val_r` are minimizer outputs relative to the
///   left/right flanks extracted for signature computation.
#[derive(Debug, Default)]
pub struct SigScratch {
    /// Reverse-complement buffer for the right flank (or other right-hand context).
    pub right_rc: Vec<u8>,

    /// Minimizer start positions on left flank.
    pub pos_l: Vec<u32>,
    /// Minimizer values on left flank.
    pub val_l: Vec<u64>,
    /// Minimizer start positions on right flank (usually on reverse-complemented right).
    pub pos_r: Vec<u32>,
    /// Minimizer values on right flank.
    pub val_r: Vec<u64>,

    /// Combined minimizer values used as the final hash input.
    pub combined: Vec<u64>,
}

/// Scratch space for breakpoint refinement.
///
/// Holds:
/// - `right_rc`: reverse complement of the right arm
/// - `prev`/`curr`: dynamic programming rows for banded edit distance
///
/// The DP buffers are sized to the refinement arm length and reused per call.
#[derive(Debug, Default)]
pub struct RefineScratch {
    pub right_rc: Vec<u8>,
    pub prev: Vec<usize>,
    pub curr: Vec<usize>,
}
