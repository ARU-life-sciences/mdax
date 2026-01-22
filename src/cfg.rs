use crate::utils::RefineMode;

/// Global configuration for mdax: shared knobs + per-detector knobs.
#[derive(Debug, Clone)]
pub struct MdaxCfg {
    pub shared: SharedCfg,
    pub fold: FoldOnlyCfg,
    pub concat: ConcatOnlyCfg,
}

impl Default for MdaxCfg {
    fn default() -> Self {
        Self {
            shared: SharedCfg::default(),
            fold: FoldOnlyCfg::default(),
            concat: ConcatOnlyCfg::default(),
        }
    }
}

/// Minimizer sampling configuration.
#[derive(Debug, Clone)]
pub struct MinimizerCfg {
    pub k: usize,
    pub w: usize,

    /// If true, use strand-specific forward-only minimizers (NtHasher<false>).
    /// If false, use canonical minimizers (strand-agnostic).
    ///
    /// Recommendation:
    /// - foldback: true
    /// - concatemer: true (for consistency with foldback)
    pub forward_only: bool,
}

impl Default for MinimizerCfg {
    fn default() -> Self {
        Self {
            k: 17,
            w: 21,
            forward_only: true,
        }
    }
}

/// Refinement configuration shared by both detectors.
#[derive(Debug, Clone)]
pub struct RefineCfg {
    pub window: usize,
    pub arm: usize,
    pub mode: RefineMode,
    pub max_ed_rate: f32,
}

impl Default for RefineCfg {
    fn default() -> Self {
        Self {
            window: 200,
            arm: 500,
            mode: RefineMode::HiFi,
            max_ed_rate: 0.25,
        }
    }
}

/// Shared knobs across foldback + concatemer detection.
#[derive(Debug, Clone)]
pub struct SharedCfg {
    pub minimizer: MinimizerCfg,

    /// Minimum number of minimizer matchpoints required for a candidate.
    pub min_matches: usize,

    /// Don’t call breakpoints too close to ends.
    pub end_guard: usize,

    /// Refinement parameters.
    pub refine: RefineCfg,

    /// Separate diagonal tolerances because foldback uses d=(p1+p2)
    /// while concat uses delta=(p2-p1).
    pub fold_diag_tol: i32,
    pub concat_diag_tol: i32,
}

impl Default for SharedCfg {
    fn default() -> Self {
        Self {
            minimizer: MinimizerCfg::default(),
            min_matches: 10,
            end_guard: 1000,
            refine: RefineCfg::default(),
            fold_diag_tol: 120,
            concat_diag_tol: 200,
        }
    }
}

/// Foldback-only gates (everything else is in SharedCfg).
#[derive(Debug, Clone)]
pub struct FoldOnlyCfg {
    /// Minimum “arm evidence span” (bp) on p1 for anti-diagonal chain.
    pub min_arm: usize,
}

impl Default for FoldOnlyCfg {
    fn default() -> Self {
        Self { min_arm: 200 }
    }
}

/// Concatemer-only gates.
#[derive(Debug, Clone)]
pub struct ConcatOnlyCfg {
    /// Minimum matched span on p1 to call a repeat.
    pub min_span: usize,

    /// Ignore short/self-local repeats (tandem repeats etc.)
    pub min_delta: usize,

    /// Require that most matchpoints cross the inferred junction.
    pub cross_frac: f32,
}

impl Default for ConcatOnlyCfg {
    fn default() -> Self {
        Self {
            min_span: 2000,
            min_delta: 2000,
            cross_frac: 0.80,
        }
    }
}
