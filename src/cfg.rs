// configuration of mdax
use crate::utils::RefineMode;

// Global configuration for mdax: shared + per-detector.
#[derive(Debug, Clone)]
pub struct MdaxCfg {
    pub shared: SharedCfg,
    pub fold: FoldOnlyCfg,
    pub fold2: FoldSecondPassCfg,
    pub sig: SigCfg,
}

impl Default for MdaxCfg {
    fn default() -> Self {
        Self {
            shared: SharedCfg::default(),
            fold: FoldOnlyCfg::default(),
            fold2: FoldSecondPassCfg::default(),
            sig: SigCfg::default(),
        }
    }
}

// Minimizer sampling configuration.
#[derive(Debug, Clone)]
pub struct MinimizerCfg {
    pub k: usize,
    pub w: usize,

    // If true, use strand-specific forward-only minimizers (NtHasher<false>).
    // If false, use canonical minimizers (strand-agnostic).
    // Should generally use forward_only=true for mdax.
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

// Refinement configuration shared by both detectors.
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

// Shared across foldback + concatemer detection.
#[derive(Debug, Clone)]
pub struct SharedCfg {
    pub minimizer: MinimizerCfg,

    // Minimum number of minimizer matchpoints required for a candidate.
    pub min_matches: usize,

    // Don’t call breakpoints too close to ends.
    pub end_guard: usize,

    // Refinement parameters.
    pub refine: RefineCfg,

    // Separate diagonal tolerances because foldback uses d=(p1+p2)
    // while concat uses delta=(p2-p1).
    pub fold_diag_tol: i32,
}

impl Default for SharedCfg {
    fn default() -> Self {
        Self {
            minimizer: MinimizerCfg::default(),
            min_matches: 10,
            end_guard: 1000,
            refine: RefineCfg::default(),
            fold_diag_tol: 120,
        }
    }
}

// Foldback-only gates.
#[derive(Debug, Clone)]
pub struct FoldOnlyCfg {
    // Minimum “arm evidence span” (bp) on p1 for anti-diagonal chain.
    pub min_arm: usize,
}

impl Default for FoldOnlyCfg {
    fn default() -> Self {
        Self { min_arm: 200 }
    }
}

#[derive(Debug, Clone)]
pub struct FoldSecondPassCfg {
    // Minimum number of reads supporting the same junction fingerprint
    // to call it "real" and therefore DO NOT correct by default.
    pub min_support: usize,

    // If supported, splits must cluster within this tolerance (bp) to be “real”.
    pub split_tol_bp: usize,

    // Minimum arm identity (refinement identity_est) to accept as a foldback candidate.
    pub min_identity: f32,
}

impl Default for FoldSecondPassCfg {
    fn default() -> Self {
        Self {
            min_support: 3,
            split_tol_bp: 100,
            min_identity: 0.60,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CallMode {
    Strict,
    Balanced,
    Permissive,
}

impl Default for CallMode {
    fn default() -> Self {
        CallMode::Balanced
    }
}

// this will be exposed on the cli, as defualts
#[derive(Debug, Clone)]
pub struct ModeTuning {
    pub min_span: usize,     // fold.min_arm + concat.min_span
    pub min_matches: usize,  // coarse detector requirement
    pub min_identity: f32,   // refinement identity gate
    pub min_support: usize,  // "real foldback" support threshold
    pub split_tol_bp: usize, // support cluster tightness
    pub max_depth: usize,    // recursion depth
    pub sig_flank_bp: usize,
    pub sig_take: usize,
    pub sig_value_shift: u8,
}

impl CallMode {
    pub fn tuning(self) -> ModeTuning {
        match self {
            // strict: fewer calls, fewer cuts, higher confidence
            CallMode::Strict => ModeTuning {
                min_matches: 30,
                min_span: 3000,
                min_identity: 0.75,
                min_support: 5,
                split_tol_bp: 30,
                max_depth: 2,
                sig_flank_bp: 800,
                sig_take: 10,
                sig_value_shift: 0,
            },
            // balanced: your current-ish defaults
            CallMode::Balanced => ModeTuning {
                min_matches: 20,
                min_span: 2000,
                min_identity: 0.60,
                min_support: 3,
                split_tol_bp: 50,
                max_depth: 5,
                sig_flank_bp: 600,
                sig_take: 8,
                sig_value_shift: 4,
            },
            // permissive: more calls/cuts (good for exploration)
            CallMode::Permissive => ModeTuning {
                min_matches: 12,
                min_span: 1200,
                min_identity: 0.50,
                min_support: 2,
                split_tol_bp: 100,
                max_depth: 8,
                sig_flank_bp: 400,
                sig_take: 6,
                sig_value_shift: 8,
            },
        }
    }
}

// for comparison to the perl pipeline from
// https://academic.oup.com/nar/article/51/15/8035/7234520#414928225
// these params approximate what that pipeline uses
// it's for testing
#[derive(Debug, Clone)]
pub struct FairnessParams {
    // detection / decision
    pub min_matches: usize,
    pub min_span: usize,
    pub min_identity: f32,
    pub min_support: usize,
    pub split_tol_bp: usize,
    pub max_depth: usize,

    // minimizers
    pub k: usize,
    pub w: usize,
    pub forward_only: bool,

    // refinement
    pub refine_mode: RefineMode,
    pub refine_window: usize,
    pub refine_arm: usize,
    pub max_ed_rate: f32,

    // misc structural guards
    pub end_guard: usize,
    pub fold_diag_tol: i32,

    pub sig_flank_bp: usize,
    pub sig_take: usize,
}

impl FairnessParams {
    pub fn baseline() -> Self {
        Self {
            // close to perl/minimap2 behaviour
            min_matches: 20,
            min_span: 1000,
            min_identity: 0.60,
            min_support: 2,
            split_tol_bp: 100,
            max_depth: 1,

            // minimizers (LOCKED)
            k: 17,
            w: 21,
            forward_only: true,

            // refinement (LOCKED)
            refine_mode: RefineMode::HiFi,
            refine_window: 200,
            refine_arm: 500,
            max_ed_rate: 0.25,

            // guards (LOCKED)
            end_guard: 1000,
            fold_diag_tol: 120,

            sig_flank_bp: 1000,
            sig_take: 12,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SigCfg {
    pub flank_bp: usize,
    pub take: usize,
    pub value_shift: u8,
}

impl Default for SigCfg {
    fn default() -> Self {
        Self {
            flank_bp: 1000,
            take: 12,
            value_shift: 0, // the default is exact
        }
    }
}
