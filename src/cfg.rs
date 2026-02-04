//! Configuration of `mdax`
use crate::utils::RefineMode;

/// Global configuration for `mdax`.
///
/// This is the top-level config object that is constructed (mostly) from CLI args.
/// It is composed of:
/// - `shared`: parameters used by *both* coarse detection and refinement steps
/// - `fold`: foldback-specific detection gates
/// - `fold2`: second-pass support-aware decision thresholds
/// - `sig`: signature/fingerprint parameters used to key junction support
#[derive(Debug, Clone)]
pub struct MdaxCfg {
    /// Parameters shared across detection modes (minimizers, guards, refinement).
    pub shared: SharedCfg,
    /// Foldback-only gates for coarse detection.
    pub fold: FoldOnlyCfg,
    /// Second-pass decision thresholds for “real vs artefact” foldbacks.
    pub fold2: FoldSecondPassCfg,
    /// Signature configuration used for support map keys.
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

/// Minimizer sampling configuration.
///
/// Minimizers are used as a cheap sketch of the sequence to find candidate self-matches
/// (e.g. foldbacks) without doing full alignment.
///
/// Tuning notes:
/// - Larger `k` reduces random matches but may miss noisy signals.
/// - Larger `w` samples fewer minimizers (faster, but lower sensitivity).
/// - `forward_only` affects strand handling and signature stability.
#[derive(Debug, Clone)]
pub struct MinimizerCfg {
    /// Minimizer k-mer length.
    pub k: usize,
    /// Minimizer window size (number of successive k-mers considered per window).
    pub w: usize,

    /// Whether minimizers are strand-specific (forward-only) or canonical.
    ///
    /// - `true`: use forward-only minimizers (strand-specific).
    ///   This is typically preferred for `mdax` detection because it keeps a consistent
    ///   orientation model and tends to reduce spurious symmetry.
    /// - `false`: use canonical minimizers (strand-agnostic); a k-mer and its reverse
    ///   complement map to the same value, which can be more stable for *signatures*.
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

/// Refinement configuration shared across detectors.
///
/// Refinement is the “expensive” local step run after coarse detection to determine
/// a more precise breakpoint (`split_pos`) and a rough identity estimate.
///
/// Intuition:
/// - `window`: how far around the coarse split we search (bp).
/// - `arm`: how much sequence from each side we compare (bp).
/// - `mode`: selects error-profile assumptions / banding.
/// - `max_ed_rate`: guardrail to reject extremely noisy/discordant candidates.
#[derive(Debug, Clone)]
pub struct RefineCfg {
    /// Half-width of the local search window around the coarse split (bp).
    pub window: usize,
    /// Length of sequence extracted from each arm during comparison (bp).
    pub arm: usize,
    /// Refinement mode (e.g. HiFi vs ONT error model / cost settings).
    pub mode: RefineMode,
    /// Maximum tolerated edit-distance rate during refinement (fraction of bases).
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

/// Configuration shared across foldback + concatemer detection.
///
/// This holds the common “coarse detector” sampling and structural safeguards.
/// Foldback-specific thresholds live in `FoldOnlyCfg`, and pass2 support thresholds
/// live in `FoldSecondPassCfg`.
#[derive(Debug, Clone)]
pub struct SharedCfg {
    /// Minimizer sketch parameters used for coarse detection.
    pub minimizer: MinimizerCfg,

    /// Minimum number of minimizer matchpoints required to consider a candidate event.
    ///
    /// Higher values reduce false positives but can miss short/noisy events.
    pub min_matches: usize,

    /// Guard region at both ends of reads (bp).
    ///
    /// Events too close to read ends are often artefactual or poorly supported, and
    /// refinement becomes unstable due to missing context.
    pub end_guard: usize,

    /// Refinement parameters (local search around the coarse split).
    pub refine: RefineCfg,

    /// Diagonal tolerance for foldback clustering (in minimizer-index space).
    ///
    /// Foldbacks cluster matchpoints around an anti-diagonal using a transform like
    /// `d = p1 + p2`. The tolerance controls how tightly we bucket those matchpoints.
    ///
    /// (Concatemers typically use `delta = p2 - p1` instead, hence “separate” tolerances.)
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

/// Foldback-only gates.
///
/// These thresholds are applied during coarse detection before refinement.
/// They capture basic geometric/evidence constraints for calling a foldback candidate.
#[derive(Debug, Clone)]
pub struct FoldOnlyCfg {
    /// Minimum span (bp) of “arm evidence” supporting the foldback on p1.
    ///
    /// This helps reject tiny, locally repetitive matches that don't represent a real
    /// foldback junction.
    pub min_arm: usize,
}

impl Default for FoldOnlyCfg {
    fn default() -> Self {
        Self { min_arm: 200 }
    }
}

/// Second-pass configuration for foldback decisions.
///
/// Pass2 uses the support map built in pass1 to decide whether a foldback junction is:
/// - “real” (genome-templated; should recur across reads), or
/// - an artefact (e.g. MDA-induced; safe to cut/correct).
#[derive(Debug, Clone)]
pub struct FoldSecondPassCfg {
    /// Minimum number of reads supporting the same junction fingerprint for “real”.
    ///
    /// If a signature reaches this support threshold, the default behavior is to
    /// treat it as *real* and therefore **do not** cut/correct (unless overridden).
    pub min_support: usize,

    /// Minimum refinement identity estimate required to accept as a foldback candidate.
    ///
    /// This is applied in pass2 before computing/looking up signatures, to avoid
    /// spending work on low-confidence calls.
    pub min_identity: f32,

    /// Minimum support identity estimate required to accept as a foldback candidate.
    pub min_support_ident: f64,
}

impl Default for FoldSecondPassCfg {
    fn default() -> Self {
        Self {
            min_support: 3,
            min_identity: 0.60,
            min_support_ident: 0.0,
        }
    }
}

/// High-level “preset” call modes.
///
/// These map to tuned thresholds that trade off:
/// - sensitivity (how many events you call / cut)
/// - specificity (false positive rate)
/// - computational cost (depth, signature params)
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CallMode {
    /// Conservative: fewer calls/cuts, higher confidence.
    Strict,
    /// Middle-ground default: practical sensitivity/specificity trade-off.
    Balanced,
    /// Aggressive: more calls/cuts; good for exploratory runs.
    Permissive,
}

impl Default for CallMode {
    fn default() -> Self {
        CallMode::Balanced
    }
}

/// Concrete parameter bundle derived from a `CallMode`.
///
/// This is convenient for CLI exposure: the user chooses a qualitative mode,
/// and you apply a consistent set of thresholds across detection, refinement,
/// recursion, and signature generation.
#[derive(Debug, Clone)]
pub struct ModeTuning {
    /// Minimum evidence span threshold (bp).
    ///
    /// In practice this may combine foldback `min_arm` and any concatemer span threshold.
    pub min_span: usize,
    /// Coarse detector requirement: minimum minimizer matchpoints.
    pub min_matches: usize,
    /// Refinement identity gate.
    pub min_identity: f32,
    /// Support threshold for calling a foldback “real”.
    pub min_support: usize,
    /// Required tightness of split clustering (bp) for “real”.
    pub split_tol_bp: usize,
    /// Maximum recursion depth for iterative cutting of nested artefacts.
    pub max_depth: usize,
    /// Signature flank size (bp) used around the split.
    pub sig_flank_bp: usize,
    /// Number of minimizers (or minimizer-derived values) retained per side in signature.
    pub sig_take: usize,
    /// Optional value shifting/quantization to stabilize signatures.
    pub sig_value_shift: u8,
}

impl CallMode {
    /// Return tuned parameter defaults for the selected `CallMode`.
    ///
    /// The intent is that these presets are internally coherent:
    /// - stricter modes demand more evidence and reduce recursion
    /// - permissive modes allow weaker evidence and increase recursion
    /// - signature params scale with strictness (bigger flanks/more takes when stricter)
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
                sig_flank_bp: 200,
                sig_take: 6,
                sig_value_shift: 8,
            },
        }
    }
}

/// Parameters chosen to approximate the legacy Perl pipeline for “fair” comparison.
///
/// This is useful for benchmarking/regression testing:
/// - keep most structural/minimizer/refinement knobs fixed (“LOCKED”)
/// - vary only the thresholds that reflect algorithmic differences
///
/// You can treat this as a reference preset and compare call rates / correction behavior.
#[derive(Debug, Clone)]
pub struct FairnessParams {
    // ---- detection / decision ----
    /// Coarse detector requirement: minimum minimizer matchpoints.
    pub min_matches: usize,
    /// Minimum evidence span threshold (bp).
    pub min_span: usize,
    /// Refinement identity gate.
    pub min_identity: f32,
    /// Support threshold for “real foldback”.
    pub min_support: usize,
    /// Required tightness of split clustering (bp) for “real”.
    pub split_tol_bp: usize,
    /// Maximum recursion depth for iterative cutting.
    pub max_depth: usize,

    // ---- minimizers ----
    /// Minimizer k-mer length (LOCKED in baseline preset).
    pub k: usize,
    /// Minimizer window size (LOCKED in baseline preset).
    pub w: usize,
    /// Strand handling for minimizers (LOCKED in baseline preset).
    pub forward_only: bool,

    // ---- refinement ----
    /// Refinement mode (LOCKED in baseline preset).
    pub refine_mode: RefineMode,
    /// Half-width of refinement search window (LOCKED in baseline preset).
    pub refine_window: usize,
    /// Arm length used during refinement (LOCKED in baseline preset).
    pub refine_arm: usize,
    /// Maximum tolerated edit-distance rate (LOCKED in baseline preset).
    pub max_ed_rate: f32,

    // ---- structural guards ----
    /// Guard region at read ends (LOCKED in baseline preset).
    pub end_guard: usize,
    /// Foldback diagonal tolerance (LOCKED in baseline preset).
    pub fold_diag_tol: i32,

    // ---- signature ----
    /// Signature flank size (bp).
    pub sig_flank_bp: usize,
    /// Number of signature values retained per side.
    pub sig_take: usize,
}

impl FairnessParams {
    /// Baseline preset approximating the legacy Perl/minimap2 pipeline.
    ///
    /// Use this for:
    /// - apples-to-apples comparisons on the same dataset
    /// - regression tests when refactoring detection/refinement
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

            // signature
            sig_flank_bp: 1000,
            sig_take: 12,
        }
    }
}

/// Signature configuration.
///
/// Signatures are compact fingerprints used as keys into the pass1 support map.
/// They are designed to be:
/// - stable across small breakpoint jitter (as much as possible)
/// - fast to compute
/// - specific enough to cluster the same junction across reads
#[derive(Debug, Clone)]
pub struct SigCfg {
    /// Half-flank length (bp) taken on either side of the split when building a signature.
    pub flank_bp: usize,
    /// How many minimizer-derived values are retained in the signature (per side).
    pub take: usize,
    /// Optional shift/quantization applied to signature values for stability.
    pub value_shift: u8,
}

impl Default for SigCfg {
    fn default() -> Self {
        Self {
            flank_bp: 1000,
            take: 12,
            value_shift: 0, // default is exact (no quantization)
        }
    }
}
