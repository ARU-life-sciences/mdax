//! Configuration of `mdax`
use crate::utils::RefineMode;

/// Global configuration for `mdax`.
///
/// This is the top-level config object that is constructed (mostly) from CLI args.
/// It is composed of:
/// - `shared`: parameters used by *both* coarse detection and refinement steps
/// - `fold`: foldback-specific detection gates
/// - `fold2`: second-pass support-aware decision thresholds
#[derive(Debug, Clone, Default)]
pub struct MdaxCfg {
    /// Parameters shared across detection modes (minimizers, guards, refinement).
    pub shared: SharedCfg,
    /// Foldback-only gates for coarse detection.
    pub fold: FoldOnlyCfg,
    /// Second-pass decision thresholds for “real vs artefact” foldbacks.
    pub fold2: FoldSecondPassCfg,
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
/// - `max_jump_clip`: maximum template-jump to probe (bp). After finding the best
///   split position, the refiner scans offsets δ = 0..max_jump_clip on the right arm
///   to find the δ that maximises identity over the `arm − δ` aligned bases.
///   A minimum effective arm of 100 bp is enforced to avoid noise at large δ.
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
    /// Maximum template-jump to probe during identity estimation (bp).
    pub max_jump_clip: usize,
}

impl Default for RefineCfg {
    fn default() -> Self {
        Self {
            window: 200,
            arm: 1200,
            mode: RefineMode::HiFi,
            max_ed_rate: 0.25,
            max_jump_clip: 1000,
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
/// Alignment-verified "rescue" detection for foldbacks the coarse detector's
/// evidence gates reject: short arms and junctions close to a read end.
///
/// The coarse gates (`min_matches`, `min_arm`, `end_guard`) use the *amount*
/// of minimizer evidence as a proxy for "this is a real inverted duplicate".
/// That proxy cannot be loosened without admitting spurious seed
/// coincidences. The rescue path instead accepts very little seed evidence
/// (`min_pairs` matched minimizer pairs on one anti-diagonal) but then
/// requires the two arms to actually align, base for base, over at least
/// `min_arm` bp at `min_identity` or better -- sensitivity from the seeds,
/// specificity from the alignment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RescueCfg {
    /// Off by default at the struct level so library callers and tests that
    /// build a bare `FoldOnlyCfg` keep the historical behaviour; the `mdax`
    /// binary switches it on.
    pub enabled: bool,
    /// Minimum matched minimizer pairs (left-arm copy + right-arm copy) with
    /// distinct values on the candidate anti-diagonal.
    pub min_pairs: usize,
    /// Minimum verified arm length (bp).
    pub min_arm: usize,
    /// Minimum identity of the verified arm alignment.
    pub min_identity: f32,
    /// Longest arm window (bp) the verifier will align.
    pub max_window: usize,
    /// How many candidate anti-diagonals to verify per read at most.
    pub max_candidates: usize,
    /// Only rescue candidates whose shorter arm runs to a read end.
    ///
    /// Measured against a reference on WGA HiFi data, 94% of the true
    /// artefacts the rescue path adds sit within one arm of a read end (the
    /// strand switch happened shortly before the fragment ended), whereas
    /// rescued candidates in the middle of a read were as often genome-encoded
    /// short inverted repeats with too little coverage to clear the support
    /// gate as they were artefacts. With this set, a candidate must have its
    /// outer arm end at the read end, and the whole of that arm must align.
    pub end_only: bool,
    /// When an arm must run to a read end, how many of its outermost bases
    /// are left out of the alignment (bp). Read ends carry adapter remnants
    /// and low-quality bases that are not part of the inverted copy; on raw
    /// ONT reads, requiring them to align cost a tenth of the identity and
    /// with it a third of otherwise verifiable candidates.
    pub end_slack: usize,
    /// Also scan both read ends with dense (every-position) short seeds.
    ///
    /// The commonest short-arm foldback is one whose junction simply falls
    /// near a read end, so its outer arm is cut short by the end of the
    /// read. Sampled minimizers are too sparse to propose those reliably on
    /// noisy reads, so the last/first `term_window` bp are matched against
    /// their own reverse complement using all `term_k`-mers instead. A
    /// terminal candidate is only accepted if the inverted copy runs all the
    /// way to the read end -- the signature of a strand switch followed by
    /// synthesis to the end of the fragment, as opposed to a small hairpin
    /// that merely sits near the end.
    pub terminal: bool,
    /// Seed length for the terminal scan.
    pub term_k: usize,
    /// How much of each read end the terminal scan looks at (bp).
    pub term_window: usize,
    /// Minimum co-diagonal dense seeds proposing a terminal candidate.
    pub term_min_seeds: usize,
}

impl RescueCfg {
    pub const fn disabled() -> Self {
        Self {
            enabled: false,
            min_pairs: 3,
            min_arm: 150,
            min_identity: 0.90,
            max_window: 500,
            max_candidates: 3,
            end_only: true,
            end_slack: 30,
            terminal: false,
            term_k: 13,
            term_window: 1200,
            term_min_seeds: 4,
        }
    }
}

impl Default for RescueCfg {
    fn default() -> Self {
        Self::disabled()
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
    /// Alignment-verified rescue of candidates the gates above reject.
    pub rescue: RescueCfg,
}

impl Default for FoldOnlyCfg {
    fn default() -> Self {
        Self {
            min_arm: 200,
            rescue: RescueCfg::disabled(),
        }
    }
}

/// Second-pass configuration for foldback decisions.
///
/// Pass2 uses the support map built in pass1 to decide whether a foldback junction is:
/// - “real” (genome-templated; should recur across reads), or
/// - an artefact (e.g. MDA-induced; safe to cut/correct).
#[derive(Debug, Clone)]
pub struct FoldSecondPassCfg {
    /// Minimum number of reads folding about the same axis for “real”.
    ///
    /// If a junction reaches this support threshold (see
    /// `fingerprint::JunctionSupport`), it is treated as *real* and therefore
    /// **not** cut/corrected.
    pub min_support: usize,

    /// Minimum refinement identity estimate required to accept as a foldback candidate.
    ///
    /// Applied in pass2 before the support lookup; what happens to reads below
    /// it is governed by `cut_low_ident`.
    pub min_identity: f32,

    /// When true, reads detected with a foldback whose arm identity is below
    /// `min_identity` still go through the support gate and are cut if
    /// unsupported, rather than labelled `low_ident` and passed through.
    /// Set from the CLI as `!--keep-low-ident` (the legacy `--cut-low-ident`
    /// flag is a no-op).
    ///
    /// Useful for ONT data where the HiFi Hamming estimator systematically
    /// underestimates arm similarity on indel-rich reads.
    pub cut_low_ident: bool,
}

impl Default for FoldSecondPassCfg {
    fn default() -> Self {
        Self {
            min_support: 3,
            min_identity: 0.60,
            cut_low_ident: false,
        }
    }
}

/// High-level “preset” call modes.
///
/// These map to tuned thresholds that trade off:
/// - sensitivity (how many events you call / cut)
/// - specificity (false positive rate)
/// - computational cost (recursion depth)
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum CallMode {
    /// Conservative: fewer calls/cuts, higher confidence.
    Strict,
    /// Middle-ground default: practical sensitivity/specificity trade-off.
    #[default]
    Balanced,
    /// Aggressive: more calls/cuts; good for exploratory runs.
    Permissive,
}

/// Concrete parameter bundle derived from a `CallMode`.
///
/// This is convenient for CLI exposure: the user chooses a qualitative mode,
/// and you apply a consistent set of thresholds across detection, refinement,
/// and recursion.
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
    /// Maximum recursion depth for iterative cutting of nested artefacts.
    pub max_depth: usize,
    /// Do not call foldbacks whose split lies within this distance of a read end (bp).
    pub end_guard: usize,
}

impl CallMode {
    /// Return tuned parameter defaults for the selected `CallMode`.
    ///
    /// The intent is that these presets are internally coherent:
    /// - stricter modes demand more evidence and reduce recursion
    /// - permissive modes allow weaker evidence and increase recursion
    pub fn tuning(self) -> ModeTuning {
        match self {
            // strict: fewer calls, fewer cuts, higher confidence
            CallMode::Strict => ModeTuning {
                min_matches: 30,
                min_span: 3000,
                min_identity: 0.75,
                min_support: 5,
                max_depth: 2,
                end_guard: 1000,
            },
            // balanced: your current-ish defaults
            // balanced.
            //
            // min_matches / min_span / end_guard were 20 / 2000 / 1000 until
            // 2026-11-18. Scored against a reference (a random 3% of a
            // ~100x C. elegans WGA HiFi library aligned to WBcel235), those
            // values found 96% of true foldbacks with arms over 4 kb but
            // only 72% at 1-2 kb, 44% at 500-1000 bp and 23-38% of those
            // within 1 kb of a read end -- 28% missed overall. 12 / 1000 /
            // 500 brings the miss rate to 21% (78% of true artefacts cut,
            // up from 72%) and assembled at least as well (hifiasm: 682
            // contigs vs 705, same N50 and BUSCO). Going further (8 / 600 /
            // 300) gained 1 more point of recall, doubled the linear reads
            // cut, and did nothing for the assembly, so stop here.
            CallMode::Balanced => ModeTuning {
                min_matches: 12,
                min_span: 1000,
                min_identity: 0.60,
                min_support: 3,
                max_depth: 5,
                end_guard: 500,
            },
            // permissive: more calls/cuts (good for exploration)
            CallMode::Permissive => ModeTuning {
                min_matches: 12,
                min_span: 1200,
                min_identity: 0.50,
                min_support: 2,
                max_depth: 8,
                end_guard: 1000,
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
            max_depth: 1,

            // minimizers (LOCKED)
            k: 17,
            w: 21,
            forward_only: true,

            // refinement (LOCKED)
            refine_mode: RefineMode::HiFi,
            refine_window: 200,
            refine_arm: 1200,
            max_ed_rate: 0.25,

            // guards (LOCKED)
            end_guard: 1000,
            fold_diag_tol: 120,
        }
    }
}
