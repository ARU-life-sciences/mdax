//! CLI definition for `mdax` chimeric read detection and clean-up.
//!
//! Here we:
//! - Declaring CLI arguments (paths, thresholds, presets)
//! - Providing clap parsing via `build_cli()`
//! - Wiring custom enums (`RefineMode`, `CallMode`) into clap as `ValueEnum`
//!
//! Design notes:
//! - Most “real knobs” live in `cfg.rs` as strongly typed config structs.
//! - The CLI exposes either:
//!   - a high-level preset (`--mode strict|balanced|permissive`), or
//!   - explicit overrides (`--min-identity`, `--min-support`, ...).
//! - Enums are implemented as `ValueEnum + FromStr` so clap can parse them cleanly.

// TODO: cut the fasta reads based on the output table

use crate::{cfg::CallMode, utils::RefineMode};
use clap::{builder::PossibleValue, value_parser, Arg, ArgMatches, Command, ValueEnum};
use std::path::PathBuf;

/// Build and parse the command-line interface, returning clap’s `ArgMatches`.
///
/// This constructs a `Command` with:
/// - positional `input` FASTA/FASTQ (potentially gzipped)
/// - output FASTA path
/// - optional TSV report destination (`-` for stdout)
/// - detection/refinement thresholds (either as a preset `--mode` or explicit overrides)
pub fn build_cli() -> ArgMatches {
    let c = Command::new("mdax")
        .version(clap::crate_version!())
        .author("Max Carter-Brown")
        .about("Remove inverted duplication chimeras from fasta files.")
        // ----------------------------
        // Inputs / outputs
        // ----------------------------
        //
        // Positional input FASTA file (optionally gzipped).
        // This is the read library to scan for foldback/palindrome artefacts.
        .arg(
            Arg::new("input")
                .help("Input file")
                .long_help(
                    "Input FASTA/FASTQ(.gz) containing reads to scan for foldback (inverted-duplication) artefacts.\n\
\n\
Foldbacks are detected *within each read* by looking for long, reverse-complement self-matches.\n\
This file is streamed: runtime is roughly linear in total bases, and memory is bounded by channel buffering.\n\
\n\
Tip: if your reader supports it, `.gz` inputs are convenient but decompression can become a bottleneck at high thread counts."
                )
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .index(1),
        )
        // Output TSV report file path.
        .arg(
            Arg::new("report")
                .help("Write a TSV report of detected concatemers")
                .long_help(
                    "Write a TSV report describing detected events and decisions.\n\
\n\
The TSV is intended for *debugging, QC, and post filtering*:\n\
- read_id: Read identifier from the FASTA/FASTQ header.\n\
- len: Read length in bases.\n\
- event: Event type detected (always foldback, may include e.g. concatemers).\n\
- called: 1 if an event was detected in pass1, else 0.\n\
- coarse_split: Approximate breakpoint index from the first-pass detection (before refinement).\n\
- refined_split: Breakpoint index after local refinement around coarse_split.\n\
- delta: refined_split - coarse_split (signed). CURRENTLY EMPTY.\n\
- matches: Number of minimizer matchpoints supporting the foldback call (pass1 evidence strength).\n\
- span_p1: Span (bp) covered by pass1 supporting matchpoints / arms. I.e. how much of the read supports the foldback geometry\n\
- p2_span: Span used/confirmed in pass2 (if you do a second-pass validation). CURRENTLY EMPTY.\n\
- cross_frac: Fraction of matchpoints crossing the split / inconsistent with a clean junction (often a “messiness” metric). CURRENTLY EMPTY.\n\
- coarse_score: Score assigned by the coarse foldback detector prior to refinement.\n\
- refined_score: Score assigned after breakpoint refinement.\n\
- identity_est: Estimated nucleotide identity between the two foldback arms.\n\
- support_n: Number of other reads supporting the same foldback junction.\n\
- support_rank_frac: Relative prominence of this junction within the dataset, computed by ranking all junctions by support_n (descending) and normalising the rank to [0,1]. Values near 1 indicate highly recurrent junctions; values near 0 indicate rare/singleton events.\n\
- support_span: Total supporting span (bp) aggregated across supporting reads.\n\
- decision: Final classification based on thresholds (real, artefact, unknown, low_ident).\n\
\n\
Use '-' to write to stdout\n\
"
                )
                .short('r')
                .long("report")
                .value_parser(value_parser!(PathBuf))
                .default_value("-"),
        )
        // Output FASTA file path.
        // TODO: by default we split the read taking the first 'half' of each read
        .arg(
            Arg::new("output")
                .help("Output FASTA file path")
                .long_help("Output FASTA containing the kept/corrected reads.\n\
\n\
For each input read, mdax writes one output record containing the sequence you decided to keep:\n\
- If no foldback is called, the read is written unchanged.\n\
- If a foldback is called and judged an artefact, mdax cuts/chops the read (possibly recursively) and writes the kept portion.\n\
- If a foldback is judged 'real' (supported across reads), mdax keeps it intact by default.\n\
")
                .short('o')
                .long("output")
                .value_parser(value_parser!(PathBuf))
                .required(true),
        )
        // Minimizer k-mer size.
        // Larger values are more specific but less tolerant of sequencing errors.
        // Typical: ONT 15–17, HiFi 19–25 (depends on read quality).
        .arg(
            Arg::new("k")
                .help("Minimizer k-mer size.")
                .short('k')
                .long_help("Minimizer k-mer size used in the *coarse detector*.\n\
\n\
How it affects foldback calling:\n\
- Larger k makes minimizers more specific, reducing random/self-repeat matches, leading to fewer false positives.\n\
- Larger k is less tolerant of sequencing errors (especially ONT), so true foldbacks may have fewer matching minimizers and therefore fewer calls.\n\
\n\
Practical guidance:\n\
- ONT: k=15–17 is typical; higher k can reduce sensitivity on noisy reads.\n\
- HiFi: k=19–25 can work well because errors are rare.\n\
\n\
Notes:\n\
- k does not directly change refinement quality; it changes which candidates reach refinement.\n\
")
                .default_value("17")
                .value_parser(value_parser!(usize)),
        )
        // Minimizer window size (in k-mers).
        // Each window chooses one minimizer; larger values sample fewer positions.
        // Typical: 10–30.
        // NOTE: must be odd!
        .arg(
            Arg::new("w")
                .help("Minimizer window size. Must be odd.")
                .long_help("Minimizer window size used in the *coarse detector*.\n\
\n\
How it affects foldback calling:\n\
- Larger w samples fewer minimizers (one per window), so faster, but fewer matchpoints.\n\
  This can reduce sensitivity, especially for short arms or noisy reads.\n\
- Smaller w samples more minimizers, so slower, but increases matchpoint density.\n\
  This generally increases sensitivity and stabilizes split estimation at the coarse stage.\n\
\n\
Rules of thumb:\n\
- If you're missing obvious foldbacks: decrease w or decrease min_matches.\n\
- If you're calling too many weak events: increase w or increase min_matches/min_span.\n\
\n\
Note: 'must be odd', see `simd_minimizers` Rust crate for more details.")
                .short('w')
                .default_value("21")
                .value_parser(value_parser!(usize)),
        )
        // The below typically over-ride preset tuning.
        // 
        // Minimum palindrome arm length (bp) required to call a foldback.
        // This is measured as the span of matched minimizer positions on one side of the junction.
        // Prevents calling short/random palindromic matches.
        .arg(
            Arg::new("min_span")
                .help("Override: minimum evidence span (bp) for calling a foldback")
                .long_help("Override: minimum evidence span (bp) required for a coarse foldback call.\n\
\n\
What 'span' means here:\n\
- It is the span of matched minimizer positions along (typically) one arm of the junction.\n\
- Conceptually: how much of the read participates in the reverse-complement self-match.\n\
\n\
How it affects foldback calling:\n\
- Increasing min-span demands longer arm evidence. Fewer calls, much fewer random palindromic hits.\n\
- Decreasing min-span allows shorter events. More calls, but more susceptibility to local repeats.\n\
\n\
Interaction with other knobs:\n\
- If you decrease w (more minimizers), you may be able to increase min-span without losing sensitivity.\n\
- If you keep min-span low, consider increasing min_identity to avoid cutting based on weak matches.")
                .long("min-span")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_identity")
              .long("min-identity")
              .help("Override: minimum refinement identity for foldback second pass")
              .long_help("Override: minimum refinement identity (identity_est) to treat a detected foldback as credible in pass2.\n\
\n\
Where it applies:\n\
- After coarse detection and refinement, `mdax` estimates similarity between the two arms.\n\
- If identity_est < min-identity, the event is treated as low-confidence (e.g. labelled 'low_ident')\n\
  and typically will not be used for support-aware 'real vs artefact' decisions.\n\
\n\
How it affects calling/cutting:\n\
- Increasing min-identity means fewer events considered; reduces over-cutting on noisy/weak alignments.\n\
- Decreasing min-identity means more events considered; increases sensitivity but may cut reads on spurious matches.\n\
\n\
Practical guidance:\n\
- ONT data often benefits from a lower identity threshold than HiFi.\n\
- If you're seeing many false artefact cuts, raise this first before touching support thresholds.")
              .value_parser(value_parser!(f32)),
        )
        .arg(
            Arg::new("min_support")
              .long("min-support")
              .help("Override: minimum support reads (coverage) for 'real' foldback")
              .long_help("Override: minimum number of reads that must share the same junction signature for the foldback to be treated as 'real'.\n\
\n\
Interpretation:\n\
- Pass1 builds a support map keyed by a fingerprint/signature of the junction.\n\
- In pass2, if a signature is observed in >= min-support reads (and split clustering is tight enough),\n\
  `mdax` treats that junction as genome-templated and does NOT cut it by default.\n\
\n\
How it affects cutting:\n\
- Increasing min-support means harder to be considered real. More foldbacks treated as artefacts therefore more cutting.\n\
- Decreasing min-support means easier to be considered real. Fewer cuts, but risk of keeping recurrent artefacts.\n\
\n\
Guidance:\n\
- For deep coverage datasets, you can raise this to be more conservative about calling 'real'.\n\
- For shallow coverage, keep it low, otherwise almost nothing will qualify as real.")
              .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max_depth")
              .long("max-depth")
              .help("Override: maximum recursion depth for chopping")
              .long_help("Override: maximum recursion depth when cutting foldback artefacts.\n\
\n\
What recursion does:\n\
- After cutting an artefact foldback, the kept fragment may still contain another foldback junction.\n\
- `mdax` can re-run detection on the kept fragment and cut again, up to max-depth times.\n\
\n\
How it affects output:\n\
- Increasing max-depth leads to more aggressive cleanup of nested/compound artefacts, but may over-trim in noisy data.\n\
- Decreasing max-depth leads to faster and safer, but may leave residual artefact structure.\n\
\n\
Guidance:\n\
- Start modest (e.g. 2–5). Raise only if TSV shows many 'artefact' decisions remaining after one cut.\n\
- Very high depth can amplify mistakes: one false cut early can cascade.
- Higher recursion, more compute!")
              .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("end_guard")
                .help("Do not call breakpoints within this distance of read ends (bp)")
                .long_help("Guard region at both ends of reads (bp).\n\
\n\
Foldback calls near read ends are often unreliable because:\n\
- there is little sequence context on one side of the split\n\
- minimizer matchpoint chains become truncated\n\
- refinement can be biased by missing arm sequence\n\
\n\
How it affects calling:\n\
- Increasing end-guard leads to fewer calls near ends (more conservative), may miss real junctions that genuinely occur close to ends.\n\
- Decreasing end-guard leads to more calls, but higher risk of unstable splits and over-cutting.\n\
\n\
Tip:\n\
- If you have many short reads, a large end-guard can suppress almost all calls.\n\
  Ensure end-guard is comfortably smaller than typical read length.")
                .long("end-guard")
                .default_value("1000")
                .value_parser(value_parser!(usize)),
        )
        // Refinement is the local, more expensive step after coarse detection.
        // These parameters generally map onto `RefineCfg`.
        .arg(
            Arg::new("refine_window")
                .help("Refinement half-window around coarse split (bp)")
                .long_help("Refinement search half-window (bp) around the coarse split position.\n\
\n\
Refinement works by testing candidate split positions around the coarse estimate.\n\
This parameter controls the maximum shift allowed during refinement.\n\
\n\
How it affects calls:\n\
- Increasing refine-window can recover splits when coarse detection is jittery, but costs more compute.\n\
- Decreasing refine-window is faster, but may 'miss' the true breakpoint if coarse split is off.\n\
\n\
If you observe large coarse to refined deltas in the TSV, increase this.\n\
If deltas are always small, you can reduce it to speed up.")
                .long("refine-window")
                .default_value("100")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("refine_arm")
                .help("Refinement arm length (bp)")
                .long_help("Arm length (bp) used during refinement to compare the two sides of the junction.\n\
\n\
Refinement typically extracts `arm` bases on each side of the candidate split and computes\n\
a similarity score (via Hamming for HiFi or banded Levenshtein for ONT).\n\
\n\
How it affects identity estimates and stability:\n\
- Increasing refine-arm gives more evidence, smoother identity estimates, fewer spurious calls; slower.\n\
- Decreasing refine-arm is faster and more local, but identity becomes noisier and can be biased by repeats.\n\
\n\
Guidance:\n\
- If you see unstable identity_est or lots of borderline calls, increase this.\n\
- If reads are short, arm too large can hit end-guard logic or reduce usable candidates.")
                .long("refine-arm")
                .default_value("200")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max_ed_rate")
                .help("ONT refinement bandwidth as fraction of arm")
                .long_help("Maximum tolerated edit-distance rate during ONT refinement.\n\
\n\
In ONT mode, refinement uses a banded edit distance. The band width is derived from the arm length\n\
and this max-ed-rate (e.g. band ≈ arm * max_ed_rate).\n\
\n\
How it affects refinement:\n\
- Increasing max-ed-rate gives wider band. More tolerant of indels/noise, but slower.\n\
- Decreasing max-ed-rate gives narrower band. Faster, but may underestimate similarity on noisy reads.\n\
\n\
This mainly impacts ONT mode; HiFi mode often uses a cheaper/Hamming-like comparison.\n\
If ONT refinement is missing obvious foldbacks, raise this slightly.\n\
If ONT refinement is slow, lower this (but watch for identity_est becoming systematically low).")
                .long("max-ed-rate")
                .default_value("0.25")
                .value_parser(value_parser!(f32)),
        )
        .arg(
            Arg::new("min_matches")
                .help("Override: minimum minimizer matches supporting a foldback candidate")
                .long_help("Override: minimum number of minimizer matchpoints required for a coarse foldback call.\n\
\n\
Matchpoints are the basic evidence units for coarse detection.\n\
They come from matching minimizers between a read segment and its reverse-complement self-match.\n\
\n\
How it affects calling:\n\
- Increasing min-matches gives fewer calls, higher specificity, better split stability.\n\
- Decreasing min-matches gives more calls (including weak ones), increased risk of calling local repeats.\n\
\n\
Interaction with k/w:\n\
- Larger k or larger w reduces matchpoint density, so you may need to lower min-matches.\n\
- Smaller k or smaller w increases matchpoint density, so you can raise min-matches to stay specific.")
                .long("min-matches")
                .value_parser(value_parser!(usize)),
        )
        // Foldbacks cluster matchpoints around an anti-diagonal of (pos1 + pos2).
        .arg(
            Arg::new("fold_diag_tol")
                .help("Foldback anti-diagonal bucketing tolerance (bp) for pos1+pos2 clustering")
                .long_help("Foldback anti-diagonal clustering tolerance.\n\
\n\
The foldback detector buckets minimizer matchpoints by an anti-diagonal coordinate, using\n\
d = pos1 + pos2. True foldbacks produce a dense cluster along an anti-diagonal;\n\
noise and repeats produce scattered points.\n\
\n\
How it affects calling:\n\
- Decreasing fold-diag-tol leads to tighter clustering and fewer calls, but can fragment evidence if reads are noisy.\n\
- Increasing fold-diag-tol leads to looser clustering and more calls, but higher chance of merging unrelated matches.\n\
\n\
If you see many borderline calls with scattered matchpoints, decrease this.\n\
If you miss foldbacks on noisy reads where matchpoints are smeared, increase this slightly.")
                .long("fold-diag-tol")
                .default_value("120")
                .value_parser(value_parser!(i32)),
        )
        // Minimizer strand handling.
        .arg(
            Arg::new("forward_only")
                .help("Use strand-specific forward-only minimizers (true) or canonical minimizers (false)")
                .long_help("Whether minimizer hashing is strand-specific.\n\
\n\
- forward-only (true): minimizers depend on the read’s forward orientation.\n\
  This tends to be good for coarse detection and avoids some symmetrical collisions.\n\
- canonical (false): a k-mer and its reverse complement map to the same minimizer value.\n\
  This can increase sensitivity in some cases, but can also increase spurious self-matches in repetitive regions.\n\
\n\
How it affects foldback calling:\n\
- canonical may increase the number of matchpoints (more hits), which can increase calls.\n\
- forward-only can be more conservative and may reduce false positives.\n\
\n\
Note: we use canonical minimizers for *signatures* even if forward-only for detection;\n\
this flag applies to the detector’s minimizer sampling.")
                .long("forward-only")
                .default_value("true")
                .value_parser(value_parser!(bool)),
        )
        // This selects a refinement “error model” and therefore affects both speed and
        // tolerance to indels/substitutions.
        .arg(
            Arg::new("refine_mode")
                .help("Refine breakpoints for reads")
                .long_help("Select refinement mode (error model) used when refining breakpoint coordinates.\n\
\n\
- hifi: fast comparison (Hamming-like) assuming low error rates.\n\
  Produces stable split positions and identity estimates when reads are accurate.\n\
- ont: banded Levenshtein-style refinement tolerant to indels/substitutions.\n\
  Slower but necessary for noisier reads.\n\
\n\
How it affects calling/cutting:\n\
- ONT mode can rescue true foldbacks that would look low-identity under HiFi assumptions.\n\
- HiFi mode can prevent over-fitting noise (calling weak alignments as foldbacks).\n\
\n\
Important:\n\
- This controls refinement; coarse detection still depends heavily on minimizer parameters and gates.\n\
- If you're using ONT reads, start with ont; if you're using HiFi reads, start with hifi.")
                .short('m')
                .long("refine-mode")
                .default_value("hifi")
                .value_parser(value_parser!(RefineMode)),
        )
        // This is the high-level knob intended for most users; fine-grained overrides exist
        // for debugging/tuning.
        .arg(
              Arg::new("mode")
                .help("Detection/correction strictness: strict reduces calls, permissive increases calls")
                .long_help("High-level preset controlling detection and correction strictness.\n\
\n\
This sets a coherent bundle of thresholds (min_matches, min_span, min_identity, min_support,\n\
max_depth, and signature parameters).\n\
\n\
- strict: fewer calls and fewer cuts. Requires strong evidence and high identity.\n\
  Best when you want to avoid over-trimming.\n\
- balanced: reasonable default trade-off.\n\
- permissive: more calls and more cutting. Useful for exploration or very messy libraries.\n\
\n\
Override precedence:\n\
- If you pass explicit overrides (e.g. --min-identity), they override the preset values (will say in help).\n\
- Use the stderr output or TSV report to confirm which effective parameters were applied.")
                .long("mode")
                .default_value("balanced")
                .value_parser(value_parser!(CallMode)),
        )
        // These impact how stable and specific the support-map signature keys are.
        .arg(
            Arg::new("sig_flank_bp")
                .long("sig-flank-bp")
                .help("Window (bp) around split for match-based fingerprinting")
                .long_help("Signature flank size (bp) used to compute the junction fingerprint.\n\
\n\
Signatures are used to build the pass1 support map and to look up support in pass2.\n\
They should be stable across small split jitter but specific enough not to collide.\n\
\n\
How it affects 'real vs artefact' decisions:\n\
- Increasing sig-flank-bp gives more context included and fewer collisions, more stable support clusters,\n\
  but slightly more compute and potentially more sensitivity to distant repeats.\n\
- Decreasing sig-flank-bp gives signature becomes more local and is faster and less affected by distant repeats,\n\
  but higher chance different junctions produce the same signature (collisions).\n\
\n\
If you see many signatures that are 'in support' but with wide split scatter, consider increasing\n\
sig_flank_bp and/or sig_take to improve specificity.")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("sig_take")
                .long("sig-take")
                .long_help("Number of minimizer-derived values retained when constructing the junction signature.\n\
\n\
This controls signature specificity:\n\
- Larger sig-take; more information in the signature, fewer collisions, better support discrimination.\n\
- Smaller sig-take; less information, more collisions, more chance of misclassifying artefacts as real (or vice versa).\n\
\n\
Cost:\n\
- Larger values are slightly slower and may fail more often if there are too few minimizers in the flank window.\n\
  (In that case you may see fewer sig_ok events.)\n\
\n\
Guidance:\n\
- If you want robust support-map behaviour, prefer increasing sig-take before relaxing split tolerance.\n\
- If you process very short reads or tiny flanks, keep take modest so signatures can still be formed.")
                .help("Number of minimizer hashes retained for foldback fingerprint")
                .value_parser(value_parser!(usize)),
        )
        // Hidden “developer mode” for reproducible benchmarking.
        .arg(
            Arg::new("fairness_baseline")
                .long("fairness-baseline")
                .hide(true) 
                .help("Developer-only: lock parameters for fair benchmarking")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of worker threads to use.")
                .long_help("Number of worker threads to parallelise foldback detection and correction.\n\
\n\
`mdax` uses a multi-threaded pipeline with a single reader and writer thread and \
multiple worker threads for detection, refinement, and classification.\n\
\n\
Notes:\n\
- Increasing threads improves performance up to the point where the workload \
becomes memory- or I/O-bound.\n\
- For compute-only runs (e.g. -o /dev/null -r /dev/null), scaling typically \
plateaus well below the total number of CPU cores.\n\
- For runs that write FASTA/TSV output, disk I/O can dominate wall-clock time, \
and increasing threads may provide little or no speedup.\n\
\n\
Default is a single core, multithreading will likely give gains up to 24-32 threads."
        )
    );


    c.get_matches()
}

// Change RefineMode to a clap ValueEnum
impl ValueEnum for RefineMode {
    fn value_variants<'a>() -> &'a [Self] {
        &[RefineMode::HiFi, RefineMode::ONT]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        Some(match self {
            RefineMode::HiFi => {
                PossibleValue::new("hifi").help("Refine breakpoints for HiFi reads (Hamming; fast)")
            }
            RefineMode::ONT => PossibleValue::new("ont")
                .help("Refine breakpoints for ONT reads (Levenshtein, slower)"),
        })
    }
}

impl std::fmt::Display for RefineMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.to_possible_value()
            .expect("no values are skipped")
            .get_name()
            .fmt(f)
    }
}

impl std::str::FromStr for RefineMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        for variant in Self::value_variants() {
            if variant.to_possible_value().unwrap().matches(s, false) {
                return Ok(*variant);
            }
        }
        Err(format!("invalid variant: {s}"))
    }
}

impl ValueEnum for CallMode {
    fn value_variants<'a>() -> &'a [Self] {
        &[CallMode::Strict, CallMode::Balanced, CallMode::Permissive]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        Some(match self {
            CallMode::Strict => PossibleValue::new("strict")
                .help("Fewer calls/cuts; higher thresholds; keeps more reads"),
            CallMode::Balanced => PossibleValue::new("balanced").help("Reasonable defaults"),
            CallMode::Permissive => PossibleValue::new("permissive")
                .help("More calls/cuts; lower thresholds; exploratory"),
        })
    }
}

impl std::fmt::Display for CallMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.to_possible_value().unwrap().get_name().fmt(f)
    }
}

impl std::str::FromStr for CallMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        for variant in Self::value_variants() {
            if variant.to_possible_value().unwrap().matches(s, false) {
                return Ok(*variant);
            }
        }
        Err(format!("invalid mode: {s}"))
    }
}
