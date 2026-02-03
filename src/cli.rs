//! CLI definition for `mdax` chimeric read detection and clean-up.
//!
//! This module is responsible for:
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
        .about("Remove inverted duplication chimeras, fast")
        // ----------------------------
        // Inputs / outputs
        // ----------------------------
        //
        // Positional input FASTA file (optionally gzipped).
        // This is the read library to scan for foldback/palindrome artefacts.
        .arg(
            Arg::new("input")
                .help("Input file")
                .value_parser(value_parser!(PathBuf))
                .required(true)
                .index(1),
        )
        // Output TSV report file path.
        .arg(
            Arg::new("report")
                .help("Write a TSV report of detected concatemers (use '-' for stdout)")
                .short('r')
                .long("report")
                .value_parser(value_parser!(PathBuf))
                .default_value("-"),
        )
        // Output FASTA file path.
        // When splitting is enabled, each input record may produce multiple output records
        // with suffixed IDs (e.g. `read123|p1`, `read123|p2`).
        // I suspect we want only the first "half" of the read.
        .arg(
            Arg::new("output")
                .help("Output FASTA (optionally .gz not supported yet)")
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
                .help("Minimizer k-mer size")
                .short('k')
                .default_value("17")
                .value_parser(value_parser!(usize)),
        )
        // Minimizer window size (in k-mers).
        // Each window chooses one minimizer; larger values sample fewer positions.
        // Typical: 10–30.
        // NOTE: must be odd!
        .arg(
            Arg::new("w")
                .help("Minimizer window size")
                .short('w')
                .default_value("21")
                .value_parser(value_parser!(usize)),
        )
        // Minimum palindrome arm length (bp) required to call a foldback.
        // This is measured as the span of matched minimizer positions on one side of the junction.
        // Prevents calling short/random palindromic matches.
        .arg(
            Arg::new("min_span")
                .help("Override: minimum evidence span (bp) for calling a foldback/concatemer")
                .long("min-span")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_identity")
              .long("min-identity")
              .help("Override: minimum refinement identity for foldback second pass")
              .value_parser(value_parser!(f32)),
        )
        .arg(
            Arg::new("min_support")
              .long("min-support")
              .help("Override: minimum support reads for 'real' foldback")
              .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("split_tol_bp")
              .long("split-tol-bp")
              .help("Override: max split span (bp) within a signature cluster")
              .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max_depth")
              .long("max-depth")
              .help("Override: maximum recursion depth for chopping")
              .value_parser(value_parser!(usize)),
        )
        // TODO: deleted min_delta, I think this was just for concatemers
        .arg(
            Arg::new("end_guard")
                .help("Do not call breakpoints within this distance of read ends (bp)")
                .long("end-guard")
                .default_value("1000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("refine_window")
                .help("Refinement half-window around coarse split (bp)")
                .long("refine-window")
                .default_value("100")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("refine_arm")
                .help("Refinement arm length (bp)")
                .long("refine-arm")
                .default_value("200")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max_ed_rate")
                .help("ONT refinement bandwidth as fraction of arm")
                .long("max-ed-rate")
                .default_value("0.25")
                .value_parser(value_parser!(f32)),
        )
        // Minimum number of minimizer matchpoints required to support a foldback candidate.
        // Higher values reduce false positives but require more coverage of the palindrome.
        .arg(
            Arg::new("min_matches")
                .help("Override: minimum minimizer matches supporting a foldback candidate")
                .long("min-matches")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("fold_diag_tol")
                .help("Foldback anti-diagonal bucketing tolerance (bp) for pos1+pos2 clustering")
                .long("fold-diag-tol")
                .default_value("120")
                .value_parser(value_parser!(i32)),
        )
        .arg(
            Arg::new("concat_diag_tol")
                .help("Concatemer delta bucketing tolerance (bp) for (pos2-pos1) clustering")
                .long("concat-diag-tol")
                .default_value("200")
                .value_parser(value_parser!(i32)),
        )
        .arg(
            Arg::new("forward_only")
                .help("Use strand-specific forward-only minimizers (true) or canonical minimizers (false)")
                .long("forward-only")
                .default_value("true")
                .value_parser(value_parser!(bool)),
        )
        // an argument for HiFi vs ONT for refinement of breakpoints
        .arg(
            Arg::new("refine_mode")
                .help("Refine breakpoints for reads")
                .short('m')
                .long("refine-mode")
                .default_value("hifi")
                .value_parser(value_parser!(RefineMode)),
        )
        .arg(
              Arg::new("mode")
                .help("Detection/correction strictness: strict reduces calls, permissive increases calls")
                .long("mode")
                .default_value("balanced")
                .value_parser(value_parser!(CallMode)),
        )
        .arg(
            Arg::new("sig_flank_bp")
                .long("sig-flank-bp")
                .help("Window (bp) around split for match-based fingerprinting")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("sig_take")
                .long("sig-take")
                .help("Number of minimizer hashes retained for foldback fingerprint")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("fairness_baseline")
                .long("fairness-baseline")
                .hide(true) 
                .help("Developer-only: lock parameters for fair benchmarking")
                .action(clap::ArgAction::SetTrue),
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
