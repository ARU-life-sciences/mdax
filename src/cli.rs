// a simple cli for fos

// use the Builder API from clap
// let's add one free positional argument

use crate::utils::RefineMode;
use clap::{Arg, ArgMatches, Command, ValueEnum, builder::PossibleValue, value_parser};
use std::path::PathBuf;

pub fn build_cli() -> ArgMatches {
    let c = Command::new("fos")
        .version(clap::crate_version!())
        .author("Max Carter-Brown")
        .about("Remove inverted duplication chimeras, fast")
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
                .help("Minimum evidence span (bp) for calling a foldback/concatemer")
                .long("min-span")
                .default_value("2000")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("min_delta")
                .help("Minimum repeat offset (bp) to ignore local tandem repeats")
                .long("min-delta")
                .default_value("2000")
                .value_parser(value_parser!(usize)),
        )
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
                .default_value("200")
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
                .help("ONT refinement bandwidth as fraction of arm (e.g. 0.2)")
                .long("max-ed-rate")
                .default_value("0.25")
                .value_parser(value_parser!(f32)),
        )
        // Minimum number of minimizer matchpoints required to support a foldback candidate.
        // Higher values reduce false positives but require more coverage of the palindrome.
        .arg(
            Arg::new("min_matches")
                .help("Minimum minimizer matches supporting a foldback candidate")
                .long("min-matches")
                .default_value("20")
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
                .value_parser(value_parser!(RefineMode)),
        );

    c.get_matches()
}

// Can also be derived with feature flag `derive`
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
