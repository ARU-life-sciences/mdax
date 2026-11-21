//! Reference-based regression test on real MDA reads.
//!
//! `data/ci_truth/<sample>.fa.gz` holds a stratified random sample of reads
//! from a public MDA-amplified ONT library, and `<sample>.labels.tsv` labels
//! each read `foldback` or `linear` from its alignment to an independent
//! reference (synth_bench/refbased_recall_eval.py; foldback = opposite-strand
//! segments within 5 kb on the same sequence). The test runs the mdax binary
//! with the defaults used for these libraries and checks the share of true
//! foldbacks cut and of linear reads cut against fixed bounds.
//!
//! The sample is a few hundred reads, so the fold-axis support gate sees
//! almost no recurrence and nearly every detected foldback is cut: this
//! guards detection and false cuts, not the support gate. The full-library
//! check that does exercise the gate is synth_bench/truth_regress.sh (farm).
//!
//! Samples (public SRA, Lee et al. 2023, PRJNA953805):
//!   celegans_wga    C. elegans, SRR24201716, truth vs WBcel235
//!   aphelenchoides  A. besseyi APTV, SRR24201688, truth vs GCA_024699835.1 (strain AORJ)

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

struct Bounds {
    sample: &'static str,
    min_foldback_cut: f64,
    max_linear_cut: f64,
}

// Measured on this fixture (deterministic across runs), 250 + 250 reads each:
//
//   build                     celegans_wga fb / lin   aphelenchoides fb / lin
//   2 Nov 2026 (pre fold-axis)    47.6% / 1.6%            55.6% / 1.6%
//   v5 (fold-axis gate)           69.2% / 4.0%            69.6% / 2.4%
//   main 2026-11-20 (rescue, cut) 79.2% / 5.6%            81.6% / 5.6%
//
// Linear cuts are higher than on the full library (2.2% / 4.6%) because the
// support gate keeps almost nothing as real on a few hundred reads. Floors sit
// ~3 points under main (one read = 0.4 points); raise them when a change lands
// that moves the numbers for good.
const BOUNDS: &[Bounds] = &[
    Bounds { sample: "celegans_wga", min_foldback_cut: 76.0, max_linear_cut: 7.0 },
    Bounds { sample: "aphelenchoides", min_foldback_cut: 78.5, max_linear_cut: 7.0 },
];

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("data/ci_truth")
}

fn read_labels(p: &Path) -> HashMap<String, String> {
    std::fs::read_to_string(p)
        .expect("labels")
        .lines()
        .filter_map(|l| {
            let mut f = l.split('\t');
            Some((f.next()?.to_string(), f.next()?.to_string()))
        })
        .collect()
}

/// read id -> decision, from mdax's report TSV (reads without a row were not detected).
fn read_decisions(p: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(p).expect("report");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap_or_else(|| panic!("no column {name}"));
    let (id_c, dec_c) = (col("read_id"), col("decision"));
    let mut out = HashMap::new();
    for l in lines {
        let f: Vec<&str> = l.split('\t').collect();
        let id = f[id_c].split_whitespace().next().unwrap_or("").to_string();
        // a read can have several rows (recursive cuts); any artefact row means it was cut
        let e = out.entry(id).or_insert_with(|| f[dec_c].to_string());
        if f[dec_c] == "artefact" {
            *e = "artefact".to_string();
        }
    }
    out
}

fn run_sample(b: &Bounds) -> Option<(f64, f64, usize, usize)> {
    let dir = data_dir();
    let reads = dir.join(format!("{}.fa.gz", b.sample));
    let labels = dir.join(format!("{}.labels.tsv", b.sample));
    if !reads.exists() || !labels.exists() {
        eprintln!("SKIP {}: fixture not present", b.sample);
        return None;
    }
    let tmp = std::env::temp_dir().join(format!("mdax_truth_{}_{}", b.sample, std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let report = tmp.join("report.tsv");
    let st = Command::new(env!("CARGO_BIN_EXE_mdax"))
        .arg(&reads)
        .args(["--mode", "permissive", "--refine-mode", "ont", "-t", "4"])
        .arg("-o")
        .arg(tmp.join("out.fa"))
        .arg("-r")
        .arg(&report)
        .output()
        .expect("run mdax");
    assert!(st.status.success(), "mdax failed: {}", String::from_utf8_lossy(&st.stderr));

    let truth = read_labels(&labels);
    let dec = read_decisions(&report);
    std::fs::remove_dir_all(&tmp).ok();

    let share = |class: &str| {
        let ids: Vec<&String> = truth.iter().filter(|(_, t)| *t == class).map(|(q, _)| q).collect();
        let cut = ids.iter().filter(|q| dec.get(q.as_str()).map(|d| d == "artefact").unwrap_or(false)).count();
        (100.0 * cut as f64 / ids.len().max(1) as f64, ids.len())
    };
    let (fb, n_fb) = share("foldback");
    let (lin, n_lin) = share("linear");
    eprintln!(
        "{}: true foldbacks cut {fb:.1}% (n={n_fb}, floor {}), linear cut {lin:.1}% (n={n_lin}, ceiling {})",
        b.sample, b.min_foldback_cut, b.max_linear_cut
    );
    Some((fb, lin, n_fb, n_lin))
}

#[test]
fn truth_regression_public_mda_ont() {
    let mut failures = Vec::new();
    for b in BOUNDS {
        if let Some((fb, lin, _, _)) = run_sample(b) {
            if fb < b.min_foldback_cut {
                failures.push(format!("{}: foldback cut {fb:.1}% < {}", b.sample, b.min_foldback_cut));
            }
            if lin > b.max_linear_cut {
                failures.push(format!("{}: linear cut {lin:.1}% > {}", b.sample, b.max_linear_cut));
            }
        }
    }
    assert!(failures.is_empty(), "truth regression failed:\n{}", failures.join("\n"));
}
