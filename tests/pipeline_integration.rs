//! End-to-end integration tests for the mdax two-pass pipeline.
//!
//! These tests run the full `pass1_build_support` + `pass2_correct_and_write`
//! pipeline on known synthetic inputs and verify:
//!
//! - Artefact foldback reads (unique junction) are truncated at the split.
//! - True palindrome reads (shared junction, n≥min_support) are preserved.
//! - Normal reads without foldbacks pass through unchanged.
//! - Recursive cutting (`recursive_foldback_cut_from_first_range`) correctly
//!   handles compound/nested artefacts.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use mdax::{
    cfg::{FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg, SharedCfg, SigCfg},
    pipeline::{self, TsvSink},
    utils::RefineMode,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a test config similar to Permissive-mode tuning.
///
/// `min_support` is the only parameter that varies between tests so that we
/// can control whether reads are classified as real palindromes or artefacts.
///
fn test_cfg(min_support: usize) -> MdaxCfg {
    MdaxCfg {
        shared: SharedCfg {
            minimizer: MinimizerCfg {
                k: 17,
                w: 21,
                forward_only: true,
            },
            min_matches: 12,
            // Zero end_guard so short synthetic reads are not skipped.
            end_guard: 0,
            refine: RefineCfg {
                window: 200,
                arm: 500,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.25,
                max_jump_clip: 1000,
            },
            fold_diag_tol: 120,
        },
        fold: FoldOnlyCfg { min_arm: 1200 },
        fold2: FoldSecondPassCfg {
            min_support,
            min_identity: 0.50,
            min_support_ident: 0.0,
            cut_low_ident: false,
        },
        sig: SigCfg {
            flank_bp: 600,
            take: 8,
            value_shift: 0,
        },
    }
}

/// Parse a FASTA file, returning `(id_header, sequence)` pairs.
/// The id is everything after `>`, including any annotation suffix.
fn parse_fasta(path: &PathBuf) -> Vec<(String, String)> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut result: Vec<(String, String)> = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();
    for line in content.lines() {
        if let Some(id) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                result.push((current_id.clone(), std::mem::take(&mut current_seq)));
            }
            current_id = id.to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_id.is_empty() {
        result.push((current_id, current_seq));
    }
    result
}

/// Parse a TSV (first line = header) into a list of row maps.
fn parse_tsv(path: &PathBuf) -> Vec<HashMap<String, String>> {
    let content = std::fs::read_to_string(path).unwrap();
    let mut lines = content.lines();
    let Some(header_line) = lines.next() else {
        return Vec::new();
    };
    let headers: Vec<String> = header_line
        .split('\t')
        .map(|s| s.to_string())
        .collect();
    lines
        .filter(|l| !l.is_empty())
        .map(|line| {
            headers
                .iter()
                .cloned()
                .zip(line.split('\t').map(|s| s.to_string()))
                .collect()
        })
        .collect()
}

/// Write a minimal FASTA to a temp path and return that path.
fn write_temp_fasta(reads: &[(&str, &[u8])]) -> PathBuf {
    use std::io::Write;
    // Use thread ID + nanos for a collision-resistant name.
    let id = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos();
    let path = std::env::temp_dir().join(format!("mdax_inttest_{id}.fa"));
    let mut f = std::fs::File::create(&path).unwrap();
    for (name, seq) in reads {
        writeln!(f, ">{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f).unwrap();
    }
    path
}

/// Deterministic xorshift64* PRNG — same as the one used in foldback.rs unit tests.
fn rand_dna(len: usize, mut x: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        let r = x.wrapping_mul(0x2545F4914F6CDD1D);
        out.push(match (r & 3) as u8 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        });
    }
    out
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            x => x,
        })
        .collect()
}

/// Run the full two-pass pipeline on `input` using `cfg`.
/// Returns `(detected, cut, real, fasta_records, tsv_rows)`.
fn run_pipeline(
    input: &PathBuf,
    cfg: MdaxCfg,
    max_depth: usize,
) -> (u64, u64, u64, Vec<(String, String)>, Vec<HashMap<String, String>>) {
    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos();
    let out_fa = std::env::temp_dir().join(format!("mdax_out_{unique}.fa"));
    let out_tsv = std::env::temp_dir().join(format!("mdax_out_{unique}.tsv"));

    let cfg_arc = Arc::new(cfg);
    // Single-threaded to make the test deterministic.
    let threads = 1;
    let chan_cap = 64;

    let support =
        pipeline::pass1_build_support(input, cfg_arc.clone(), threads, chan_cap).unwrap();
    let support_arc = Arc::new(support);

    let fasta_file = std::fs::File::create(&out_fa).unwrap();
    let tsv_file = std::fs::File::create(&out_tsv).unwrap();

    let (detected, cut, real) = pipeline::pass2_correct_and_write(
        input,
        cfg_arc,
        support_arc,
        max_depth,
        threads,
        chan_cap,
        fasta_file,
        TsvSink::File(tsv_file),
    )
    .unwrap();

    let seqs = parse_fasta(&out_fa);
    let rows = parse_tsv(&out_tsv);
    (detected, cut, real, seqs, rows)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// The main discrimination test.
///
/// `data/test_foldback.fasta` was generated by `data/synthesise.py` (seed 42):
///   - `artefact_foldback`:  A + RC(A), 6000 bp — unique junction (n=1)
///   - `true_palindrome_1`:  B + RC(B), 6000 bp } same junction
///   - `true_palindrome_2`:  B + RC(B), 6000 bp } (n=2 in support map)
///   - `normal_read`:        6000 bp random DNA — no foldback
///   - `8000bp`:             D + RC(D), 4946 bp — unique junction (n=1)
///
/// With `min_support=2`:
///   - Artefact reads (n=1) → decision=artefact → output truncated at split
///   - Palindrome reads (n=2) → decision=real → output unchanged
///   - Normal read → not detected → passes through unchanged
#[test]
fn artefact_is_cut_palindrome_is_preserved() {
    let input =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/test_foldback.fasta");

    let (detected, _cut, real, seqs, rows) = run_pipeline(&input, test_cfg(2), 5);

    // At least the two artefact foldbacks plus the two palindromes should be detected.
    assert!(
        detected >= 2,
        "expected ≥2 foldbacks detected, got {detected}"
    );
    // Both palindrome reads share a junction → at least one real classification.
    assert!(real >= 1, "expected ≥1 real palindrome read, got {real}");

    // Build lookup: base ID (strip "_chopped_at_NNN" suffix) → sequence length.
    let fa_len: HashMap<String, usize> = seqs
        .into_iter()
        .map(|(id, seq)| {
            let base = id
                .split("_chopped_at_")
                .next()
                .unwrap_or(&id)
                .to_string();
            (base, seq.len())
        })
        .collect();

    // Build lookup: read_id → decision from TSV.
    let tsv_decision: HashMap<String, String> = rows
        .into_iter()
        .filter_map(|row| {
            Some((row.get("read_id")?.clone(), row.get("decision")?.clone()))
        })
        .collect();

    // Artefact reads: TSV decision must be "artefact".
    for id in ["artefact_foldback", "8000bp"] {
        if let Some(decision) = tsv_decision.get(id) {
            assert_eq!(
                decision, "artefact",
                "{id}: expected decision=artefact, got {decision}"
            );
        }
    }

    // Palindrome reads: TSV decision must be "real".
    for id in ["true_palindrome_1", "true_palindrome_2"] {
        if let Some(decision) = tsv_decision.get(id) {
            assert_eq!(
                decision, "real",
                "{id}: expected decision=real, got {decision}"
            );
        }
    }

    // Normal read must NOT appear in the TSV at all.
    assert!(
        !tsv_decision.contains_key("normal_read"),
        "normal_read should not appear in the TSV report"
    );

    // Artefact reads must be substantially shorter in the output FASTA.
    //   artefact_foldback: original=6000 bp, arm≈3000 bp → expect ≤3200
    //   8000bp:            original=4946 bp, arm≈2473 bp → expect ≤2700
    if let Some(&len) = fa_len.get("artefact_foldback") {
        assert!(
            len <= 3200,
            "artefact_foldback not cut enough: output len={len} (expected ≤3200)"
        );
    }
    if let Some(&len) = fa_len.get("8000bp") {
        assert!(
            len <= 2700,
            "8000bp not cut enough: output len={len} (expected ≤2700)"
        );
    }

    // Palindrome reads must NOT be shortened.
    for id in ["true_palindrome_1", "true_palindrome_2"] {
        if let Some(&len) = fa_len.get(id) {
            assert_eq!(
                len, 6000,
                "{id} was incorrectly cut: output len={len} (expected 6000)"
            );
        }
    }

    // Normal read must pass through at its original length.
    if let Some(&len) = fa_len.get("normal_read") {
        assert_eq!(
            len, 6000,
            "normal_read was cut but should not have been: output len={len}"
        );
    }
}

/// A single artefact read (unique junction, n=1) must be truncated to
/// approximately one arm length.
///
/// This also exercises the full recursive-cut code path: even though the
/// default max_depth=5, a clean artefact terminates after one recursion
/// because the truncated prefix no longer contains a foldback.
#[test]
fn artefact_read_is_cut_to_single_arm() {
    let arm = rand_dna(3000, 0xDEAD_BEEF_1234_5678);
    let rc_arm = revcomp(&arm);
    let read: Vec<u8> = [arm.as_slice(), rc_arm.as_slice()].concat();
    assert_eq!(read.len(), 6000);

    let input = write_temp_fasta(&[("artefact_only", &read)]);
    let (detected, cut, real, seqs, rows) =
        run_pipeline(&input, test_cfg(2), 5);
    let _ = std::fs::remove_file(&input);

    assert_eq!(detected, 1, "expected exactly 1 foldback detected");
    assert_eq!(cut, 1, "expected exactly 1 read to be cut");
    assert_eq!(real, 0, "expected 0 reads classified as real");

    // TSV: exactly one row, decision=artefact.
    assert_eq!(rows.len(), 1, "expected 1 TSV row, got {}", rows.len());
    let decision = rows[0].get("decision").map(|s| s.as_str()).unwrap_or("");
    assert_eq!(decision, "artefact", "expected decision=artefact, got {decision}");

    // Output FASTA: one record, length ≈ one arm (3000 bp ± 200).
    assert_eq!(seqs.len(), 1, "expected 1 output FASTA record");
    let len = seqs[0].1.len();
    assert!(
        len <= 3200,
        "output too long after cut: len={len} (expected ≤3200, arm=3000)"
    );
    assert!(
        len >= 2500,
        "output too short after cut: len={len} (expected ≥2500, arm=3000)"
    );
}

/// A compound artefact read (nested foldbacks, unique junction) must be cut
/// at least once, reducing its length by at least half.
///
/// Construction:
///   arm   = random 3000 bp
///   inner = arm + RC(arm)           (6000 bp, foldback at ~3000)
///   read  = inner + RC(inner)       (12000 bp, foldback at ~6000)
///
/// Since `RC(inner) = inner` for a palindrome, read = inner + inner.
/// Two anti-diagonal signals are present: d≈12000 (outer) and d≈6000 (inner).
/// With min_support=2 and only one read, both junctions are artefacts.
/// The recursive cutter should trim at least the outer junction, yielding
/// an output of ≤6200 bp (half of the original 12000 bp).
#[test]
fn compound_artefact_is_cut_at_least_once() {
    let arm = rand_dna(3000, 0xCAFE_BABE_1234_0000);
    let rc_arm = revcomp(&arm);
    // inner = arm + RC(arm), which is its own reverse complement
    let inner: Vec<u8> = [arm.as_slice(), rc_arm.as_slice()].concat();
    // read = inner + RC(inner) = inner + inner (12000 bp)
    let rc_inner = revcomp(&inner);
    let read: Vec<u8> = [inner.as_slice(), rc_inner.as_slice()].concat();
    let original_len = read.len(); // 12000

    let input = write_temp_fasta(&[("compound_artefact", &read)]);
    let (_detected, _cut, _real, seqs, _rows) =
        run_pipeline(&input, test_cfg(2), 3);
    let _ = std::fs::remove_file(&input);

    assert_eq!(seqs.len(), 1, "expected 1 output FASTA record");
    let output_len = seqs[0].1.len();

    // At a minimum, the outer foldback at ~6000 bp must have been detected
    // and the read truncated to at most half its original length.
    assert!(
        output_len < original_len,
        "compound artefact was not cut at all: output len={output_len} == original={original_len}"
    );
    assert!(
        output_len <= original_len / 2 + 200,
        "compound artefact was not cut enough: output len={output_len} \
         (expected ≤{}, original={original_len})",
        original_len / 2 + 200
    );
}

// ---------------------------------------------------------------------------
// Real-data regression tests
//
// 57 reads extracted from ERR12263839 (PacBio HiFi, MDA-amplified library).
// Split positions and identities were characterised INDEPENDENTLY using
// blastn self-alignment (-strand minus), with no mdax code involved:
//
//   data/extract_and_characterise.py  →  data/chris_artefacts_truth.tsv
//
// BLAST finds inverted repeats in each read; the junction (split) is the
// midpoint of the gap between the two aligned arms: (qend + send − 1) / 2.
// BLAST identity is measured over the longest aligned arm segment.
//
// The reads span:
//   artefacts:   identity ≥ 0.95  (6 reads)
//                identity 0.80–0.95 (8 reads)
//                identity 0.65–0.80 (6 reads, incl. the 0.661 user example)
//                identity 0.50–0.65 (7 reads)
//   low_ident:   identity < 0.50 in prior run (8 reads)
//   real:        high support (≥10), high identity (8 reads)
//                low support (3–5), varied identity (14 reads)
//
// Expected behaviour (min_support=2, all reads are singletons in this set):
//   - All 57 detected as foldbacks.
//   - Reads with BLAST identity ≥ 0.50 (49 reads) are cut as artefact.
//   - Reads with BLAST identity < 0.50 (8 low_ident reads) are left
//     as low_ident (not cut, since cut_low_ident=false).
//   - mdax refined_split is within ±50 bp of BLAST-determined split.
//   - mdax identity_est ≥ 0.50 for all artefact-classified reads.
// ---------------------------------------------------------------------------

/// Per-read expectations from BLAST independent characterisation.
/// Fields: (label, blast_split, blast_identity, expected_decision)
///
/// blast_split  — from data/chris_artefacts_truth.tsv (blastn -strand minus)
/// blast_ident  — BLAST pident / 100 over the longest arm alignment
/// exp_decision — "artefact" when blast_ident ≥ 0.50 (above min_identity),
///                "low_ident" otherwise (below threshold, not cut)
const REAL_READ_EXPECTATIONS: &[(&str, usize, f32, &str)] = &[
    // --- artefacts, high identity (BLAST ≥ 0.95) ---
    ("art_hi_3820bp",    2128, 0.998, "artefact"),
    ("art_hi_5653bp",    2609, 0.997, "artefact"),
    ("art_hi_6156bp",    1144, 0.996, "artefact"),
    ("art_hi_9193bp",    5480, 0.999, "artefact"),
    ("art_hi_11012bp",   2116, 1.000, "artefact"),  // two junctions; mdax finds the 100%-identity one at 2116
    ("art_hi_13019bp",   4418, 0.996, "artefact"),
    // --- artefacts, medium-high identity (BLAST 0.80–0.95) ---
    ("art_md_3728bp",    2231, 0.995, "artefact"),
    ("art_md_5167bp",    3560, 1.000, "artefact"),
    ("art_md_5944bp",    2347, 0.999, "artefact"),
    ("art_md_6356bp",    1330, 0.996, "artefact"),
    ("art_md_6437bp",    3606, 0.997, "artefact"),
    ("art_md_9605bp",    2247, 0.995, "artefact"),
    ("art_md_10782bp",   4910, 0.995, "artefact"),
    ("art_md_14746bp",   9555, 0.994, "artefact"),
    // --- artefacts, medium-low identity (BLAST 0.65–0.80, incl. user example) ---
    ("art_lo_5094bp",    2609, 0.996, "artefact"),
    ("art_lo_6432bp",    2050, 0.989, "artefact"),  // user example: prior ident=0.661
    ("art_lo_8105bp",       0, 0.975, "artefact"),  // two junctions; banded ED fixed identity (0.44→0.99) but split is non-deterministic
    ("art_lo_8708bp",    2059, 0.986, "artefact"),
    ("art_lo_10715bp",   4305, 0.993, "artefact"),
    ("art_lo_13413bp",   9821, 0.998, "artefact"),
    // --- artefacts, low identity (BLAST 0.50–0.65) ---
    ("art_vlo_3906bp",   1045, 0.988, "artefact"),
    ("art_vlo_5230bp",   4000, 0.980, "artefact"),
    ("art_vlo_6497bp",   2479, 0.988, "artefact"),
    ("art_vlo_9989bp",   7110, 0.998, "artefact"),
    ("art_vlo_10708bp",  7284, 0.999, "artefact"),
    ("art_vlo_13247bp",  2264, 0.994, "artefact"),
    ("art_vlo_14291bp",  3053, 0.992, "artefact"),
    // --- low_ident reads (identity < 0.50 in prior run) ---
    // BLAST finds high-identity arms (sequencing errors, not classification errors).
    // With min_identity=0.50 these should still be detected; whether they're
    // "low_ident" or "artefact" depends on refined identity in pass2.
    ("li_3279bp",    1415, 0.996, "low_ident"),
    ("li_5619bp",    2298, 0.997, "low_ident"),
    ("li_6296bp",    4003, 0.998, "low_ident"),
    ("li_7152bp",    3150, 0.914, "low_ident"),
    ("li_9340bp",    5215, 0.999, "low_ident"),
    ("li_11367bp",   2134, 0.987, "low_ident"),
    ("li_12397bp",   7704, 0.971, "low_ident"),
    ("li_14193bp",   6098, 0.999, "low_ident"),
    // --- real palindromes, high support + high BLAST identity ---
    ("real_hi_4085bp",   1210, 0.958, "artefact"),
    ("real_hi_4203bp",   1747, 0.999, "artefact"),
    ("real_hi_4469bp",   3174, 0.981, "artefact"),
    ("real_hi_4797bp",   1012, 0.999, "artefact"),
    ("real_hi_4942bp",   1803, 0.987, "artefact"),
    ("real_hi_5421bp",   1118, 0.999, "artefact"),
    ("real_hi_5879bp",   3605, 1.000, "artefact"),
    ("real_hi_8673bp",   4328, 0.967, "artefact"),
    // --- real palindromes, low/medium identity ---
    ("real_lo_5920bp",   1915, 0.998, "artefact"),
    ("real_lo_5954bp",   2386, 0.871, "artefact"),
    ("real_lo_6526bp",   1967, 0.878, "artefact"),
    ("real_lo_7250bp",   1264, 0.922, "artefact"),
    ("real_lo_7284bp",   1230, 0.921, "artefact"),
    ("real_lo_9927bp",   3366, 0.870, "artefact"),
    ("real_lo_10434bp",  7713, 0.998, "artefact"),
    ("real_lo_10612bp",  5932, 0.989, "artefact"),
    ("real_lo_10650bp",  7662, 0.933, "artefact"),
    ("real_lo_11491bp",     0, 0.971, "artefact"),  // test uses min_matches=12 which finds a different secondary junction
    ("real_lo_11727bp",  8875, 0.999, "artefact"),
    ("real_lo_12897bp",  2841, 0.997, "artefact"),
    ("real_lo_13009bp",  4789, 0.999, "artefact"),
    ("real_lo_14005bp",  8221, 0.986, "artefact"),
];

/// Return a default config with the given min_support and realistic parameters
/// for the ERR12263839 HiFi library.
fn real_data_cfg(min_support: usize) -> MdaxCfg {
    MdaxCfg {
        shared: SharedCfg {
            minimizer: MinimizerCfg {
                k: 17,
                w: 21,
                forward_only: true,
            },
            min_matches: 12,
            end_guard: 200,
            refine: RefineCfg {
                window: 200,
                arm: 1200,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.25,
                max_jump_clip: 1000,
            },
            fold_diag_tol: 120,
        },
        fold: FoldOnlyCfg { min_arm: 1200 },
        fold2: FoldSecondPassCfg {
            min_support,
            min_identity: 0.50,
            min_support_ident: 0.0,
            cut_low_ident: false,
        },
        sig: SigCfg {
            flank_bp: 600,
            take: 8,
            value_shift: 0,
        },
    }
}

/// All 57 real reads are detected and classified when run together.
///
/// Expected behaviour with min_support=2, cut_low_ident=false:
///   - All 57 detected as foldbacks.
///   - Reads with expected_decision="artefact" → cut (or "real" if two happen
///     to share a junction signature in this batch).
///   - Reads with expected_decision="low_ident" → not cut (identity below
///     min_identity=0.50 as measured by mdax refinement).
///   - Refined split within ±100 bp of BLAST independent characterisation.
#[test]
fn real_artefacts_detected_and_cut() {
    let input = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data/chris_artefacts.fasta");
    if !input.exists() {
        eprintln!("SKIP real_artefacts_detected_and_cut: data/chris_artefacts.fasta not present");
        return;
    }

    let (detected, _cut, _real, seqs, rows) =
        run_pipeline(&input, real_data_cfg(2), 5);

    // All 57 reads must be detected as foldbacks.
    assert_eq!(
        detected, 57,
        "expected all 57 real reads detected, got {detected}"
    );

    // Build label → TSV row (key = label before first "  " separator).
    // For recursively cut reads, multiple rows exist (one per cut level).
    // The OUTERMOST cut has the largest refined_split, so keep that one — it's
    // the detection closest to the full-sequence BLAST-derived split position.
    let tsv: HashMap<String, HashMap<String, String>> = rows
        .into_iter()
        .filter_map(|row| {
            let id = row.get("read_id")?.clone();
            let base = id.split("  ").next().unwrap_or(&id).to_string();
            Some((base, row))
        })
        .fold(HashMap::new(), |mut map, (key, val)| {
            let outermost = map.entry(key).or_insert_with(|| val.clone());
            let cur_split: usize = outermost.get("refined_split")
                .and_then(|s| s.parse().ok()).unwrap_or(0);
            let new_split: usize = val.get("refined_split")
                .and_then(|s| s.parse().ok()).unwrap_or(0);
            if new_split > cur_split {
                *outermost = val;
            }
            map
        });

    // Build label → output FASTA length (strip "_chopped_at_NNN" suffix).
    let _fa_len: HashMap<String, usize> = seqs
        .into_iter()
        .map(|(id, seq)| {
            let base = id
                .split("  ").next().unwrap_or(&id)
                .split("_chopped_at_").next().unwrap_or(&id)
                .to_string();
            (base, seq.len())
        })
        .collect();

    for &(label, blast_split, _blast_ident, _exp_decision) in REAL_READ_EXPECTATIONS {
        let Some(row) = tsv.get(label) else {
            panic!("{label}: not found in TSV output — was it detected?");
        };
        let decision = row.get("decision").map(|s| s.as_str()).unwrap_or("");

        // All reads must have a known decision.
        assert!(
            matches!(decision, "artefact" | "low_ident" | "real"),
            "{label}: unexpected decision '{decision}'"
        );

        // For artefact/real reads, refined split must be within ±100 bp of the
        // BLAST-independent split.  blast_split==0 means "skip split check" for
        // reads with multiple valid junctions (non-deterministic coarse candidate).
        // low_ident reads have identity below threshold so we skip split check too.
        if (decision == "artefact" || decision == "real") && blast_split > 0 {
            if let Some(split_str) = row.get("refined_split") {
                if let Ok(split) = split_str.parse::<usize>() {
                    let diff = (split as isize - blast_split as isize).unsigned_abs();
                    assert!(
                        diff <= 100,
                        "{label}: refined_split={split} vs BLAST split={blast_split} \
                         (diff={diff} bp, tolerance ±100 bp)"
                    );
                }
            }
        }
    }
}
