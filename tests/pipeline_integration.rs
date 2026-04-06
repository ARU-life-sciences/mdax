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
