//! `irx --min-arm` must bound the length of *each* arm, not the IR footprint.
//!
//! Regression test: mdax's coarse bin span is the whole footprint
//! (arm + spacer + arm), so a short-armed IR with a long spacer used to pass a
//! 2 kb "min_arm" gate with ~150 bp arms.

use std::path::PathBuf;
use std::process::Command;

/// Deterministic pseudo-random ACGT (simple LCG; no extra dev-dependencies).
fn random_dna(n: usize, seed: &mut u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(*seed >> 62) as usize]
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// Coordinates (left arm start, right arm end) of each planted IR.
struct Planted {
    contig_len: usize,
    short_ir: (usize, usize),
    long_ir: (usize, usize),
}

/// Background with two planted perfect IRs:
/// - short: 300 bp arms, 20 kb spacer (footprint 20.6 kb)
/// - long: 3 kb arms, 5 kb spacer
fn build_contig() -> (Vec<u8>, Planted) {
    let mut seed = 42u64;
    let mut s = random_dna(5_000, &mut seed);

    let short_start = s.len();
    let arm = random_dna(300, &mut seed);
    s.extend_from_slice(&arm);
    s.extend(random_dna(20_000, &mut seed));
    s.extend(revcomp(&arm));
    let short_end = s.len();

    s.extend(random_dna(5_000, &mut seed));

    let long_start = s.len();
    let arm = random_dna(3_000, &mut seed);
    s.extend_from_slice(&arm);
    s.extend(random_dna(5_000, &mut seed));
    s.extend(revcomp(&arm));
    let long_end = s.len();

    s.extend(random_dna(5_000, &mut seed));

    let contig_len = s.len();
    (
        s,
        Planted {
            contig_len,
            short_ir: (short_start, short_end),
            long_ir: (long_start, long_end),
        },
    )
}

struct Row {
    start: usize,
    end: usize,
    la0: usize,
    ra0: usize,
    arm_len: usize,
}

fn run_irx(fasta: &PathBuf, extra: &[&str]) -> Vec<Row> {
    let out = Command::new(env!("CARGO_BIN_EXE_irx"))
        .arg(fasta)
        .args(["--threads", "1"])
        .args(extra)
        .output()
        .expect("failed to run irx");
    assert!(out.status.success(), "irx failed: {}", String::from_utf8_lossy(&out.stderr));

    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Row {
                start: f[1].parse().unwrap(),
                end: f[2].parse().unwrap(),
                la0: f[11].parse().unwrap(),
                ra0: f[13].parse().unwrap(),
                arm_len: f[20].parse().unwrap(),
            }
        })
        .collect()
}

fn overlaps(r: &Row, (s, e): (usize, usize)) -> bool {
    r.start < e && s < r.end
}

/// Each test gets its own file: the tests in this binary run in parallel,
/// and `std::fs::write` truncates before writing, so a shared path lets one
/// test's `irx` open the file while another test has just emptied it.
fn write_fasta(tag: &str) -> (PathBuf, Planted) {
    let (seq, planted) = build_contig();
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("irx_min_arm.{tag}.fa"));
    let mut rec = b">ctg1\n".to_vec();
    rec.extend_from_slice(&seq);
    rec.push(b'\n');
    std::fs::write(&path, rec).unwrap();
    (path, planted)
}

#[test]
fn default_min_arm_rejects_short_arms_with_long_spacer() {
    let (fa, planted) = write_fasta("default");
    let rows = run_irx(&fa, &[]);

    for r in &rows {
        assert!(r.arm_len >= 2000, "row with arm_len {} < default --min-arm 2000", r.arm_len);
    }
    assert!(
        rows.iter().any(|r| overlaps(r, planted.long_ir)),
        "long-armed IR (3 kb arms) was not reported"
    );
    assert!(
        !rows.iter().any(|r| overlaps(r, planted.short_ir)),
        "short-armed IR (300 bp arms, 20 kb spacer) passed --min-arm 2000"
    );
}

#[test]
fn gff_matches_tsv_and_is_well_formed() {
    let (fa, planted) = write_fasta("gff");
    let gff_path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("irx_min_arm.gff3");
    let rows = run_irx(&fa, &["--min-arm", "200", "--gff", gff_path.to_str().unwrap()]);
    assert!(rows.len() >= 2, "expected both planted IRs, got {} rows", rows.len());

    let gff = std::fs::read_to_string(&gff_path).unwrap();
    let mut lines = gff.lines();
    assert_eq!(lines.next(), Some("##gff-version 3"));
    assert!(gff.contains(&format!("##sequence-region ctg1 1 {}", planted.contig_len)));

    let feats: Vec<Vec<&str>> = gff
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    assert!(feats.iter().all(|f| f.len() == 9), "every feature line has 9 columns");

    let parents: Vec<&Vec<&str>> = feats.iter().filter(|f| f[2] == "inverted_repeat").collect();
    assert_eq!(parents.len(), rows.len(), "one GFF parent per TSV row");

    // Sorted by start (tabix requirement).
    let starts: Vec<usize> = feats.iter().map(|f| f[3].parse().unwrap()).collect();
    assert!(starts.windows(2).all(|w| w[0] <= w[1]), "features not sorted by start");

    for p in &parents {
        let id = p[8].split(';').next().unwrap().strip_prefix("ID=").unwrap();
        let (ps, pe): (usize, usize) = (p[3].parse().unwrap(), p[4].parse().unwrap());
        let arms: Vec<&Vec<&str>> = feats
            .iter()
            .filter(|f| f[8].contains(&format!("Parent={id};")))
            .collect();
        assert_eq!(arms.len(), 2, "{id}: expected two arms");
        let (l, r) = (arms.iter().find(|f| f[6] == "+").unwrap(), arms.iter().find(|f| f[6] == "-").unwrap());
        assert_eq!(l[2], "repeat_component");
        let (ls, le): (usize, usize) = (l[3].parse().unwrap(), l[4].parse().unwrap());
        let (rs, re): (usize, usize) = (r[3].parse().unwrap(), r[4].parse().unwrap());
        assert!(ps == ls && pe == le.max(re), "{id}: parent must span the arms");
        assert!(ls <= le && rs <= re && ls <= rs);

        // GFF is 1-based closed: start = TSV la0 + 1 for the same IR.
        let tsv = rows
            .iter()
            .find(|t| t.la0.min(t.ra0) + 1 == ls)
            .unwrap_or_else(|| panic!("{id}: no TSV row with la0 = {}", ls - 1));
        assert!(p[8].contains(&format!("arm_len={};", tsv.arm_len)));
    }
}

#[test]
fn lower_min_arm_admits_short_arms() {
    let (fa, planted) = write_fasta("lower");
    let rows = run_irx(&fa, &["--min-arm", "200"]);

    for r in &rows {
        assert!(r.arm_len >= 200, "row with arm_len {} < --min-arm 200", r.arm_len);
    }
    assert!(
        rows.iter().any(|r| overlaps(r, planted.short_ir)),
        "short-armed IR not reported with --min-arm 200"
    );
}
