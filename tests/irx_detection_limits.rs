//! Regression tests for irx detection limits found by stress-testing
//! (validation/stress.py in the irx survey repo):
//! - an N-gap inside an IR's spacer or arm must not hide or halve the IR;
//! - nested IRs (inner IR inside the outer one's spacer) must both be reported;
//! - more than 4 IRs within one window must all be reported (top-K cap).

use std::path::PathBuf;
use std::process::Command;

struct Rng(u64);
impl Rng {
    fn base(&mut self) -> u8 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        b"ACGT"[(self.0 >> 62) as usize]
    }
    fn dna(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.base()).collect()
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Planted IR: (left arm start, left arm end, right arm start, right arm end).
type Truth = (usize, usize, usize, usize);

struct Builder {
    rng: Rng,
    seq: Vec<u8>,
    truth: Vec<Truth>,
}

impl Builder {
    fn new(seed: u64) -> Self {
        Self { rng: Rng(seed), seq: Vec::new(), truth: Vec::new() }
    }
    fn random(&mut self, n: usize) {
        let d = self.rng.dna(n);
        self.seq.extend(d);
    }
    fn gap(&mut self, n: usize) {
        self.seq.extend(std::iter::repeat(b'N').take(n));
    }
    /// Left arm, then `spacer` (closure appends it), then revcomp(left arm).
    fn ir(&mut self, arm: usize, spacer: impl FnOnce(&mut Self)) {
        let left = self.rng.dna(arm);
        let la0 = self.seq.len();
        self.seq.extend(&left);
        let la1 = self.seq.len();
        spacer(self);
        let ra0 = self.seq.len();
        self.seq.extend(revcomp(&left));
        self.truth.push((la0, la1, ra0, self.seq.len()));
    }
}

struct Call {
    arms: Truth,
    arm_len: usize,
    locus_id: String,
    locus_n_irs: usize,
    pass: String,
}

fn run(name: &str, seq: &[u8]) -> Vec<Call> {
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("{name}.fa"));
    let mut rec = format!(">{name}\n").into_bytes();
    rec.extend_from_slice(seq);
    rec.push(b'\n');
    std::fs::write(&path, rec).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_irx"))
        .arg(&path)
        .args(["--threads", "1"])
        .output()
        .expect("failed to run irx");
    assert!(out.status.success(), "irx failed: {}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let header: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
    let col = |name: &str| header.iter().position(|h| h.trim_start_matches('#') == name).unwrap();
    let (c_arm, c_locus, c_n, c_pass) =
        (col("arm_len"), col("locus_id"), col("locus_n_irs"), col("pass"));
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), header.len(), "row/header column mismatch: {l}");
            let n = |i: usize| f[i].parse::<usize>().unwrap();
            Call {
                arms: (n(11), n(12), n(13), n(14)),
                arm_len: n(c_arm),
                locus_id: f[c_locus].to_string(),
                locus_n_irs: n(c_n),
                pass: f[c_pass].to_string(),
            }
        })
        .collect()
}

fn found<'a>(calls: &'a [Call], t: &Truth) -> Option<&'a Call> {
    let ov = |a0: usize, a1: usize, b0: usize, b1: usize| a0 < b1 && b0 < a1;
    calls.iter().find(|c| ov(c.arms.0, c.arms.1, t.0, t.1) && ov(c.arms.2, c.arms.3, t.2, t.3))
}

#[test]
fn n_gap_in_spacer_does_not_hide_ir() {
    let mut b = Builder::new(1);
    b.random(60_000);
    b.ir(3_000, |b| {
        b.random(2_000);
        b.gap(1_000);
        b.random(2_000);
    });
    b.random(60_000);
    let calls = run("n_spacer", &b.seq);
    assert!(found(&calls, &b.truth[0]).is_some(), "IR with an N-gap in its spacer was missed");
}

#[test]
fn n_gap_in_arm_does_not_halve_arm() {
    let mut b = Builder::new(2);
    b.random(60_000);
    let left = b.rng.dna(5_000);
    let la0 = b.seq.len();
    b.seq.extend(&left[..2_000]);
    b.gap(500);
    b.seq.extend(&left[2_500..]);
    let la1 = b.seq.len();
    b.random(5_000);
    let ra0 = b.seq.len();
    b.seq.extend(revcomp(&left));
    let truth = (la0, la1, ra0, b.seq.len());
    b.random(60_000);
    let calls = run("n_arm", &b.seq);
    let c = found(&calls, &truth).expect("IR with an N-gap in its arm was missed");
    assert!(c.arm_len >= 4_500, "arm split by the gap: {} bp of 5,000", c.arm_len);
}

#[test]
fn nested_irs_with_different_centres_both_reported() {
    let mut b = Builder::new(3);
    b.random(60_000);
    b.ir(5_000, |b| {
        b.random(2_000);
        b.ir(3_000, |b| b.random(2_000));
        b.random(8_000);
    });
    b.random(60_000);
    let calls = run("nested", &b.seq);
    for (i, t) in b.truth.iter().enumerate() {
        assert!(found(&calls, t).is_some(), "nested IR {i} of 2 missed");
    }
}

#[test]
fn more_than_four_irs_in_one_window_all_reported() {
    let mut b = Builder::new(4);
    b.random(60_000);
    for _ in 0..8 {
        b.ir(3_000, |b| b.random(5_000));
        b.random(10_000);
    }
    b.random(60_000);
    let calls = run("dense", &b.seq);
    let n = b.truth.iter().filter(|t| found(&calls, t).is_some()).count();
    assert_eq!(n, 8, "only {n} of 8 densely packed IRs reported");
}

#[test]
fn inverted_copies_of_one_sequence_form_one_locus() {
    // X ... rc(X) ... rc(X): two IRs sharing X as their left arm -> one locus of 2.
    let mut b = Builder::new(5);
    b.random(60_000);
    let x = b.rng.dna(3_000);
    let x0 = b.seq.len();
    b.seq.extend(&x);
    let x1 = b.seq.len();
    b.random(8_000);
    let c1 = b.seq.len();
    b.seq.extend(revcomp(&x));
    let c1e = b.seq.len();
    b.random(20_000);
    let c2 = b.seq.len();
    b.seq.extend(revcomp(&x));
    let c2e = b.seq.len();
    b.random(60_000);
    let calls = run("copies", &b.seq);
    let a = found(&calls, &(x0, x1, c1, c1e)).expect("X / first copy not reported");
    let c = found(&calls, &(x0, x1, c2, c2e)).expect("X / second copy not reported");
    assert_eq!(a.locus_id, c.locus_id, "IRs sharing an arm got different loci");
    assert!(a.locus_n_irs >= 2);
}

#[test]
fn large_footprint_ir_found_by_large_pass() {
    // 20 kb arms, 360 kb spacer: footprint 400 kb, beyond the local pass (complete
    // to 200 kb, blind past 250 kb) but within the large pass (complete to 1 Mb).
    let mut b = Builder::new(6);
    b.random(100_000);
    b.ir(20_000, |b| b.random(360_000));
    b.random(100_000);
    let calls = run("large_fp", &b.seq);
    let c = found(&calls, &b.truth[0]).expect("400 kb-footprint IR missed");
    assert_eq!(c.pass, "large");
    assert!(c.arm_len >= 19_500, "arm {} of 20,000", c.arm_len);
}
