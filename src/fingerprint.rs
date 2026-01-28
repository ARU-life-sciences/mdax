// a junction fingerprint module for foldback structures

use gxhash::{GxHasher, HashMap, HashMapExt};
use std::hash::{Hash, Hasher};

use crate::{
    cfg::{MdaxCfg, SharedCfg},
    foldback, io,
    minimizer::{self, sampled_minimizers_into},
    scratch::{self, SigScratch},
    utils::revcomp_in_place,
};

// Compute a reproducible fingerprint for a foldback junction at `split`.
// Returns None if flanks exceed bounds or not enough minimizers.
pub fn foldback_signature(
    seq: &[u8],
    split: usize,
    shared: &SharedCfg,
    flank_bp: usize,
    take: usize,
    scratch: &mut SigScratch,
) -> Option<u64> {
    if split < flank_bp || split + flank_bp > seq.len() {
        return None;
    }

    let left = &seq[split - flank_bp..split];
    let right = &seq[split..split + flank_bp];

    // canonical minimizers for signature stability
    // (avoid cloning the whole config each call)
    let mcfg = crate::cfg::MinimizerCfg {
        k: shared.minimizer.k,
        w: shared.minimizer.w,
        forward_only: false,
    };

    // left minimizers -> scratch.val_l
    sampled_minimizers_into(left, &mcfg, &mut scratch.pos_l, &mut scratch.val_l);

    // right minimizers, but RC orientation
    scratch.right_rc.clear();
    scratch.right_rc.extend_from_slice(right);
    crate::utils::revcomp_in_place(&mut scratch.right_rc);

    sampled_minimizers_into(
        &scratch.right_rc,
        &mcfg,
        &mut scratch.pos_r,
        &mut scratch.val_r,
    );

    if scratch.val_l.len() < take || scratch.val_r.len() < take {
        return None;
    }

    // sort values and take smallest `take`
    scratch.val_l.sort_unstable();
    scratch.val_r.sort_unstable();

    // Build combined list without allocating a new vec each time
    scratch.combined.clear();
    scratch.combined.extend_from_slice(&scratch.val_l[..take]);
    scratch.combined.extend_from_slice(&scratch.val_r[..take]);
    scratch.combined.sort_unstable();

    let mut h = GxHasher::default();
    scratch.combined.hash(&mut h);
    Some(h.finish())
}

// gather the support statistics for a given junction fingerprint

#[derive(Debug, Clone, Default)]
pub struct SupportStats {
    pub n: usize,
    pub min_split: usize,
    pub max_split: usize,
    pub mean_ident: f32,
}

impl SupportStats {
    pub fn new(split: usize, ident: f32) -> Self {
        Self {
            n: 1,
            min_split: split,
            max_split: split,
            mean_ident: ident,
        }
    }

    pub fn update(&mut self, split: usize, ident: f32) {
        self.n += 1;
        self.min_split = self.min_split.min(split);
        self.max_split = self.max_split.max(split);
        // incremental mean
        self.mean_ident += (ident - self.mean_ident) / self.n as f32;
    }

    pub fn split_span(&self) -> usize {
        self.max_split.saturating_sub(self.min_split)
    }
}

// TODO: should this be in foldback.rs?
pub fn foldback_pass1_support<P: AsRef<std::path::Path>>(
    fasta_path: P,
    cfg: &MdaxCfg,
    scratch: &mut scratch::FoldScratch,
) -> anyhow::Result<HashMap<u64, SupportStats>> {
    let mut reader = io::open_fasta_reader(fasta_path)?;
    let mut support: HashMap<u64, SupportStats> = HashMap::new();
    let mut sig_scratch = SigScratch::default();

    while let Some(rec) = reader.next() {
        let r = rec?;
        let seq = r.seq();

        let fb = match foldback::detect_foldback(&seq, &cfg.shared, &cfg.fold, scratch) {
            Some(x) => x,
            None => continue,
        };

        let refined = match foldback::refine_breakpoint(&seq, fb.split_pos, &cfg.shared) {
            Some(rf) => rf,
            None => continue,
        };

        if refined.identity_est < cfg.fold2.min_identity {
            continue;
        }

        // TODO: maybe make these parameters configurable, but probably should not need to
        const DEFAULT_FLANK_BP: usize = 1000;
        const DEFAULT_TAKE: usize = 12;

        let sig = match foldback_signature(
            &seq,
            refined.split_pos,
            &cfg.shared,
            DEFAULT_FLANK_BP,
            DEFAULT_TAKE,
            &mut sig_scratch,
        ) {
            Some(s) => s,
            None => continue,
        };

        // enter into the map
        support
            .entry(sig)
            .and_modify(|st| st.update(refined.split_pos, refined.identity_est))
            .or_insert_with(|| SupportStats::new(refined.split_pos, refined.identity_est));
    }

    Ok(support)
}

pub fn is_real_foldback(sig: u64, support: &HashMap<u64, SupportStats>, cfg: &MdaxCfg) -> bool {
    let Some(st) = support.get(&sig) else {
        return false;
    };
    st.n >= cfg.fold2.min_support && st.split_span() <= cfg.fold2.split_tol_bp
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cfg::{
        ConcatOnlyCfg, FoldOnlyCfg, FoldSecondPassCfg, MinimizerCfg, RefineCfg, SharedCfg,
    };
    use crate::utils::RefineMode;

    use noodles::fasta::{
        self as fasta,
        record::{Definition, Sequence},
    };
    use std::fs::File;
    use std::io::BufWriter;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn pseudo_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let b = match (x >> 62) & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            out.push(b);
        }
        out
    }

    fn repeat_bytes(pat: &[u8], times: usize) -> Vec<u8> {
        let mut v = Vec::with_capacity(pat.len() * times);
        for _ in 0..times {
            v.extend_from_slice(pat);
        }
        v
    }

    fn make_clean_foldback(left: &[u8], join: &[u8]) -> (Vec<u8>, usize) {
        let mut rc = left.to_vec();
        crate::utils::revcomp_in_place(&mut rc);

        let mut s = Vec::with_capacity(left.len() + join.len() + rc.len());
        s.extend_from_slice(left);
        s.extend_from_slice(join);
        let split = left.len() + join.len();
        s.extend_from_slice(&rc);
        (s, split)
    }

    fn shared_for_sig(k: usize, w: usize, min_matches: usize) -> SharedCfg {
        SharedCfg {
            minimizer: MinimizerCfg {
                k,
                w,
                forward_only: true,
            },
            min_matches,
            end_guard: 0,
            refine: RefineCfg {
                window: 0,
                arm: 0,
                mode: RefineMode::HiFi,
                max_ed_rate: 0.0,
            },
            fold_diag_tol: 100,
            concat_diag_tol: 100,
        }
    }

    fn tmp_fasta_path(prefix: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("{prefix}_{now}.fa"));
        p
    }

    fn write_fasta(path: &std::path::Path, records: &[(&str, &[u8])]) {
        let f = File::create(path).expect("create temp fasta");
        let mut w = fasta::io::Writer::new(BufWriter::new(f));
        for (id, seq) in records {
            let definition = Definition::new(*id, None);
            let sequence = Sequence::from(seq.to_vec());
            let record = fasta::Record::new(definition, sequence);
            w.write_record(&record).expect("write record");
        }
    }

    fn cfg_for_pass1() -> MdaxCfg {
        // Make it easy to detect + refine + fingerprint.
        MdaxCfg {
            shared: SharedCfg {
                minimizer: MinimizerCfg {
                    k: 9,
                    w: 11,
                    forward_only: true,
                },
                min_matches: 5,
                end_guard: 0,
                refine: RefineCfg {
                    window: 50,
                    arm: 300, // refinement arm; independent of fingerprint flank=1000
                    mode: RefineMode::HiFi,
                    max_ed_rate: 0.25,
                },
                fold_diag_tol: 50,
                concat_diag_tol: 200,
            },
            fold: FoldOnlyCfg {
                min_arm: 200, // span evidence for foldback detector
            },
            concat: ConcatOnlyCfg::default(),
            fold2: FoldSecondPassCfg {
                min_support: 3,
                split_tol_bp: 100,
                min_identity: 0.60,
            },
        }
    }

    #[test]
    fn foldback_pass1_support_finds_one_cluster() {
        let cfg = cfg_for_pass1();

        // Need >= 1000bp on each side of junction for foldback_signature().
        // So left must be >= 1000, and overall read long enough.
        let left = pseudo_dna(2500, 42);
        let (seq, _true_split) = make_clean_foldback(&left, b"");

        // 3 identical foldback reads -> should become one signature with n=3
        // plus a decoy random read that should be ignored
        let decoy = pseudo_dna(6000, 999);

        let path = tmp_fasta_path("mdax_fold2_pass1");
        write_fasta(
            &path,
            &[("r1", &seq), ("r2", &seq), ("r3", &seq), ("decoy", &decoy)],
        );
        let mut scratch = scratch::FoldScratch::new();

        let support = foldback_pass1_support(&path, &cfg, &mut scratch).expect("pass1 support");

        // remove file
        let _ = std::fs::remove_file(&path);

        assert_eq!(
            support.len(),
            1,
            "expected exactly one foldback fingerprint cluster"
        );

        let (_sig, st) = support.iter().next().unwrap();
        assert_eq!(st.n, 3, "should have 3 supporting reads");
        assert_eq!(st.split_span(), 0, "identical reads should cluster tightly");
        assert!(st.mean_ident >= cfg.fold2.min_identity);
    }

    #[test]
    fn signature_is_deterministic_same_input() {
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 40); // 640bp
        let (seq, split) = make_clean_foldback(&left, b"");

        let shared = shared_for_sig(9, 11, 5);

        let sig1 = foldback_signature(&seq, split, &shared, 200, 8).expect("sig1");
        let sig2 = foldback_signature(&seq, split, &shared, 200, 8).expect("sig2");

        // should be the same signature for same input
        assert_eq!(sig1, sig2);
    }

    #[test]
    fn signature_matches_for_same_junction_across_reads() {
        // Two reads that share same left flank content around the junction
        // but have extra prefix/suffix elsewhere should still match as long as
        // we fingerprint using the same split and flank region content.
        let left = repeat_bytes(b"ACGTTGCAACGTTGCA", 50); // 800bp
        let (seq1, split1) = make_clean_foldback(&left, b"");
        let (mut seq2, split2) = make_clean_foldback(&left, b"");

        // add a prefix (doesn't affect the junction *content* if we also shift split)
        let prefix = repeat_bytes(b"GATTACA", 20);
        let mut with_prefix = Vec::with_capacity(prefix.len() + seq2.len());
        with_prefix.extend_from_slice(&prefix);
        with_prefix.extend_from_slice(&seq2);
        seq2 = with_prefix;

        let shared = shared_for_sig(9, 11, 5);

        // fingerprint around the true junction positions in each read
        let sig1 = foldback_signature(&seq1, split1, &shared, 200, 8).expect("sig1");
        let sig2 = foldback_signature(&seq2, split2 + prefix.len(), &shared, 200, 8).expect("sig2");

        assert_eq!(sig1, sig2);
    }

    #[test]
    fn signature_changes_if_junction_flanks_change() {
        let left1 = repeat_bytes(b"ACGTTGCAACGTTGCA", 40);
        let left2 = repeat_bytes(b"TTTTGGGGCCCCAAAA", 40);

        let (seq1, split1) = make_clean_foldback(&left1, b"");
        let (seq2, split2) = make_clean_foldback(&left2, b"");

        let shared = shared_for_sig(9, 11, 5);

        let sig1 = foldback_signature(&seq1, split1, &shared, 200, 8).expect("sig1");
        let sig2 = foldback_signature(&seq2, split2, &shared, 200, 8).expect("sig2");

        assert_ne!(sig1, sig2);
    }

    #[test]
    fn supportstats_updates_and_span() {
        let mut st = SupportStats::new(1000, 0.8);
        st.update(1010, 0.6);
        st.update(990, 0.9);

        assert_eq!(st.n, 3);
        assert_eq!(st.min_split, 990);
        assert_eq!(st.max_split, 1010);
        assert_eq!(st.split_span(), 20);

        // mean should be in range
        assert!(st.mean_ident >= 0.6 && st.mean_ident <= 0.9);
    }

    #[test]
    fn is_real_foldback_respects_threshold_and_tolerance() {
        use crate::cfg::{ConcatOnlyCfg, FoldOnlyCfg, MdaxCfg};

        // minimal cfg stub
        let cfg = MdaxCfg {
            shared: shared_for_sig(9, 11, 5),
            fold: FoldOnlyCfg { min_arm: 20 },
            concat: ConcatOnlyCfg {
                min_span: 2000,
                min_delta: 2000,
                cross_frac: 0.8,
            },
            fold2: crate::cfg::FoldSecondPassCfg {
                min_support: 3,
                split_tol_bp: 50,
                min_identity: 0.6,
            },
        };

        let mut support: HashMap<u64, SupportStats> = HashMap::new();
        let sig = 123456_u64;

        // only 2 supports => not real
        support.insert(sig, {
            let mut st = SupportStats::new(1000, 0.8);
            st.update(1005, 0.8);
            st
        });
        assert!(!is_real_foldback(sig, &support, &cfg));

        // 3 supports but too spread => not real
        support.insert(sig, {
            let mut st = SupportStats::new(1000, 0.8);
            st.update(1100, 0.8);
            st.update(1200, 0.8);
            st
        });
        assert!(!is_real_foldback(sig, &support, &cfg));

        // 3 supports and tight cluster => real
        support.insert(sig, {
            let mut st = SupportStats::new(1000, 0.8);
            st.update(1010, 0.8);
            st.update(1020, 0.8);
            st
        });
        assert!(is_real_foldback(sig, &support, &cfg));
    }
}
