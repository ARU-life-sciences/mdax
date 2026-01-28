mod cfg;
mod cli;
mod concatemer;
mod fingerprint;
mod foldback;
mod io;
mod minimizer;
mod scratch;
mod utils;

use anyhow::Result;
use calm_io::stderrln;
use std::io::Write;
use std::path::PathBuf;

use crate::{
    cfg::{
        ConcatOnlyCfg, FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg, SharedCfg,
    },
    fingerprint::{foldback_pass1_support, is_real_foldback},
    foldback::recursive_foldback_cut,
    scratch::SigScratch,
    utils::RefineMode,
};

fn main() -> Result<()> {
    let args = cli::build_cli();

    let fasta = args
        .get_one::<PathBuf>("input")
        .expect("not allowed no input");
    let output = args
        .get_one::<PathBuf>("output")
        .expect("not allowed no output");
    let report = args.get_one::<PathBuf>("report").unwrap();
    let k = *args.get_one::<usize>("k").expect("defaulted k of 17");
    let w = *args.get_one::<usize>("w").expect("defaulted w of 21");
    let min_span = *args
        .get_one::<usize>("min_span")
        .expect("defaulted min_span of 2000");
    let min_matches = *args
        .get_one::<usize>("min_matches")
        .expect("defaulted min_matches of 40");
    let refine_mode = *args
        .get_one::<RefineMode>("refine_mode")
        .unwrap_or(&RefineMode::HiFi);
    // also all defaulted
    let min_delta = *args.get_one::<usize>("min_delta").unwrap();
    let end_guard = *args.get_one::<usize>("end_guard").unwrap();
    let refine_window = *args.get_one::<usize>("refine_window").unwrap();
    let refine_arm = *args.get_one::<usize>("refine_arm").unwrap();
    let max_ed_rate = *args.get_one::<f32>("max_ed_rate").unwrap();
    let fold_diag_tol = *args.get_one::<i32>("fold_diag_tol").unwrap();
    let concat_diag_tol = *args.get_one::<i32>("concat_diag_tol").unwrap();
    // forward-only minimizers toggle
    let forward_only = *args.get_one::<bool>("forward_only").unwrap();

    // cfg here:
    let cfg = MdaxCfg {
        shared: SharedCfg {
            minimizer: MinimizerCfg { k, w, forward_only },
            min_matches,
            end_guard,
            refine: RefineCfg {
                window: refine_window,
                arm: refine_arm,
                mode: refine_mode,
                max_ed_rate,
            },
            fold_diag_tol,
            concat_diag_tol,
        },
        fold: FoldOnlyCfg {
            // This is the foldback “arm evidence span” threshold.
            // If you want a dedicated fold option, add a separate CLI flag.
            min_arm: min_span,
        },
        concat: ConcatOnlyCfg {
            min_span,
            min_delta,
            cross_frac: 0.80, // optionally CLI-control this later
        },
        fold2: FoldSecondPassCfg {
            // TODO: sensible defaults + add to cli
            min_support: 3,
            split_tol_bp: 50,
            min_identity: 0.6,
        },
    };
    let mut fold_scratch = scratch::FoldScratch::new();

    // first pass
    let support = foldback_pass1_support(&fasta, &cfg, &mut fold_scratch)?;
    stderrln!("Foldback support map entries: {}", support.len())?;

    let mut fasta_reader = io::open_fasta_reader(fasta)?;
    let mut fasta_writer = io::open_fasta_writer(output)?;

    // output tsv
    let mut tsv = io::open_tsv_writer(report)?;
    writeln!(
        tsv,
        "read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est\tsupport_n\tsupport_span\tdecision"
    )?;

    // tweakables (eventually CLI or cfg.fold2)
    let max_depth = 5usize;

    let mut num_foldbacks = 0;
    let mut num_foldbacks_cut = 0usize;
    let mut num_foldbacks_kept_real = 0usize;

    let mut sig_scratch = SigScratch::default();

    while let Some(record) = fasta_reader.next() {
        let r = record?;
        let seq = r.seq();
        let id = std::str::from_utf8(r.id())?;
        let len = seq.len();

        let fold = foldback::detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut fold_scratch);

        // ---- FASTA output (foldback correction) ----
        if let Some(_) = fold.as_ref() {
            // Decide whether to cut recursively
            let kept = recursive_foldback_cut(
                &seq,
                &cfg,
                &support,
                max_depth,
                &mut fold_scratch,
                &mut sig_scratch,
            )?;

            // i.e. if we got a recursive cut
            if kept.len() < seq.len() {
                // we actually cut at least once
                num_foldbacks_cut += 1;
                let out_id = io::split_id(id, 1);
                io::write_fasta_record(&mut fasta_writer, &out_id, kept)?;
            } else {
                io::write_fasta_record(&mut fasta_writer, id, &seq)?;
            }
        } else {
            // no foldback => passthrough
            io::write_fasta_record(&mut fasta_writer, id, &seq)?;
        }

        // ---- Foldback TSV row ----
        if let Some(fb) = fold {
            num_foldbacks += 1;

            let refined = foldback::refine_breakpoint(&seq, fb.split_pos, &cfg.shared);

            if let Some(rf) = refined {
                // signature + decision if possible
                let (decision, support_n, support_span) =
                    if rf.identity_est >= cfg.fold2.min_identity {
                        if let Some(sig) = fingerprint::foldback_signature(
                            &seq,
                            rf.split_pos,
                            &cfg.shared,
                            1000,
                            12,
                            &mut sig_scratch,
                        ) {
                            if let Some(st) = support.get(&sig) {
                                let real = is_real_foldback(sig, &support, &cfg);
                                if real {
                                    num_foldbacks_kept_real += 1;
                                }
                                (
                                    if real { "real" } else { "artefact" },
                                    st.n,
                                    st.split_span(),
                                )
                            } else {
                                ("unknown", 0, 0)
                            }
                        } else {
                            ("unknown", 0, 0)
                        }
                    } else {
                        ("low_ident", 0, 0)
                    };

                writeln!(
                    tsv,
                    "{id}\t{len}\tfoldback\t1\t{}\t{}\t\t{}\t{}\t\t\t{}\t{}\t{:.3}\t{}\t{}\t{}",
                    fb.split_pos,
                    rf.split_pos,
                    fb.matches,
                    fb.span,
                    fb.score,
                    rf.score,
                    rf.identity_est,
                    support_n,
                    support_span,
                    decision
                )?;
            } else {
                writeln!(
                    tsv,
                    "{id}\t{len}\tfoldback\t1\t{}\t\t\t{}\t{}\t\t\t{}\t\t\t0\t0\tunknown",
                    fb.split_pos, fb.matches, fb.span, fb.score
                )?;
            }
        }
    }

    stderrln!("Total foldbacks detected: {}", num_foldbacks)?;
    stderrln!("Total foldbacks cut (>=1 recursion): {}", num_foldbacks_cut)?;
    stderrln!(
        "Total foldbacks classified real (if enabled): {}",
        num_foldbacks_kept_real
    )?;

    Ok(())
}
