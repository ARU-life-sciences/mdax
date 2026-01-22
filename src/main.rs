mod cfg;
mod cli;
mod concatemer;
mod foldback;
mod io;
mod minimizer;
mod utils;

use anyhow::Result;
use calm_io::stderrln;
use std::io::Write;
use std::path::PathBuf;

use crate::{
    cfg::{ConcatOnlyCfg, FoldOnlyCfg, MdaxCfg, MinimizerCfg, RefineCfg, SharedCfg},
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

    let fasta_reader = io::open_fasta_reader(fasta)?;
    let mut fasta_writer = io::open_fasta_writer(output)?;

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
    };

    // output tsv
    let mut tsv = io::open_tsv_writer(report)?;
    writeln!(
        tsv,
        "read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est"
    )?;

    let mut num_foldbacks = 0;
    let mut num_concatemers = 0;

    for record in fasta_reader.records() {
        let r = record?;
        let seq = r.seq();
        let id = r.id();
        let len = seq.len();

        let fold = foldback::detect_foldback(seq, &cfg.shared, &cfg.fold);
        let concat = concatemer::detect_concatemer(seq, &cfg.shared, &cfg.concat);

        // ---- FASTA output (foldback correction) ----
        if let Some(fb) = fold.as_ref() {
            // refine split if possible; otherwise fall back to coarse
            let split = match foldback::refine_breakpoint(seq, fb.split_pos, &cfg.shared) {
                Some(rf) => rf.split_pos,
                None => fb.split_pos,
            };

            // guard against pathological split values
            let split = split.min(seq.len());

            // suffix the ID
            let out_id = io::split_id(id, 1);
            io::write_fasta_record(&mut fasta_writer, &out_id, &seq[..split])?;
        } else {
            // no foldback => passthrough
            io::write_fasta_record(&mut fasta_writer, id, seq)?;
        }

        // ---- Foldback row ----
        if let Some(fb) = fold {
            let foldback::FoldBreakpoint {
                split_pos,
                score,
                matches,
                span,
            } = fb;

            let refined = foldback::refine_breakpoint(seq, split_pos, &cfg.shared);

            num_foldbacks += 1;

            if let Some(rf) = refined {
                writeln!(
                    tsv,
                    "{id}\t{len}\tfoldback\t1\t{}\t{}\t\t{}\t{}\t\t\t{}\t{}\t{:.3}",
                    split_pos, rf.split_pos, matches, span, score, rf.score, rf.identity_est
                )?;
            } else {
                writeln!(
                    tsv,
                    "{id}\t{len}\tfoldback\t1\t{}\t\t\t{}\t{}\t\t\t{}\t\t",
                    split_pos, matches, span, score
                )?;
            }
        }

        // ---- Concatemer row ----
        if let Some(cb) = concat {
            let concatemer::ConcatBreakpoint {
                split_pos,
                score,
                delta,
                matches,
                span,
                cross_frac,
                p2_span,
            } = cb;

            let refined =
                concatemer::refine_concatemer_breakpoint(seq, split_pos, delta, &cfg.shared);

            num_concatemers += 1;

            if let Some(rf) = refined {
                writeln!(
                    tsv,
                    "{id}\t{len}\tconcatemer\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{:.3}",
                    split_pos,
                    rf.split_pos,
                    delta,
                    matches,
                    span,
                    p2_span,
                    cross_frac,
                    score,
                    rf.score,
                    rf.identity_est
                )?;
            } else {
                writeln!(
                    tsv,
                    "{id}\t{len}\tconcatemer\t1\t{}\t\t{}\t{}\t{}\t{}\t{:.3}\t{}\t\t",
                    split_pos, delta, matches, span, p2_span, cross_frac, score
                )?;
            }
        }
    }

    stderrln!("Total foldbacks detected: {}", num_foldbacks)?;
    stderrln!("Total concatemers detected: {}", num_concatemers)?;

    Ok(())
}
