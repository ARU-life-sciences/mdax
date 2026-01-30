mod cfg;
mod cli;
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
        CallMode, FairnessParams, FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg,
        SharedCfg, SigCfg,
    },
    fingerprint::{foldback_pass1_support, is_real_foldback},
    foldback::recursive_foldback_cut,
    scratch::{RefineScratch, SigScratch},
    utils::RefineMode,
};

fn fmt_param<T: std::fmt::Display + PartialEq>(
    name: &str,
    effective: T,
    mode_default: T,
) -> String {
    if effective == mode_default {
        format!("{name}={effective}")
    } else {
        format!("{name}={effective} (override; mode={mode_default})")
    }
}

fn main() -> Result<()> {
    let args = cli::build_cli();

    // for testing
    let fairness_baseline = args.get_flag("fairness_baseline")
        || std::env::var("MDAX_FAIRNESS")
            .map(|v| v == "1")
            .unwrap_or(false);

    // inputs and outputs
    let fasta = args
        .get_one::<PathBuf>("input")
        .expect("not allowed no input");
    let output = args
        .get_one::<PathBuf>("output")
        .expect("not allowed no output");
    let report = args.get_one::<PathBuf>("report").unwrap();

    // deal with the baseline mode
    let mode = *args
        .get_one::<CallMode>("mode")
        .unwrap_or(&CallMode::Balanced);
    let t = mode.tuning();

    // params overriding tuning
    let mut max_depth = *args.get_one::<usize>("max_depth").unwrap_or(&t.max_depth);
    let mut min_matches = *args
        .get_one::<usize>("min_matches")
        .unwrap_or(&t.min_matches);
    let mut min_span = *args.get_one::<usize>("min_span").unwrap_or(&t.min_span);
    let mut min_identity = *args
        .get_one::<f32>("min_identity")
        .unwrap_or(&t.min_identity);
    let mut min_support = *args
        .get_one::<usize>("min_support")
        .unwrap_or(&t.min_support);
    let mut split_tol_bp = *args
        .get_one::<usize>("split_tol_bp")
        .unwrap_or(&t.split_tol_bp);

    let mut k = *args.get_one::<usize>("k").expect("defaulted k of 17");
    let mut w = *args.get_one::<usize>("w").expect("defaulted w of 21");
    // w must be odd, so error out here if it's not
    if w % 2 == 0 {
        anyhow::bail!("Window size `w` must be odd, got {}", w);
    }
    let mut forward_only = *args.get_one::<bool>("forward_only").unwrap();

    let mut refine_mode = *args.get_one::<RefineMode>("refine_mode").unwrap();
    let mut refine_window = *args.get_one::<usize>("refine_window").unwrap();
    let mut refine_arm = *args.get_one::<usize>("refine_arm").unwrap();
    let mut max_ed_rate = *args.get_one::<f32>("max_ed_rate").unwrap();

    let mut end_guard = *args.get_one::<usize>("end_guard").unwrap();
    let mut fold_diag_tol = *args.get_one::<i32>("fold_diag_tol").unwrap();

    // TODO: override on cli later?
    let mut sig_flank_bp = t.sig_flank_bp;
    let mut sig_take = t.sig_take;

    let msg = [
        fmt_param("min_matches", min_matches, t.min_matches),
        fmt_param("min_span", min_span, t.min_span),
        fmt_param("min_identity", min_identity, t.min_identity),
        fmt_param("min_support", min_support, t.min_support),
        fmt_param("split_tol_bp", split_tol_bp, t.split_tol_bp),
        fmt_param("max_depth", max_depth, t.max_depth),
    ]
    .join("\n  ");

    let fairness = if fairness_baseline {
        Some(FairnessParams::baseline())
    } else {
        None
    };

    if let Some(f) = fairness.as_ref() {
        // detection thresholds
        min_matches = f.min_matches;
        min_span = f.min_span;
        min_identity = f.min_identity;
        min_support = f.min_support;
        split_tol_bp = f.split_tol_bp;
        max_depth = f.max_depth;

        // minimizers
        k = f.k;
        w = f.w;
        forward_only = f.forward_only;

        // refinement
        refine_mode = f.refine_mode;
        refine_window = f.refine_window;
        refine_arm = f.refine_arm;
        max_ed_rate = f.max_ed_rate;

        // guards
        end_guard = f.end_guard;
        fold_diag_tol = f.fold_diag_tol;
    }

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
        },
        fold: FoldOnlyCfg {
            // This is the foldback “arm evidence span” threshold.
            // If you want a dedicated fold option, add a separate CLI flag.
            min_arm: min_span,
        },
        fold2: FoldSecondPassCfg {
            // TODO: sensible defaults + add to cli
            min_support,
            split_tol_bp,
            min_identity,
        },
        sig: SigCfg {
            flank_bp: sig_flank_bp,
            take: sig_take,
        },
    };

    if fairness_baseline {
        stderrln!("FAIRNESS BASELINE ENABLED (developer mode)")?;
        stderrln!(
            "Locked parameters: k={} w={} forward_only={} refine={:?}",
            k,
            w,
            forward_only,
            refine_mode
        )?;
    }

    stderrln!("Mode={mode:?} effective parameters:\n  {msg}")?;

    let mut fold_scratch = scratch::FoldScratch::new();

    // first pass to store support for foldbacks
    let support = foldback_pass1_support(&fasta, &cfg, &mut fold_scratch)?;
    stderrln!("Foldback support map entries: {}", support.len())?;

    let mut n1 = 0usize;
    let mut n2 = 0usize;
    let mut n3 = 0usize;
    let mut nge4 = 0usize;
    let mut spans: Vec<usize> = Vec::new();

    for st in support.values() {
        match st.n {
            0 | 1 => n1 += 1,
            2 => n2 += 1,
            3 => n3 += 1,
            _ => nge4 += 1,
        }
        spans.push(st.split_span());
    }
    spans.sort_unstable();

    stderrln!(
        "Support clusters: total={}  n=1:{} n=2:{} n=3:{} n>=4:{}  span_median:{} span_p95:{}",
        support.len(),
        n1,
        n2,
        n3,
        nge4,
        spans.get(spans.len() / 2).copied().unwrap_or(0),
        spans.get((spans.len() * 95) / 100).copied().unwrap_or(0),
    )?;

    let mut fasta_reader = io::open_fasta_reader(fasta)?;
    let mut fasta_writer = io::open_fasta_writer(output)?;

    // output tsv
    let mut tsv = io::open_tsv_writer(report)?;
    writeln!(
        tsv,
        "read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est\tsupport_n\tsupport_span\tdecision"
    )?;

    let mut num_foldbacks = 0;
    let mut num_foldbacks_cut = 0usize;
    let mut num_foldbacks_kept_real = 0usize;

    let mut sig_scratch = SigScratch::default();
    let mut refine_scratch = &mut RefineScratch::default();

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

            let refined = foldback::refine_breakpoint(
                &seq,
                fb.split_pos,
                &cfg.shared,
                &mut fold_scratch.refine,
            );

            if let Some(rf) = refined {
                // signature + decision if possible
                let (decision, support_n, support_span) =
                    if rf.identity_est >= cfg.fold2.min_identity {
                        if let Some(sig) = fingerprint::foldback_signature(
                            &seq,
                            rf.split_pos,
                            &cfg.shared,
                            cfg.sig.flank_bp,
                            cfg.sig.take,
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
