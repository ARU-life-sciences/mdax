use anyhow::{Result, bail};
use calm_io::stderrln;
use std::path::PathBuf;
use std::sync::Arc;

use mdax::{
    cfg::{
        CallMode, FairnessParams, FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg,
        SharedCfg, SigCfg,
    },
    cli, elog,
    pipeline,
    utils::{RefineMode, compact_histogram},
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

    elog!(
        "WELCOME",
        "This is `mdax`, a foldback detector v{}",
        env!("CARGO_PKG_VERSION")
    );
    elog!(
        "WELCOME",
        "Developed by Max Brown at Anglia Ruskin University <max.carter-brown@aru.ac.uk>"
    );

    // for testing
    let fairness_baseline = args.get_flag("fairness_baseline")
        || std::env::var("MDAX_FAIRNESS")
            .map(|v| v == "1")
            .unwrap_or(false);

    // pass1-only diagnostic mode: report signature-cluster stats, skip pass2
    // entirely (no cutting, no FASTA/TSV written). See --qc-only long help.
    let qc_only = args.get_flag("qc_only");

    // inputs and outputs
    let fasta = args
        .get_one::<PathBuf>("input")
        .expect("not allowed no input");
    // `output` is only required when not --qc-only (enforced by clap's
    // required_unless_present); safe to unwrap once we're past the
    // qc_only early-return below.
    let output = args.get_one::<PathBuf>("output");
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
    // TODO: include these in fairness?
    let sig_flank_bp = *args
        .get_one::<usize>("sig_flank_bp")
        .unwrap_or(&t.sig_flank_bp);
    let sig_take = *args.get_one::<usize>("sig_take").unwrap_or(&t.sig_take);
    let sig_value_shift = *args
        .get_one::<u8>("sig_value_shift")
        .unwrap_or(&t.sig_value_shift);

    let mut k = *args.get_one::<usize>("k").expect("defaulted k of 17");
    let mut w = *args.get_one::<usize>("w").expect("defaulted w of 21");
    // w must be odd, so error out here if it's not
    if w % 2 == 0 {
        anyhow::bail!("Window size `w` must be odd, got {}", w);
    }
    let mut forward_only = *args.get_one::<bool>("forward_only").unwrap();

    let mut refine_mode = *args.get_one::<RefineMode>("refine_mode").unwrap();

    // For ONT mode, automatically use smaller k/w when the user hasn't
    // overridden them explicitly.  With k=17 and ~10% ONT error, each 17-mer
    // has only ~2% probability of surviving in both arms of the foldback
    // (cross-match probability ≈ (0.9)^34), yielding ~2 expected matches per
    // 2 kb arm — well below any practical min_matches threshold.
    // k=11, w=11 raises this to ~14 expected matches (cross-match ≈ (0.9)^22),
    // recovering most of the missed reads at negligible precision cost.
    if refine_mode == RefineMode::ONT {
        let k_explicit = args.value_source("k")
            == Some(clap::parser::ValueSource::CommandLine);
        let w_explicit = args.value_source("w")
            == Some(clap::parser::ValueSource::CommandLine);
        if !k_explicit {
            k = 11;
        }
        if !w_explicit {
            w = 11;
        }
    }
    let mut refine_window = *args.get_one::<usize>("refine_window").unwrap();
    let mut refine_arm = *args.get_one::<usize>("refine_arm").unwrap();
    let mut max_ed_rate = *args.get_one::<f32>("max_ed_rate").unwrap();
    let max_jump_clip = *args.get_one::<usize>("max_jump_clip").unwrap();

    let mut end_guard = *args.get_one::<usize>("end_guard").unwrap_or(&t.end_guard);
    let mut fold_diag_tol = *args.get_one::<i32>("fold_diag_tol").unwrap();

    let ont_k_note = if refine_mode == RefineMode::ONT {
        let k_explicit = args.value_source("k")
            == Some(clap::parser::ValueSource::CommandLine);
        let w_explicit = args.value_source("w")
            == Some(clap::parser::ValueSource::CommandLine);
        match (k_explicit, w_explicit) {
            (false, false) => format!(" [ONT auto: k={k}, w={w}]"),
            (false, _)     => format!(" [ONT auto: k={k}]"),
            (_, false)     => format!(" [ONT auto: w={w}]"),
            (true, true)   => String::new(),
        }
    } else {
        String::new()
    };

    for legacy in ["sig_flank_bp", "sig_take", "sig_value_shift"] {
        if args.value_source(legacy) == Some(clap::parser::ValueSource::CommandLine) {
            elog!(
                "INIT",
                "note: --{} no longer has any effect (support is now counted per fold axis, not per locus signature)",
                legacy.replace('_', "-")
            );
        }
    }

    let msg = [
        fmt_param("min_matches", min_matches, t.min_matches),
        fmt_param("min_span", min_span, t.min_span),
        fmt_param("min_identity", min_identity, t.min_identity),
        fmt_param("min_support", min_support, t.min_support),
        fmt_param("end_guard", end_guard, t.end_guard),
        fmt_param("max_depth", max_depth, t.max_depth),
        format!("k={k}  w={w}{ont_k_note}"),
        format!("refine_mode={refine_mode:?}  refine_arm={refine_arm}  refine_window={refine_window}"),
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
    let mut cfg = MdaxCfg {
        shared: SharedCfg {
            minimizer: MinimizerCfg { k, w, forward_only },
            min_matches,
            end_guard,
            refine: RefineCfg {
                window: refine_window,
                arm: refine_arm,
                mode: refine_mode,
                max_ed_rate,
                max_jump_clip,
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
            min_identity,
            min_support_ident: 0.0,
            cut_low_ident: !args.get_flag("keep_low_ident"),
        },
        sig: SigCfg {
            flank_bp: sig_flank_bp,
            take: sig_take,
            value_shift: sig_value_shift,
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

    // choose threads + channel capacity
    let available_threads = std::thread::available_parallelism().map(|n| n.get())?;

    let threads = *args.get_one::<usize>("threads").unwrap();
    if threads > available_threads {
        bail!(
            "Requested {} threads, but only {} are available",
            threads,
            available_threads
        );
    }

    // keep this modest to avoid buffering too much sequence
    let chan_cap = 512usize;

    elog!("INIT", "Starting with mode: {mode:?}");
    elog!(
        "INIT",
        "Threads used/available: {threads}/{available_threads}; Specified chunk size: {}",
        chan_cap
    );
    elog!("INIT", "Parameter Summary:\n  {msg}");

    let cfg_arc = Arc::new(cfg.clone());

    // PASS 1
    let (support, total_reads) = pipeline::pass1_build_support(&fasta, cfg_arc.clone(), threads, chan_cap)?;
    // Per-read junction support: how many reads (including itself) fold
    // about the same axis as each detected foldback.
    elog!(
        "SUPPORT",
        "foldbacks with a junction key: {} (of {} input reads); shared keys: {}",
        support.reads,
        total_reads,
        support.len()
    );
    elog!("SUPPORT HISTOGRAM", "reads by number of reads sharing their fold axis:");
    let histogram: std::collections::HashMap<usize, usize> = support
        .n_hist
        .iter()
        .enumerate()
        .filter(|&(n, &c)| n > 0 && c > 0)
        .map(|(n, &c)| (n, c as usize))
        .collect();
    compact_histogram(&histogram, 30)?;

    // Density-aware threshold. The number of *unrelated* artefacts that
    // share a fold axis by chance grows with how densely junctions cover the
    // genome, so a fixed min_support is too lax for deep libraries (measured
    // against a reference: at ~140x, 15% of true artefacts were kept as
    // "real"; at ~1500x, half). Estimate the chance rate from the data and
    // raise the bar until under 1% of chance co-locations clear it. This
    // only ever raises the threshold, and an explicit --min-support (or the
    // frozen fairness baseline) disables it.
    let min_support_explicit =
        args.value_source("min_support") == Some(clap::parser::ValueSource::CommandLine);
    let lambda = support.chance_lambda();
    if !fairness_baseline && !min_support_explicit {
        let adaptive = support.adaptive_min_support(min_support, 0.01);
        if adaptive != min_support {
            elog!(
                "SUPPORT CALIBRATION",
                "min_support {} -> {} (chance co-location rate lambda={:.3}, plus measured tail)",
                min_support,
                adaptive,
                lambda
            );
            min_support = adaptive;
            cfg.fold2.min_support = min_support;
        }
    }

    let real_reads_est: u64 = support.n_hist.iter().skip(min_support.min(support.n_hist.len())).sum();
    elog!(
        "QC SUMMARY",
        "mode={mode:?} lambda={:.3} min_support={} max_support_n={} reads_at_or_above_min_support={}",
        lambda,
        min_support,
        support.max_n,
        real_reads_est,
    );

    if qc_only {
        elog!(
            "QC ONLY",
            "--qc-only set: stopping after pass1, no reads were cut or written"
        );
        return Ok(());
    }

    // PASS 2
    // PASS 2 writers
    let output = output.expect("required_unless_present(\"qc_only\") guarantees this when qc_only is false");
    let fasta_out = output.clone();
    let tsv_out = report.clone();

    let fasta_file = std::fs::File::create(&fasta_out)
        .map_err(|e| anyhow::anyhow!("failed to create {}: {e}", fasta_out.display()))?;

    // report can be "-" (stdout)
    let tsv_sink: pipeline::TsvSink = if tsv_out.to_string_lossy() == "-" {
        pipeline::TsvSink::Stdout
    } else {
        let f = std::fs::File::create(&tsv_out)
            .map_err(|e| anyhow::anyhow!("failed to create {}: {e}", tsv_out.display()))?;
        pipeline::TsvSink::File(f)
    };

    let support_arc = Arc::new(support);

    // Pass1 didn't need fold2.min_support (only pass2's real/artefact gate
    // does), so it's safe that `cfg_arc` above was built before calibration --
    // but pass2 must see the calibrated value, so build a fresh Arc here.
    let cfg_arc_pass2 = Arc::new(cfg.clone());

    let (foldback_reads, cut_reads, real_reads) = pipeline::pass2_correct_and_write(
        &fasta,
        cfg_arc_pass2,
        support_arc.clone(),
        max_depth,
        threads,
        chan_cap,
        fasta_file,
        tsv_sink,
    )?;

    elog!("SUMMARY", "Foldback reads detected: {}", foldback_reads);
    elog!(
        "SUMMARY",
        "Foldback reads cut (>=1 recursion): {}",
        cut_reads
    );
    elog!("SUMMARY", "Foldback reads classified real: {}", real_reads);

    Ok(())
}
