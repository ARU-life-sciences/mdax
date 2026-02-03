mod cfg;
mod cli;
mod fingerprint;
mod foldback;
mod io;
mod minimizer;
mod pipeline;
mod scratch;
mod utils;

use anyhow::Result;
use calm_io::stderrln;
use std::path::PathBuf;
use std::sync::Arc;

use crate::{
    cfg::{
        CallMode, FairnessParams, FoldOnlyCfg, FoldSecondPassCfg, MdaxCfg, MinimizerCfg, RefineCfg,
        SharedCfg, SigCfg,
    },
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
    elog!(
        "WELCOME",
        "This is `mdax`, a foldback detector v{}",
        env!("CARGO_PKG_VERSION")
    );
    elog!(
        "WELCOME",
        "Developed by Max Brown at Anglia Ruskin University <max.carter-brown@aru.ac.uk>"
    );

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
    // TODO: include these in fairness?
    let sig_flank_bp = *args
        .get_one::<usize>("sig_flank_bp")
        .unwrap_or(&t.sig_flank_bp);
    let sig_take = *args.get_one::<usize>("sig_take").unwrap_or(&t.sig_take);

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

    let msg = [
        fmt_param("min_matches", min_matches, t.min_matches),
        fmt_param("min_span", min_span, t.min_span),
        fmt_param("min_identity", min_identity, t.min_identity),
        fmt_param("min_support", min_support, t.min_support),
        fmt_param("split_tol_bp", split_tol_bp, t.split_tol_bp),
        fmt_param("max_depth", max_depth, t.max_depth),
        fmt_param("sig_flank_bp", sig_flank_bp, t.sig_flank_bp),
        fmt_param("sig_take", sig_take, t.sig_take),
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
            // TODO: maybe use in cli later
            value_shift: t.sig_value_shift,
        },
    };

    // so we use hifi initial, for fast ont
    let mut shared_fast = cfg.shared.clone();
    shared_fast.refine.mode = RefineMode::HiFi;
    // shared_fast.refine.window = shared_fast.refine.window.min(50);
    // shared_fast.refine.arm    = shared_fast.refine.arm.min(120);

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
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .max(1);

    // keep this modest to avoid buffering too much sequence
    let chan_cap = 512usize;

    elog!("INIT", "Starting with mode: {mode:?}");
    elog!(
        "INIT",
        "Threads available: {threads}; Specified chunk size: {}",
        chan_cap
    );
    elog!("INIT", "Parameter Summary:\n  {msg}");

    let cfg_arc = Arc::new(cfg);

    // PASS 1
    let support = pipeline::pass1_build_support(&fasta, cfg_arc.clone(), threads, chan_cap)?;
    elog!("SUPPORT MAP", "Number of entries: {}", support.len());

    // Collect histogram of clusters into buckets based on their size
    let mut histogram = std::collections::HashMap::new();
    for st in support.values() {
        *histogram.entry(st.n).or_insert(0) += 1;
    }

    // Print histogram
    elog!("SUPPORT HISTOGRAM", "Cluster Size Distribution:");
    compact_histogram(&histogram, 30)?;

    let total_n: u64 = support.values().map(|s| s.n as u64).sum();
    elog!(
        "SUPPORT UNIQUENESS",
        "clusters={} total_n={} avg_n={:.3}",
        support.len(),
        total_n,
        (total_n as f64 / support.len().max(1) as f64),
    );

    elog!(
        "SUPPORT SANITY",
        "total_clusters={} total_support_reads={} max_n={}",
        support.len(),
        support.values().map(|s| s.n as u64).sum::<u64>(),
        support.values().map(|s| s.n).max().unwrap_or(0),
    );

    // PASS 2
    // PASS 2 writers
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

    let (num_foldbacks, num_cut, num_real) = pipeline::pass2_correct_and_write(
        &fasta,
        cfg_arc.clone(),
        support_arc.clone(),
        max_depth,
        threads,
        chan_cap,
        fasta_file,
        tsv_sink,
    )?;

    elog!("SUMMARY", "Total foldbacks detected: {}", num_foldbacks);
    elog!(
        "SUMMARY",
        "Total foldbacks cut (>=1 recursion): {}",
        num_cut
    );
    elog!(
        "SUMMARY",
        "Total foldbacks classified real (if enabled): {}",
        num_real
    );

    Ok(())
}
