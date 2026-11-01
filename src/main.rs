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
    fingerprint::{MAX_PLAUSIBLE_CLUSTER_N, is_real_foldback_stats},
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

/// Weighted median of histogram keys `n`, restricted to entries with `n >= min_n`.
///
/// `histogram` maps support count -> number of support-map entries with that
/// count (as built from `SupportStats::n` after pass1). Used to estimate the
/// dataset's own "typical" recurrence depth among loci that recur at all, so
/// `min_support` can be calibrated to it instead of staying a fixed constant
/// blind to how much coverage this particular dataset actually achieved.
fn weighted_median_at_least(
    histogram: &std::collections::HashMap<usize, usize>,
    min_n: usize,
) -> Option<usize> {
    let mut entries: Vec<(usize, i64)> = histogram
        .iter()
        .filter(|&(&n, _)| n >= min_n)
        .map(|(&n, &c)| (n, c as i64))
        .collect();
    if entries.is_empty() {
        return None;
    }
    entries.sort_unstable_by_key(|&(n, _)| n);
    let total: i64 = entries.iter().map(|&(_, c)| c).sum();
    let half = (total + 1) / 2; // 1-indexed median position
    let mut cum = 0i64;
    for (n, c) in entries {
        cum += c;
        if cum >= half {
            return Some(n);
        }
    }
    None
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

    let mut end_guard = *args.get_one::<usize>("end_guard").unwrap();
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

    let msg = [
        fmt_param("min_matches", min_matches, t.min_matches),
        fmt_param("min_span", min_span, t.min_span),
        fmt_param("min_identity", min_identity, t.min_identity),
        fmt_param("min_support", min_support, t.min_support),
        fmt_param("max_depth", max_depth, t.max_depth),
        fmt_param("sig_flank_bp", sig_flank_bp, t.sig_flank_bp),
        fmt_param("sig_take", sig_take, t.sig_take),
        fmt_param("sig_value_shift", sig_value_shift, t.sig_value_shift),
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
            cut_low_ident: *args.get_one::<bool>("cut_low_ident").unwrap_or(&false),
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

    // Self-calibrate min_support to this dataset's own observed recurrence
    // level, unless the user explicitly passed --min-support (or this is the
    // fairness baseline, which must stay behaviorally frozen for legacy
    // comparison).
    //
    // min_support alone is a *fixed absolute* recurrence threshold, blind to
    // local sequencing depth. MDA coverage is notoriously uneven, so a real
    // inverted repeat sitting in a low-coverage patch of the library can
    // easily be spanned by fewer reads than a hardcoded preset value purely
    // by chance -- in which case mdax's default behaviour is to cut it in
    // every read that covers it, systematically damaging assembly
    // contiguity at exactly the loci where long-read spanning matters most.
    //
    // We only ever LOWER the threshold here, clamped to the preset/override
    // value the user already chose, so this can't make a well-covered
    // dataset more permissive than requested -- it only rescues datasets
    // where the achieved recurrence is genuinely lower than the preset
    // assumed.
    let min_support_explicit =
        args.value_source("min_support") == Some(clap::parser::ValueSource::CommandLine);

    if !fairness_baseline && !min_support_explicit {
        if let Some(typical_recurrence) = weighted_median_at_least(&histogram, 2) {
            let adaptive = typical_recurrence.clamp(2, min_support);
            if adaptive < min_support {
                elog!(
                    "SUPPORT CALIBRATION",
                    "min_support {} -> {} (dataset's typical recurrence among recurring loci is {})",
                    min_support,
                    adaptive,
                    typical_recurrence
                );
                min_support = adaptive;
                cfg.fold2.min_support = min_support;
            }
        }
    }

    elog!(
        "SUPPORT SANITY",
        "total_clusters={} total_support_reads={} max_n={}",
        support.len(),
        support.values().map(|s| s.n as u64).sum::<u64>(),
        support.values().map(|s| s.n).max().unwrap_or(0),
    );

    // Delegate to the same `is_real_foldback_stats()` pass2 actually gates on,
    // rather than re-deriving the logic inline -- this diagnostic had drifted
    // out of sync with the real gate before (it hardcoded the current gate
    // conditions rather than calling the shared function), which would have
    // silently stopped reflecting any future change to the real gate.
    let mut ge_support_clusters = 0usize;
    let mut ge_support_ident_clusters = 0usize;
    let mut real_clusters = 0usize;

    for st in support.values() {
        if st.n >= cfg.fold2.min_support {
            ge_support_clusters += 1;

            let ident_ok = cfg.fold2.min_support_ident <= 0.0
                || (st.mean_ident as f64) >= cfg.fold2.min_support_ident;

            if ident_ok {
                ge_support_ident_clusters += 1;
            }
        }
        if is_real_foldback_stats(st, &cfg) {
            real_clusters += 1;
        }
    }

    elog!(
        "SUPPORT GATE",
        "clusters: total={} ge_support={} ge_support_ident={} real={}",
        support.len(),
        ge_support_clusters,
        ge_support_ident_clusters,
        real_clusters
    );

    // Machine-parseable one-line summary for the --qc-only workflow (and
    // useful in a normal run too): the max signature-cluster size is the
    // statistic that flagged a real permissive-mode over-fragmentation case
    // during benchmarking (one dataset's largest cluster was 3-6x any other
    // tested dataset's). A large cluster here means many reads collapsed
    // onto the same signature -- consistent with a repetitive genomic
    // region being conflated into one "recurring locus" rather than a
    // single true recurring artefact, which --mode permissive's coarser
    // signature (smaller sig_flank_bp/sig_take, non-zero sig_value_shift)
    // and deeper cut recursion (max_depth) are more prone to.
    let large_cluster_threshold = 50usize;
    let max_cluster_n = support.values().map(|s| s.n).max().unwrap_or(0);
    let large_clusters = support
        .values()
        .filter(|s| s.n > large_cluster_threshold)
        .count();
    let reads_in_large_clusters: u64 = support
        .values()
        .filter(|s| s.n > large_cluster_threshold)
        .map(|s| s.n as u64)
        .sum();
    elog!(
        "QC SUMMARY",
        "mode={mode:?} max_cluster_n={} clusters_gt_{}={} reads_in_clusters_gt_{}={}",
        max_cluster_n,
        large_cluster_threshold,
        large_clusters,
        large_cluster_threshold,
        reads_in_large_clusters,
    );

    // Advisory only (not a gate change): --mode permissive's coarser signature
    // is the one mode shown, on real data, to occasionally collapse an
    // unrelated repeat family into one inflated "recurring locus" and
    // fragment the assembly as a result (see test_single_nematode/README.md,
    // "microlaimidae anomaly"). Thresholds below are set from that one
    // confirmed case (max_cluster_n 407, 0.26% of reads in large clusters)
    // vs. three unaffected species (max_cluster_n 70-130, <=0.05%) -- a
    // small n, so treat this as a heuristic flag worth a second look, not a
    // verdict. --mode strict/balanced are not flagged: not observed on this
    // dataset and their tighter signature/shallower recursion make the same
    // collision mechanism less likely (unvalidated beyond that one case).
    if mode == CallMode::Permissive {
        // Denominator must match what the thresholds above were calibrated
        // against in test_single_nematode/README.md: fractions there (e.g.
        // epsilonema 70/368,903 = 0.02%) are against *all* input reads, not
        // just the ones that passed detection/refinement/identity/signature
        // gates (total_n, used above for avg_n -- a much smaller, dataset-
        // dependent denominator that would make this fraction read higher
        // than the calibration and could false-fire on a healthy dataset).
        let frac_in_large_clusters = reads_in_large_clusters as f64 / total_reads.max(1) as f64;
        if max_cluster_n > MAX_PLAUSIBLE_CLUSTER_N || frac_in_large_clusters > 0.001 {
            elog!(
                "QC MODE ADVICE",
                "max_cluster_n={max_cluster_n} ({:.2}% of all input reads in clusters_gt_{large_cluster_threshold}) \
is high for --mode permissive on this dataset -- consider rerunning with --mode balanced or \
--mode strict and comparing assembly/BUSCO outcome before committing to permissive. \
(Heuristic from a single confirmed case; see test_single_nematode/README.md.)",
                frac_in_large_clusters * 100.0,
            );
        }
    }

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
