//! The `mdax` processing pipeline called from main
//! We implement a two step pipeline:
//!  1) pass1: build support map from input reads (parallel)
//!  2) pass2: correct reads using support map (parallel), write output FASTA + TSV
//! Note: output is unordered, which is much faster

use anyhow::Result;
use crossbeam_channel as chan;
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;
use std::thread;

use crate::{
    cfg::MdaxCfg,
    fingerprint::{JunctionSupport, axis_tol, junction_keys},
    foldback, io,
    scratch::{FoldScratch, SigScratch},
    utils::RefineMode,
};

/// Lightweight counters for periodic debug/progress reporting.
///
/// This struct is intentionally tiny and `Clone` so it can be:
/// - kept per-thread (cheap to copy/reset)
/// - merged/summarised when emitting progress lines
#[derive(Default, Debug, Clone)]
struct DebugCounts {
    /// How many records were scheduled/seen
    jobs: u64,
    /// How many had a breakpoint refined to a final coordinate
    detected: u64,
    /// How many had a breakpoint refined to a final coordinate
    refined: u64,
    /// How many passed an identity/quality threshold
    ident_ok: u64,
    /// How many produced a valid signature
    sig_ok: u64,
    /// (pass2 only) how many signatures were found in the support map
    sig_in_support: u64,
    /// (pass2 only) how many signatures were NOT found in the support map
    sig_not_in_support: u64,
}

impl DebugCounts {
    /// Increment one of the counters by name.
    ///
    /// This is a pragmatic helper for call-sites that want to do something like
    /// `counts.bump("detected")` without carrying around a pile of `counts.foo += 1`.
    ///
    /// TODO: probably change this to an enum
    fn bump(&mut self, field: &str) {
        match field {
            "jobs" => self.jobs += 1,
            "detected" => self.detected += 1,
            "refined" => self.refined += 1,
            "ident_ok" => self.ident_ok += 1,
            "sig_ok" => self.sig_ok += 1,
            "sig_in_support" => self.sig_in_support += 1,
            "sig_not_in_support" => self.sig_not_in_support += 1,
            _ => {}
        }
    }
}

/// Print an error message with a tag prefix to stderr.
#[macro_export]
macro_rules! elog {
    ($tag:expr, $($arg:tt)*) => {{
        let prefix = format!("[{tag}] {message}", tag = $tag, message = format!($($arg)*));
        calm_io::stderrln!("{}", prefix)?;
    }}
}

/// Print debug output every N iterations.
///
/// This is a tiny throttling macro to avoid spamming stderr in tight loops.
/// It behaves like `eprintln!`, but only fires when `i` is a multiple of `every`.
///
/// Example:
/// ```ignore
/// dbg_every!(idx, 10_000, "processed {} reads", idx);
/// ```
///
/// Safety note:
/// - `every == 0` would panic due to modulo-by-zero; ensure callers pass a non-zero value.
macro_rules! dbg_every {
    ($i:expr, $every:expr, $tag:expr, $($arg:tt)*) => {
        if ($i % $every) == 0 {
            elog!($tag, $($arg)*);
        }
    }
}

/// Destination for TSV output.
///
/// This enum makes it easy to support either:
/// - streaming TSV to stdout
/// - writing TSV to a file handle.
pub enum TsvSink {
    /// Write TSV rows to standard output.
    Stdout,
    /// Write TSV rows to the given file handle.
    File(std::fs::File),
}

impl TsvSink {
    /// Convert the sink into a buffered, `Write + Send` trait object.
    fn into_writer(self) -> Box<dyn std::io::Write + Send> {
        match self {
            TsvSink::Stdout => Box::new(std::io::BufWriter::new(std::io::stdout())),
            TsvSink::File(f) => Box::new(std::io::BufWriter::new(f)),
        }
    }
}

/// A unit of work sent to worker threads.
///
/// Workers typically operate on "one input record at a time"
#[derive(Debug)]
struct Job {
    /// raw FASTA/FASTQ identifier bytes
    id: Vec<u8>,
    /// nucleotide sequence as bytes (shared, no pass2 copying)
    seq: Arc<Vec<u8>>,
}

#[derive(Debug)]
struct BatchJob {
    jobs: Vec<Job>,
}

/// Output produced by workers during pass2.
///
/// This struct carries everything the downstream "writer/collector"
#[derive(Debug)]
struct OutRec {
    /// final FASTA header bytes (base id + optional annotation)
    id: Vec<u8>,
    /// shared backing storage
    seq: Arc<Vec<u8>>,
    /// slice window into seq (no copy)
    keep: Range<usize>,
    /// optional TSV line
    tsv_row: Option<String>,
}

/// Pass 1 (parallel): scan the input once and count, for every detected
/// foldback, how many reads fold about the same axis (see
/// `fingerprint::JunctionSupport`).
/// Returns that support plus the total number of input records scanned
/// (not just the ones that produced a signature) -- callers computing a
/// fraction against "how collision-prone is this dataset" should divide by
/// this, not by the support map's own read count, which only covers reads
/// that passed detection/refinement/identity/signature gates and can be a
/// small, dataset-dependent fraction of the input.
pub fn pass1_build_support<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    threads: usize,
    chan_cap: usize,
) -> Result<(JunctionSupport, u64)> {
    let (tx, rx) = chan::bounded::<BatchJob>(chan_cap);
    const BATCH_SIZE: usize = 256;

    // Reader thread
    let input_path = input.as_ref().to_owned();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut r = io::open_fasta_reader(input_path)?;

        let mut batch: Vec<Job> = Vec::with_capacity(BATCH_SIZE);

        while let Some(rec) = r.next() {
            let rec = rec?;
            let id = rec.id().to_vec();
            let seq = Arc::new(rec.seq().to_vec());

            batch.push(Job { id, seq });

            if batch.len() == BATCH_SIZE {
                tx.send(BatchJob {
                    jobs: std::mem::take(&mut batch),
                })
                .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
                batch = Vec::with_capacity(BATCH_SIZE);
            }
        }

        if !batch.is_empty() {
            tx.send(BatchJob { jobs: batch })
                .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
        }

        drop(tx);
        Ok(())
    });

    // Worker threads
    let mut worker_handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let cfg = cfg.clone();
        let rx = rx.clone();

        worker_handles.push(thread::spawn(
            move || -> Result<(Vec<Vec<(u64, u64)>>, u64)> {
                // One junction-key set per read with a detected foldback.
                let mut groups: Vec<Vec<(u64, u64)>> = Vec::new();
                let tol = axis_tol(&cfg);

                let mut fold_scratch = FoldScratch::new();
                let mut dbg = DebugCounts::default();

                while let Ok(batch) = rx.recv() {
                    for job in batch.jobs {
                        dbg.bump("jobs");

                        dbg_every!(
                            dbg.jobs,
                            50_000,
                            "PASS1",
                            "jobs={} detected={} keyed={}",
                            dbg.jobs,
                            dbg.detected,
                            dbg.sig_ok
                        );

                        let seq = &job.seq;

                        if foldback::detect_foldback(
                            seq,
                            &cfg.shared,
                            &cfg.fold,
                            &mut fold_scratch,
                        )
                        .is_none()
                        {
                            continue;
                        }
                        dbg.bump("detected");

                        // Every detected foldback registers its fold axis,
                        // whatever its estimated arm identity: the keys
                        // need no refinement, and a diverged genome-encoded
                        // repeat needs its support counted just as much as
                        // a pristine one does.
                        junction_keys(
                            &fold_scratch.best_trip,
                            tol,
                            &mut fold_scratch.jk_tmp,
                            &mut fold_scratch.jk_out,
                        );
                        if fold_scratch.jk_out.is_empty() {
                            continue;
                        }
                        dbg.bump("sig_ok");

                        groups.push(fold_scratch.jk_out.clone());
                    }
                }

                elog!(
                    "PASS1 WORKER DONE",
                    "jobs={}, detected={}, keyed={}",
                    dbg.jobs,
                    dbg.detected,
                    groups.len(),
                );

                Ok((groups, dbg.jobs))
            },
        ));
    }

    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;

    let mut all_groups: Vec<Vec<(u64, u64)>> = Vec::new();
    let mut total_reads: u64 = 0;
    for h in worker_handles {
        let (groups, jobs) = h
            .join()
            .map_err(|_| anyhow::anyhow!("worker thread panicked"))??;
        all_groups.extend(groups);
        total_reads += jobs;
    }

    let support = JunctionSupport::build(&all_groups);

    Ok((support, total_reads))
}

/// Pass 2 (parallel): detect foldbacks, decide whether they're "real", and (optionally) correct reads.
pub fn pass2_correct_and_write<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    support: Arc<JunctionSupport>,
    max_depth: usize,
    threads: usize,
    chan_cap: usize,
    fasta_file: std::fs::File,
    tsv_sink: TsvSink,
) -> Result<(u64, u64, u64)> {
    let (tx, rx) = chan::bounded::<BatchJob>(chan_cap);
    const BATCH_SIZE: usize = 256;
    let (out_tx, out_rx) = chan::bounded::<OutRec>(chan_cap * 8);

    // Reader thread
    let input_path = input.as_ref().to_owned();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut r = io::open_fasta_reader(input_path)?;

        let mut batch: Vec<Job> = Vec::with_capacity(BATCH_SIZE);

        while let Some(rec) = r.next() {
            let rec = rec?;
            let id = rec.id().to_vec();
            let seq = Arc::new(rec.seq().to_vec());

            batch.push(Job { id, seq });

            if batch.len() == BATCH_SIZE {
                tx.send(BatchJob {
                    jobs: std::mem::take(&mut batch),
                })
                .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
                batch = Vec::with_capacity(BATCH_SIZE);
            }
        }

        if !batch.is_empty() {
            tx.send(BatchJob { jobs: batch })
                .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
        }

        drop(tx);
        Ok(())
    });

    // Writer thread
    let writer_handle = thread::spawn(move || -> Result<()> {
        use std::io::Write;

        // Bigger BufWriter than default can help when doing large sequential writes.
        // (default is usually 8 KiB; we want something much larger)
        const FASTA_BUFW_CAP: usize = 8 * 1024 * 1024; // 8 MiB
        const TSV_BUFW_CAP: usize = 2 * 1024 * 1024; // 2 MiB

        // How big our user-space batching buffers grow before we flush into BufWriter.
        const FASTA_FLUSH_BYTES: usize = 32 * 1024 * 1024; // 32 MiB
        const TSV_FLUSH_BYTES: usize = 4 * 1024 * 1024; // 4 MiB

        let mut fasta_out = std::io::BufWriter::with_capacity(FASTA_BUFW_CAP, fasta_file);
        let mut tsv_out = {
            // keep your existing sink behavior, but wrap it in a BufWriter with capacity
            let w = tsv_sink.into_writer();
            std::io::BufWriter::with_capacity(TSV_BUFW_CAP, w)
        };

        // Write TSV header once (via bytes, no formatting)
        tsv_out.write_all(
        b"read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est\tgap_est\tsupport_n\tsupport_rank_frac\tsupport_span\tdecision\tsig_band\n"
    )?;

        // User-space batching buffers
        let mut fasta_buf: Vec<u8> = Vec::with_capacity(FASTA_FLUSH_BYTES.min(2 * 1024 * 1024));
        let mut tsv_buf: Vec<u8> = Vec::with_capacity(TSV_FLUSH_BYTES.min(512 * 1024));

        while let Ok(rec) = out_rx.recv() {
            // ---- FASTA ----
            // Append: '>' + id + '\n' + seq + '\n'
            fasta_buf.push(b'>');
            fasta_buf.extend_from_slice(&rec.id);
            fasta_buf.push(b'\n');

            let seq = &rec.seq[rec.keep.start..rec.keep.end];
            fasta_buf.extend_from_slice(seq);
            fasta_buf.push(b'\n');

            if fasta_buf.len() >= FASTA_FLUSH_BYTES {
                fasta_out.write_all(&fasta_buf)?;
                fasta_buf.clear();
            }

            // ---- TSV ----
            if let Some(line) = rec.tsv_row {
                tsv_buf.extend_from_slice(line.as_bytes());
                tsv_buf.push(b'\n');

                if tsv_buf.len() >= TSV_FLUSH_BYTES {
                    tsv_out.write_all(&tsv_buf)?;
                    tsv_buf.clear();
                }
            }
        }

        // Final flush
        if !fasta_buf.is_empty() {
            fasta_out.write_all(&fasta_buf)?;
        }
        if !tsv_buf.is_empty() {
            tsv_out.write_all(&tsv_buf)?;
        }

        fasta_out.flush()?;
        tsv_out.flush()?;
        Ok(())
    });

    // Worker threads
    let mut worker_handles = Vec::with_capacity(threads);

    for _ in 0..threads {
        let cfg = cfg.clone();
        let support = support.clone();
        let rx = rx.clone();
        let out_tx = out_tx.clone();

        worker_handles.push(thread::spawn(move || -> Result<(u64, u64, u64)> {
            let mut num_foldbacks: u64 = 0;
            let mut num_cut: u64 = 0;
            let mut num_real: u64 = 0;

            let mut fold_scratch = FoldScratch::new();
            let mut sig_scratch = SigScratch::default();

            let mut shared_fast = cfg.shared.clone();
            shared_fast.refine.mode = RefineMode::HiFi;

            let mut dbg = DebugCounts::default();
            let tol = axis_tol(&cfg);

            while let Ok(batch) = rx.recv() {
                for job in batch.jobs {
                    dbg.bump("jobs");
                    dbg_every!(
                        dbg.jobs,
                        50_000,
                        "PASS2",
                        "jobs={} detected={} refined={} ident_ok={} sig_ok={} in_support={}",
                        dbg.jobs,
                        dbg.detected,
                        dbg.refined,
                        dbg.ident_ok,
                        dbg.sig_ok,
                        dbg.sig_in_support
                    );

                    let len = job.seq.len();
                    let mut keep: Range<usize> = 0..len;
                    let mut annot: Option<Vec<u8>> = None;
                    let mut tsv_row: Option<String> = None;

                    // Detect once
                    let fold = foldback::detect_foldback(
                        &job.seq,
                        &cfg.shared,
                        &cfg.fold,
                        &mut fold_scratch,
                    );

                    if fold.is_some() {
                        dbg.bump("detected");
                    }

                    if let Some(fb) = fold {
                        num_foldbacks += 1;

                        // Refine once for main call, optionally once more for stable sig coordinates.
                        // When mode is already HiFi, shared_fast == cfg.shared so we skip the
                        // second call and reuse the first result directly.
                        let refined_call = foldback::refine_breakpoint(
                            &job.seq,
                            fb.split_pos,
                            &cfg.shared,
                            &mut fold_scratch.refine,
                            fb.span,
                            fb.arm_end,
                            fb.arm_start_right,
                        )?;

                        let refined_sig = if cfg.shared.refine.mode == RefineMode::HiFi {
                            refined_call.clone()
                        } else {
                            foldback::refine_breakpoint(
                                &job.seq,
                                fb.split_pos,
                                &shared_fast,
                                &mut fold_scratch.refine,
                                fb.span,
                                fb.arm_end,
                                fb.arm_start_right,
                            )?
                        };

                        // When the user-mode refiner (e.g. ONT banded Levenshtein) returns
                        // None — typically because every candidate in the window exceeds the
                        // band — fall back to the HiFi result used for signature coherence.
                        // Without this fallback the read is silently dropped from the TSV
                        // and never cut, which is the dominant cause of low recall on ONT.
                        let rf_call_opt = refined_call.or_else(|| refined_sig.clone());

                        if let Some(rf_call) = rf_call_opt {
                            dbg.bump("refined");

                            let sig_ident = refined_sig
                                .as_ref()
                                .map_or(rf_call.identity_est, |rf_sig| rf_sig.identity_est);

                            let (decision, support_n, support_rank_frac, support_span, sig_band) =
                                if sig_ident >= cfg.fold2.min_identity {
                                    dbg.bump("ident_ok");

                                    junction_keys(
                                        &fold_scratch.best_trip,
                                        tol,
                                        &mut fold_scratch.jk_tmp,
                                        &mut fold_scratch.jk_out,
                                    );
                                    let n_elems = fold_scratch.jk_out.len();
                                    let key0 = fold_scratch.jk_out.first().map(|k| k.0);
                                    if n_elems > 0 {
                                        dbg.bump("sig_ok");
                                    }
                                    let n = support.support_n(&fold_scratch.jk_out);
                                    if n >= 2 {
                                        dbg.bump("sig_in_support");
                                    } else {
                                        dbg.bump("sig_not_in_support");
                                    }
                                    let real = n >= cfg.fold2.min_support;

                                    if !real {
                                        annot = Some(
                                            format!("_chopped_at_{}", rf_call.split_pos)
                                                .into_bytes(),
                                        );

                                        let keep2 = foldback::recursive_foldback_cut_from_first_range(
                                            &job.seq,
                                            keep.clone(),
                                            fb.clone(),
                                            rf_call.clone(),
                                            &cfg,
                                            &support,
                                            max_depth,
                                            &mut fold_scratch,
                                            &mut sig_scratch,
                                        )?;

                                        if keep2.start != 0 || keep2.end != len {
                                            num_cut += 1;
                                        }
                                        keep = keep2;
                                    } else {
                                        num_real += 1;
                                    }

                                    (if real { "real" } else { "artefact" }, n, 0.0f32, n_elems, key0)
                                } else if cfg.fold2.cut_low_ident {
                                    // Identity gate failed but the caller opted in to
                                    // cutting anyway. The support gate still applies:
                                    // a low-identity foldback that many reads share
                                    // an axis with is a diverged genomic repeat, not
                                    // an artefact.
                                    junction_keys(
                                        &fold_scratch.best_trip,
                                        tol,
                                        &mut fold_scratch.jk_tmp,
                                        &mut fold_scratch.jk_out,
                                    );
                                    let n_elems = fold_scratch.jk_out.len();
                                    let key0 = fold_scratch.jk_out.first().map(|k| k.0);
                                    let n = support.support_n(&fold_scratch.jk_out);
                                    if n >= cfg.fold2.min_support {
                                        num_real += 1;
                                        ("real", n, 0.0, n_elems, key0)
                                    } else {
                                    annot = Some(
                                        format!("_chopped_at_{}", rf_call.split_pos)
                                            .into_bytes(),
                                    );
                                    // Same convention as the main cutter: drop
                                    // one arm, keep any spacer, and keep
                                    // whichever side is longer rather than
                                    // always the left.
                                    let half_gap = rf_call.gap_est / 2;
                                    let arm_end = rf_call.split_pos.saturating_sub(half_gap);
                                    let arm_start = (rf_call.split_pos + half_gap).min(len);
                                    keep = if len - arm_end > arm_start {
                                        arm_end..len
                                    } else {
                                        0..arm_start
                                    };
                                    num_cut += 1;
                                    ("artefact", n, 0.0, n_elems, key0)
                                    }
                                } else {
                                    ("low_ident", 0, 0.0, 0, None)
                                };

                            // TSV row: only now do we build id_str
                            let id_str =
                                std::str::from_utf8(&job.id).unwrap_or("<nonutf8_id>");
                            let sig_band_str = sig_band
                                .map(|b| format!("{b:x}"))
                                .unwrap_or_else(|| "NA".to_string());
                            tsv_row = Some(format!(
                                "{id}\t{len}\tfoldback\t1\t{}\t{}\t\t{}\t{}\t\t\t{}\t{}\t{:.3}\t{}\t{}\t{:.6}\t{}\t{}\t{}",
                                fb.split_pos,
                                rf_call.split_pos,
                                fb.matches,
                                fb.span,
                                fb.score,
                                rf_call.score,
                                rf_call.identity_est,
                                rf_call.gap_est,
                                support_n,
                                support_rank_frac,
                                support_span,
                                decision,
                                sig_band_str,
                                id = id_str,
                                len = len
                            ));

                        }
                    }

                    // Final FASTA header bytes (base id + optional annotation bytes)
                    let mut out_id = job.id.clone();
                    if let Some(a) = annot {
                        out_id.extend_from_slice(&a);
                    }

                    out_tx
                        .send(OutRec {
                            id: out_id,
                            seq: job.seq.clone(),
                            keep,
                            tsv_row,
                        })
                        .map_err(|e| anyhow::anyhow!("out channel send failed: {e}"))?;
                }
            }

            elog!(
                "PASS2 WORKER DONE",
                "reads={} detected={} refined={} cut_reads={} real_reads={} reads_with_support_sig={} reads_without_support_sig={}",
                dbg.jobs,
                dbg.detected,
                dbg.refined,
                num_cut,
                num_real,
                dbg.sig_in_support,
                dbg.sig_not_in_support
            );

            Ok((num_foldbacks, num_cut, num_real))
        }));
    }

    // Let writer exit when workers are done
    drop(out_tx);
    drop(rx);

    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;

    let mut tot_fold = 0u64;
    let mut tot_cut = 0u64;
    let mut tot_real = 0u64;

    for h in worker_handles {
        let (a, b, c) = h.join().map_err(|_| anyhow::anyhow!("worker panicked"))??;
        tot_fold += a;
        tot_cut += b;
        tot_real += c;
    }

    writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    Ok((tot_fold, tot_cut, tot_real))
}
