//! The `mdax` processing pipeline called from main
//! We implement a two step pipeline:
//!  1) pass1: build support map from input reads (parallel)
//!  2) pass2: correct reads using support map (parallel), write output FASTA + TSV
//! Note: output is unordered, which is much faster

use anyhow::Result;
use crossbeam_channel as chan;
use gxhash::{HashMap, HashMapExt};
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;
use std::thread;

use crate::{
    cfg::MdaxCfg,
    fingerprint::{
        SupportStats, foldback_signature, foldback_signature_from_local_matches, is_real_foldback,
    },
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

/// Pass 1 (parallel): scan the input once and build the **support map**.
pub fn pass1_build_support<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    threads: usize,
    chan_cap: usize,
) -> Result<HashMap<u64, SupportStats>> {
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
            move || -> Result<HashMap<u64, SupportStats>> {
                let mut local: HashMap<u64, SupportStats> = HashMap::new();

                let mut fold_scratch = FoldScratch::new();
                let mut sig_scratch = SigScratch::default();

                // pass1 uses cheap refinement for signature coherence
                let mut shared_fast = cfg.shared.clone();
                shared_fast.refine.mode = RefineMode::HiFi;

                let mut dbg = DebugCounts::default();

                while let Ok(batch) = rx.recv() {
                    for job in batch.jobs {
                        dbg.bump("jobs");

                        dbg_every!(
                            dbg.jobs,
                            50_000,
                            "PASS1",
                            "jobs={} detected={} refined={} ident_ok={} sig_ok={}",
                            dbg.jobs,
                            dbg.detected,
                            dbg.refined,
                            dbg.ident_ok,
                            dbg.sig_ok
                        );

                        let seq = &job.seq;

                        let Some(fb) = foldback::detect_foldback(
                            seq,
                            &cfg.shared,
                            &cfg.fold,
                            &mut fold_scratch,
                        ) else {
                            continue;
                        };
                        dbg.bump("detected");

                        let Some(rf) = foldback::refine_breakpoint(
                            seq,
                            fb.split_pos,
                            &shared_fast,
                            &mut fold_scratch.refine,
                        )?
                        else {
                            continue;
                        };
                        dbg.bump("refined");

                        if rf.identity_est < cfg.fold2.min_identity {
                            continue;
                        }
                        dbg.bump("ident_ok");

                        let sig = foldback_signature_from_local_matches(
                            &mut fold_scratch.best_matches,
                            rf.split_pos,
                            cfg.sig.flank_bp,
                            cfg.sig.take,
                            cfg.sig.value_shift,
                        )
                        .or_else(|| {
                            let q = 150usize;
                            let split_q = (rf.split_pos / q) * q;
                            foldback_signature(
                                seq,
                                split_q,
                                &cfg.shared,
                                cfg.sig.flank_bp,
                                cfg.sig.take,
                                &mut sig_scratch,
                                cfg.sig.value_shift,
                            )
                        });

                        let Some(sig) = sig else { continue };
                        dbg.bump("sig_ok");

                        local
                            .entry(sig)
                            .and_modify(|st| st.update(rf.split_pos, rf.identity_est))
                            .or_insert_with(|| SupportStats::new(rf.split_pos, rf.identity_est));
                    }
                }

                elog!(
                    "PASS1 WORKER DONE",
                    "jobs={}, detected={}, refined={}, ident_ok={}, sig_ok={}, unique_sigs={}",
                    dbg.jobs,
                    dbg.detected,
                    dbg.refined,
                    dbg.ident_ok,
                    dbg.sig_ok,
                    local.len(),
                );

                Ok(local)
            },
        ));
    }

    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;

    // Reduce locals into global
    let mut support: HashMap<u64, SupportStats> = HashMap::new();
    for h in worker_handles {
        let local = h
            .join()
            .map_err(|_| anyhow::anyhow!("worker thread panicked"))??;

        for (sig, st) in local {
            support
                .entry(sig)
                .and_modify(|g| {
                    let total_n = g.n + st.n;
                    let mean = if total_n > 0 {
                        (g.mean_ident * g.n as f32 + st.mean_ident * st.n as f32) / total_n as f32
                    } else {
                        0.0
                    };
                    g.n = total_n;
                    g.min_split = g.min_split.min(st.min_split);
                    g.max_split = g.max_split.max(st.max_split);
                    g.mean_ident = mean;
                })
                .or_insert(st);
        }
    }

    Ok(support)
}

/// Pass 2 (parallel): detect foldbacks, decide whether they're "real", and (optionally) correct reads.
pub fn pass2_correct_and_write<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    support: Arc<HashMap<u64, SupportStats>>,
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
        b"read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est\tsupport_n\tsupport_span\tdecision\n"
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
            let mut miss_sig_sample: Vec<u64> = Vec::new();
            let mut hit_sig_sample: Vec<(u64, usize)> = Vec::new();

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

                    // Precompute match-derived signature (only if fold exists)
                    let sig_from_matches_pre: Option<u64> = fold.as_ref().and_then(|fb| {
                        foldback_signature_from_local_matches(
                            &mut fold_scratch.best_matches,
                            fb.split_pos,
                            cfg.sig.flank_bp,
                            cfg.sig.take,
                            cfg.sig.value_shift,
                        )
                    });

                    if let Some(fb) = fold {
                        num_foldbacks += 1;

                        // Refine once for main call, optionally once more for stable sig coordinates
                        let refined_call = foldback::refine_breakpoint(
                            &job.seq,
                            fb.split_pos,
                            &cfg.shared,
                            &mut fold_scratch.refine,
                        )?;

                        let refined_sig = foldback::refine_breakpoint(
                            &job.seq,
                            fb.split_pos,
                            &shared_fast,
                            &mut fold_scratch.refine,
                        )?;

                        if let Some(rf_call) = refined_call {
                            dbg.bump("refined");

                            let (sig_split, sig_ident) = if let Some(rf_sig) = refined_sig.as_ref()
                            {
                                (rf_sig.split_pos, rf_sig.identity_est)
                            } else {
                                (rf_call.split_pos, rf_call.identity_est)
                            };

                            let (decision, support_n, support_span) =
                                if sig_ident >= cfg.fold2.min_identity {
                                    dbg.bump("ident_ok");

                                    let sig = sig_from_matches_pre.or_else(|| {
                                        let q = 150usize;
                                        let split_q = (sig_split / q) * q;
                                        foldback_signature(
                                            &job.seq,
                                            split_q,
                                            &cfg.shared,
                                            cfg.sig.flank_bp,
                                            cfg.sig.take,
                                            &mut sig_scratch,
                                            cfg.sig.value_shift,
                                        )
                                    });

                                    if let Some(sig) = sig {
                                        dbg.bump("sig_ok");

                                        let (real, n, span) = if let Some(st) = support.get(&sig)
                                        {
                                            dbg.bump("sig_in_support");
                                            if hit_sig_sample.len() < 5 {
                                                hit_sig_sample.push((sig, st.n));
                                            }
                                            (is_real_foldback(sig, &support, &cfg), st.n, st.split_span())
                                        } else {
                                            dbg.bump("sig_not_in_support");
                                            if miss_sig_sample.len() < 5 {
                                                miss_sig_sample.push(sig);
                                            }
                                            (false, 0, 0)
                                        };

                                        if !real {
                                            annot = Some(
                                                format!("_chopped_at_{}", rf_call.split_pos)
                                                    .into_bytes(),
                                            );

                                            let keep2 = foldback::recursive_foldback_cut_from_first_range(
                                                &job.seq,
                                                keep.clone(),   // <— use current keep window
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

                                        (if real { "real" } else { "artefact" }, n, span)
                                    } else {
                                        ("unknown", 0, 0)
                                    }
                                } else {
                                    ("low_ident", 0, 0)
                                };

                            // TSV row: only now do we build id_str
                            let id_str =
                                std::str::from_utf8(&job.id).unwrap_or("<nonutf8_id>");
                            tsv_row = Some(format!(
                                "{id}\t{len}\tfoldback\t1\t{}\t{}\t\t{}\t{}\t\t\t{}\t{}\t{:.3}\t{}\t{}\t{}",
                                fb.split_pos,
                                rf_call.split_pos,
                                fb.matches,
                                fb.span,
                                fb.score,
                                rf_call.score,
                                rf_call.identity_est,
                                support_n,
                                support_span,
                                decision,
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
