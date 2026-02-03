//! The `mdax` processing pipeline called from main
//! We implement a two step pipeline:
//!  1) pass1: build support map from input reads (parallel)
//!  2) pass2: correct reads using support map (parallel), write output FASTA + TSV
//! Note: output is unordered, which is much faster

use anyhow::Result;
use crossbeam_channel as chan;
use gxhash::{HashMap, HashMapExt};
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
    utils::{RefineMode, write_fasta},
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
    /// monotonically increasing index
    idx: u64,
    /// raw FASTA/FASTQ identifier bytes
    id: Vec<u8>,
    /// nucleotide sequence as bytes
    seq: Vec<u8>,
}

/// Output produced by workers during pass2.
///
/// This struct carries everything the downstream "writer/collector"
#[derive(Debug)]
struct OutRec {
    /// same index as `Job`
    idx: u64,
    /// record identifier as String for writing fasta headers/tsv rows
    id: String,
    /// the sequence to emit
    kept_seq: Vec<u8>,
    /// optional TSV line
    tsv_row: Option<String>,
}

/// Pass 1 (parallel): scan the input once and build the **support map**.
///
/// The support map is keyed by a foldback "signature" (e.g., minimizer-derived fingerprint)
/// and stores aggregated evidence (`SupportStats`) about where/how strongly that signature
/// appears across reads.
///
/// Concurrency / architecture:
/// - A single **reader thread** streams FASTA records and pushes `Job`s into a bounded channel.
/// - `threads` **worker threads** pop `Job`s, run detect -> refine -> signature, and accumulate
///   a **thread-local** `HashMap<u64, SupportStats>`.
/// - The main thread **reduces** worker-local maps into a single global `HashMap`.
///
/// Memory behavior:
/// - Input is streamed; the dominant buffering is bounded by `chan_cap` jobs.
/// - Each worker’s local support map grows with the number of unique signatures it observes.
///
/// Implementation note:
/// - We intentionally use a *cheap/robust refinement mode* (HiFi) when producing signatures
///   during pass1, to make signature keys less sensitive to small split jitter.
pub fn pass1_build_support<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    threads: usize,
    chan_cap: usize,
) -> Result<HashMap<u64, SupportStats>> {
    // bounded job channel backpressures the reader if workers fall behind.
    let (tx, rx) = chan::bounded::<Job>(chan_cap);

    // Reader thread
    // Reads FASTA records and sends owned `Job`s to workers.
    // needletail yields borrowed slices, so we clone `id` and `seq` to move them across threads.
    let input_path = input.as_ref().to_owned();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut r = io::open_fasta_reader(input_path)?;
        let mut idx: u64 = 0;

        while let Some(rec) = r.next() {
            let rec = rec?;
            // needletail returns borrowed data; clone into owned Vec for crossing threads
            let id = rec.id().to_vec();
            let seq = rec.seq().to_vec();

            tx.send(Job { idx, id, seq })
                .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
            idx += 1;
        }
        drop(tx);
        Ok(())
    });

    // Worker threads: each returns a local map
    let mut worker_handles = Vec::with_capacity(threads);
    for _ in 0..threads {
        let cfg = cfg.clone();
        let rx = rx.clone();

        worker_handles.push(thread::spawn(
            move || -> Result<HashMap<u64, SupportStats>> {
                let mut local: HashMap<u64, SupportStats> = HashMap::new();

                // Per-thread scratch buffers reused across reads to avoid repeated allocations.
                let mut fold_scratch = FoldScratch::new();
                let mut sig_scratch = SigScratch::default();

                // IMPORTANT: pass1 uses cheap refinement (HiFi) for support/signature coherence.
                // We do this by cloning shared and forcing mode=HiFi.
                let mut shared_fast = cfg.shared.clone();
                shared_fast.refine.mode = RefineMode::HiFi;

                let mut dbg = DebugCounts::default();
                let mut seen_sig_sample: Vec<u64> = Vec::new(); // tiny sample for inspection

                while let Ok(job) = rx.recv() {
                    dbg.bump("jobs");

                    // progress heartbeat per worker: every 50k reads
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

                    let seq = job.seq;

                    // 1) Coarse detection: find a foldback candidate and seed matchpoints.
                    // detect_foldback typically populates `fold_scratch.best_matches` as a side effect.
                    let Some(fb) =
                        foldback::detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut fold_scratch)
                    else {
                        continue;
                    };
                    dbg.bump("detected");

                    // 2) Refine breakpoint using *cheap* mode for stability across reads.
                    let Some(rf) = foldback::refine_breakpoint(
                        &seq,
                        fb.split_pos,
                        &shared_fast,
                        &mut fold_scratch.refine,
                    ) else {
                        continue;
                    };
                    dbg.bump("refined");

                    // 3) Identity filter: only keep evidence strong enough to be meaningful in support.
                    if rf.identity_est < cfg.fold2.min_identity {
                        continue;
                    }
                    dbg.bump("ident_ok");

                    // 4) Compute signature.
                    //
                    // Prefer a signature derived from the matchpoints that supported the best bin:
                    // - More tolerant to split jitter (since matchpoints encode the event pattern)
                    // - Avoids depending too heavily on the exact refined split coordinate
                    //
                    // Fallback: flank-based signature, but quantize split to increase coherence.
                    // TODO: make quantization configurable?
                    let sig = foldback_signature_from_local_matches(
                        &mut fold_scratch.best_matches,
                        rf.split_pos,
                        cfg.sig.flank_bp,
                        cfg.sig.take,
                        cfg.sig.value_shift,
                    )
                    .or_else(|| {
                        // Fallback to flank-based signature (still useful, but more split-sensitive).
                        // Quantize split so pass1/pass2 behave consistently.
                        let q = 150usize;
                        let split_q = (rf.split_pos / q) * q;

                        foldback_signature(
                            &seq,
                            split_q,
                            &cfg.shared,
                            cfg.sig.flank_bp,
                            cfg.sig.take,
                            &mut sig_scratch,
                            cfg.sig.value_shift,
                        )
                    });

                    let Some(sig) = sig else {
                        continue;
                    };

                    dbg.bump("sig_ok");

                    // Keep a tiny sample for debugging so you can eyeball signature stability.
                    // TODO: improve output
                    if seen_sig_sample.len() < 5 {
                        seen_sig_sample.push(sig);
                    }

                    // 5) Update local support stats for this signature.
                    // We store split + identity estimates so pass2 can judge "real vs artefact".
                    local
                        .entry(sig)
                        .and_modify(|st| st.update(rf.split_pos, rf.identity_est))
                        .or_insert_with(|| SupportStats::new(rf.split_pos, rf.identity_est));
                }

                // At worker exit: print final counts + small signature sample
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

    // Wait for reader
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
                    // merge stats
                    // merge via repeated update is okay; we can do exact merge too
                    // We'll do exact merge:
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
///
/// Output:
/// - Always writes a FASTA record for each input read (possibly corrected/trimmed).
/// - Optionally writes a TSV row for reads where an event was called.
///
/// Concurrency / architecture:
/// - Reader thread streams input records into a bounded `Job` channel.
/// - Worker threads:
///   1) run detection/refinement once,
///   2) compute a signature,
///   3) consult the pass1 `support` map,
///   4) decide real vs artefact,
///   5) optionally perform recursive cutting for artefacts,
///   6) send an `OutRec` to the writer.
/// - A single writer thread serializes FASTA + TSV output.
///
/// Returns:
/// `(num_foldbacks, num_cut, num_real)` aggregated across workers.
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
    use std::io::Write;

    // Jobs are bounded to apply backpressure.
    let (tx, rx) = chan::bounded::<Job>(chan_cap);

    // Output channel is larger because writing can be slower than detection.
    // Buffering here prevents workers from stalling on transient writer slowness.
    let (out_tx, out_rx) = chan::bounded::<OutRec>(chan_cap * 8);

    // Reader thread
    // Reads FASTA records and pushes owned jobs to workers.
    // Closing `tx` is essential so worker loops terminate.
    let input_path = input.as_ref().to_owned();
    let reader_handle = thread::spawn(move || -> Result<()> {
        let mut r = io::open_fasta_reader(input_path)?;
        let mut idx: u64 = 0;

        while let Some(rec) = r.next() {
            let rec = rec?;
            let id = rec.id().to_vec();
            let seq = rec.seq().to_vec();

            tx.send(Job { idx, id, seq })
                .map_err(|e| anyhow::anyhow!("channel send failed: {e}"))?;
            idx += 1;
        }

        // critical: close the jobs channel so workers stop
        drop(tx);
        Ok(())
    });

    // Writer thread (single)
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut fasta_out = std::io::BufWriter::new(fasta_file);
        let mut tsv_writer = tsv_sink.into_writer();

        // TSV header once
        writeln!(
            tsv_writer,
            "read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est\tsupport_n\tsupport_span\tdecision"
        )?;

        // This loop ends when ALL out_tx senders are dropped (after workers exit)
        while let Ok(rec) = out_rx.recv() {
            // FASTA
            write_fasta(&mut fasta_out, &rec.id, &rec.kept_seq)?;

            // TSV
            if let Some(line) = rec.tsv_row {
                writeln!(tsv_writer, "{line}")?;
            }
        }

        tsv_writer.flush()?;
        Ok(())
    });

    // Worker threads
    // Each worker processes jobs independently and reports summary counts.
    let mut worker_handles = Vec::with_capacity(threads);

    for _ in 0..threads {
        let cfg = cfg.clone();
        let support = support.clone();
        let rx = rx.clone();
        let out_tx = out_tx.clone();

        worker_handles.push(thread::spawn(move || -> Result<(u64, u64, u64)> {
            // Worker-local tallies (returned to main for aggregation).
            let mut num_foldbacks: u64 = 0;
            let mut num_cut: u64 = 0;
            let mut num_real: u64 = 0;

            let mut fold_scratch = FoldScratch::new();
            let mut sig_scratch = SigScratch::default();

            // cheap refine for signature lookup coherence
            let mut shared_fast = cfg.shared.clone();
            // using HiFi mode tends to reduce signature instability due to small split noise.
            shared_fast.refine.mode = RefineMode::HiFi;

            // Debug counters and small samples to inspect signature hits/misses.
            let mut dbg = DebugCounts::default();
            let mut miss_sig_sample: Vec<u64> = Vec::new();
            let mut hit_sig_sample: Vec<(u64, usize)> = Vec::new(); // (sig, n)

            while let Ok(job) = rx.recv() {
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

                let id_str = std::str::from_utf8(&job.id)?.to_string();
                // to annotate the header
                let mut header_annotation = String::new();
                let len = job.seq.len();

                // 1) Coarse detect ONCE
                // detect_foldback may also populate `fold_scratch.best_matches`, which we reuse
                // to compute a match-based signature without re-running detection.
                let fold = foldback::detect_foldback(&job.seq, &cfg.shared, &cfg.fold, &mut fold_scratch);

                if fold.is_some() {
                    dbg.bump("detected");
                }

                // Precompute signature from the match evidence produced by detect_foldback.
                // Only valid if `fold.is_some()`, since best_matches is a side-effect.
                let sig_from_matches_pre: Option<u64> = fold.as_ref().and_then(|fb| {
                    foldback_signature_from_local_matches(
                        &mut fold_scratch.best_matches,
                        fb.split_pos,          // coarse split is fine for match-based signature
                        cfg.sig.flank_bp,
                        cfg.sig.take,
                        cfg.sig.value_shift,
                    )
                });

                // We'll decide what sequence we keep as a slice, then allocate once at end.
                let mut kept_slice: &[u8] = &job.seq;

                // TSV row (optional)
                let mut tsv_row: Option<String> = None;

                if let Some(fb) = fold {
                    num_foldbacks += 1;

                    // 2) Refine ONCE
                    let refined_call = foldback::refine_breakpoint(
                        &job.seq,
                        fb.split_pos,
                        &cfg.shared,
                        &mut fold_scratch.refine,
                    );

                    // optional: HiFi refine to get stable sig_split/sig_ident if you want
                    let refined_sig = foldback::refine_breakpoint(
                        &job.seq,
                        fb.split_pos,
                        &shared_fast,
                        &mut fold_scratch.refine,
                    );

                    if let Some(rf_call) = refined_call {
                        dbg.bump("refined");

                        // Prefer the "stable" HiFi-derived split/identity if available for signature
                        // coherence; otherwise fall back to the main refined call.
                        let (sig_split, sig_ident) = if let Some(rf_sig) = refined_sig.as_ref() {
                            (rf_sig.split_pos, rf_sig.identity_est)
                        } else {
                            (rf_call.split_pos, rf_call.identity_est)
                        };

                        // 3) Decide + recurse WITHOUT double detection
                        // We gate on identity first; low-identity events are not trusted.
                        let (decision, support_n, support_span) = if sig_ident >= cfg.fold2.min_identity {
                            dbg.bump("ident_ok");

                            // Prefer evidence-based signature for support lookup
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

                                // Support-aware classification:
                                // - If the signature is present in support, it *may* represent a real,
                                //   genome-templated foldback.
                                // - If absent/weak, treat as likely MDA artefact.
                                let (real, n, span) = if let Some(st) = support.get(&sig) {
                                    dbg.bump("sig_in_support");
                                    if hit_sig_sample.len() < 5 {
                                        hit_sig_sample.push((sig, st.n));
                                    }
                                    (is_real_foldback(sig, &support, &cfg), st.n, st.split_span())
                                } else {
                                    if miss_sig_sample.len() < 5 {
                                        miss_sig_sample.push(sig);
                                    }
                                    (false, 0, 0)
                                };

                                // IMPORTANT: only cut if NOT real
                                if !real {
                                    header_annotation = format!(
                                        "_chopped_at_{}", rf_call.split_pos
                                    );
                                    // Seed recursion with the already-computed detection + refinement
                                    // to avoid redundant work.
                                    let kept = foldback::recursive_foldback_cut_from_first(
                                        &job.seq,
                                        fb.clone(),
                                        rf_call.clone(),
                                        &cfg,
                                        &support,
                                        max_depth,
                                        &mut fold_scratch,
                                        &mut sig_scratch,
                                    )?;

                                    if kept.len() < job.seq.len() {
                                        num_cut += 1;
                                    }
                                    kept_slice = kept;
                                } else {
                                    num_real += 1;
                                }

                                (if real { "real" } else { "artefact" }, n, span)
                            } else {
                                // Signature failed: can't consult support, so label as unknown.
                                ("unknown", 0, 0)
                            }
                        } else {
                            // Refined identity too low to trust the call.
                            ("low_ident", 0, 0)
                        };

                        // TSV output (no extra detect/refine)
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

                // Allocate once for output record
                let kept_seq: Vec<u8> = kept_slice.to_vec();

                // Construct the FASTA header
                let complete_id = format!("{}{}", id_str, header_annotation);

                out_tx
                    .send(OutRec {
                        idx: job.idx,
                        id: complete_id,
                        kept_seq,
                        tsv_row,
                    })
                    .map_err(|e| anyhow::anyhow!("out channel send failed: {e}"))?;
            }


            elog!(
                "PASS2 WORKER DONE",
                "jobs={}, detected={}, refined={}, cut={}, real={}, support_hits={}, support_misses={}",
                dbg.jobs,
                dbg.detected,
                dbg.refined,
                num_cut,
                num_real,
                dbg.sig_in_support,
                miss_sig_sample.len()
            );

            Ok((num_foldbacks, num_cut, num_real))
        }));
    }

    // IMPORTANT:
    // drop the main copies so the writer can terminate once workers are done
    drop(out_tx);
    drop(rx); // drop main rx clone (workers still have theirs)

    // Join reader
    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;

    // Join workers and sum stats
    let mut tot_fold = 0u64;
    let mut tot_cut = 0u64;
    let mut tot_real = 0u64;

    for h in worker_handles {
        let (a, b, c) = h.join().map_err(|_| anyhow::anyhow!("worker panicked"))??;
        tot_fold += a;
        tot_cut += b;
        tot_real += c;
    }

    // After workers exit, all their out_tx clones drop → out_rx closes -> writer exits
    writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    Ok((tot_fold, tot_cut, tot_real))
}
