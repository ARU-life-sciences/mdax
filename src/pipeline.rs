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
    fingerprint::{
        SupportMap, best_support, foldback_signature_bands,
        foldback_signature_from_local_matches_bands, is_real_foldback_stats,
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
) -> Result<SupportMap> {
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
            move || -> Result<(Vec<(u64, usize, f32)>, Vec<(u64, u64)>)> {
                // One (band, split, ident) record per sig_ok read (any one of
                // that read's own bands serves as the lookup key into the
                // final union-find -- see SupportMap::build), plus the edges
                // connecting that read's own bands together so the merge
                // step can resolve them all to one canonical cluster.
                let mut records: Vec<(u64, usize, f32)> = Vec::new();
                let mut edges: Vec<(u64, u64)> = Vec::new();

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
                            fb.span,
                            fb.arm_end,
                            fb.arm_start_right,
                        )?
                        else {
                            continue;
                        };
                        dbg.bump("refined");

                        if rf.identity_est < cfg.fold2.min_identity {
                            continue;
                        }
                        dbg.bump("ident_ok");

                        let bands = foldback_signature_from_local_matches_bands(
                            &mut fold_scratch.best_matches,
                            rf.split_pos,
                            cfg.sig.flank_bp,
                            cfg.sig.take,
                            cfg.sig.value_shift,
                        )
                        .or_else(|| {
                            let q = 150usize;
                            let split_q = (rf.split_pos / q) * q;
                            foldback_signature_bands(
                                seq,
                                split_q,
                                &cfg.shared,
                                cfg.sig.flank_bp,
                                cfg.sig.take,
                                &mut sig_scratch,
                                cfg.sig.value_shift,
                            )
                        });

                        let Some(bands) = bands else { continue };
                        dbg.bump("sig_ok");

                        // This read's own bands are, by construction, all evidence for
                        // the same locus -- union them so the merge step treats them as
                        // one cluster even if another read only ever shares one of them
                        // (see SupportMap::build / fingerprint::UnionFind).
                        for i in 1..bands.len() {
                            edges.push((bands[0], bands[i]));
                        }
                        if let Some(&rep) = bands.first() {
                            records.push((rep, rf.split_pos, rf.identity_est));
                        }
                    }
                }

                elog!(
                    "PASS1 WORKER DONE",
                    "jobs={}, detected={}, refined={}, ident_ok={}, sig_ok={}, records={}",
                    dbg.jobs,
                    dbg.detected,
                    dbg.refined,
                    dbg.ident_ok,
                    dbg.sig_ok,
                    records.len(),
                );

                Ok((records, edges))
            },
        ));
    }

    reader_handle
        .join()
        .map_err(|_| anyhow::anyhow!("reader thread panicked"))??;

    // Merge every worker's per-read records and within-read band edges, then
    // resolve the whole dataset's bands to canonical clusters in one shot.
    // See `SupportMap::build` for why this can't be done as an incremental
    // per-band reduce: a band's final cluster membership can depend on an
    // edge discovered by a read processed by a different thread.
    let mut all_records: Vec<(u64, usize, f32)> = Vec::new();
    let mut all_edges: Vec<(u64, u64)> = Vec::new();
    for h in worker_handles {
        let (records, edges) = h
            .join()
            .map_err(|_| anyhow::anyhow!("worker thread panicked"))??;
        all_records.extend(records);
        all_edges.extend(edges);
    }

    let support = SupportMap::build(&all_records, &all_edges);

    Ok(support)
}

/// Pass 2 (parallel): detect foldbacks, decide whether they're "real", and (optionally) correct reads.
pub fn pass2_correct_and_write<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    support: Arc<SupportMap>,
    max_depth: usize,
    threads: usize,
    chan_cap: usize,
    fasta_file: std::fs::File,
    tsv_sink: TsvSink,
) -> Result<(u64, u64, u64)> {
    let (tx, rx) = chan::bounded::<BatchJob>(chan_cap);
    const BATCH_SIZE: usize = 256;
    let (out_tx, out_rx) = chan::bounded::<OutRec>(chan_cap * 8);

    // Precompute support_rank_frac for every signature (once)
    // rank_frac: 1.0 = top-ranked (largest n), 0.0 = bottom-ranked.
    let rank_map: std::collections::HashMap<u64, f32> = {
        let mut v: Vec<(u64, usize)> = support.iter().map(|(sig, st)| (*sig, st.n)).collect();

        // Sort by descending n, then ascending sig for deterministic tie-break.
        v.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let m = v.len();
        let mut rm = std::collections::HashMap::with_capacity(m);

        if m <= 1 {
            // If 0 or 1 signatures, everything is trivially "top".
            for (sig, _) in v {
                rm.insert(sig, 1.0);
            }
        } else {
            for (i, (sig, _n)) in v.into_iter().enumerate() {
                let rank = i + 1; // 1..=m
                let frac = 1.0 - ((rank - 1) as f32 / (m - 1) as f32);
                rm.insert(sig, frac);
            }
        }
        rm
    };
    let rank_map = std::sync::Arc::new(rank_map);

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
        let rank_map = rank_map.clone();

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

                            let (sig_split, sig_ident) = if let Some(rf_sig) = refined_sig.as_ref()
                            {
                                (rf_sig.split_pos, rf_sig.identity_est)
                            } else {
                                (rf_call.split_pos, rf_call.identity_est)
                            };

                            let (decision, support_n, support_rank_frac, support_span, sig_band) =
                                if sig_ident >= cfg.fold2.min_identity {
                                    dbg.bump("ident_ok");

                                    // Compute the local-match signature using the *refined* split
                                    // (sig_split), not the coarse split from detect_foldback.
                                    // Pass 1 also uses the refined split, so this keeps the two
                                    // passes consistent and ensures support-map lookups work for
                                    // genuine palindromes.
                                    //
                                    // Note: this destructively filters fold_scratch.best_matches.
                                    // The recursive cutter's first iteration will find an empty
                                    // window and fall back to the flank signature, which also
                                    // uses the refined split and therefore still matches pass 1.
                                    let bands_from_matches = foldback_signature_from_local_matches_bands(
                                        &mut fold_scratch.best_matches,
                                        sig_split,
                                        cfg.sig.flank_bp,
                                        cfg.sig.take,
                                        cfg.sig.value_shift,
                                    );

                                    let bands = bands_from_matches.or_else(|| {
                                        let q = 150usize;
                                        let split_q = (sig_split / q) * q;
                                        foldback_signature_bands(
                                            &job.seq,
                                            split_q,
                                            &cfg.shared,
                                            cfg.sig.flank_bp,
                                            cfg.sig.take,
                                            &mut sig_scratch,
                                            cfg.sig.value_shift,
                                        )
                                    });

                                    if let Some(bands) = bands {
                                        dbg.bump("sig_ok");

                                        let (real, n, rank_frac, span, band_val) = if let Some((band, st)) = best_support(&bands, &support) {
                                            dbg.bump("sig_in_support");
                                            if hit_sig_sample.len() < 5 {
                                                hit_sig_sample.push((band, st.n));
                                            }
                                            let rf = *rank_map.get(&band).unwrap_or(&0.0);
                                            (
                                                is_real_foldback_stats(st, &cfg),
                                                st.n,
                                                rf,
                                                st.split_span(),
                                                Some(band),
                                            )
                                        } else {
                                            dbg.bump("sig_not_in_support");
                                            if miss_sig_sample.len() < 5 && !bands.is_empty() {
                                                miss_sig_sample.push(bands[0]);
                                            }
                                            (false, 0, 0.0, 0, None)
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

                                        (if real { "real" } else { "artefact" }, n, rank_frac, span, band_val)
                                    } else {
                                        ("unknown", 0, 0.0, 0, None)
                                    }
                                } else if cfg.fold2.cut_low_ident {
                                    // Identity gate failed but the caller opted-in to
                                    // cutting anyway (useful for ONT where HiFi Hamming
                                    // systematically underestimates arm similarity).
                                    annot = Some(
                                        format!("_chopped_at_{}", rf_call.split_pos)
                                            .into_bytes(),
                                    );
                                    keep = 0..rf_call.split_pos;
                                    num_cut += 1;
                                    ("artefact", 0, 0.0, 0, None)
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
