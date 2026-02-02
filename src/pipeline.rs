use anyhow::Result;
use crossbeam_channel as chan;
use gxhash::{HashMap, HashMapExt};
use noodles::fasta::{
    Record,
    record::{Definition, Sequence},
};
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;
use std::thread;

use crate::{
    cfg::MdaxCfg,
    fingerprint::{
        SupportStats, foldback_signature, foldback_signature_from_matches, is_real_foldback,
    },
    foldback, io,
    scratch::{FoldScratch, SigScratch},
    utils::RefineMode,
};

#[derive(Default, Debug, Clone)]
struct DebugCounts {
    jobs: u64,
    detected: u64,
    refined: u64,
    ident_ok: u64,
    sig_ok: u64,
    sig_in_support: u64, // pass2 only
}

impl DebugCounts {
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

macro_rules! dbg_every {
    ($i:expr, $every:expr, $($arg:tt)*) => {
        if ($i % $every) == 0 {
            eprintln!($($arg)*);
        }
    }
}

pub enum TsvSink {
    Stdout,
    File(std::fs::File),
}

impl TsvSink {
    fn into_writer(self) -> Box<dyn std::io::Write + Send> {
        match self {
            TsvSink::Stdout => Box::new(std::io::BufWriter::new(std::io::stdout())),
            TsvSink::File(f) => Box::new(std::io::BufWriter::new(f)),
        }
    }
}

/// A job sent to workers.
#[derive(Debug)]
struct Job {
    idx: u64,
    id: Vec<u8>,
    seq: Vec<u8>,
}

/// Output from workers for pass2.
#[derive(Debug)]
struct OutRec {
    idx: u64,
    id: String,
    kept_seq: Vec<u8>,
    tsv_row: Option<String>,
}

/// Pass1 (parallel): build support map.
/// This streams the input; memory bounded by channel capacity.
/// Each worker accumulates a local HashMap and we reduce at the end.
pub fn pass1_build_support<P: AsRef<Path>>(
    input: P,
    cfg: Arc<MdaxCfg>,
    threads: usize,
    chan_cap: usize,
) -> Result<HashMap<u64, SupportStats>> {
    let (tx, rx) = chan::bounded::<Job>(chan_cap);

    // Reader thread
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
                    dbg_every!(dbg.jobs, 50_000, "[pass1 worker] jobs={} detected={} refined={} ident_ok={} sig_ok={}",
                        dbg.jobs, dbg.detected, dbg.refined, dbg.ident_ok, dbg.sig_ok);

                    let seq = job.seq;

                    let Some(fb) =
                        foldback::detect_foldback(&seq, &cfg.shared, &cfg.fold, &mut fold_scratch)
                    else {
                        continue;
                    };
                    dbg.bump("detected");

                    let Some(rf) = foldback::refine_breakpoint(
                        &seq,
                        fb.split_pos,
                        &shared_fast,
                        &mut fold_scratch.refine,
                    ) else {
                        continue;
                    };
                    dbg.bump("refined");

                    if rf.identity_est < cfg.fold2.min_identity {
                        continue;
                    }
                    dbg.bump("ident_ok");

                    // Prefer evidence-based signature from the matched minimizer seeds
                    // that supported the best anti-diagonal bin (more tolerant to split jitter).
                    let sig = foldback_signature_from_matches(&mut fold_scratch.best_vals, cfg.sig.take, cfg.sig.value_shift)
                        .or_else(|| {
                            // Fallback to flank-based signature (still useful, but more split-sensitive).
                            // Quantize split so pass1/pass2 behave consistently.
                            let q = 50usize;
                            let split_q = (rf.split_pos / q) * q;

                            foldback_signature(
                                &seq,
                                split_q,
                                &cfg.shared,
                                cfg.sig.flank_bp,
                                cfg.sig.take,
                                &mut sig_scratch,
                                cfg.sig.value_shift
                            )
                        });

                    let Some(sig) = sig else { continue; };

                    dbg.bump("sig_ok");

                    if seen_sig_sample.len() < 5 {
                        seen_sig_sample.push(sig);
                    }

                    local
                        .entry(sig)
                        .and_modify(|st| st.update(rf.split_pos, rf.identity_est))
                        .or_insert_with(|| SupportStats::new(rf.split_pos, rf.identity_est));
                }

                // At worker exit: print final counts + small signature sample
                eprintln!(
                    "[pass1 worker done] jobs={} detected={} refined={} ident_ok={} sig_ok={} local_sigs={} sample={:?}",
                    dbg.jobs, dbg.detected, dbg.refined, dbg.ident_ok, dbg.sig_ok, local.len(), seen_sig_sample
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

    let (tx, rx) = chan::bounded::<Job>(chan_cap);
    let (out_tx, out_rx) = chan::bounded::<OutRec>(chan_cap);

    // ----------------
    // Reader thread
    // ----------------
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

    // ----------------
    // Writer thread (single)
    // ----------------
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut fasta_writer = noodles::fasta::io::Writer::new(fasta_file);
        let mut tsv_writer = tsv_sink.into_writer();

        // TSV header once
        writeln!(
            tsv_writer,
            "read_id\tlen\tevent\tcalled\tcoarse_split\trefined_split\tdelta\tmatches\tspan_p1\tp2_span\tcross_frac\tcoarse_score\trefined_score\tidentity_est\tsupport_n\tsupport_span\tdecision"
        )?;

        let mut next: u64 = 0;
        let mut buf: BTreeMap<u64, OutRec> = BTreeMap::new();

        // This loop ends when ALL out_tx senders are dropped (after workers exit)
        while let Ok(rec) = out_rx.recv() {
            buf.insert(rec.idx, rec);

            while let Some(r) = buf.remove(&next) {
                // FASTA
                let def = Definition::new(r.id.as_str(), None);
                let seq = Sequence::from(r.kept_seq);
                let record = Record::new(def, seq);
                fasta_writer.write_record(&record)?;

                // TSV
                if let Some(line) = r.tsv_row {
                    writeln!(tsv_writer, "{line}")?;
                }

                next += 1;
            }
        }

        // Flush any remaining in-order records (should normally be none)
        while let Some(r) = buf.remove(&next) {
            let def = Definition::new(r.id.as_str(), None);
            let seq = Sequence::from(r.kept_seq);
            let record = Record::new(def, seq);
            fasta_writer.write_record(&record)?;

            if let Some(line) = r.tsv_row {
                writeln!(tsv_writer, "{line}")?;
            }
            next += 1;
        }

        tsv_writer.flush()?;
        Ok(())
    });

    // ----------------
    // Worker threads
    // ----------------
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

            // cheap refine for signature lookup coherence
            let mut shared_fast = cfg.shared.clone();
            shared_fast.refine.mode = RefineMode::HiFi;

            let mut dbg = DebugCounts::default();
            let mut miss_sig_sample: Vec<u64> = Vec::new();
            let mut hit_sig_sample: Vec<(u64, usize)> = Vec::new(); // (sig, n)

            while let Ok(job) = rx.recv() {
                dbg.bump("jobs");
                dbg_every!(dbg.jobs, 50_000, "[pass2 worker] jobs={} detected={} refined={} ident_ok={} sig_ok={} in_support={}",
                    dbg.jobs, dbg.detected, dbg.refined, dbg.ident_ok, dbg.sig_ok, dbg.sig_in_support);

                let id_str = std::str::from_utf8(&job.id)?.to_string();
                let len = job.seq.len();

                let fold = foldback::detect_foldback(&job.seq, &cfg.shared, &cfg.fold, &mut fold_scratch);
                let sig_from_matches_pre = fold.as_ref().and_then(|_fb| {
                    foldback_signature_from_matches(
                        &mut fold_scratch.best_vals,
                        cfg.sig.take,
                        cfg.sig.value_shift,
                    )
                });

                if fold.is_some() {
                    dbg.bump("detected");
                }

                // kept sequence
                let kept_seq: Vec<u8> = if fold.is_some() {
                    let kept = foldback::recursive_foldback_cut(
                        &job.seq,
                        &cfg,
                        &support,
                        max_depth,
                        &mut fold_scratch,
                        &mut sig_scratch,
                    )?;
                    if kept.len() < job.seq.len() {
                        num_cut += 1;
                    }
                    kept.to_vec()
                } else {
                    job.seq.clone()
                };

                // TSV row (optional)
                let mut tsv_row: Option<String> = None;

                if let Some(fb) = fold {
                    num_foldbacks += 1;

                    let refined_call = foldback::refine_breakpoint(
                        &job.seq,
                        fb.split_pos,
                        &cfg.shared,
                        &mut fold_scratch.refine,
                    );

                    let refined_sig = foldback::refine_breakpoint(
                        &job.seq,
                        fb.split_pos,
                        &shared_fast,
                        &mut fold_scratch.refine,
                    );

                    if let Some(rf_call) = refined_call {
                        dbg.bump("refined");
                        let (sig_split, sig_ident) = if let Some(rf_sig) = refined_sig.as_ref() {
                            (rf_sig.split_pos, rf_sig.identity_est)
                        } else {
                            (rf_call.split_pos, rf_call.identity_est)
                        };

                        let (decision, support_n, support_span) = if sig_ident >= cfg.fold2.min_identity {
                            dbg.bump("ident_ok");
                            
                            // Prefer evidence-based signature for support lookup.
                        let sig = sig_from_matches_pre
                            .or_else(|| {
                                // fallback: flank-based, using quantized split for consistency
                                let q = 50usize;
                                let split_q = (sig_split / q) * q;

                                foldback_signature(
                                    &job.seq,
                                    split_q,
                                    &cfg.shared,
                                    cfg.sig.flank_bp,
                                    cfg.sig.take,
                                    &mut sig_scratch,
                                    cfg.sig.value_shift
                                )
                            });
                            if let Some(sig) = sig {
                                dbg.bump("sig_ok");
                                if let Some(st) = support.get(&sig) {
                                    dbg.bump("sig_in_support");
                                    if hit_sig_sample.len() < 5 {
                                        hit_sig_sample.push((sig, st.n));
                                    }
                                    let real = is_real_foldback(sig, &support, &cfg);
                                    if real { num_real += 1; }
                                    (if real { "real" } else { "artefact" }, st.n, st.split_span())
                                } else {
                                    if miss_sig_sample.len() < 5 {
                                        miss_sig_sample.push(sig);
                                    }
                                    ("unknown", 0, 0)
                                }
                            } else {
                                ("unknown", 0, 0)
                            }
                        } else {
                            ("low_ident", 0, 0)
                        };

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

                out_tx.send(OutRec {
                    idx: job.idx,
                    id: id_str,
                    kept_seq,
                    tsv_row,
                }).map_err(|e| anyhow::anyhow!("out channel send failed: {e}"))?;
            }

            eprintln!(
                "[pass2 worker done] jobs={} detected={} refined={} ident_ok={} sig_ok={} in_support={} miss_sigs={:?} hit_sigs={:?}",
                dbg.jobs, dbg.detected, dbg.refined, dbg.ident_ok, dbg.sig_ok, dbg.sig_in_support,
                miss_sig_sample, hit_sig_sample
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

    // After workers exit, all their out_tx clones drop → out_rx closes → writer exits
    writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("writer thread panicked"))??;

    Ok((tot_fold, tot_cut, tot_real))
}
