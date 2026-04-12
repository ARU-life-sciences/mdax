# test_celegans/

Benchmark `mdax` on two real single-worm *C. elegans* MDA libraries with contrasting assembly quality and artefact rates. No external ground truth — the primary comparison is artefact rate between the two samples.

## Samples

Data is **not tracked in git** (stored in `data/`, gitignored):

| Sample | File | Size | BUSCO | Contigs | Protocol | Notes |
|--------|------|------|-------|---------|----------|-------|
| `nrCaeEleg92` | `data/nrCaeEleg92_ccs.trim.fasta.gz` | 2.2 GB | 82.3% | 1,917 | Repli-G DLB rep 2, 120 min | moderate artefact rate |
| `nrCaeEleg95` | `data/nrCaeEleg95_ccs.trim.fasta.gz` | 8.6 GB | 64.5% | 7,944 | Repli-G DLB 120 min | high artefact rate |

Both are PacBio CCS (HiFi) libraries. The lower BUSCO score in nrCaeEleg95 is consistent with a higher proportion of MDA foldback artefacts displacing true genomic content.

## Running

```bash
./benchmark.sh                  # process up to 200,000 reads per sample (fast, representative)
MAX_READS=0 ./benchmark.sh      # process all reads (slow for nrCaeEleg95)
THREADS=16  ./benchmark.sh      # override thread count (default: 8)
```

## Scripts

| Script | Description |
|--------|-------------|
| `benchmark.sh` | Runs `mdax --cut-low-ident` on each sample, produces per-sample `.report.tsv` and `.out.fasta`, then calls `report.py`. |
| `report.py` | Summarises and compares mdax outputs across samples. Reports artefact rates, identity distributions, read length stats before/after correction. |
| `dotplots.sh` | Samples detected foldback reads by identity bin, extracts sequences, and generates self-dotplots using `flexidot`. |
| `probe_identity.py` | Per-position arm match-rate analysis — plots base-by-base palindromic identity decaying from the junction outward. |

## Key results (200k reads per sample, `--cut-low-ident`)

| | nrCaeEleg92 (82.3% BUSCO) | nrCaeEleg95 (64.5% BUSCO) |
|---|---|---|
| Foldback detected | 85,508 (42.8%) | 77,249 (38.6%) |
| → artefact cut | 85,137 (42.6%) | 76,714 (38.4%) |
| → real palindrome | 371 (0.2%) | 535 (0.3%) |
| Bases removed | ~422 Mb | ~396 Mb |
| Fraction of artefact read retained | ~49.8% | ~50.0% |

The higher BUSCO sample (nrCaeEleg92) shows a slightly higher artefact detection rate, consistent with more uniform coverage enabling better support-map resolution.

## Output files (produced by benchmark.sh)

| File | Description |
|------|-------------|
| `out/<sample>.report.tsv` | Full mdax TSV report per sample |
| `out/<sample>.out.fasta` | Corrected FASTA per sample |
| `out/summary.tsv` | Machine-readable per-sample summary (written by `report.py --out-tsv`) |
| `dotplots/` | Self-dotplot PNGs organised by identity bin |

---

## Interpretation

This benchmark provides a real-world scale check on two *C. elegans* MDA libraries with contrasting quality, without external ground truth. The primary signal is the difference in artefact rates between samples.

**Artefact rates are high and expected.** Both samples show ~40% artefact rates, consistent with published MDA artefact estimates and with the BUSCO scores (lower BUSCO = more artefact reads displacing real genomic content). The nrCaeEleg95 sample (64.5% BUSCO) has a somewhat lower detection rate than nrCaeEleg92 (82.3% BUSCO), which is counter-intuitive at first glance — higher artefact burden should yield more detections. The explanation is that uniform coverage (reflected in the higher BUSCO score) produces denser support maps, allowing weaker evidence to be resolved correctly. In highly contaminated libraries, many reads compete for the same support positions.

**Bases removed at scale.** Approximately 400–420 Mb of sequence is removed per 200 k-read subsample, corresponding to ~50% of each detected artefact read. This matches the theoretical expectation for a foldback where the cut retains the left arm (roughly half the read length).

**Dotplots confirm structural classes.** The `dotplots/` directory contains self-dotplots of sampled foldback reads stratified by identity bin. The expected diagonal + anti-diagonal structure is visible at all identity levels down to 0.30–0.40, providing visual confirmation that detections at low identity are genuine foldback artefacts rather than false positives.

**Practical guidance for real MDA datasets.** Use `--cut-low-ident` to recover all foldback classes including diverged ones. Process in parallel (`--threads N`); at default settings `mdax` completes 200 k HiFi reads in ~2 minutes on an 8-core machine. The `out.fasta` file is the corrected input for downstream assembly or alignment.
