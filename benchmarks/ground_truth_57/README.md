# ground_truth_57/

mdax vs. CADECT scored against real, independently-curated ground truth
(`data/chris_artefacts.fasta` / `chris_artefacts_truth.tsv`, 57 HiFi CCS
reads: 35 true foldback artefacts, 22 genuine reads, from `ERR12263839`).
Unlike `benchmarks/cadect_vs_mdax/` (60k-read microlaimidae subsample,
tool-vs-tool agreement only, no ground truth), this set has independently
assigned true/false labels, so precision/recall are real, not proxies.

## Commands

```
mdax data/chris_artefacts.fasta -r mdax_report.tsv -o mdax_corrected.fa
conda run -n cadect python3 CADECT/CADECTv1.2.py -R data/chris_artefacts_raw.fastq -w 500 -o cadect_out
python3 compare_truth.py
```

mdax run with defaults (`--mode balanced`, min_support=3, auto-calibrated
down to 2 by `SUPPORT CALIBRATION` for this dataset's observed recurrence).

## Results

| | precision | recall | F1 | FPR on genuine | FNR on artefact |
|---|---|---|---|---|---|
| mdax | 0.630 | **1.000** | 0.773 | **90.9%** (20/22) | 0.0% (0/35) |
| CADECT | 0.600 | 0.257 | 0.360 | 27.3% (6/22) | 74.3% (26/35) |

mdax catches every true artefact (recall 1.0) but cuts 20 of the 22
genuine reads too. CADECT misses three-quarters of true artefacts but is
much less likely to cut a genuine read.

## Why mdax's precision is this low here — read before citing this number

This is the mechanism `synth_bench/README.md` warned about, caught live:
mdax's only evidence that a foldback is *real* genomic structure is
cross-read recurrence of the same junction signature (`min_support`).
`[SUPPORT UNIQNESS]` for this run: 486 clusters, 494 supporting reads,
**avg_n=1.016, max_n=2** — i.e. almost every junction in this dataset is
a singleton. That's not a flaw in mdax's logic being exposed by real
data; it's a property of how `chris_artefacts.fasta` was built: one
curated example per category (`art_hi_*`, `real_lo_*`, etc.), not deep
multi-read coverage of the same genomic locus. The cross-read-support
gate *cannot* work by construction on a set with ~1 read per true locus,
genuine or artefact. On real whole-genome WGA/MDA coverage, a true
inverted repeat is typically sampled by many reads if coverage is
adequate — this 57-read set doesn't test that regime at all.

So: this result is a legitimate, real demonstration that mdax's
precision degrades badly when per-locus coverage is low/singleton
(consistent with the colleague-reported assembly regressions motivating
`synth_bench/`), but it is **not** a general "mdax has 91% false-positive
rate" claim — report it as "on a low-coverage/singleton-replicate
dataset" specifically, and corroborate with `synth_bench`'s controlled
coverage-tier sweep (`cov_tier` column in `per_locus.tsv`) before using
it as a headline number anywhere.

## Caveats

- n=57 is small; CIs are wide.
- CADECT run with default `-w 500`; not tuned per-dataset.
- `chris_artefacts.fasta` was used to help build/tune `is_real_foldback`
  historically (see `data/README.md`) — not a fully held-out test set for
  mdax, though it was not used for CADECT.
- 1 mdax read classified `low_ident` (excluded from the confusion matrix,
  counted under "uncalled/other"); counting it as artefact doesn't change
  recall (already 1.0) or precision materially (34/35 -> 34/36 if forced
  artefact).
