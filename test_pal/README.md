# test_pal/

Comparison of `mdax` against the Strobl et al. (2023) Perl-based foldback baseline, using a 20,000-read ONT subsample from a real WGA library.

**Reference**: Strobl et al. (2023), *Nucleic Acids Research* — https://doi.org/10.1093/nar/gkad573

## Data

Input data is **not tracked in git** (stored in `out/`, gitignored):

| File | Description |
|------|-------------|
| `out/subsample/SRR24201687.first20k.fastq.gz` | First 20,000 ONT reads from SRR24201687 |
| `out/baseline/SRR24201687_first20k/` | Baseline outputs from Strobl et al. Perl script (iterations 1–2, flagged IDs, breakpoints) |

**Accession**: SRR24201687 (NCBI SRA)

## Running

```bash
# single run with default flags
./test_pal.sh

# full parameter sweep (cached; delete sweep_cache/ to re-run)
./sweep.sh

# specific sweep groups
./sweep.sh modes
./sweep.sh refine
./sweep.sh identity
```

## How it works

1. `test_pal.sh` canonicalises baseline IDs (strips `-dup…` suffixes), runs `mdax`, and normalises output into `mdax.calls.norm.tsv`.
2. `test_pal.py` computes confusion-matrix metrics against the baseline:
   - **Strict**: `decision == "artefact"` only.
   - **Generous**: `decision == "artefact"` or `"low_ident"`.
3. Breakpoint accuracy: absolute and signed delta between `refined_split` and the baseline iteration-2 breakpoints.

## Results

Metrics are generous (artefact + low_ident both counted as true positive) unless noted.

| Configuration | Precision | Recall | F1 |
|---|---|---|---|
| `balanced --refine-mode hifi` (default, `test_pal.sh`) | 0.998 | 0.647 | 0.785 |
| `permissive --refine-mode ont --cut-low-ident` (best ONT, from sweep) | **0.995** | **0.803** | **0.889** |

The default run uses balanced/HiFi mode with no additional flags — this is what `test_pal.sh` executes. The sweep explores all combinations; see `sweep_summary.tsv` for machine-readable per-configuration metrics.

## Output files (produced by test_pal.sh)

| File | Description |
|------|-------------|
| `mdax.report.tsv` | Full mdax TSV report |
| `mdax.out.fa` | Corrected FASTA output |
| `mdax_calls.tsv` | Normalised mdax calls (id, refined_split, decision) |
| `baseline.flagged.canon.ids` | Canonicalised baseline positive IDs |
| `baseline.bp.tsv` | Baseline breakpoints (id, position) |
| `mdax_only.ids` | Reads called by mdax but not baseline |
| `baseline_only.ids` | Reads called by baseline but not mdax |
| `sweep_cache/` | Cached per-configuration metrics (one TSV per run) |
| `sweep_summary.tsv` | Machine-readable table of all sweep configurations (written by `sweep.sh`) |

---

## Interpretation

This benchmark validates `mdax` on ONT data against an independent published tool (Strobl et al. 2023) and shows how configuration choices affect recall vs precision.

**ONT requires different settings to HiFi.** The default HiFi configuration (`--mode balanced --refine-mode hifi`) achieves ~65% recall on this ONT dataset. This is substantially lower than the ~80% recall of ONT mode — the Hamming-based refinement is tuned for the low error rate of HiFi CCS reads and places split positions less accurately in noisier ONT sequence, causing some reads to fall below the identity threshold. Switching to `--refine-mode ont` (banded Levenshtein) raises generous recall from 0.647 to 0.803.

**Recommended ONT flags.** `--mode permissive --refine-mode ont --cut-low-ident` achieves the best balance (F1 ≈ 0.889, precision 0.995, recall 0.803). The `--mode permissive` flag lowers the support-map threshold, accepting weaker evidence that is typical in noisier ONT data. `--cut-low-ident` ensures reads with a confirmed foldback but lower arm identity are still corrected rather than flagged and kept.

**Precision is consistently high.** Across all tested configurations precision is ≥ 0.995 — the tool almost never calls a false positive. The main lever for improving performance on ONT is recall, not false-positive suppression.

**Breakpoint accuracy.** With HiFi mode, median absolute breakpoint delta vs Strobl et al. iteration-2 is ~823 bp (mean 1716 bp) — acceptable for gross artefact removal but noisier than ONT mode. ONT refinement places splits closer to the structural boundary.
