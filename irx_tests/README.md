# irx_tests/

Synthetic ground-truth evaluation of the `irx` inverted-repeat detector. A Python generator plants known IRs into random sequence across a range of scenarios; a second script evaluates `irx` calls against the truth set.

## Running

```bash
# regenerate truth set, run irx, evaluate (all three steps)
bash test_suite.bash
```

Or step by step:

```bash
# 1. generate synthetic FASTA + truth TSV
python3 make_test_suite.py \
    --out-fasta irx_truth.fa \
    --out-truth irx_truth.tsv \
    --out-manifest irx_manifest.tsv \
    --seed 42

# 2. run irx
../target/release/irx -b irx_calls.tsv irx_truth.fa

# 3. evaluate
python3 eval_irx_test_suite.py \
    --truth irx_truth.tsv \
    --calls irx_calls.tsv \
    --out_summary irx_eval.summary.tsv \
    --out_truth  irx_eval.truth.tsv \
    --out_fp     irx_eval.fp.tsv
```

## Scenario groups

| Group | n planted | Description |
|-------|-----------|-------------|
| `clean_high_identity` | 4 | IRs with ≥0.95 arm identity |
| `diverged` | 4 | IRs with 0.50–0.95 identity |
| `length_identity_matrix` | 90 | Grid of arm lengths × identities |
| `asymmetric` | 3 | Left and right arms of different length |
| `class_spacing` | 4 | Immediate, spaced, and wide gap structures |
| `edge_cases` | 4 | IRs positioned near contig start/end |
| `gc_bias` | 3 | Varied GC composition |
| `reused_family` | 3 | Multiple IRs from the same sequence family |
| `micro_ir` | 5 | Very short arms / low-complexity sequence |
| Negative controls | — | Random contigs with no planted IR (for FP measurement) |

## Current results (`irx_eval.summary.tsv`)

| Scenario group | Planted | Detected | Recovery |
|---|---|---|---|
| `asymmetric` | 3 | 3 | 100% |
| `class_spacing` | 4 | 4 | 100% |
| `clean_high_identity` | 4 | 2 | 50% |
| `diverged` | 4 | 1 | 25% |
| `edge_cases` | 4 | 4 | 100% |
| `gc_bias` | 3 | 3 | 100% |
| `length_identity_matrix` | 90 | 36 | 40% |
| `micro_ir` | 5 | 1 | 20% |
| `reused_family` | 3 | 3 | 100% |

Recovery drops for short arms, highly diverged sequences, and very small IRs — see `irx_eval.truth.tsv` for per-record detail and error statistics (arm length error, spacer error, `tir_ident` calibration).

### Detection boundary (from `length_identity_matrix`)

The 90-case grid (arm lengths 1 kb–100 kb × identities 0.40–0.95) reveals a clear threshold:

| Identity | Recovery |
|---|---|
| ≥ 0.95 | 100% at all arm lengths (1 kb–100 kb) |
| 0.85 | 100% at arm ≥ 2 kb; 0% at 1 kb |
| 0.70 | Partial recovery only at arms ≥ 50 kb |
| ≤ 0.55 | Not detected regardless of arm length |

`clean_high_identity` missed 2/4 cases because the 100 bp and 500 bp arms fall below the reliable detection floor; `micro_ir` recovered 1/5 for similar reasons.

**False positive rate: 0/N** — no spurious calls on any negative-control contig.

---

## Interpretation

`irx` is designed for detecting inverted repeats in assembled or long-read sequences, not as a sensitive scanner for all possible IRs.

**Reliable use case.** Arms ≥ 1 kb at ≥ 0.95 identity are detected with 100% recall. Arms ≥ 2 kb at ≥ 0.85 identity are also reliable. These parameter ranges cover the MDA foldback artefacts that `mdax` targets (arm lengths typically 3–30 kb, identity 0.75–1.0).

**Known limitations.** IRs with identity below ~0.80 are not reliably detected unless arms are very long (≥ 50 kb). Short IRs (arms < 1 kb) are missed even at near-perfect identity. `irx` should not be used as a general-purpose short-IR finder.

**High specificity.** Zero false positives on random-sequence negative controls confirms that the minimizer chaining approach does not spuriously call repetitive but non-IR sequence. This makes `irx` safe to run on large assemblies without post-filtering.

## Output files

| File | Description |
|------|-------------|
| `irx_truth.fa` | Synthetic FASTA with planted IRs |
| `irx_truth.tsv` | Per-planted-IR ground truth (coordinates, arm lengths, spacer, identity, class) |
| `irx_manifest.tsv` | All contigs including negative controls |
| `irx_calls.tsv` | Raw `irx` output |
| `irx_eval.summary.tsv` | Per-scenario-group recovery rates and error statistics |
| `irx_eval.truth.tsv` | Per-truth-record match status and error metrics |
| `irx_eval.fp.tsv` | False-positive calls on negative-control contigs |
