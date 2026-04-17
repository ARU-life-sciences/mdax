# tests/chris_crossvalidation/

Cross-validation of mdax against Chris Laumer's independent soft-clip pipeline
for ERR12263839 (PacBio HiFi / CCS, MDA-amplified, *C. elegans*).

Chris's pipeline aligns reads to the reference (minimap2), extracts soft-clipped
sequences (≥ 200 bp), re-maps them to the genome, and classifies each clip by
where it lands relative to the primary alignment.  This is a completely
independent approach to mdax — reference-based vs. reference-free — making it
ideal for cross-validation.

See `../../chris_data/README.md` for dataset description and the key finding that
the modal template-switch gap is 660–1320 bp, with ~61% of reads having gap > 660 bp.

---

## Tests in this directory

| Script | Test # | What it checks |
|---|---|---|
| `crossmatch.py` | 1 + 3 | Per-class comparison and tandem concatemer deep-dive |
| `jump_clip_sweep.sh` | 2 | How `--max-jump-clip` affects artefact/low_ident counts |
| `gap_scatter.py` | 4 | Correlation of mdax `gap_est` vs Chris's supplementary distance |

Test 5 (Rust regression test for large-gap detection) lives in
`../pipeline_integration.rs` as `large_gap_artefact_wide_jump_clip` and
`large_gap_artefact_narrow_jump_clip`.

---

## Prerequisites

```bash
# Build the release binary (required for jump_clip_sweep.sh)
cargo build --release

# Python packages (required for crossmatch.py and gap_scatter.py)
pip install matplotlib numpy    # only for gap_scatter.py --plot
# No extra packages needed for crossmatch.py
```

The data files live in `../../chris_data/` (gitignored — large files).
All Python scripts exit cleanly with a SKIP message if the files are absent.

---

## Test 1 + 3: Read-level crossmatch and tandem concatemer analysis

**Script:** `crossmatch.py`

**What it does:**

Joins the mdax report (`ERR12263839.tsv`) with Chris's per-clip classification
(`artifact_classifications_head.tsv`) on read ID.  For each of Chris's artifact
classes it reports:

- How many reads appear in both tables
- Distribution of mdax decisions (artefact / low_ident / real)
- Mean/median `gap_est` and `identity_est` from mdax
- Tandem concatemer unit sizes and copy counts from Chris's `detail` field

**Expected behaviour:**

| Chris class | Expected mdax decision | Reasoning |
|---|---|---|
| `intrachromosomal_chimera` | mostly `artefact` | gap > 0, mdax δ-scan should find it |
| `local_inversion` | mostly `artefact`, small `gap_est` | clip maps near split, gap ≈ 0 |
| `palindromic_hairpin` | mostly `real` | recurring junction → support map → real |
| `tandem_rc_concatemer` | mostly `artefact` | mdax sees the overall foldback |

**Run:**

```bash
python3 crossmatch.py
python3 crossmatch.py --verbose             # show per-read detail for small classes
python3 crossmatch.py --out joined.tsv      # save joined table
```

---

## Test 2: max-jump-clip sensitivity sweep

**Script:** `jump_clip_sweep.sh`

**What it does:**

Runs mdax on the 57-read ground-truth set (`data/chris_artefacts.fasta`) and,
if available, the full ERR12263839 dataset, at four `--max-jump-clip` values:
500, 1000 (default), 2000, 4000.

**Hypothesis (from gap distribution analysis):**

The modal gap bin is 660–1320 bp.  The default `--max-jump-clip = 1000` means
`max_delta = 500`, which can only handle gaps up to ~1000 bp.  For the full
dataset we expect:

- artefact count to increase from `max_jump_clip=1000` → `2000` (reads in the
  modal bin rescued)
- diminishing returns above `3000 bp`

For the 57-read test set, most `li_*` reads are expected to STAY `low_ident`
across all values — their low identity is caused by a very long arm
(arm_len 5000–9000 bp) being compared over a short refine window (1200 bp),
**not** by a template-switch gap.  The sweep confirms whether max_jump_clip is
the cause.

**Run:**

```bash
chmod +x jump_clip_sweep.sh
./jump_clip_sweep.sh
./jump_clip_sweep.sh                   # use default binary (../../target/release/mdax)
THREADS=8 ./jump_clip_sweep.sh         # more threads
MDAX_BIN=/custom/path/mdax ./jump_clip_sweep.sh
```

---

## Test 4: gap_est vs supplementary distance scatter

**Script:** `gap_scatter.py`

**What it does:**

For reads present in both the mdax report and Chris's supplementary-patterns
table, plots `gap_est` (mdax, read-coordinate gap) against `distance` (Chris,
genomic-coordinate gap between primary and supplementary alignments).

If both methods measure the same template-switch gap, the scatter should lie
close to the y=x line with high Pearson correlation.  Systematic divergence
indicates:

- Reads with multiple supplementary alignments (Chris picks the closest one,
  which may not correspond to the primary foldback junction)
- Tandem concatemers (gap structure is non-trivial; neither measure is correct)
- Reads where `gap_est > max_jump_clip` (mdax underestimates, Chris does not)

**Run:**

```bash
python3 gap_scatter.py
python3 gap_scatter.py --plot gap_scatter.png      # save a scatter plot
python3 gap_scatter.py --plot gap_scatter.png --max-gap 5000
```

---

## Test 5: Rust regression test for large-gap detection

**Location:** `../pipeline_integration.rs`
**Functions:** `large_gap_artefact_wide_jump_clip`, `large_gap_artefact_narrow_jump_clip`

**What it does:**

Creates a deterministic synthetic read:

```
arm (1500 bp random) | gap (2000 bp random) | RC(arm) (1500 bp)   →  5000 bp total
```

The true template-switch gap is exactly 2000 bp.  The delta-scan needs
`best_delta ≈ 1000 = gap/2`, so requires `max_jump_clip ≥ 2000`.

**`large_gap_artefact_wide_jump_clip`** (`max_jump_clip = 4000`):
- Asserts `decision = artefact`
- Asserts `gap_est ∈ [1600, 2400]` bp (within ±400 bp of truth)
- Asserts `identity_est ≥ 0.90` (arms are identical — should be near 1.0)

**`large_gap_artefact_hint_fallback`** (`max_jump_clip = 800`):
- Verifies the `use_hint` fallback: when `gap_from_hint > max_delta * 2`, mdax
  detects that the δ-scan can't reach the arm boundaries, and places the
  comparison window directly on arm boundaries from the minimizer chain.
- Asserts `gap_est ≈ 2000 bp` (hint path reports the true gap from arm boundaries,
  not the capped δ-scan value)
- Asserts `identity_est ≥ 0.90` (window correctly placed at arm edges)
- Asserts `decision = artefact`

**Run:**

```bash
cargo test large_gap_artefact
```

**Key finding — the hint fallback:**

The `use_hint` condition fires when `gap_from_hint > max_delta * 2`, skipping
the δ-scan and placing the comparison window directly on minimizer-chain arm
boundaries.  This means mdax is robust to small `max_jump_clip` **when the
coarse detector produces reliable arm boundaries**.

`max_jump_clip` still matters for reads where the minimizer chain is sparse
(noisy reads, short arms, low coverage) and `arm_end = 0` — in those cases
`use_hint` does not trigger and the δ-scan is used with its cap.  The
`jump_clip_sweep.sh` test probes this on real data.

---

## Class correspondence table

| Chris's `artifact_class` | mdax equivalent | Notes |
|---|---|---|
| `intrachromosomal_chimera` | `artefact` (gap_est > 0) | Main foldback class |
| `local_inversion` | `artefact` (gap_est ≈ 0) | Clean hairpin |
| `palindromic_hairpin` | `real` | Recurring junction, not cut |
| `tandem_rc_concatemer` | `artefact` | mdax detects the outer foldback |
| `interchromosomal_chimera` | not detected | mdax works within a single read |
| `ambiguous` / `unmapped_clip_unknown` | varies | Weak signal in both tools |
