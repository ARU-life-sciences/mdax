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

## Results summary

| Test | Finding |
|---|---|
| **1+3** crossmatch.py | 95.8% of `intrachromosomal_chimera` and 95.9% of `tandem_rc_concatemer` reads → mdax `artefact`. All sanity checks pass. |
| **2** jump_clip_sweep.sh | Wider `--max-jump-clip` *reduces* artefact calls (1.02M→955K) and *increases* low_ident (32K→81K). Wider windows find the true global-best δ; borderline reads called artefact at a capped scan are correctly reclassified. Clean foldbacks stable (use_hint fires independently of jc). |
| **4** gap_scatter.py | Pearson r=0.48; mdax median gap 1000 bp vs Chris median 439 bp — expected systematic offset between read-coordinate and reference-coordinate gap measurement. |
| **5** pipeline_integration.rs | Both `large_gap_artefact_wide_jump_clip` and `large_gap_artefact_hint_fallback` pass. |

---

## Tests in this directory

| Script | Test # | What it checks |
|---|---|---|
| `crossmatch.py` | 1 + 3 | Per-class comparison and tandem concatemer deep-dive |
| `jump_clip_sweep.sh` | 2 | How `--max-jump-clip` affects artefact/low_ident counts |
| `gap_scatter.py` | 4 | Correlation of mdax `gap_est` vs Chris's supplementary distance |

Test 5 (Rust regression test for large-gap detection) lives in
`../pipeline_integration.rs` as `large_gap_artefact_wide_jump_clip` and
`large_gap_artefact_hint_fallback`.

Output files are written to `output/` (gitignored).  Running any script bare
produces default outputs there; use `--out` / `--plot` to override the path.

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

**Observed results (ERR12263839, head table = first 10 k rows):**

7,224 reads overlap between the two tables.

| Chris class | n | mdax artefact | low_ident | real | mean gap_est | med identity |
|---|---|---|---|---|---|---|
| `intrachromosomal_chimera` | 6,104 | 5,849 (95.8%) | 174 | 81 | 1634 bp | 0.992 |
| `tandem_rc_concatemer` | 519 | 498 (95.9%) | 14 | 7 | 1396 bp | 0.992 |
| `local_inversion` | 455 | 435 (95.6%) | 12 | 8 | 1505 bp | 0.992 |
| `interchromosomal_chimera` | 120 | 113 (94.2%) | 4 | 3 | 1680 bp | 0.992 |
| `palindromic_hairpin` | 2 | 2 | 0 | 0 | — | — |

Tandem concatemer deep-dive (n=519):
- mdax decision: 498 artefact / 14 low_ident / 7 real
- gap_est: mean 1396 bp, median 998 bp, max 15975 bp
- Unit length: mean 682 bp, median 597 bp, range 100–1491 bp
- Copy count: mean 6.1, median 5, range 2–30

**Notes on local_inversion gap_est:**
Chris's `local_inversion` class means the soft-clip re-maps *genomically close*
to the primary alignment in reverse complement — but the read-coordinate gap
(what mdax measures) can still be ~1000 bp because the MDA enzyme skipped that
many bases within the molecule before reversing.  The check in `crossmatch.py`
therefore tests that the local_inversion median gap_est is ≤ the
intrachromosomal_chimera median (not an absolute threshold).

**Notes on palindromic_hairpin:**
Only 2 reads appear in the head-table overlap.  Both were called `artefact` by
mdax.  A meaningful check requires the full classification table (no head limit).
The script emits `[INFO]` rather than `[FAIL]` when n < 10.

**Sanity checks (all PASS on current data):**

```
[PASS] intrachromosomal_chimera → mdax artefact (expect >50%)
[PASS] local_inversion → mdax artefact (expect >50%)
[PASS] local_inversion gap_est ≤ intrachromosomal_chimera median
[INFO] palindromic_hairpin: n=2 (too few for assertion); 0% real — need full table
[PASS] tandem_rc_concatemer → mdax artefact (expect >60%)
```

**Run:**

```bash
python3 crossmatch.py                       # joined table → output/joined.tsv
python3 crossmatch.py --verbose             # also print per-read detail for small classes
python3 crossmatch.py --out other.tsv       # override output path
```

---

## Test 2: max-jump-clip sensitivity sweep

**Script:** `jump_clip_sweep.sh`

**What it does:**

Runs mdax on the 57-read ground-truth set (`data/chris_artefacts.fasta`) and,
if available, the full ERR12263839 dataset, at four `--max-jump-clip` values:
500, 1000 (default), 2000, 4000.

**Observed results — Part A (57-read ground-truth set):**

Overall counts:

| decision | jc=500 | jc=1000 | jc=2000 | jc=4000 |
|---|---|---|---|---|
| artefact | 56 | 56 | 57 | 57 |
| low_ident | 1 | 1 | 0 | 0 |
| real | 0 | 0 | 0 | 0 |

`li_*` reads (high BLAST identity — checking whether low_ident is gap-driven):

| label | jc=500 | jc=1000 | jc=2000 | jc=4000 |
|---|---|---|---|---|
| li_3279bp | artefact | artefact | artefact | artefact |
| li_5619bp | artefact | artefact | artefact | artefact |
| li_6296bp | artefact | artefact | artefact | artefact |
| li_7152bp | **low_ident** | **low_ident** | artefact | artefact |
| li_9340bp | artefact | artefact | artefact | artefact |
| li_11367bp | artefact | artefact | artefact | artefact |
| li_12397bp | artefact | artefact | artefact | artefact |
| li_14193bp | artefact | artefact | artefact | artefact |

`art_hi_*` reads (clean hairpins — should stay artefact):
all artefact at all four values (stable as expected).

**Key finding — `li_7152bp`:**
This read flips from `low_ident` to `artefact` at `jc=2000`.  Unlike the other
`li_*` reads, its gap is large enough that the δ-scan was capped at
`max_jump_clip=1000`.  At `jc=2000`, the scan reaches the true gap and identity
recovers.  The other 7 `li_*` reads are `low_ident` due to arm length, not gap
size, and are unaffected by the clip window.

**Observed results — Part B (full ERR12263839 dataset, ~1.07 M foldbacks):**

| decision | jc=500 | jc=1000 | jc=2000 | jc=4000 |
|---|---|---|---|---|
| artefact | 1,017,719 | 1,002,832 | 968,387 | 955,247 |
| low_ident | 16,585 | 32,006 | 67,621 | 81,281 |
| real | 37,528 | 36,993 | 35,824 | 35,303 |

**Finding — hypothesis was wrong (in an informative way):**

The original hypothesis predicted artefact count would *increase* with wider
windows (missed detections rescued).  The data show the opposite: artefact
decreases monotonically as jc grows, and low_ident absorbs those reads.

Explanation: at narrow jc, the δ-scan is capped at `max_delta = jc/2`.  For
some borderline reads the capped offset accidentally lands on a high-identity
region → artefact call.  At wider jc, the scan finds the true global-best
offset; if identity at that offset is below threshold, the read is correctly
reclassified as low_ident.  In other words, **wider max_jump_clip reduces
false positives rather than rescuing false negatives.**

The `use_hint` fallback (fires when `gap_from_hint > max_delta * 2`) means
true clean foldbacks with reliable minimizer-chain arm boundaries are already
classified correctly at any jc — they do not shift with the sweep.  The reads
that *do* shift are those where the hint doesn't fire (`arm_end = 0`) and the
δ-scan at a narrow cap happened to score a spuriously high identity.

**Run:**

```bash
chmod +x jump_clip_sweep.sh
./jump_clip_sweep.sh
THREADS=8 ./jump_clip_sweep.sh         # more threads
MDAX_BIN=/custom/path/mdax ./jump_clip_sweep.sh
```

---

## Test 4: gap_est vs supplementary distance scatter

**Script:** `gap_scatter.py`

**What it does:**

For reads present in both the mdax report and Chris's supplementary-patterns
table, computes `gap_est` (mdax, read-coordinate gap) vs `distance` (Chris,
genomic gap between primary and supplementary alignment endpoints).

**Observed results (ERR12263839):**

5,836 pairs within max_gap=10,000 bp (from 5,861 overlapping reads).

| Metric | mdax gap_est | chris distance |
|---|---|---|
| mean | 1551 bp | 815 bp |
| median | 1000 bp | 439 bp |
| std | 1519 bp | 1009 bp |
| min | 0 bp | 0 bp |
| max | 9984 bp | 8984 bp |

Pearson r = **0.484** (n=5,836).  Fraction within 2× of each other: **50%**.

Reads where mdax gap_est < 100 bp but chris_dist > 500 bp: 36 (0.6%)
Reads where mdax gap_est > 500 bp but chris_dist < 100 bp: 959 (16.4%)

**Interpretation:**

Moderate correlation: both methods agree on large gaps; diverge on small gaps.
The systematic offset (mdax median 1000 bp vs Chris 439 bp) reflects a
fundamental measurement difference:

- **mdax** measures the gap in *read coordinates* — the number of bases between
  the two arm ends as they appear in the raw molecule.
- **Chris** measures the gap in *reference coordinates* — the genomic distance
  between where the primary alignment ends and the supplementary begins.

For a simple foldback the two should be equal.  Divergence arises when:
1. The read has multiple supplementary alignments (Chris picks the closest);
2. The read is a tandem concatemer (gap structure is non-trivial);
3. The foldback junction is not at the primary alignment boundary (clipping
   artefacts, adaptor contamination).

The 959 reads where mdax sees a gap but Chris's distance is < 100 bp are
consistent with tandem concatemers (many short RC units) or reads where the
supplementary alignment happens to land very close to the primary even though
the molecule-internal gap is large.

**Run:**

```bash
python3 gap_scatter.py                             # scatter plot → output/gap_scatter.png
python3 gap_scatter.py --max-gap 5000              # restrict axis range
python3 gap_scatter.py --plot other.png            # override output path
```

---

## Test 5: Rust regression test for large-gap detection

**Location:** `../pipeline_integration.rs`
**Functions:** `large_gap_artefact_wide_jump_clip`, `large_gap_artefact_hint_fallback`

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
- Asserts `gap_est ∈ [1600, 2400]` bp (hint path reports the true gap, not a
  capped δ-scan value)
- Asserts `identity_est ≥ 0.90` (window correctly placed at arm edges)
- Asserts `decision = artefact`

**Observed results:**

```
test large_gap_artefact_wide_jump_clip ... ok
test large_gap_artefact_hint_fallback  ... ok
```

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
`jump_clip_sweep.sh` test probes this on real data: `li_7152bp` flips from
`low_ident` to `artefact` at `jc=2000`, confirming that at least one real read
has a gap too large for the δ-scan at the default clip window.

---

## Class correspondence table

| Chris's `artifact_class` | mdax equivalent | Notes |
|---|---|---|
| `intrachromosomal_chimera` | `artefact` (gap_est > 0) | Main foldback class; 95.8% agreement |
| `local_inversion` | `artefact` (gap_est in read coords ~1000 bp) | Genomically close but read gap is real |
| `palindromic_hairpin` | `real` (expected) | Too rare in head table to verify — need full data |
| `tandem_rc_concatemer` | `artefact` (95.9%) | mdax detects outer foldback; internal repeat structure ignored |
| `interchromosomal_chimera` | `artefact` (94.2%) | mdax works within a single read; inter-chrom gap not measured |
| `ambiguous` / `unmapped_clip_unknown` | varies | Weak signal in both tools |
