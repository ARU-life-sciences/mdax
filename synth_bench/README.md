# synth_bench/

Synthetic-genome, ground-truthed benchmark for testing mdax's effect on
**downstream assembly contiguity** -- not just per-read classification
accuracy. `test_pal/` and `test_celegans/` validate whether mdax's
artefact *calls* agree with a baseline; nothing in the repo previously
checked whether *acting on those calls* helps or hurts assembly. This
directory closes that gap.

> The motivation below describes the pre-v0.2.0 locus-signature gate; the fold-axis gate has since replaced it, and `--cut-low-ident` is now a no-op.

## Other scripts

| Script | Purpose |
|---|---|
| `truth_regress.sh submit BIN NAME` / `check NAME` | Full-library truth regression on four sets (farm); exits 1 if below floors |
| `refbased_recall_eval.py`, `celegans_refbased_eval.py` | Reference-based foldback truth and scoring |
| `measure_residual.sh`, `residual_location.py` | Foldbacks left in cut reads |
| `noir_*.py` | Foldbacks with no inverted copy: junction geometry, k-mer dip, oracle cuts |
| `retest_*_arm.sh`, `retest_summarise.sh` | hifiasm / Flye + BUSCO re-assembly arms |

## Motivation

Colleagues observed that mdax-purged HiFi datasets produced **worse**
assemblies (lower N50, less contiguous) than either the raw uncorrected
reads or a reference-based chimera-removal method. Code review
(`src/fingerprint.rs::is_real_foldback`) suggests why: mdax's only
evidence for "this foldback is real genomic structure, don't cut it" is
whether the same junction signature recurs in `>= min_support` reads
(default 3). Real inverted repeats in low/uneven-coverage regions of an
MDA library can easily fail to reach that bar by chance, in which case
mdax cuts them in *every* read that covers them -- deterministically, not
occasionally -- specifically at the repeat boundaries where long-read
contiguity matters most.

This harness builds a synthetic genome with real inverted repeats at
**controlled, known** coverage tiers, arm lengths and identities, so we
can test that mechanism directly and measure its effect on an actual
hifiasm assembly, with full ground truth.

## Pipeline

```
make_genome.py       -> genome.fa + repeats_truth.bed + coverage_profile.tsv
simulate_reads.py     -> reads.fa + reads_truth.tsv
  (clean reads sampled with per-locus coverage-tier bias;
   ~20% converted to synthetic MDA foldback chimeras, independent of the
   real repeat loci, with identity/gap classes modelled on data/README.md)
mdax (several configs)  -> corrected/<arm>.fa
refbased_clip.py         -> refbased/reads.refbased.fa
  (reference-based proxy: minimap2 to the truth genome, clip reads whose
   alignment blocks overlap in reference space on opposite strands --
   the foldback signature -- down to the primary block. Never touches
   real repeat-spanning reads by construction, since those map as one
   contiguous linear alignment.)
hifiasm                  -> asm/<arm>.p_ctg.fa   (per arm)
evaluate.py               -> eval/summary.tsv, eval/per_locus.tsv
```

`run_all.sh` drives the whole thing end to end.

## Running

```bash
# quick sanity run (~1-2 min, small genome)
./run_all.sh run_small 300000 20 4

# larger, more statistically meaningful run
./run_all.sh run_main 3000000 30 4
```

Args: `outdir genome_len target_coverage threads`. Each stage is skipped
if its output already exists, so you can tweak one script and re-run.
`MDAX_BIN`, `MINIMAP2`, `HIFIASM`, `GFATOOLS`, `SEED` are overridable via
env vars (defaults point at `/software/team301/{minimap2,hifiasm,gfatools}`
and `target/release/mdax`).

## Arms compared

| Arm | Meaning |
|---|---|
| `raw` | no correction |
| `mdax_balanced` / `mdax_balanced_cutlow` | mdax defaults |
| `mdax_strict` | mdax `--mode strict` (min_support=5) |
| `mdax_permissive` | mdax `--mode permissive --cut-low-ident` (min_support=2, lower evidence thresholds) |
| `mdax_support1` | mdax with `--min-support 1` (accept any recurrence-of-1 as real) |
| `mdax_support8` | mdax with `--min-support 8` (demand heavy recurrence) |
| `refbased` | reference-based clipping proxy |

## Reading the output

`eval/summary.tsv` / stdout table: per-arm N50, contig count, genome
fraction covered, and IR-locus spanned/broken/missing counts.

`eval/per_locus.tsv`: one row per real inverted-repeat locus (with its
true arm length / gap / identity / coverage tier) x one column per arm,
value `spanned` / `broken` / `missing`. **This is the mechanistic
evidence**: filter to rows where `raw` and `refbased` say `spanned` but an
`mdax_*` arm says `broken`/`missing` -- those are loci where mdax
specifically destroyed real genomic structure that every other method
preserved. Cross-reference against `cov_tier` and `min_support` used by
each arm to confirm the coverage-dependence mechanism.

## Caveats / what this harness does NOT model

- Coverage bias is a simple per-locus multiplier, not a realistic MDA
  amplification-bias model (real bias is more extreme and spatially
  correlated). Treat `cov_tier` as a controlled knob, not a claim about
  real MDA statistics.
- The reference-based comparator is a cheap proxy (minimap2 split-alignment
  clipping), not a specific published tool. It is deliberately simple so
  its *false positive rate on real repeats is provably ~0 by construction*
  -- which is the property that matters for this comparison -- but its
  recall on chimeras is a rough approximation.
- Single haploid genome, uniform background GC, no long-range structural
  repeats beyond what's explicitly inserted. This isolates the foldback/
  inverted-repeat mechanism cleanly but isn't a realism claim about any
  particular organism.
- `hifiasm` is run with default settings (no purge-dups tuning etc.);
  absolute N50 numbers are only meaningful *relative to other arms on the
  same simulated dataset*, not as absolute assembly-quality benchmarks.
