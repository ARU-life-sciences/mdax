# test_single_nematode/

Benchmarks `mdax` against the palindrome-correction method from Lee et al.'s
single-worm nematode genome project, on real ONT whole-genome-amplification
(WGA) data. Compares three arms per species — **raw** (no correction),
**theirs** (`remove_palindrome.sh`, the paper's own method), **mdax** — all
assembled identically downstream, so the only variable is the correction
method itself.

## Source data

- Paper/pipeline: [Single-nematode-genome](https://zenodo.org/records/8144614)
  (Yi-Chien Lee et al., Academia Sinica) — a Zenodo-archived academic
  website containing their `remove_palindrome.sh` pipeline (copied into
  `scripts/`, see [Their method](#their-method-remove_palindromesh) below)
- Raw ONT reads: NCBI BioProject
  [PRJNA953805](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA953805)
  ("Free-living nematodes Genome sequencing and assembly"), 13 species total
- Published assemblies (ground truth): NCBI GenBank, same BioProject,
  `assembly_method: Flye v.2.9.1`, pipeline = Flye → Racon/Medaka polish →
  NextPolish (Illumina) → `purge_dups` → contamination screen. All at
  `assembly_level: Contig` (no chromosome-scale scaffolding) — directly
  comparable in kind to what this benchmark produces, just with additional
  polishing/dedup/decontam steps this benchmark deliberately omits (see
  [Caveats](#caveats)).

4 of the 13 species were run (picked for small ONT read size, plus one
different clade for spread): Epsilonema sp. ZAB3_2, Microlaimidae sp. YZB2_3,
Linhomoeus sp. GSCO2_2 (all Chromadorea), Trileptium ribeirensis (Enoplida).

## Method

Per species: `remove_palindrome.sh` (2 iterations, self-vs-self
`minimap2 -x ava-ont` + midpoint chop) and `mdax --mode permissive
--refine-mode ont --cut-low-ident` run independently on the same raw reads;
each of the three resulting read sets (raw/theirs/mdax) assembled with
`Flye --nano-raw --meta` (matching the paper's assembler choice), scored
with BUSCO `nematoda_odb10` (genome mode). See `scripts/run_species_pipeline.sh`
for the full LSF-orchestrated pipeline.

## Results

BUSCO nematoda_odb10, genome mode. Bold = best value per metric within a
species (excluding the mdax-permissive anomaly on microlaimidae — see below).

| Species | Arm | Complete | Single | Dup | Frag | Missing | Contigs | N50 | Total length |
|---|---|---|---|---|---|---|---|---|---|
| **Epsilonema sp.** | raw | 46.8% | 35.4% | 11.4% | 5.8% | 47.4% | 14,179 | 25 kb | 204.9 Mb |
| | theirs | **49.3%** | **38.1%** | 11.2% | 5.6% | **45.1%** | 9,934 | **37 kb** | 182.1 Mb |
| | mdax | 48.3% | 38.0% | **10.3%** | 5.7% | 46.0% | **9,813** | 35 kb | **180.7 Mb** |
| **Microlaimidae sp.** | raw | 16.1% | 15.7% | 0.4% | 5.2% | 78.7% | 45,291 | 20 kb | 630.7 Mb |
| | theirs | **20.1%** | **19.6%** | **0.5%** | 5.0% | **74.9%** | **42,179** | **21 kb** | **605.1 Mb** |
| | mdax (permissive) ⚠️ | 11.9% | 6.4% | 5.5% | 3.2% | 84.9% | 40,978 | 34 kb | 1,084.4 Mb |
| **Linhomoeus sp.** | raw | 35.8% | 28.4% | 7.4% | 3.9% | 60.3% | 21,289 | 20 kb | 254.3 Mb |
| | theirs | **36.2%** | **29.1%** | 7.1% | 3.8% | **60.0%** | 20,795 | 19 kb | 246.4 Mb |
| | mdax | 35.9% | 29.0% | **6.9%** | **3.6%** | 60.5% | **20,421** | 19 kb | **242.0 Mb** |
| **Trileptium ribeirensis** | raw | 21.6% | 20.6% | 1.0% | **3.4%** | 75.0% | 42,276 | 34 kb | 813.0 Mb |
| | theirs | **22.9%** | **21.8%** | 1.1% | 3.8% | **73.3%** | 39,932 | **37 kb** | 799.6 Mb |
| | mdax | 22.0% | 21.0% | 1.0% | **3.4%** | 74.6% | **38,771** | 36 kb | **786.1 Mb** |

### Ground truth (published GenBank assemblies, same BioProject)

| Species | Accession | Total length | Contigs | N50 |
|---|---|---|---|---|
| Trileptium ribeirensis | GCA_030248225 | 746.0 Mb | 30,783 | 42.0 kb |
| Linhomoeus sp. | GCA_030248115 | 240.5 Mb | 14,458 | 28.9 kb |
| Microlaimidae sp. | GCA_030248175 | 575.8 Mb | 35,431 | 24.1 kb |

Our raw/theirs/mdax total lengths land within 1-10% of the published final
genome size on every species (these are genuinely large genomes, not
inflated ones) — total assembly length was never really the gap. The
persistent gap is **contiguity**: our best arms still run ~20-45% more
contigs and 12-34% lower N50 than published. That gap is most plausibly
`purge_dups` (haplotig/duplicate removal), a step outside any foldback
corrector's scope by design — see write-up in the session this benchmark
came from for the reasoning.

## The microlaimidae anomaly, and mode comparison

`mdax --mode permissive` produced a dramatically *worse* assembly than doing
nothing at all on microlaimidae specifically: total length nearly doubled
(1,084 Mb vs. 631 Mb raw / 576 Mb published), duplicated BUSCOs jumped 10x
(0.4% → 5.5%), complete BUSCO dropped below raw (11.9% vs. 16.1%). This did
not happen on the other three species.

Tested `--mode balanced` and `--mode strict` on the same raw reads:

| Mode | Complete | Single | Dup | Contigs | N50 | Total length |
|---|---|---|---|---|---|---|
| permissive ⚠️ | 11.9% | 6.4% | 5.5% | 40,978 | 34 kb | 1,084.4 Mb |
| balanced | 17.4% | 17.0% | 0.4% | 40,143 | 21 kb | 582.3 Mb |
| **strict** | **18.9%** | **18.4%** | 0.5% | **39,946** | 22 kb | **595.0 Mb** |

Both `balanced` and `strict` fully resolve it — total length lands right at
the published 575.8 Mb, duplication returns to normal. `strict` is the best
of the three mdax modes on this species, just shy of `theirs` (20.1%
complete) but clearly ahead of raw and permissive.

### Root cause (investigation only, no code changes)

`--mode` changes three things at once, not just a confidence threshold
(`src/cfg.rs`, `CallMode::tuning()`): the detection threshold
(`min_matches`/`min_identity`), the recursion cap on nested cutting
(`max_depth`: strict=2, permissive=8), and the **signature coarseness**
(`sig_flank_bp`/`sig_take`/`sig_value_shift`: permissive uses a much
smaller flank window, fewer minimizers, and heavier value-quantization). A
coarser signature is more collision-prone — reads from genuinely different
loci that happen to be similar (e.g. members of a repeat family) can get
folded into the same signature bucket, inflating the cross-read support
count that `is_real_foldback_stats()` (`src/fingerprint.rs:369`) uses as its
sole "real" gate.

### Attempted fix: signature-collision guard via split-position spread (ruled out, 2026-10-29)

`ModeTuning` already had a `split_tol_bp` field (strict=30, balanced=50,
permissive=100 bp) defined but never wired to anything (`#[allow(dead_code)]`)
— its doc comment suggested gating "real" on how tightly a signature's
observed split positions (`SupportStats::split_span()`) cluster across reads,
on the theory that a true recurring locus's reads localize to a consistent
split position while reads from unrelated loci colliding under a coarsened
signature would not.

Wired this into `is_real_foldback_stats()` and tested directly against the
microlaimidae raw reads with `--qc-only` (fast, no reassembly needed). Result:
it **does not discriminate** — it rejects "real" uniformly across every
species, not just the collision-affected one:

| Species | `real` clusters, guard off | `real` clusters, guard on (100bp) |
|---|---|---|
| Epsilonema (fine) | 3,215 | 105 (-97%) |
| Linhomoeus (fine) | 1,335 | 23 (-98%) |
| Microlaimidae (blowup) | 2,266 | 86 (-96%) |

The reduction is ~96-98% on *every* species regardless of whether it has the
collision problem. Root cause of the failure: `split` is a read-space
position (where the breakpoint falls within that particular read), which
varies with each read's own fragment boundaries regardless of whether the
locus is a true recurring one or a collision — it carries no signal about
inter-read relatedness. **Reverted** (`is_real_foldback_stats()`, `cfg.rs`,
`cli.rs` back to pre-change state) rather than ship a change that would
silently over-cut real palindromes in permissive mode. Left here so a future
session doesn't re-attempt the same dead end.

### Shipped instead: advisory-only mode suggestion (2026-10-29)

Since the classifier itself couldn't be safely fixed with the signal
available, `mdax` now prints a `QC MODE ADVICE` line (in both normal and
`--qc-only` runs) when `--mode permissive` and `max_cluster_n > 200` (or
>0.1% of support-map observations sit in clusters >50) — the threshold is
set from this directory's one confirmed case (microlaimidae: 407) vs. three
unaffected species (70-130), so treat it as a heuristic flag worth a second
look, not a verdict. It recommends comparing against `--mode balanced`/
`--mode strict` before committing to permissive; it does not change cutting
behavior itself. Verified it fires only on microlaimidae across all 4 species
tested here.

Also fixed in the same pass: the `--qc-only` `SUPPORT GATE` diagnostic had
its own inline reimplementation of the real/artefact decision instead of
calling `is_real_foldback_stats()` directly — i.e. it was a second copy of
the gate logic that could silently drift from the one pass2 actually uses.
Now calls the shared function; output is unchanged (verified byte-identical
on microlaimidae) but can't drift again.

Directly measurable from the existing `report.tsv` per-read output (no new
runs needed) — the `support_n` column, i.e. how many reads collapsed onto
the same signature:

| Species | Total reads | Max signature cluster | Reads in clusters >50 | Outcome (permissive) |
|---|---|---|---|---|
| Epsilonema | 368,903 | 70 | 70 (0.02%) | fine |
| Linhomoeus | 150,055 | 73 | 73 (0.05%) | fine |
| Trileptium | 286,229 | 130 | 130 (0.05%) | fine |
| **Microlaimidae** | 376,892 | **407** | **970 (0.26%)** | **blowup** |

*Read-count caveat (2026-10-30): the total-reads column above is stale —
re-checked against the actual files on disk and all four are substantially
larger (microlaimidae alone is ~1.45M reads, not 376,892; likely true of
all four). Doesn't change the qualitative picture (microlaimidae's cluster
is still 3-6x the others', scaled or not), just flagging so a future
session doesn't trust these counts for anything quantitative.*

### Fixed: sig_flank_bp 200 -> 400 (2026-10-30)

The "candidate follow-up" flagged above — check max signature-cluster size
under permissive-mode parameters — was carried out. Swept `--sig-flank-bp`
200-800bp (and `--sig-value-shift` 0 vs. 8, which turned out to have no
measurable effect at either setting — it's a bare bit-shift on a
well-distributed hash, not a "fuzzy-match near-duplicates" operation)
against all 4 species, qc-only, on matched 200k-read subsamples (first
confirmed microlaimidae's default-vs-subsample scaling is near-linear, so
subsample numbers are representative):

| Species | fb=200 (old default) | fb=400 (new default) |
|---|---|---|
| Microlaimidae | 63 | 31 |
| Epsilonema | 45 | 24 |
| Linhomoeus | 8 | 6 |
| Trileptium | 9 | 5 |

400bp is strictly better than 200bp on every species tested — shipped as
the new permissive-mode default (`src/cfg.rs`). Important caveat: this is
**not** a "bigger is always safer" relationship — 800bp regresses
microlaimidae back to 456 (as bad as the unfixed default), almost
certainly by pulling a second, more distant repeat family into the
window. Don't push this knob further without re-sweeping all 4 species.

Also worth noting: the union-find support-merging fix shipped earlier the
same session (commit `d993618`, aimed at the *opposite* problem — genuine
inverted repeats losing real support to signature fragmentation) measurably
*worsened* this specific collision case before the flank_bp fix landed
(max_cluster_n 407 -> 454) — expected, since that fix exists to merge more
aggressively, which is exactly the wrong direction when the underlying
signature is already too collision-prone. The two fixes pull in opposite
directions on purpose, for opposite failure modes; this is a reminder that
either one in isolation can look worse on the other's test case.

**Closed, confirmed at the assembly level (2026-10-31)**: mdax (permissive,
fixed binary) -> Flye -> BUSCO on microlaimidae, same method as the
original benchmark (job chain 374170 -> 374177 -> 374182, output under
`.../mdax_single_nematode_benchmark/{corrected,asm,busco}/microlaimidae/mdax_permissive_fixed/`):

| Arm | Complete | Single | Dup | Contigs | N50 | Total length |
|---|---|---|---|---|---|---|
| raw | 16.1% | 15.7% | 0.4% | 45,291 | 20 kb | 630.7 Mb |
| theirs | 20.1% | 19.6% | 0.5% | 42,179 | 21 kb | 605.1 Mb |
| balanced | 17.4% | 17.0% | 0.4% | 40,143 | 21 kb | 582.3 Mb |
| strict | 18.9% | 18.4% | 0.5% | 39,946 | 22 kb | 595.0 Mb |
| permissive (old, buggy) ⚠️ | 11.9% | 6.4% | **5.5%** | 40,978 | 34 kb | **1,084.4 Mb** |
| **permissive (fixed)** | **16.4%** | 16.0% | **0.4%** | 40,922 | 19 kb | **566.3 Mb** |
| published (ground truth) | — | — | — | 35,431 | 24.1 kb | 575.8 Mb |

The fix fully resolves the original bug: total length collapses back from
1,084 Mb to 566 Mb -- now the closest of any arm to the published 575.8
Mb -- and duplicated BUSCOs drop from 5.5% to 0.4%, in line with every
other arm. Complete BUSCO also improves (11.9% -> 16.4%, now ahead of raw
and competitive with balanced/strict). The old 34kb N50 was an artifact
of the duplication, not real contiguity; 19kb now that it's gone is
consistent with the other arms. `sig_flank_bp` 200 -> 400 is the whole
fix; nothing else needed changing for this case.

Microlaimidae's largest signature cluster is 3-6x bigger than any other
species', and it has 7-14x more reads embedded in large clusters. Consistent
with an unusually large/repetitive family in this genome that permissive
mode's coarse signature + deep recursion (`max_depth=8`) conflates into one
massively "recurring" locus and fragments accordingly — while strict/balanced's
tighter signature and shallower recursion don't let the same collision
happen. Not yet fixed in code; a plausible cheap QC check (max
signature-cluster size under permissive-mode parameters) is flagged as a
candidate follow-up, not implemented.

## The nrCaeEleg94 anomaly: balanced mode's sig_take=3 was unsafe on real repeat-rich genomes (2026-11-02)

An ultrareview of this session's union-find merge fix (see `src/fingerprint.rs`'s
`MAX_PLAUSIBLE_CLUSTER_N` doc comment) flagged that it has no bound on
transitive merging: a single band value shared between two otherwise-
unrelated reads unions *all* of each read's own bands together, and that
chains through the dataset via any other read bridging two components.
Checked the farm batch for the separate C. elegans rerun and 11-species
suite (both submitted 2026-10-30, before this was caught) and found it
had already happened, badly:

- Every C. elegans sample/mode combination except the two PiMmS samples
  showed `max_cluster_n` well above the 70-130 "fine" range documented
  above — several in the thousands (nrCaeEleg94 balanced: 9766).
- 9 of 10 new nematode-suite species (permissive mode) also exceeded it.
- Flye + BUSCO had already completed on top of these corrections for
  most of the nematode suite, and hifiasm had completed for several
  C. elegans samples across all 3 modes.

Investigated on nrCaeEleg94 (balanced mode, the worst offender) with the
exact detected/cut/"real" read counts at each step:

| Fix attempted | max cluster | real clusters | real reads | real % |
|---|---|---|---|---|
| none (as-submitted) | 9,766 | 285,737 | 2,623,972 | **75.0%** |
| hard cap (`MAX_PLAUSIBLE_CLUSTER_N=200`) | 9,766 (unenforced, just gated) | 285,659 | — | ~75% |
| stop-value filter (band values shared by >200 reads can't bridge) | 6,321 | 285,721 | — | ~75% |
| + revert `sig_take` 3→8 | 3,593 | 277,590 | — | ~74% (barely moved) |
| + also raise `min_support` 3→8 | 3,593 | 36,536 | 425,350 | **12.2%** |

Two things stood out. First, capping or filtering the *merge* step barely
helped — the dominant failure mode wasn't a few catastrophic cascades,
it was pervasive: 823,095 of ~906,000 clusters sat in the 0-9 size
bucket, and huge numbers of those still cleared `min_support=3` by
chance on a genome this repetitive. A merge-time fix can't touch that.
Second, `sig_take=3` turns out to have been the real regression:
`sig_take` 8→3 (committed 2026-10-30, reasoning in `src/cfg.rs` git
history) was validated only against `run_main` (a low-repeat-content
*synthetic* genome, where it was a clean win: hard-case FPR 4.4%→0.0%)
and a collision-risk spot check against microlaimidae (which also isn't
especially repeat-dense, so the check passed cleanly). Neither caught
that take=3's weaker 2-minimizer leave-one-out bands collide far too
easily on an actually repeat-rich genome like C. elegans.

**Fixed**: `sig_take` reverted 3→8 in `src/cfg.rs` (commit `a9c010a`).

**`min_support`: investigated, deliberately left at its default of 3
(2026-11-05) — raising it is not a fix.** Reverting `sig_take` alone only
dropped nrCaeEleg94's real-cluster count 3% (285,737→277,590); raising
`min_support` to 8 as well dropped the *aggregate* "real%" to a
reassuring-looking 12.2%. But that aggregate was the wrong number to
optimize: built real per-read ground truth for nrCaeEleg94 by aligning a
51,423-read sample to the actual C. elegans reference genome (GCF_000002985.6
/ WBcel235, 100,286,401 bp — downloaded fresh, a vastly better reference
than anything available for microlaimidae) and applying `synth_bench/
refbased_clip.py`'s overlapping-opposite-strand-alignment foldback test
unmodified. That gave 48,619 confirmed true artefacts and 1,948 confirmed
true genuine reads (856 unmapped reads excluded as uninformative) in the
sample. Joining this against mdax's own decisions at min_support 5/8/10:

| min_support | artefact recall (correctly cut) | genuine-read FPR (wrongly cut) |
|---|---|---|
| 5  | 76.2% | 71.7% |
| 8  | 88.0% | 81.5% |
| 10 | 92.1% | 84.9% |

Raising `min_support` makes genuine-read FPR *worse*, not better — it
only looked like an improvement in aggregate because artefacts
outnumber genuine reads ~25:1 in the detected-foldback population here,
so correctly reclassifying more artefacts as cut dominates the aggregate
number and hides what's actually happening to the much smaller genuine
population underneath it. This was a real methodological miss caught
only by building actual ground truth instead of trusting the aggregate;
see `synth_bench/celegans_refbased_eval.py` and
`synth_bench/extract_reads_by_id.py` for the method, reusable on any
other sample with a real reference genome available.

**Decision (2026-11-05): leave `min_support` at 3.** No value tested in
5-10 actually solves "don't cut real reads" on this genome — the
genuine-read damage (72-85% FPR) looks structural (support fragmentation
for genuinely low-recurrence real loci, the same class of problem
`synth_bench/README.md` originally motivated this whole investigation
with), not something this threshold, in this range, can fix. Raising it
trades a small amount of additional artefact recall for more genuine-read
damage, which is the wrong direction to default to. Documented here as a
known, unresolved limitation rather than solved: a real fix needs either
a smarter per-cluster signal (not just a count) or a different mechanism
entirely (e.g. the learned-classifier idea discussed earlier in this
project's history), not a `min_support` retune.

**The union-find merge fix itself (`src/fingerprint.rs::SupportMap`) was
not reverted** — it fixes a real, separately-validated problem (under-
counting support for genuine recurring loci) and the evidence above
shows the dominant cause here was signature specificity (`sig_take`),
not the merge logic. The stop-value filter it gained during this
investigation stays too: harmless where not needed, and it demonstrably
reduces (if not eliminates) the worst individual cascades.

**Lesson for future signature-specificity changes**: validate against a
genuinely repeat-rich real genome, not just a low-repeat synthetic one
or a single real species that happens not to be especially repetitive.
`run_main` and microlaimidae both said take=3 was safe; C. elegans said
otherwise, decisively.

**Farm impact**: the whole batch (both suites) needs rerunning against
the reverted binary once min_support is decided. Separately, most of the
already-submitted downstream Flye/BUSCO jobs are now permanently stuck
anyway — `bjobs -p` shows "Dependency condition invalid or never
satisfied" on essentially every one, because their LSF dependency
referenced an upstream job ID old enough to have aged out of the
scheduler's history. Not a congestion issue; they will never run as-is.

**Farm batch results (2026-11-09/10)**: after the full remediation, the
11-species nematode suite and most of the C. elegans rerun samples
confirm mdax performs comparably to or better than both raw and the
"theirs" alternative, no catastrophic failures. The nrCaeEleg93/94/95
samples confirm the documented min_support limitation above at the real
assembly level (50-55% duplicated BUSCOs). A separate, unrelated finding
-- a reproducible `hifiasm` (third-party assembler, not mdax) crash on
nrCaeEleg33's permissive-mode assembly -- is documented in full in
`mdax_celegans_rerun/README.md` on lustre (not duplicated here since it's
a tooling issue, not an mdax one); decided to document and not chase
further.

## Caveats

- **Single-pass only.** No Illumina polishing, no `purge_dups`, no
  decontamination — deliberately, to isolate the correction-method
  comparison. Absolute BUSCO/N50 numbers here are not comparable to the
  paper's own published headline figures (which are protein-mode BUSCO
  against the `eukaryota` lineage on a fully-finished genome, not
  genome-mode `nematoda_odb10` on a bare Flye assembly).
- 4 of 13 species tested, picked mostly for small ONT read size. Not a
  claim that this generalizes to the full species set.
- microlaimidae was run through the full mode sweep; the other three
  species were only tested under the default `--mode permissive`.

## Scripts

| Script | Description |
|---|---|
| `run_species_pipeline.sh` | Full driver: submits the 9-job LSF chain (2 corrections → 3 Flye → 3 BUSCO) for one species. |
| `01_palindrome_removal.sh` | Wraps the authors' `remove_palindrome.sh` (2 iterations). |
| `02_run_mdax.sh` | Runs `mdax --refine-mode ont --cut-low-ident`. |
| `03_flye.sh` | `Flye --nano-raw --meta`. |
| `04_busco.sh` | BUSCO `nematoda_odb10`, genome mode. |
| `remove_palindrome.sh`, `pafIdentifyPalimdrom.pl`, `fastq_partition.and.chop.palindrome.pl` | Lee et al.'s original scripts, unmodified, copied from the Zenodo archive for provenance. |

## Data

Raw ONT reads are symlinked in `data/` (gitignored — points to
`/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_single_nematode_benchmark/reads/`).
Full working directory (corrected reads, Flye assemblies, LSF logs) lives
alongside at the same lustre path, outside git.

## Results

Full BUSCO `short_summary.txt` for every run is under `results/<species>_<arm>/`.

## 2026-11-18: the gate was counting loci, not junctions — replaced with fold-axis support

**How it was found.** Every earlier ground-truth check only scored reads mdax
had already detected. Taking a *random* 3% of all nrCaeEleg92 reads (41,675),
aligning them to WBcel235, and calling a read a true foldback when it splits
into opposite-strand segments within 5 kb on one chromosome
(`synth_bench/refbased_recall_eval.py`) showed half the library is foldback
artefact, and that mdax was cutting far fewer of them than the logs implied:

| build / settings (balanced, nrCaeEleg92) | true foldbacks cut | kept as "real" | not detected | linear reads cut |
|---|---|---|---|---|
| 30 Oct (`sig_take=3`; what every `mdax_celegans_rerun` assembly used) | 37.5% | 34.7% | 27.8% | 0.8% |
| 2 Nov (`sig_take=8`) | 54.3% | 15.1% | 27.8% | 1.0% |
| 2 Nov, gate off (`--min-support 1000000`) | 69.4% | 0 | 27.8% | 1.4% |
| fold-axis support, old detection thresholds | 71.7% | 0.5% | 27.8% | 1.2% |
| fold-axis support, new balanced defaults | 78.1% | 0.7% | 21.2% | 2.1% |

On nrCaeEleg94 (~720x) the 30 Oct build cut 15.6% and kept 51.5% as "real";
fold-axis support cuts 73.1% and keeps 1.0%.

**Why the old gate failed.** The signature was the bottom-k minimizers within
±600 bp of the split: it names a locus. nrCaeEleg92 has ~534k detected
junctions over 100 Mb, one per ~190 bp, so unrelated artefacts share a window
by chance; nrCaeEleg94 has one per ~29 bp. The wrongly-"real" reads were not
PCR duplicates (no cluster had fewer than two distinct (length, split)
pairs). The old adaptive calibration made it worse by *lowering*
`min_support` to the dataset's "typical recurrence" of 2, which is the
chance-collision level; it fired on 5 of 12 C. elegans samples. This also
supersedes the "clonal amplification" explanation for nrCaeEleg94 above as
the main cause, though the measured support tail there is heavier than
Poisson (decaying ~0.7x per step), consistent with clonal copies or coverage
hotspots contributing.

**The replacement** (`fingerprint::JunctionSupport`). For each foldback, key
the matched minimizers nearest the axis by `(minimizer value, distance to its
own inverted copy)`. That distance is twice the minimizer's distance to the
symmetry axis, is read-coordinate-free, and needs no refinement, so two reads
share a key only if they fold about the same axis. A read's support is the
median, over its keys, of the number of reads carrying each. The threshold is
calibrated upward from the data: a Poisson bound from the fraction of
foldbacks sharing an axis with nobody, then extended along the measured
histogram tail (nrCaeEleg92: 3 → 7-8; nrCaeEleg94: 3 → 15). `--sig-flank-bp`,
`--sig-take`, `--sig-value-shift`, `MAX_PLAUSIBLE_CLUSTER_N` and the
stop-value filter no longer take part in the gate.

**Does read-level recall turn into assembly gains? Yes** (nrCaeEleg92,
hifiasm, `nematoda_odb10`):

| arm | contigs | N50 (Mb) | largest (Mb) | BUSCO C | D |
|---|---|---|---|---|---|
| raw | 1248 | 0.14 | 0.89 | 88.4% | 7.2% |
| 30 Oct balanced | 894 | 0.25 | 1.32 | 93.9% | 5.4% |
| 2 Nov balanced | 792 | 0.32 | 1.44 | 95.6% | 3.7% |
| 2 Nov, gate off | 691 | 0.36 | 1.97 | 96.3% | 2.7% |
| 2 Nov, gate off, loose detection (8/600/300) | 823 | 0.31 | 1.68 | 94.4% | 3.1% |
| fold-axis, old detection (20/2000/1000) | 705 | 0.39 | 1.96 | 96.6% | 3.4% |
| fold-axis, new balanced defaults (12/1000/500) | 682 | 0.39 | 1.82 | 96.5% | 2.9% |
| fold-axis, loose detection (8/600/300) | 702 | 0.39 | 1.68 | 96.0% | 2.9% |

Two things to take from it. Cutting more true artefacts helps, but false cuts
hurt: loose detection with no gate is worse than tight detection with no
gate, and with the gate in place loose detection buys nothing over the new
defaults. And the three fold-axis arms sit within run-to-run noise of each
other on N50 (an earlier build of the same settings gave 0.41), so the
detection-threshold change is a small gain, the gate change the large one.

**ONT** (celegans_wga, same method with `map-ont`, permissive mode): recall is
unchanged at 65.4% because the old gate was already cutting ~99% of what it
detected there; linear reads cut fall 2.4% → 1.7% and linear reads protected
as "real" rise 0.1% → 0.7%. The ONT limit is detection: 20% of foldbacks with
arms under 500 bp and 37% at 500-1000 bp are found.

**Still open.** ~21% of true HiFi foldbacks (34% on ONT) go undetected, nearly
all short-armed or within 500 bp of a read end. Non-foldback chimeras (1.4% of
reads) are out of scope by design. The arm-identity estimate is noisy (only
64% of true HiFi artefacts score >= 0.995). The numbers above come from two
C. elegans libraries and one ONT library.
