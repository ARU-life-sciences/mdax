# `mdax`

Read in HiFi or ONT raw data generated from Multiple Displacement Amplification methods, and correct reads. This tool currently only corrects 'foldback' or 'invert' chimeras (`mdax`). We also provide a binary to search for long Inverted Repeats (IR's) in assemblies (`irx`).

**Contents:** [Install](#install) · [Quick start & recommended settings](#quick-start--recommended-settings) · [Full CLI reference](#usage) · [Validation](#validation) · [Tests and benchmarks](#tests-and-benchmarks) · [Performance](#performance) · [Algorithm overview](#algorithm-overview) · [`irx`](#irx-detecting-natural-irs-from-assemblies)

## Install

A Linux x86-64 binary is attached to each [release](https://github.com/ARU-life-sciences/mdax/releases) (needs AVX2/AES-NI, glibc ≥ 2.35). Or build from source:

```bash
# install Rust first (https://rust-lang.org/tools/install/)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/ARU-life-sciences/mdax
cd mdax
cargo build --release      # SIMD flags are set in .cargo/config.toml
# or, to put mdax and irx in PATH
cargo install --path=.
```

## Quick start & recommended settings

```bash
mdax --output <output.fa> <in.fa> > <out.table>
```

Default outputs are the fasta file of corrected reads, and a table of all detected foldback chimeras in the readset. The full flag reference is below in [Usage](#usage) — this section covers the settings you actually need for the common cases.

### Recommended by data type

| Data type | Recommended command |
|---|---|
| PacBio HiFi / CCS | `mdax --mode balanced --refine-mode hifi -o out.fa -r out.tsv in.fa` |
| ONT | `mdax --mode permissive --refine-mode ont -o out.fa -r out.tsv in.fa` |

A detected foldback is cut unless enough other reads fold about the *same axis* (see [Algorithm overview](#algorithm-overview)); that applies to foldbacks with a low identity estimate too, which are cut by default (`--keep-low-ident` passes them through instead; the old `--cut-low-ident` flag is still accepted and does nothing, as are the retired locus-signature flags `--sig-flank-bp`, `--sig-take` and `--sig-value-shift`; mdax prints a note when one is given). For ONT, `--refine-mode ont` automatically switches to smaller k-mers (k=11, w=11) for better sensitivity on error-prone reads, unless you set `-k`/`-w` explicitly — the effective values are echoed in the parameter summary printed at startup.

**On `--min-support`:** leave it alone unless you have a reason. The number of unrelated artefacts that share a fold axis purely by chance grows with library depth, so after pass 1 mdax estimates that chance rate from the data and raises the threshold to match (never below the mode default) — check stderr for a `SUPPORT CALIBRATION` line to see what it picked (e.g. 3 → 7 on a ~100x library, 3 → 14 on a ~720x one). An explicit `--min-support` disables this and applies one value whatever the depth.

### Looking before you commit to assembly (`--qc-only`)

`--qc-only` runs pass1 (detection + fold-axis support counting) and then exits — no cutting, no FASTA/TSV written, `-o` not required. Pass1 is far cheaper than pass2 plus a downstream assembly.

The `QC SUMMARY` line reports `lambda` (the mean number of other reads sharing a foldback's axis by chance — a direct read-out of how saturated the library is with artefact junctions), the calibrated `min_support`, and how many foldback reads sit at or above it and would therefore be kept as genome-templated.

```bash
mdax --mode permissive --refine-mode ont --qc-only reads.fastq.gz
```

### Examples

Basic HiFi correction, writing the report to a file instead of stdout:
```bash
mdax --mode balanced --refine-mode hifi \
    -o corrected.fa -r report.tsv reads.fasta.gz
```

ONT correction, using more threads:
```bash
mdax --mode permissive --refine-mode ont -t 16 \
    -o corrected.fa -r report.tsv reads.fastq.gz
```

Conservative run for a dataset you don't fully trust yet — fewer, higher-confidence cuts:
```bash
mdax --mode strict --refine-mode hifi -o corrected.fa -r report.tsv reads.fasta.gz
```

Report only, no correction — inspect calls before deciding whether to trust them (nothing is written to `/dev/null` that you need):
```bash
mdax --mode balanced -o /dev/null -r report.tsv reads.fasta.gz
```

## Usage

Current usage is:

```txt
Remove inverted duplication chimeras from fasta files.

Usage: mdax [OPTIONS] <input>

Arguments:
  <input>
          Input FASTA/FASTQ(.gz) containing reads to scan for foldback (inverted-duplication) artefacts.
          
          Foldbacks are detected *within each read* by looking for long, reverse-complement self-matches.
          This file is streamed: runtime is roughly linear in total bases, and memory is bounded by channel buffering.
          
          Tip: if your reader supports it, `.gz` inputs are convenient but decompression can become a bottleneck at high thread counts.

Options:
  -r, --report <report>
          Write a TSV report describing detected events and decisions.
          
          The TSV is intended for *debugging, QC, and post filtering*:
          - read_id: Read identifier from the FASTA/FASTQ header.
          - len: Read length in bases.
          - event: Event type detected (always foldback, may include e.g. concatemers).
          - called: 1 if an event was detected in pass1, else 0.
          - coarse_split: Approximate breakpoint index from the first-pass detection (before refinement).
          - refined_split: Breakpoint index after local refinement around coarse_split.
          - delta: refined_split - coarse_split (signed). CURRENTLY EMPTY.
          - matches: Number of minimizer matchpoints supporting the foldback call (pass1 evidence strength).
          - span_p1: Span (bp) covered by pass1 supporting matchpoints / arms. I.e. how much of the read supports the foldback geometry
          - p2_span: Span used/confirmed in pass2 (if you do a second-pass validation). CURRENTLY EMPTY.
          - cross_frac: Fraction of matchpoints crossing the split / inconsistent with a clean junction (often a “messiness” metric). CURRENTLY EMPTY.
          - coarse_score: Score assigned by the coarse foldback detector prior to refinement.
          - refined_score: Score assigned after breakpoint refinement.
          - identity_est: Estimated nucleotide identity between the two foldback arms.
          - gap_est: Estimated template-switch gap size (bp) between the two arms; 0 for clean hairpin folds.
          - support_n: Number of reads (including this one) folding about the same axis; compared against min_support.
          - support_rank_frac: DEPRECATED, always 0. Kept so column positions do not change.
          - support_span: Number of junction-key elements this foldback had (0 if it had too few to be keyed).
          - decision: Final classification (real, artefact, low_ident).
          - sig_band: First junction key of this foldback, as hex (NA if it was not keyed). For grouping rows only.
          
          Use '-' to write to stdout
          
          
          [default: -]

  -o, --output <output>
          Output FASTA containing the kept/corrected reads.
          
          For each input read, mdax writes one output record containing the sequence you decided to keep:
          - If no foldback is called, the read is written unchanged.
          - If a foldback is called and judged an artefact, mdax cuts/chops the read (possibly recursively) and writes the kept portion.
          The kept portion is whichever side of the split is longer (left arm, or
          everything from the right arm onward), not always the left arm -- a
          foldback near the start of a long read would otherwise discard almost the
          whole read even though the tail past it is ordinary, assemblable sequence.
          Nested/compound artefacts are still fully corrected regardless of which
          side was kept first.
          - If a foldback is judged 'real' (supported across reads), mdax keeps it intact by default.
          

  -k <k>
          Minimizer k-mer size used in the *coarse detector*.
          
          How it affects foldback calling:
          - Larger k makes minimizers more specific, reducing random/self-repeat matches, leading to fewer false positives.
          - Larger k is less tolerant of sequencing errors (especially ONT), so true foldbacks may have fewer matching minimizers and therefore fewer calls.
          
          Practical guidance:
          - ONT: k=15–17 is typical; higher k can reduce sensitivity on noisy reads.
          - HiFi: k=19–25 can work well because errors are rare.
          
          Notes:
          - k does not directly change refinement quality; it changes which candidates reach refinement.
          
          
          [default: 17]

  -w <w>
          Minimizer window size used in the *coarse detector*.
          
          How it affects foldback calling:
          - Larger w samples fewer minimizers (one per window), so faster, but fewer matchpoints.
          This can reduce sensitivity, especially for short arms or noisy reads.
          - Smaller w samples more minimizers, so slower, but increases matchpoint density.
          This generally increases sensitivity and stabilizes split estimation at the coarse stage.
          
          Rules of thumb:
          - If you're missing obvious foldbacks: decrease w or decrease min_matches.
          - If you're calling too many weak events: increase w or increase min_matches/min_span.
          
          Note: 'must be odd', see `simd_minimizers` Rust crate for more details.
          
          [default: 21]

      --min-span <min_span>
          Override: minimum evidence span (bp) required for a coarse foldback call.
          
          What 'span' means here:
          - It is the span of matched minimizer positions along (typically) one arm of the junction.
          - Conceptually: how much of the read participates in the reverse-complement self-match.
          
          How it affects foldback calling:
          - Increasing min-span demands longer arm evidence. Fewer calls, much fewer random palindromic hits.
          - Decreasing min-span allows shorter events. More calls, but more susceptibility to local repeats.
          
          Interaction with other knobs:
          - If you decrease w (more minimizers), you may be able to increase min-span without losing sensitivity.
          - If you keep min-span low, consider increasing min_identity to avoid cutting based on weak matches.

      --min-identity <min_identity>
          Override: minimum refinement identity (identity_est) to treat a detected foldback as credible in pass2.
          
          Where it applies:
          - After coarse detection and refinement, `mdax` estimates similarity between the two arms.
          - If identity_est < min-identity, the event is treated as low-confidence (e.g. labelled 'low_ident')
          and typically will not be used for support-aware 'real vs artefact' decisions.
          
          How it affects calling/cutting:
          - Increasing min-identity means fewer events considered; reduces over-cutting on noisy/weak alignments.
          - Decreasing min-identity means more events considered; increases sensitivity but may cut reads on spurious matches.
          
          Practical guidance:
          - ONT data often benefits from a lower identity threshold than HiFi.
          - If you're seeing many false artefact cuts, raise this first before touching support thresholds.

      --min-support <min_support>
          Override: minimum number of reads (including the read itself) that must fold about the same axis for a foldback to be treated as 'real'.
          
          Interpretation:
          - Pass1 records, for every detected foldback, where its symmetry axis sits
          relative to the minimizers next to it.
          - A genome-encoded inverted repeat folds about the same axis in every read
          that spans it; independent amplification artefacts at the same locus fold
          about different ones. In pass2 a foldback whose axis is shared by
          >= min-support reads is treated as genome-templated and is NOT cut.
          
          How it affects cutting:
          - Increasing min-support means harder to be considered real. More foldbacks treated as artefacts therefore more cutting.
          - Decreasing min-support means easier to be considered real. Fewer cuts, but risk of keeping recurrent artefacts.
          
          Depth dependence and auto-calibration:
          - The deeper the library, the more unrelated artefacts share an axis purely
          by chance, so a fixed threshold that is right at 30x is far too lax at
          1000x.
          - Unless you pass this flag explicitly, `mdax` calibrates after pass1: it
          estimates the chance co-location rate from the data and raises the
          threshold until chance alone rarely reaches it, then follows the measured
          tail of the support histogram further out if the library is over-dispersed
          (heavily amplified libraries are). It never goes below the mode value.
          Check stderr for a `SUPPORT CALIBRATION` line showing what it picked.
          
          Guidance:
          - Leave this alone unless you have a reason: an explicit value disables the
          calibration and applies uniformly whatever the library's depth.

      --max-depth <max_depth>
          Override: maximum recursion depth when cutting foldback artefacts.
          
          What recursion does:
          - After cutting an artefact foldback, the kept fragment may still contain another foldback junction.
          - `mdax` can re-run detection on the kept fragment and cut again, up to max-depth times.
          
          How it affects output:
          - Increasing max-depth leads to more aggressive cleanup of nested/compound artefacts, but may over-trim in noisy data.
          - Decreasing max-depth leads to faster and safer, but may leave residual artefact structure.
          
          Guidance:
          - Start modest (e.g. 2–5). Raise only if TSV shows many 'artefact' decisions remaining after one cut.
          - Very high depth can amplify mistakes: one false cut early can cascade.
          - Higher recursion, more compute!

      --end-guard <end_guard>
          Guard region at both ends of reads (bp).
          
          Foldback calls near read ends are often unreliable because:
          - there is little sequence context on one side of the split
          - minimizer matchpoint chains become truncated
          - refinement can be biased by missing arm sequence
          
          How it affects calling:
          - Increasing end-guard leads to fewer calls near ends (more conservative), may miss real junctions that genuinely occur close to ends.
          - Decreasing end-guard leads to more calls, but higher risk of unstable splits and over-cutting.
          
          Tip:
          - If you have many short reads, a large end-guard can suppress almost all calls.
          Ensure end-guard is comfortably smaller than typical read length.

      --refine-window <refine_window>
          Refinement search half-window (bp) around the coarse split position.
          
          Refinement works by testing candidate split positions around the coarse estimate.
          This parameter controls the maximum shift allowed during refinement.
          
          How it affects calls:
          - Increasing refine-window can recover splits when coarse detection is jittery, but costs more compute.
          - Decreasing refine-window is faster, but may 'miss' the true breakpoint if coarse split is off.
          
          If you observe large coarse to refined deltas in the TSV, increase this.
          If deltas are always small, you can reduce it to speed up.
          
          [default: 100]

      --refine-arm <refine_arm>
          Arm length (bp) used during refinement to compare the two sides of the junction.
          
          Refinement typically extracts `arm` bases on each side of the candidate split and computes
          a similarity score (via Hamming for HiFi or banded Levenshtein for ONT).
          
          How it affects identity estimates and stability:
          - Increasing refine-arm gives more evidence, smoother identity estimates, fewer spurious calls; slower.
          - Decreasing refine-arm is faster and more local, but identity becomes noisier and can be biased by repeats.
          
          Guidance:
          - If you see unstable identity_est or lots of borderline calls, increase this.
          - If reads are short, arm too large can hit end-guard logic or reduce usable candidates.
          
          [default: 1200]

      --max-jump-clip <max_jump_clip>
          After finding the best split position, the refiner scans offsets δ = 0..max-jump-clip
          on the right arm to find the δ that maximises identity over `arm − δ` bases.
          
          In MDA foldback artefacts the enzyme may slip 0–N bp before beginning the reverse copy,
          leaving a short uncopied region immediately right of the junction ('template jump').
          Without this correction those bases inflate the apparent mismatch count and drive
          identity_est down even when the arms are near-identical.
          
          A minimum effective arm of 100 bp is always enforced regardless of this setting,
          so the actual ceiling is min(max-jump-clip, refine-arm - 100).
          
          Guidance:
          - Increase if you observe reads with very low identity_est but clear anti-diagonal
          foldback structure in dotplots, especially if gap_est is hitting the current ceiling.
          - Requires --refine-arm >= max-jump-clip + 100 to take full effect.
          
          [default: 1000]

      --max-ed-rate <max_ed_rate>
          Maximum tolerated edit-distance rate during ONT refinement.
          
          In ONT mode, refinement uses a banded edit distance. The band width is derived from the arm length
          and this max-ed-rate (e.g. band ≈ arm * max_ed_rate).
          
          How it affects refinement:
          - Increasing max-ed-rate gives wider band. More tolerant of indels/noise, but slower.
          - Decreasing max-ed-rate gives narrower band. Faster, but may underestimate similarity on noisy reads.
          
          This mainly impacts ONT mode; HiFi mode often uses a cheaper/Hamming-like comparison.
          If ONT refinement is missing obvious foldbacks, raise this slightly.
          If ONT refinement is slow, lower this (but watch for identity_est becoming systematically low).
          
          [default: 0.25]

      --min-matches <min_matches>
          Override: minimum number of minimizer matchpoints required for a coarse foldback call.
          
          Matchpoints are the basic evidence units for coarse detection.
          They come from matching minimizers between a read segment and its reverse-complement self-match.
          
          How it affects calling:
          - Increasing min-matches gives fewer calls, higher specificity, better split stability.
          - Decreasing min-matches gives more calls (including weak ones), increased risk of calling local repeats.
          
          Interaction with k/w:
          - Larger k or larger w reduces matchpoint density, so you may need to lower min-matches.
          - Smaller k or smaller w increases matchpoint density, so you can raise min-matches to stay specific.

      --fold-diag-tol <fold_diag_tol>
          Foldback anti-diagonal clustering tolerance.
          
          The foldback detector buckets minimizer matchpoints by an anti-diagonal coordinate, using
          d = pos1 + pos2. True foldbacks produce a dense cluster along an anti-diagonal;
          noise and repeats produce scattered points.
          
          How it affects calling:
          - Decreasing fold-diag-tol leads to tighter clustering and fewer calls, but can fragment evidence if reads are noisy.
          - Increasing fold-diag-tol leads to looser clustering and more calls, but higher chance of merging unrelated matches.
          
          If you see many borderline calls with scattered matchpoints, decrease this.
          If you miss foldbacks on noisy reads where matchpoints are smeared, increase this slightly.
          
          [default: 120]

      --forward-only <forward_only>
          Whether minimizer hashing is strand-specific.
          
          - forward-only (true): minimizers depend on the read’s forward orientation.
          This tends to be good for coarse detection and avoids some symmetrical collisions.
          - canonical (false): a k-mer and its reverse complement map to the same minimizer value.
          This can increase sensitivity in some cases, but can also increase spurious self-matches in repetitive regions.
          
          How it affects foldback calling:
          - canonical may increase the number of matchpoints (more hits), which can increase calls.
          - forward-only can be more conservative and may reduce false positives.
          
          Note: this flag applies to the detector’s minimizer sampling.
          
          [default: true]
          [possible values: true, false]

  -m, --refine-mode <refine_mode>
          Select refinement mode (error model) used when refining breakpoint coordinates.
          
          - hifi: fast comparison (Hamming-like) assuming low error rates.
          Produces stable split positions and identity estimates when reads are accurate.
          - ont: banded Levenshtein-style refinement tolerant to indels/substitutions.
          Slower but necessary for noisier reads.
          
          How it affects calling/cutting:
          - ONT mode can rescue true foldbacks that would look low-identity under HiFi assumptions.
          - HiFi mode can prevent over-fitting noise (calling weak alignments as foldbacks).
          
          Important:
          - This controls refinement; coarse detection still depends heavily on minimizer parameters and gates.
          - If you're using ONT reads, start with ont; if you're using HiFi reads, start with hifi.

          Possible values:
          - hifi: Refine breakpoints for HiFi reads (Hamming; fast)
          - ont:  Refine breakpoints for ONT reads (Levenshtein, slower)
          
          [default: hifi]

      --mode <mode>
          High-level preset controlling detection and correction strictness.
          
          This sets a coherent bundle of thresholds (min_matches, min_span, min_identity, min_support,
          max_depth, end_guard).
          
          - strict: fewer calls and fewer cuts. Requires strong evidence and high identity.
          Best when you want to avoid over-trimming.
          - balanced: reasonable default trade-off.
          - permissive: more calls and more cutting. Useful for exploration or very messy libraries.
          
          Override precedence:
          - If you pass explicit overrides (e.g. --min-identity), they override the preset values (will say in help).
          - Use the stderr output or TSV report to confirm which effective parameters were applied.

          Possible values:
          - strict:     Fewer calls/cuts; higher thresholds; keeps more reads
          - balanced:   Reasonable defaults
          - permissive: More calls/cuts; lower thresholds; exploratory
          
          [default: balanced]

      --no-rescue
          By default, when a read's foldback evidence is too thin for the coarse gates (--min-matches, --min-span, --end-guard), mdax still considers any anti-diagonal with a few matched minimizer pairs, and accepts it only if the two arms then align base-for-base (--rescue-min-arm bp at --rescue-min-identity or better). This is what finds foldbacks with short arms or a junction close to a read end. With this flag only the coarse gates are used.

      --rescue-anywhere
          By default the rescue path only accepts a candidate whose shorter arm runs to a read end. With this flag, verified candidates anywhere in the read are accepted. On reference-scored WGA HiFi data the mid-read candidates were as often genuine short inverted repeats as artefacts, so this is off by default.

      --rescue-terminal <rescue_terminal>
          Whether the rescue path also matches the first and last ~1.2 kb of each read against its own reverse complement using every k-mer (--rescue-term-k) rather than sampled minimizers, accepting a candidate only if the inverted copy runs all the way to the read end and aligns. Off by default: on reference-scored HiFi and ONT WGA data the minimizer-seeded rescue already proposed nearly every candidate this finds, and the few extra calls were false nearly as often as true.
          
          [possible values: on, off]

      --rescue-term-k <rescue_term_k>
          Rescue path: seed length for the read-end scan [default: 13 hifi, 9 ont]

      --rescue-min-arm <rescue_min_arm>
          Rescue path: minimum verified arm length (bp) [default: 150 hifi, 200 ont]

      --rescue-min-identity <rescue_min_identity>
          Rescue path: minimum identity of the verified arm alignment [default: 0.90 hifi, 0.65 ont]

      --rescue-min-pairs <rescue_min_pairs>
          Rescue path: minimum matched minimizer pairs proposing a candidate [default: 3]

      --keep-low-ident
          By default a detected foldback whose estimated arm identity falls below --min-identity is still cut, unless enough other reads fold about the same axis (the same support gate every other foldback goes through). With this flag such reads are passed through unchanged and labelled low_ident.
          
          Why cutting is the default: measured against a reference genome on WGA HiFi data, ~97% of detected foldbacks with a low identity estimate were true amplification artefacts whose identity had simply been under-estimated, and the estimate is systematically low on indel-rich (ONT / CLR) reads.

      --qc-only
          Stop after pass1 (detection + fold-axis support counting) and report the SUPPORT HISTOGRAM / SUPPORT CALIBRATION / QC SUMMARY diagnostics that a normal run already logs, then exit -- no pass2, no cutting, no FASTA/TSV written. --output/-o is not required in this mode.
          
          Intended as a cheap look at a new dataset before committing to correction + assembly: the QC SUMMARY line reports the chance co-location rate (lambda), the calibrated min_support, and how many foldback reads would be kept as genome-templated.
          
          Note --mode changes pass1's own detection parameters, so comparing modes means rerunning this once per --mode, not one shared pass.

  -t, --threads <threads>
          Number of worker threads to parallelise foldback detection and correction.
          
          `mdax` uses a multi-threaded pipeline with a single reader and writer thread and multiple worker threads for detection, refinement, and classification.
          
          Notes:
          - Increasing threads improves performance up to the point where the workload becomes memory- or I/O-bound.
          - For compute-only runs (e.g. -o /dev/null -r /dev/null), scaling typically plateaus well below the total number of CPU cores.
          - For runs that write FASTA/TSV output, disk I/O can dominate wall-clock time, and increasing threads may provide little or no speedup.
          
          Default is a single core, multithreading will likely give gains up to 24-32 threads.
          
          [default: 1]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Validation

Full evaluation (read-level truth on four MDA read sets, hifiasm/Flye + BUSCO on 27 libraries): [`reports/mdax_eval/mdax_eval.pdf`](reports/mdax_eval/mdax_eval.pdf), also attached to the v0.2.0 release. In short: worth using on MDA HiFi (40–62% fewer contigs, 2–4× N50 vs raw), smaller gains on MDA ONT/CLR, no use on unamplified libraries.

## Tests and benchmarks

| Directory | Tool | Type | Data | README |
|---|---|---|---|---|
| `tests/` | both | Rust integration tests, incl. reference-truth regression | 57 real HiFi reads, 2 × 500 public ONT reads, synthetic | [tests/README.md](tests/README.md) |
| `synth_bench/` | `mdax` | Synthetic harness; farm-side truth scorers | — | [synth_bench/README.md](synth_bench/README.md) |
| `test_pal/` | `mdax` | ONT comparison vs Strobl 2023 Perl baseline | 20k ONT reads (SRR24201687) | [test_pal/README.md](test_pal/README.md) |
| `test_celegans/` | `mdax` | Real HiFi benchmark, two MDA samples | 200k reads × 2 *C. elegans* samples | [test_celegans/README.md](test_celegans/README.md) |
| `irx_tests/` | `irx` | Synthetic ground-truth evaluation | ~120 planted IRs, seed 42 | [irx_tests/README.md](irx_tests/README.md) |
| `data/` | both | Test data and generators | — | [data/README.md](data/README.md) |

Test data for `test_pal/` and `test_celegans/` is not tracked in git (large files). The `test_pal/` and `test_celegans/` numbers below predate v0.2.0; current numbers are in the [evaluation report](#validation).

### Rust integration tests (`tests/`)

```bash
cargo test
```

Tests the full two-pass pipeline on synthetic and 57 real HiFi reads with BLAST-derived ground truth, and (`truth_regression.rs`) checks foldbacks cut / linear reads cut on 1,000 public MDA ONT reads labelled against a reference (`data/ci_truth/`). Deterministic, a few seconds; run in CI on every push. `synth_bench/truth_regress.sh` runs the full-library version on the farm.

### ONT baseline comparison (`test_pal/`)

20,000 ONT reads from SRR24201687; ground truth from Strobl et al. 2023.

Metrics are generous (artefact + low_ident both counted as true positive).

| Configuration | Precision | Recall | F1 |
|---|---|---|---|
| `balanced --refine-mode hifi` (default) | 0.998 | 0.647 | 0.785 |
| `permissive --refine-mode ont --cut-low-ident`¹ | **0.995** | **0.803** | **0.888** |

¹ `--cut-low-ident` is now the default and the flag is a no-op.

```bash
cd test_pal && ./test_pal.sh       # single run
cd test_pal && ./sweep.sh          # parameter sweep (cached)
```

### *C. elegans* MDA benchmark (`test_celegans/`)

Two real PacBio CCS single-worm MDA samples with contrasting BUSCO scores. No external ground truth; primary metric is artefact detection rate.

| | nrCaeEleg92 (82.3% BUSCO) — pre-fix, not refreshed | nrCaeEleg95 (64.5% BUSCO) — post-fix |
|---|---|---|
| Foldback detected | 85,508 (42.8%) | 77,249 (38.6%) |
| → artefact cut | 85,137 (42.6%) | 70,445 (35.2%) |
| → real palindrome | 371 (0.2%) | 6,804 (3.4%) |
| Bases removed | ~422 Mb | ~366 Mb |
| Fraction of artefact read retained | ~49.8% | ~50.0% |

Detection is unchanged by the current build (nrCaeEleg95: 77,249 both before and after); what changed is the real/artefact split — a 12.7x increase in junctions correctly classified "real" (535 → 6,804), because `min_support` now self-calibrates to this dataset's actual recurrence instead of using a fixed threshold. Confirmed at the assembly level too: reassembling a corrected nrCaeEleg95 subsample with hifiasm gives N50 216kb → 370kb (+71%) and BUSCO complete 92.2% → 96.6% versus no correction. See [test_celegans/README.md](test_celegans/README.md) for the full breakdown; nrCaeEleg92's raw reads weren't available to refresh.

```bash
cd test_celegans && ./benchmark.sh           # 200k reads per sample
MAX_READS=0 THREADS=16 ./benchmark.sh       # all reads, 16 threads
```

### Single-worm nematode benchmark (`test_single_nematode/`)

Compares mdax against the palindrome-correction method from a published
single-nematode-genome study (Lee et al., PRJNA953805), on real ONT WGA data
across 4 species, with published GenBank assemblies as ground truth. Its
README also records the history of the real/artefact gate, including the
locus-signature collisions that the fold-axis support count replaced.

```bash
cd test_single_nematode && cat README.md   # full results tables + writeup
```

### `irx` synthetic evaluation (`irx_tests/`)

~120 synthetic IRs planted across scenario groups (asymmetric arms, diverged, edge cases, GC bias, etc.). Recovery is high for typical cases; drops for very short or highly-diverged IRs.

```bash
cd irx_tests && bash test_suite.bash
```

## Performance

A full MDA library (HiFi or ONT) takes 3–20 minutes on 16 threads and under 4 GB of memory. `ont` mode uses Myers' bit-parallel edit distance for refinement and is comparable in speed.

## Algorithm overview

It's a two step pipeline:

- Pass 1: count fold-axis support
    - Reads are scanned for foldbacks using coarse minimizer-based self-similarity (forward minimizers matched against the read's own reverse complement, clustered by anti-diagonal).
    - For each foldback, the matched minimizers nearest the symmetry axis are recorded together with the distance between each minimizer and its inverted copy. That distance is twice the minimizer's distance from the axis, so two reads share a `(minimizer, distance)` key only if they fold about the same axis.
    - Reads whose evidence is too thin for those gates get a second look ("rescue"): any anti-diagonal with a few matched minimizer pairs is taken as a candidate, and accepted only if its shorter arm runs to a read end and aligns base-for-base against the reverse complement of the other (`--rescue-min-arm`, `--rescue-min-identity`). This is what finds a strand switch that happened shortly before the fragment ended. Candidates in the middle of a read are not rescued (`--rescue-anywhere` allows them): on real data those were as often genuine short inverted repeats as artefacts. An optional dense every-k-mer scan of the read ends (`--rescue-terminal on`) is off by default because it found almost nothing extra.
    - A genome-encoded inverted repeat folds about the same axis in every read that spans it. Independent amplification artefacts at the same locus fold about different ones. Counting reads per axis (rather than per locus) is what separates the two in deep libraries, where artefact junctions are dense enough that every locus "recurs".
    - The support threshold is then calibrated to the library: the chance co-location rate is estimated from the fraction of foldbacks nobody shares an axis with, and the threshold is raised until chance alone rarely reaches it, following the measured tail of the support histogram for over-dispersed (heavily amplified) libraries.

- Pass 2: correct reads
    - Each detected foldback is refined to a breakpoint and classified: at or above the support threshold it is genomic (`real`) and left alone; below it, it is an artefact.
    - Artefacts are cut at the strand switch, keeping the longer side, and the kept side is re-examined recursively for nested artefacts. When the arms are separated by an un-mirrored spacer, the spacer is dropped on ONT; on HiFi mdax keeps whichever side of it library k-mer counts say is genomic.

Measured against a reference (random read samples from *C. elegans* WGA HiFi libraries aligned to WBcel235; `synth_bench/refbased_recall_eval.py`), this cuts 79% of true foldback artefacts at ~100x and 74% at ~720x (77.5% on MDA ONT), cutting about 2% of normal reads, with under 1% of true artefacts kept as `real`. Almost all of what remains is not a detection failure: about a fifth of reads with a nearby strand switch contain no inverted copy at all (the read ended before the new strand had been copied back as far as the old one), so no self-comparison can see them. Among foldbacks that do contain at least 100 bp of inverted copy, 95–99% are cut. Cutting the rest at reference-detected junctions gave no assembly gain. In the TSV, `support_n` is the number of reads (including the read itself) sharing the fold axis, `support_span` is the number of key elements the junction had, and `support_rank_frac` is no longer computed (always 0).

## `irx`: detecting 'natural' IR's from assemblies

We provide the program `irx` which uses the core underlying logic from `mdax` to identify IR's from genome assembly data. It's very fast and parallelises over windows.

irx is best suited to detection of moderate-to-large inverted repeats with strong minimizer support. In synthetic tests, recovery was high for kilobase-scale and larger IRs, including asymmetric, edge-positioned, and wide-spacer cases, while recovery dropped for very short (~ <500bp) and/or strongly diverged repeats (~ <75% identity). The alignment-based tir_ident metric was generally accurate for recovered calls. Reported arm lengths tend to be slightly conservative, and reported spacers slightly inflated, relative to truth. These observations are based on the defualt general parameters.

### Usage

This outputs a table of putative IR's, and extracts with `-f` the repeats into a fasta file.

```bash
# most basic usage
irx -b out.bed in.fasta 
# save fasta output
irx -b out.bed -f out.fasta in.fasta
# save output diagnostic HTML, save bed through redirection instead
irx --html out.html in.fa > out.bed
```

More detailed usage.

```txt
`irx` scans an assembly FASTA(.gz) for inverted-repeat (IR) candidates using the `mdax` foldback detector.

Workflow per contig:
- generate overlapping windows
- in each window, find up to K candidate foldback bins (top-k anti-diagonals)
- refine breakpoint position for each candidate
- derive left/right arm bounds from minimizer matchpoints
- best-first deduplicate across windows
- emit BED-like TSV rows

Output is 0-based, half-open coordinates (BED convention).

TSV columns (header line begins with '#'):
#contig	start	end	name	score	strand	break_pos	identity_est	tir_ident	matches	span	la0	la1	ra0	ra1	win_start	win_end	kept_pts	bin

Column meanings:
contig         First token of FASTA record id
start,end      Candidate interval covering both arms (contig coords)
name           Constant 'IR'
score          identity_est scaled to 0..1000 (rounded)
strand         '.' (IR is not strand-specific here)
break_pos      Refined split position in contig coords
identity_est   Refinement identity estimate (0..1). Computed locally around breakpoint during refinement.
tir_ident      Estimated sequence identity of the IR arms (la0..la1 vs revcomp(ra0..ra1))
matches        Coarse supporting minimizer matches (fb.matches)
span           Coarse span between arms (fb.span, bp)
la0,la1        Left arm bounds (contig coords)
ra0,ra1        Right arm bounds (contig coords)
contig_len     Length of the contig/chromosome (bp)
win_start,end  Window bounds that produced this hit (contig coords)
kept_pts       Number of matchpoints retained near refined anti-diagonal (0 if unfiltered fallback)
bin            Anti-diagonal bin id (debug / dedup key component)



Usage: irx [OPTIONS] <input>

Arguments:
  <input>
          Input assembly in FASTA format. If the filename ends with .gz, `irx` will transparently decompress it.
          
          Note: `irx` will trim to the first whitespace-delimited token for output.

Options:
  -b, --bed <bed>
          Output path for the BED-like TSV. Use '-' to write to stdout.
          
          The first line is a header prefixed with '#'.
          
          [default: -]

  -f, --fasta <fasta>
          Write the sequences corresponding to accepted inverted-repeat (IR)
          intervals to a FASTA file.

          Each emitted IR interval is extracted from the input assembly using the 
          reported `start`–`end` coordinates and written as a FASTA record. The FASTA 
          header contains the contig name and interval coordinates so sequences can be 
          easily traced back to their source region.

          If '-' is supplied as the filename, FASTA output is written to stdout.
          This is useful for piping IR sequences directly into downstream tools
          (e.g. alignment, repeat analysis, or clustering).

          If this option is not provided, no FASTA sequences are written and only the
          BED-like interval table is produced.

      --fasta-wrap <fasta_wrap>
          Set the line width used when writing FASTA sequences.

          Sequences are wrapped to the specified number of characters per line in the
          output FASTA file. This improves readability and compatibility with tools
          that expect fixed-width FASTA formatting.
          
          A value of 0 disables wrapping and writes each sequence on a single line.
          
          The default value of 60 follows a common FASTA formatting convention.
          
          [default: 60]

      --max-interval-bp <max_interval_bp>
          Maximum length of an emitted IR interval (bp). Helps prevent window-scale inflation.
          
          [default: 250000]

      --window-len <window_len>
          Length of the sliding window used to scan contigs.
          
          - If 0, the entire contig is processed as a single window.
          - Otherwise, windows are [start, start+window_len) clamped to contig length.
          
          Smaller windows can increase sensitivity to local IRs but may increase runtime.
          
          [default: 250000]

      --step <step>
          Step size between consecutive window starts.
          
          Smaller step => more overlap => higher sensitivity but more duplicate discoveries.
          Deduplication is applied downstream, but runtime still increases with overlap.
          
          [default: 50000]

      --hits-per-window <hits_per_window>
          Maximum number of coarse foldback candidates per window.
          
          Internally, `mdax` bins minimizer matchpoints by quantised anti-diagonal (d = p1 + p2). `--hits-per-window` selects the top-K bins by a cheap proxy rank (count + span).
          
          Higher values may find multiple IRs within the same window, at increased cost.
          
          [default: 4]

      --dedup-bp <dedup_bp>
          Quantisation size (in bp) used in deduplication.
          
          We do two-stage dedup:
          A) breakpoint/bin dedup on (break_pos/q, diag_bin)
          B) overlap dedup on intervals quantised to (start/q, end/q)
          
          Larger values merge more aggressively (fewer outputs) and speed up overlap checks.
          
          [default: 200]

      --diag-widen <diag_widen>
          When deriving arm bounds near the refined split, we keep matchpoints where:
          |(p1+p2) - 2*split| <= fold_diag_tol * diag_widen
          
          Smaller values are stricter (tighter arms, fewer cap hits).
          Larger values are more permissive (can inflate intervals).
          
          [default: 4]

      --trim-lo <trim_lo>
          After filtering matchpoints near the split, we derive arm bounds using
          trimmed quantiles rather than min/max to avoid single-point inflation.
          
          Example: --trim-lo 0.10 keeps the 10th percentile as the lower bound.
          
          [default: 0.10]

      --trim-hi <trim_hi>
          Upper quantile used for trimmed arm bounds.
          Example: --trim-hi 0.90 uses the 90th percentile as the upper bound.
          
          [default: 0.90]

      --classify-ir
          Disable IR class columns (arm_len, spacer, ir_class) in TSV.

      --immediate-bp <immediate_bp>
          Max spacer (bp) to call an IR 'immediate'.
          
          [default: 2000]

      --wide-bp <wide_bp>
          Min spacer (bp) to call an IR 'wide'.
          
          [default: 100000]

      --threads <threads>
          Number of worker threads.
          
          - 0: use Rayon default (typically = number of logical CPUs)
          - N>0: use exactly N threads via a *local* Rayon thread pool
          
          Using a local pool ensures the requested thread count is respected even if another crate initialises the Rayon global pool.
          
          [default: 0]

      --cap-policy <cap_policy>
          When the inferred IR interval length exceeds --max-interval-bp:
          clamp     Clamp to the capdrop      Discard the candidate entirely
          penalize  Clamp, but strongly reduce its rank so non-capped hits win dedup
          
          Recommended: penalize (prevents huge piles of exactly-cap intervals).
          
          [default: penalize]
          [possible values: clamp, drop, penalize]

      --html <html>
          Write an interactive HTML dashboard containing embedded IR data and diagnostic plots.
          
          The report is self-contained except for any external JS libraries referenced by the template.

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

#### TSV output

Core bits:

- contig: Contig / chromosome name (first FASTA header token)
- start: IR interval start (0-based, inclusive)
- end: IR interval end (0-based, exclusive)
- name: Always IR
- score: Identity estimate ×1000 (rounded; 0–1000 scale)
- strand: Always . (IR is strand-symmetric)

Breakpoints & similarity:

- break_pos: Refined split position (putative IR center)
- identity_est: Estimated arm identity (0–1 float)
- tir_ident: sequence identity between the inferred inverted-repeat arms (la0..la1 vs revcomp(ra0..ra1)).
- matches: Number of supporting minimizer matches
- span: Coarse arm span estimate from foldback geometry

Arm coordinates:

- la0, la1: Left arm start/end
- ra0, ra1: Right arm start/end

Windows:

- win_start, win_end: Sliding window that produced this hit
- kept_pts: Number of diagonal-consistent matchpoints used for arm bounds
- bin: Anti-diagonal bin ID (used internally for ranking & dedup)

Annotation:

- arm_len: Length of the shorter arm (bp)
- spacer: Distance between arms (ra0 - la1)
- ir_class: Structure class of IR (immediate, spaced, or wide). Configurable in CLI.

#### Fasta header format

```txt
>contig:start-end|break=BREAK|ident=I|tir_ident=T|matches=M|span=S|bin=B|class=C|spacer=SP
```

- contig:start-end — contig name (first FASTA token) and emitted interval in contig coordinates, 0-based half-open.
- break=BREAK — refined breakpoint (foldback split) position on the contig, in bp.
- ident=I — refinement identity estimate between the two arms near the breakpoint (approximate; not a full-length alignment), formatted as 0.XXX.
- tir_ident=T — an accurate identity between the two terminal inverted repeats, derived from `block-aligner`.
- matches=M — number of supporting minimizer matchpoints in the coarse detector.
- span=S — coarse span estimate between arms from the detector (bp).
- bin=B — anti-diagonal bin id used internally for candidate discovery/dedup.
- class=C — IR class inferred from arm spacing: immediate, spaced, wide, or unknown (if arm bounds unavailable).
- spacer=SP — estimated spacer length between arms (ra0 - la1 after arm normalisation); can be <=0 if arms abut/overlap.

#### The stderr output

For example:

```console 
[irx] OW119596.1: wins=662 cands=634 emitted=92 clamped=0 (0.0%) near_cap=0 (0.0%) arm_missing=0 precap_too_big=0 stageA_dups=0 stageB_dups=542
```

- Record name of the fasta
- wins: Number of sliding windows scanned on the contig
- cands: Total candidates produced by the window scanning after coarse and refine filters, before per contig de-dup
- emitted: Number of candidates emitted to TSV/FASTA after de-dup
- clamped (X%): Amongst emitted candidates, how many were clamped due to `--max-interval-bp`. This is policy dependent. % (clamped / emitted)
- near_cap: Number of emitted IR intervals whose final length is ≥ 95% of `--max-interval-bp`.
- arm_missing: Number of emitted IRs where no reliable arm bounds were derived from minimizer matchpoints.
- precap_too_big: Number of emitted IRs whose raw inferred interval exceeded --max-interval-bp before any clamping/drop/penalization.
- stageA_dups: duplicates removed by Stage A dedup (breakpoint/bin key).
- stageB_dups: duplicates removed by Stage B dedup (interval overlap on quantized coordinates).

