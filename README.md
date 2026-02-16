# `mdax`

Read in HiFi or ONT raw data generated from Multiple Displacement Amplification methods, and correct reads. This tool currently only corrects 'foldback' or 'invert' chimeras (`mdax`). We also provide a binary to search for long Inverted Repeats (IR's) in assemblies (`irx`).

## Install

Currently you'll have to clone this repo, have Rust installed, and compile:

```bash
# install Rust first (https://rust-lang.org/tools/install/)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# clone and cd
git clone github.com/ARU-life-sciences/mdax
cd mdax
# use SIMD instructions
RUSTFLAGS='-C target-feature=+aes,+sse2,+avx2' cargo build --release
# or if you want in PATH
RUSTFLAGS='-C target-feature=+aes,+sse2,+avx2' cargo install --path=.
```

## Usage

There are lots of options (see below), but to get started:

```bash
mdax --output <output.fa> <in.fa> > <out.table>
```

Default outputs are the fasta file of corrected reads, and a table of all detected foldback chimeras in the readset.

Current usage is:

```bash
Remove inverted duplication chimeras from fasta files.

Usage: mdax [OPTIONS] --output <output> <input>

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
          - support_n: Number of other reads supporting the same foldback junction.
          - support_span: Total supporting span (bp) aggregated across supporting reads.
          - decision: Final classification based on thresholds (real, artefact, unknown, low_ident).
          
          Use '-' to write to stdout
          
          
          [default: -]

  -o, --output <output>
          Output FASTA containing the kept/corrected reads.
          
          For each input read, mdax writes one output record containing the sequence you decided to keep:
          - If no foldback is called, the read is written unchanged.
          - If a foldback is called and judged an artefact, mdax cuts/chops the read (possibly recursively) and writes the kept portion.
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
          Override: minimum number of reads that must share the same junction signature for the foldback to be treated as 'real'.
          
          Interpretation:
          - Pass1 builds a support map keyed by a fingerprint/signature of the junction.
          - In pass2, if a signature is observed in >= min-support reads (and split clustering is tight enough),
          `mdax` treats that junction as genome-templated and does NOT cut it by default.
          
          How it affects cutting:
          - Increasing min-support means harder to be considered real. More foldbacks treated as artefacts therefore more cutting.
          - Decreasing min-support means easier to be considered real. Fewer cuts, but risk of keeping recurrent artefacts.
          
          Guidance:
          - For deep coverage datasets, you can raise this to be more conservative about calling 'real'.
          - For shallow coverage, keep it low, otherwise almost nothing will qualify as real.

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
          
          [default: 1000]

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
          
          [default: 200]

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
          
          Note: we use canonical minimizers for *signatures* even if forward-only for detection;
          this flag applies to the detector’s minimizer sampling.
          
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
          max_depth, and signature parameters).
          
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

      --sig-flank-bp <sig_flank_bp>
          Signature flank size (bp) used to compute the junction fingerprint.
          
          Signatures are used to build the pass1 support map and to look up support in pass2.
          They should be stable across small split jitter but specific enough not to collide.
          
          How it affects 'real vs artefact' decisions:
          - Increasing sig-flank-bp gives more context included and fewer collisions, more stable support clusters,
          but slightly more compute and potentially more sensitivity to distant repeats.
          - Decreasing sig-flank-bp gives signature becomes more local and is faster and less affected by distant repeats,
          but higher chance different junctions produce the same signature (collisions).
          
          If you see many signatures that are 'in support' but with wide split scatter, consider increasing
          sig_flank_bp and/or sig_take to improve specificity.

      --sig-take <sig_take>
          Number of minimizer-derived values retained when constructing the junction signature.
          
          This controls signature specificity:
          - Larger sig-take; more information in the signature, fewer collisions, better support discrimination.
          - Smaller sig-take; less information, more collisions, more chance of misclassifying artefacts as real (or vice versa).
          
          Cost:
          - Larger values are slightly slower and may fail more often if there are too few minimizers in the flank window.
          (In that case you may see fewer sig_ok events.)
          
          Guidance:
          - If you want robust support-map behaviour, prefer increasing sig-take before relaxing split tolerance.
          - If you process very short reads or tiny flanks, keep take modest so signatures can still be formed.

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

## Performance

It's reasonably fast, processing ~100Gb fasta data in 5/6 minutes in `hifi` mode, slightly longer in `ont` mode on my local Mac with 8 threads. More performance to follow.

## Algorithm overview

It's a two step pipeline:

- Pass 1: Build Support Map
    - Reads are scanned to detect potential foldback events using coarse minimizer-based self-similarity.
    - Breakpoints are refined, and events are filtered using identity thresholds and stability metrics.
    - A support map is generated, clustering junction signatures observed across reads.

- Pass 2: Correct Reads
    - Reads are reprocessed, and each foldback event is classified as a genomic event (real) or an artefact based on the support map.
    - Artefacts are recursively trimmed or corrected, ensuring high specificity against real, genome-templated events.

## `irx`: detecting 'natural' IR's from assemblies

We provide the program `irx` which uses the core underlying logic from `mdax` to identify IR's from assembly data. It's very fast and parallelises over windows. I'll expand on this but for now...

### Usage

This outputs a table of putative IR's, and extracts with `-f` the repeats into a fasta file.

```bash
irx -b ir_repeats.tsv -f ir_repeats.fa assembly.fasta.gz
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
>contig:start-end|break=BREAK|ident=I|matches=M|span=S|bin=B|class=C|spacer=SP
```

- contig:start-end — contig name (first FASTA token) and emitted interval in contig coordinates, 0-based half-open.
- break=BREAK — refined breakpoint (foldback split) position on the contig, in bp.
- ident=I — refinement identity estimate between the two arms near the breakpoint (approximate; not a full-length alignment), formatted as 0.XXX.
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

