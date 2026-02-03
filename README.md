# `mdax`

Read in HiFi or ONT raw data generated from Multiple Displacement Amplification methods, and correct reads. This tool currently only corrects 'foldback' or 'invert' chimeras.

More information to follow.

There are lots of options (see below), but to get started:

```bash
mdax --output <output.fa> <in.fa> > <out.table>
```

Default outputs are the fasta file of corrected reads, and a table of all detected foldback chimeras in the readset.

Current usage is:

```bash
Remove inverted duplication chimeras, fast

Usage: mdax [OPTIONS] --output <output> <input>

Arguments:
  <input>  Input file

Options:
  -r, --report <report>
          Write a TSV report of detected concatemers (use '-' for stdout) [default: -]
  -o, --output <output>
          Output FASTA (optionally .gz not supported yet)
  -k <k>
          Minimizer k-mer size [default: 17]
  -w <w>
          Minimizer window size [default: 21]
      --min-span <min_span>
          Override: minimum evidence span (bp) for calling a foldback/concatemer
      --min-identity <min_identity>
          Override: minimum refinement identity for foldback second pass
      --min-support <min_support>
          Override: minimum support reads for 'real' foldback
      --split-tol-bp <split_tol_bp>
          Override: max split span (bp) within a signature cluster
      --max-depth <max_depth>
          Override: maximum recursion depth for chopping
      --end-guard <end_guard>
          Do not call breakpoints within this distance of read ends (bp) [default: 1000]
      --refine-window <refine_window>
          Refinement half-window around coarse split (bp) [default: 100]
      --refine-arm <refine_arm>
          Refinement arm length (bp) [default: 200]
      --max-ed-rate <max_ed_rate>
          ONT refinement bandwidth as fraction of arm [default: 0.25]
      --min-matches <min_matches>
          Override: minimum minimizer matches supporting a foldback candidate
      --fold-diag-tol <fold_diag_tol>
          Foldback anti-diagonal bucketing tolerance (bp) for pos1+pos2 clustering [default: 120]
      --concat-diag-tol <concat_diag_tol>
          Concatemer delta bucketing tolerance (bp) for (pos2-pos1) clustering [default: 200]
      --forward-only <forward_only>
          Use strand-specific forward-only minimizers (true) or canonical minimizers (false) [default: true] [possible values: true, false]
  -m, --refine-mode <refine_mode>
          Refine breakpoints for reads [default: hifi] [possible values: hifi, ont]
      --mode <mode>
          Detection/correction strictness: strict reduces calls, permissive increases calls [default: balanced] [possible values: strict, balanced, permissive]
      --sig-flank-bp <sig_flank_bp>
          Window (bp) around split for match-based fingerprinting
      --sig-take <sig_take>
          Number of minimizer hashes retained for foldback fingerprint
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```

## Performance

It's reasonably fast, processing ~10Gb of gzipped fasta data in a few minutes in `hifi` mode, slightly longer in `ont` mode on my local Mac with 8 threads. More performance to follow.

## Algorithm overview

It's a two step pipeline:

- Pass 1: Build Support Map
    - Reads are scanned to detect potential foldback events using coarse minimizer-based self-similarity.
    - Breakpoints are refined, and events are filtered using identity thresholds and stability metrics.
    - A support map is generated, clustering junction signatures observed across reads.

- Pass 2: Correct Reads
    - Reads are reprocessed, and each foldback event is classified as a genomic event (real) or an artefact based on the support map.
    - Artefacts are recursively trimmed or corrected, ensuring high specificity against real, genome-templated events.
