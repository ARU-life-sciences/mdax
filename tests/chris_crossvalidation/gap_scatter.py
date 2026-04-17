#!/usr/bin/env python3
"""
gap_scatter.py — Test 4: mdax gap_est vs Chris's supplementary-alignment distance

Background
----------
Both mdax and Chris's pipeline produce a measure of the template-switch gap —
the number of bases the MDA enzyme skipped before reversing:

  mdax:    gap_est = 2 * best_delta from the delta-scan (bp)
  Chris:   distance = min genomic distance between primary and supplementary
           alignment endpoints (supplementary_patterns table, bp)

These measure slightly different things:
  - mdax measures the gap within the read (read-local coordinates)
  - Chris measures the gap in reference coordinates (genomic)

For a simple intrachromosomal foldback, the two should correlate tightly.
Divergence indicates:
  - Complex reads (multiple supplementary alignments)
  - Reads where the foldback junction is not the primary alignment boundary
  - Tandem concatemers, where the gap structure is non-trivial

Input files
-----------
  ../../chris_data/ERR12263839.tsv              mdax report
  ../../chris_data/supplementary_patterns_head.tsv  Chris's split-read events

Output
------
  Prints summary statistics.
  If matplotlib is available and --plot is given, saves a scatter PNG.

Usage
-----
  python3 gap_scatter.py [--plot gap_scatter.png] [--max-gap 10000]
"""

import argparse
import csv
import os
import statistics
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "..", "chris_data")

MDAX_TSV  = os.path.join(DATA_DIR, "ERR12263839.tsv")
SUPP_TSV  = os.path.join(DATA_DIR, "supplementary_patterns_head.tsv")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--plot",    default=None, metavar="PATH.png",
                help="Save scatter plot to this PNG path (requires matplotlib)")
ap.add_argument("--max-gap", type=int, default=10000, metavar="BP",
                help="Exclude pairs where either gap measure exceeds this (default: 10000)")
args = ap.parse_args()

# ---------------------------------------------------------------------------
# Guard: check data files exist
# ---------------------------------------------------------------------------

MISSING = [p for p in (MDAX_TSV, SUPP_TSV) if not os.path.exists(p)]
if MISSING:
    print("SKIP: required data files not found:", file=sys.stderr)
    for p in MISSING:
        print(f"  {p}", file=sys.stderr)
    sys.exit(0)

# ---------------------------------------------------------------------------
# Load mdax: short_id → gap_est
# ---------------------------------------------------------------------------

print("Loading mdax report …", flush=True)
mdax_gap: dict[str, float] = {}
with open(MDAX_TSV, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        short_id = row["read_id"].split()[0]
        try:
            mdax_gap[short_id] = float(row["gap_est"])
        except (KeyError, ValueError):
            pass
print(f"  {len(mdax_gap):,} reads with gap_est in mdax report")

# ---------------------------------------------------------------------------
# Load Chris's supplementary patterns
# Columns: read_id, prim_chrom, prim_start, prim_end, prim_strand,
#          prim_aln_len, supp_chrom, supp_start, supp_end, supp_strand,
#          supp_aln_len, distance, event_type
#
# Keep the MINIMUM distance per read (some reads have multiple supplementaries;
# the closest one is most likely to correspond to the primary foldback gap).
# ---------------------------------------------------------------------------

print("Loading supplementary patterns …", flush=True)
chris_dist: dict[str, float] = {}
with open(SUPP_TSV, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rid = row.get("read_id", "")
        raw_dist = row.get("distance", "inter")
        if raw_dist == "inter":
            continue  # interchromosomal: no comparable distance
        try:
            d = float(raw_dist)
        except ValueError:
            continue
        # Keep minimum distance per read
        if rid not in chris_dist or d < chris_dist[rid]:
            chris_dist[rid] = d

print(f"  {len(chris_dist):,} reads with intra-chromosomal supplementary distance")

# ---------------------------------------------------------------------------
# Overlap and paired values
# ---------------------------------------------------------------------------

overlap = sorted(set(mdax_gap) & set(chris_dist))
print(f"  Overlap: {len(overlap):,} reads in both tables\n")

if not overlap:
    print("No overlap found — the head tables (10k rows) may not share reads.")
    print("Run with the full tables for a meaningful comparison.")
    sys.exit(0)

pairs = [
    (mdax_gap[rid], chris_dist[rid])
    for rid in overlap
    if mdax_gap[rid] <= args.max_gap and chris_dist[rid] <= args.max_gap
]
print(f"  Pairs within max_gap={args.max_gap:,} bp: {len(pairs):,}")

if len(pairs) < 5:
    print("  Too few pairs for meaningful statistics.")
    sys.exit(0)

mdax_vals  = [p[0] for p in pairs]
chris_vals = [p[1] for p in pairs]

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

print("=" * 60)
print("  Summary statistics")
print("=" * 60)
print(f"  {'Metric':<25}  {'mdax gap_est':>13}  {'chris distance':>14}")
print(f"  {'─'*25}  {'─'*13}  {'─'*14}")
for label, vals in [("mean (bp)", (statistics.mean(mdax_vals),
                                   statistics.mean(chris_vals))),
                    ("median (bp)", (statistics.median(mdax_vals),
                                    statistics.median(chris_vals))),
                    ("std (bp)", (statistics.stdev(mdax_vals) if len(mdax_vals) > 1 else float("nan"),
                                  statistics.stdev(chris_vals) if len(chris_vals) > 1 else float("nan"))),
                    ("min (bp)", (min(mdax_vals), min(chris_vals))),
                    ("max (bp)", (max(mdax_vals), max(chris_vals)))]:
    print(f"  {label:<25}  {vals[0]:>13.1f}  {vals[1]:>14.1f}")

# Pearson correlation
def pearson(xs, ys):
    if len(xs) < 3:
        return float("nan")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denom = (sum((x - mx)**2 for x in xs) * sum((y - my)**2 for y in ys)) ** 0.5
    return num / denom if denom else float("nan")

r = pearson(mdax_vals, chris_vals)
print(f"\n  Pearson r = {r:.4f}  (n={len(pairs):,})")

# Fraction within 2× of each other
within_2x = sum(1 for m, c in pairs if c > 0 and 0.5 <= m/c <= 2.0) / len(pairs)
print(f"  Fraction within 2× of each other: {within_2x:.1%}")

# Reads where mdax gap_est ≈ 0 but chris distance is large
low_mdax_high_chris = [(m, c) for m, c in pairs if m < 100 and c > 500]
high_mdax_low_chris = [(m, c) for m, c in pairs if m > 500 and c < 100]
print(f"\n  mdax_gap<100 but chris_dist>500: {len(low_mdax_high_chris):,} reads")
print(f"  mdax_gap>500 but chris_dist<100: {len(high_mdax_low_chris):,} reads")
print()

# ---------------------------------------------------------------------------
# Interpretation hints
# ---------------------------------------------------------------------------

print("  Interpretation:")
if r > 0.7:
    print("  Strong positive correlation — both methods measuring the same gap.")
elif r > 0.4:
    print("  Moderate correlation — agreement on large gaps; divergence on small gaps.")
else:
    print("  Weak correlation — consider whether head tables have enough overlap,")
    print("  or whether multi-supplementary reads are skewing Chris's distance.")

if len(low_mdax_high_chris) > len(pairs) * 0.1:
    print("  Many reads where mdax sees no gap but Chris sees a large genomic distance.")
    print("  Possible causes: reads where the gap exceeds --max-jump-clip (1000 bp),")
    print("  or reads with complex multi-supplementary structure.")

if len(high_mdax_low_chris) > len(pairs) * 0.1:
    print("  Many reads where mdax sees a gap but Chris's supplementary distance is small.")
    print("  Possible causes: tandem concatemers, or the closest supplementary alignment")
    print("  happens to land near the primary even though the true gap is large.")

# ---------------------------------------------------------------------------
# Optional scatter plot
# ---------------------------------------------------------------------------

if args.plot:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=(7, 6))

        # Scatter, log-scale to handle the wide range
        x = np.array(chris_vals)
        y = np.array(mdax_vals)
        ax.scatter(x, y, s=10, alpha=0.4, color="#3B8BD4", rasterized=True)

        # y=x reference line
        lim = max(args.max_gap, max(x.max(), y.max()) * 1.1)
        ax.plot([0, lim], [0, lim], "--", color="#888", lw=1, label="y = x")

        ax.set_xlabel("Chris's supplementary-alignment distance (bp)")
        ax.set_ylabel("mdax gap_est (bp)")
        ax.set_title(f"Template-switch gap: mdax vs Chris's pipeline\n"
                     f"n={len(pairs):,}  Pearson r={r:.3f}  max_gap={args.max_gap:,} bp")
        ax.set_xlim(0, args.max_gap)
        ax.set_ylim(0, args.max_gap)
        ax.legend(fontsize=8, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(args.plot, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Scatter plot saved to: {args.plot}")
    except ImportError:
        print("  matplotlib not available — skipping scatter plot.")
        print("  Install with:  pip install matplotlib numpy")
