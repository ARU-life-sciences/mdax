#!/usr/bin/env python3
"""Diagnostic A: for each real inverted-repeat locus, do its covering reads
share ONE signature cluster (support just under-sampled -- the original
hypothesis) or get split across MANY distinct sig_band clusters (support
fragmentation -- the alternative hypothesis raised after cov_tier showed no
FPR gradient)?

Requires a report.tsv built with the temporary `sig_band` diagnostic column
(pipeline.rs, reverted after use).
"""
import csv
import sys
from collections import defaultdict

truth_path, report_path, bed_path = sys.argv[1], sys.argv[2], sys.argv[3]

overlap_loci = {}
fully_spans = {}
true_class = {}
with open(truth_path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rid = row["read_id"]
        true_class[rid] = "artefact" if row["type"] == "chimera" else "genuine"
        fully_spans[rid] = row["fully_spans_locus"] == "1"
        overlap_loci[rid] = [x for x in row["overlap_loci"].split(",") if x and x != "NA"]

cov_tier = {}
with open(bed_path) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        cov_tier[p[0]] = float(p[7])

pred = {}
band = {}
with open(report_path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rid = row["read_id"].split()[0]
        d = row["decision"]
        pred[rid] = "genuine" if d == "real" else ("artefact" if d == "artefact" else "other")
        band[rid] = row["sig_band"]

# hard case: genuine reads fully spanning >=1 real locus
hard = [r for r in overlap_loci if true_class.get(r) == "genuine" and fully_spans.get(r) and overlap_loci[r]]
print(f"Hard-case reads (genuine, fully spans a real IR locus): {len(hard)}")

by_locus = defaultdict(list)
for r in hard:
    # fully_spans_locus=1 with possibly >1 overlap_loci entries; attribute to
    # all listed loci (rare multi-locus reads), consistent with classify_eval.py.
    for locus in overlap_loci[r]:
        by_locus[locus].append(r)

print(f"Distinct real IR loci covered by >=1 hard-case read: {len(by_locus)}\n")

rows = []
for locus, reads in by_locus.items():
    bands_here = [band[r] for r in reads if band.get(r) not in (None, "NA")]
    distinct_bands = set(bands_here)
    no_sig = sum(1 for r in reads if band.get(r) in (None, "NA"))
    fp = sum(1 for r in reads if pred.get(r) == "artefact")
    tn = sum(1 for r in reads if pred.get(r) == "genuine")
    other = sum(1 for r in reads if pred.get(r) == "other")
    rows.append(dict(
        locus=locus, cov_tier=cov_tier.get(locus, float("nan")), n_reads=len(reads),
        n_with_sig=len(bands_here), n_distinct_bands=len(distinct_bands),
        n_no_sig=no_sig, fp=fp, tn=tn, other=other,
    ))

rows.sort(key=lambda r: -r["n_distinct_bands"])

print(f"{'locus':<8}{'cov_tier':>9}{'n_reads':>8}{'n_w_sig':>8}{'distinct_bands':>15}{'no_sig':>7}{'FP':>4}{'TN':>4}{'other':>6}")
for r in rows:
    print(f"{r['locus']:<8}{r['cov_tier']:>9.3f}{r['n_reads']:>8}{r['n_with_sig']:>8}"
          f"{r['n_distinct_bands']:>15}{r['n_no_sig']:>7}{r['fp']:>4}{r['tn']:>4}{r['other']:>6}")

n_frag = sum(1 for r in rows if r["n_distinct_bands"] > 1)
n_single = sum(1 for r in rows if r["n_distinct_bands"] == 1)
n_zero = sum(1 for r in rows if r["n_distinct_bands"] == 0)
print(f"\nSummary over {len(rows)} loci with >=1 hard-case read:")
print(f"  loci where covering reads share exactly 1 signature cluster: {n_single}")
print(f"  loci where covering reads FRAGMENT into >1 distinct signature clusters: {n_frag}")
print(f"  loci where no covering read got any signature at all (sig build failed): {n_zero}")

multi_read_loci = [r for r in rows if r["n_with_sig"] >= 2]
print(f"\nOf {len(multi_read_loci)} loci with >=2 signature-bearing reads "
      f"(where fragmentation is even possible):")
frag_among_multi = sum(1 for r in multi_read_loci if r["n_distinct_bands"] > 1)
print(f"  fragmented (>1 distinct band): {frag_among_multi} "
      f"({100*frag_among_multi/len(multi_read_loci):.1f}%)")
print(f"  single shared band despite >=2 reads: {len(multi_read_loci) - frag_among_multi}")

# Does fragmentation correlate with FP rate at the locus level?
frag_loci = [r for r in rows if r["n_distinct_bands"] > 1]
single_loci = [r for r in rows if r["n_distinct_bands"] == 1]


def fp_rate(subset):
    tot_fp = sum(r["fp"] for r in subset)
    tot_n = sum(r["n_reads"] for r in subset)
    return tot_fp / tot_n if tot_n else float("nan")


print(f"\nFP rate (reads wrongly cut) at fragmented loci:  {fp_rate(frag_loci):.3f} (n_loci={len(frag_loci)})")
print(f"FP rate (reads wrongly cut) at single-cluster loci: {fp_rate(single_loci):.3f} (n_loci={len(single_loci)})")
