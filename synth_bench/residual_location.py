#!/usr/bin/env python3
"""For reads mdax cut that still carry a foldback after correction: is the
surviving junction sitting at the new read end (i.e. the cut left part of
the same fold behind, typically the un-mirrored spacer) or in the interior
(an independent junction mdax never examined)? Also: how much of the
surviving junction's shorter arm is explained by the retained spacer.

usage: residual_location.py TRUTH_SAM CORR_SAM REPORT RAW_FA CORR_FA
"""
import collections
import importlib.util
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mr", os.path.join(here, "measure_residual.py"))
# reuse the parsing helpers without running measure_residual's report
src = open(os.path.join(here, "measure_residual.py")).read()
helpers = src.split("# ---------------------------------------------------------------- load")[0]
helpers = helpers.replace("truth_sam, corr_sam, report, raw_fa, corr_fa = sys.argv[1:6]", "")
ns = {}
exec(compile(helpers, "measure_residual_helpers", "exec"), ns)
load_sam, fold_junctions, load_fa = ns["load_sam"], ns["fold_junctions"], ns["load_fa"]

import csv

truth_sam, corr_sam, report, raw_fa, corr_fa = sys.argv[1:6]
tsegs, _ = load_sam(truth_sam)
csegs, crlen = load_sam(corr_sam)
rows = {}
for r in csv.DictReader(open(report), delimiter="\t"):
    rows[r["read_id"].split()[0]] = r
raw = load_fa(raw_fa)
corr = load_fa(corr_fa)

END = 300  # a surviving junction within this many bp of a cut end counts as "at the cut"
res = collections.Counter()
by_jn = collections.defaultdict(collections.Counter)
for q, ss in csegs.items():
    js = fold_junctions(ss)
    if not js or not fold_junctions(tsegs.get(q, [])):
        continue
    r = rows.get(q)
    if r is None or r["decision"] != "artefact":
        continue
    rr, cc = raw.get(q), corr.get(q)
    if rr is None or cc is None or len(cc) == len(rr):
        continue
    ks = rr.find(cc)
    if ks < 0:
        continue
    ke = ks + len(cc)
    L = len(cc)
    cut_left = ks > 0      # left end of the corrected read is a cut end
    cut_right = ke < len(rr)
    nj = len(fold_junctions(tsegs[q]))
    for a, b in js:
        jpos = (a[1] + b[0]) / 2
        short = min(a[1] - a[0], b[1] - b[0])
        # which arm touches the cut end, and is that arm short (i.e. a stub)?
        at_cut = (cut_left and a[0] < END and jpos < L / 2) or (cut_right and L - b[1] < END and jpos > L / 2)
        stub = (cut_left and a[0] < END and (a[1] - a[0]) == short) or (cut_right and L - b[1] < END and (b[1] - b[0]) == short)
        key = "stub arm at the cut end (fold partly left behind)" if (at_cut and stub) else \
              "junction near cut end, long arm at end" if at_cut else "interior junction (not at a cut end)"
        res[key] += 1
        by_jn["single-junction raw read" if nj == 1 else "multi-junction raw read"][key] += 1
tot = sum(res.values())
print(f"surviving foldback junctions in reads mdax cut: {tot}")
for k, v in res.most_common():
    print(f"  {k:52s} {v:7d} {100 * v / tot:5.1f}%")
for g, c in by_jn.items():
    t = sum(c.values())
    print(f"  {g}: n={t}  " + "  ".join(f"{k.split(' (')[0]}={100 * v / t:.1f}%" for k, v in c.most_common()))
