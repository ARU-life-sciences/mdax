#!/usr/bin/env python3
"""Join refbased_clip.py's reference-based ground truth against mdax's
own real/cut decisions at several min_support candidates, for the same
read sample."""
import csv
import sys

refbased_tsv, *report_paths = sys.argv[1:]

truth = {}  # read_id -> "artefact" | "genuine" (unmapped reads excluded: no evidence either way)
with open(refbased_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["action"] == "clipped":
            truth[row["read_id"]] = "artefact"
        elif row["reason"].split(";")[0] in ("no_supplementary", "supplementary_not_foldback_like"):
            truth[row["read_id"]] = "genuine"
        # else: unmapped / not_in_sam / degenerate_clip -> no label, excluded

print(f"Ground truth (reference-based, {len(truth)} reads): "
      f"artefact={sum(1 for v in truth.values() if v=='artefact')} "
      f"genuine={sum(1 for v in truth.values() if v=='genuine')}\n")

for path in report_paths:
    pred = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rid = row["read_id"].split()[0]
            if rid not in truth:
                continue
            d = row["decision"]
            pred[rid] = "genuine" if d == "real" else ("artefact" if d == "artefact" else "other")

    tp = fp = tn = fn = other = 0
    for rid, t in truth.items():
        p = pred.get(rid, "missing")
        if p not in ("artefact", "genuine"):
            other += 1
            continue
        if t == "artefact" and p == "artefact":
            tp += 1
        elif t == "genuine" and p == "artefact":
            fp += 1
        elif t == "genuine" and p == "genuine":
            tn += 1
        elif t == "artefact" and p == "genuine":
            fn += 1

    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    rec = tp / (tp + fn) if (tp + fn) else float("nan")
    fpr = fp / (fp + tn) if (fp + tn) else float("nan")
    print(f"=== {path.split('/')[-1]} ===")
    print(f"  TP={tp} FP={fp} TN={tn} FN={fn} other(no-sig/missing)={other}")
    print(f"  precision={prec:.3f}  recall={rec:.3f}  FPR(wrongly-kept-genuine-cut-as-real)={fpr:.3f}")
    print(f"  Of {tp+fn} true artefacts: recall (correctly cut) = {rec*100:.1f}%")
    print(f"  Of {fp+tn} true genuine reads: {fp} wrongly cut (FPR={fpr*100:.1f}%)\n")
