#!/usr/bin/env python3
"""Score mdax and CADECT against the 57-read chris_artefacts ground truth."""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parent.parent

truth = {}  # label -> true_class ("artefact" | "genuine")
accession = {}  # label -> accession (e.g. ERR12263839.4269)
with open(REPO / "data/chris_artefacts_truth.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        label = row["label"]
        truth[label] = "genuine" if label.startswith("real_") else "artefact"
        accession[label] = row["read_id"].split()[0]

mdax_pred = {}  # label -> "artefact" | "genuine" | "low_ident"
with open(HERE / "mdax_report.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        label = row["read_id"].split()[0]
        d = row["decision"]
        mdax_pred[label] = "genuine" if d == "real" else d  # artefact | low_ident

cadect_by_acc = {}  # accession -> "artefact" | "genuine" | "short"
with open(HERE / "cadect_out/classification_table.txt") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["Read ID"]
        cls = row["Classification"]
        cadect_by_acc[acc] = {
            "putative_concatemers": "artefact",
            "non_concatemers": "genuine",
            "short": "short",
        }[cls]

labels = sorted(truth)
print(f"Total ground-truth reads: {len(labels)}  "
      f"(artefact={sum(1 for l in labels if truth[l]=='artefact')}, "
      f"genuine={sum(1 for l in labels if truth[l]=='genuine')})\n")


def score(name, pred_fn):
    tp = fp = tn = fn = other = 0
    rows = []
    for label in labels:
        t = truth[label]
        p = pred_fn(label)
        rows.append((label, t, p))
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
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else float("nan")
    print(f"=== {name} ===")
    print(f"  TP={tp} FP={fp} TN={tn} FN={fn}  (uncalled/other={other})")
    print(f"  precision={prec:.3f}  recall={rec:.3f}  F1={f1:.3f}")
    print(f"  false-positive rate on genuine reads: {fp}/{fp+tn} = "
          f"{fp/(fp+tn)*100 if (fp+tn) else float('nan'):.1f}%")
    print(f"  false-negative rate on artefact reads: {fn}/{tp+fn} = "
          f"{fn/(tp+fn)*100 if (tp+fn) else float('nan'):.1f}%\n")
    return rows


mdax_rows = score("mdax", lambda label: mdax_pred.get(label, "missing"))
cadect_rows = score(
    "CADECT",
    lambda label: cadect_by_acc.get(accession[label], "missing"),
)

print("=== Per-read detail (label, truth, mdax, CADECT) ===")
for label in labels:
    m = mdax_pred.get(label, "missing")
    c = cadect_by_acc.get(accession[label], "missing")
    flag = "" if (m == truth[label] and c == truth[label]) else "  <-- disagreement with truth"
    print(f"  {label:20s} truth={truth[label]:8s} mdax={m:10s} cadect={c:10s}{flag}")

print("\n=== Genuine reads only (the hard case) ===")
for label in labels:
    if truth[label] != "genuine":
        continue
    m = mdax_pred.get(label, "missing")
    c = cadect_by_acc.get(accession[label], "missing")
    print(f"  {label:20s} mdax={m:10s} cadect={c:10s}")
