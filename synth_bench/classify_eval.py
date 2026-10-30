#!/usr/bin/env python3
"""Per-read classification scoring against synth_bench ground truth.

Complements evaluate.py (which scores downstream *assemblies*). This scores
each arm's *per-read call* directly against reads_truth.tsv -- the thing
mdax's report.tsv and CADECT's classification_table.txt both produce before
any assembly is built.

True class:
  type == "chimera" -> artefact
  type == "clean"   -> genuine

The hard case is genuine reads with fully_spans_locus == 1: a clean read
that fully contains a real inverted-repeat locus. This is the synthetic
analogue of the "real_*" reads in data/chris_artefacts_truth.tsv, now at
scale and with known identity/gap/cov_tier per locus, so false-positive
rate can be stratified by cov_tier (the specific mechanism synth_bench's
README hypothesizes: low-coverage real IRs get cut because their support
count never reaches min_support).

Usage:
  classify_eval.py --reads-truth reads_truth.tsv --repeats-bed repeats_truth.bed \
      --mdax-report mdax_balanced=mdax_balanced.report.tsv \
      --mdax-report mdax_strict=mdax_strict.report.tsv \
      --cadect-table cadect/classification_table.txt \
      --outdir eval_out
"""
import argparse
import csv
import os


def load_truth(path):
    truth = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            true_class = "artefact" if row["type"] == "chimera" else "genuine"
            truth[row["read_id"]] = {
                "true_class": true_class,
                "chimera_class": row["chimera_class"],
                "fully_spans_locus": row["fully_spans_locus"] == "1",
                "overlap_loci": [x for x in row["overlap_loci"].split(",") if x and x != "NA"],
            }
    return truth


def load_bed_cov_tier(path):
    cov = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            cov[parts[0]] = float(parts[7])
    return cov


def load_mdax_pred(path):
    pred = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            read_id = row["read_id"].split()[0]
            d = row["decision"]
            pred[read_id] = "genuine" if d == "real" else ("artefact" if d == "artefact" else "other")
    return pred


def load_cadect_pred(path):
    pred = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cls = row["Classification"]
            pred[row["Read ID"]] = {
                "putative_concatemers": "artefact",
                "non_concatemers": "genuine",
                "short": "other",
            }[cls]
    return pred


def confusion(truth, pred, read_ids):
    tp = fp = tn = fn = other = 0
    for rid in read_ids:
        t = truth[rid]["true_class"]
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
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else float("nan")
    fpr = fp / (fp + tn) if (fp + tn) else float("nan")
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, other=other, precision=prec, recall=rec, f1=f1, fpr=fpr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads-truth", required=True)
    ap.add_argument("--repeats-bed", required=True)
    ap.add_argument("--mdax-report", action="append", default=[], help="NAME=PATH, repeatable")
    ap.add_argument("--cadect-table", help="CADECT classification_table.txt")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    truth = load_truth(args.reads_truth)
    cov_tier = load_bed_cov_tier(args.repeats_bed)
    all_ids = list(truth)
    chimera_ids = [r for r in all_ids if truth[r]["true_class"] == "artefact"]
    genuine_ids = [r for r in all_ids if truth[r]["true_class"] == "genuine"]
    hard_ids = [r for r in genuine_ids if truth[r]["fully_spans_locus"]]
    easy_genuine_ids = [r for r in genuine_ids if not truth[r]["fully_spans_locus"]]

    arms = []
    for spec in args.mdax_report:
        name, path = spec.split("=", 1)
        arms.append((name, load_mdax_pred(path)))
    if args.cadect_table:
        arms.append(("CADECT", load_cadect_pred(args.cadect_table)))

    summary_lines = []
    summary_rows = []

    def report(title, ids, pred_by_arm):
        summary_lines.append(f"\n=== {title} (n={len(ids)}) ===")
        for name, pred in pred_by_arm:
            c = confusion(truth, pred, ids)
            summary_rows.append(dict(subset=title, arm=name, **c))
            summary_lines.append(
                f"  {name:<20} TP={c['tp']:>5} FP={c['fp']:>5} TN={c['tn']:>5} FN={c['fn']:>5} "
                f"other={c['other']:>4}  precision={c['precision']:.3f} recall={c['recall']:.3f} "
                f"F1={c['f1']:.3f} FPR={c['fpr']:.3f}"
            )

    report("OVERALL", all_ids, arms)
    report("genuine, fully-spans-real-IR-locus (hard case)", hard_ids, arms)
    report("genuine, no real-IR-locus overlap (easy case)", easy_genuine_ids, arms)

    for cls in sorted(set(truth[r]["chimera_class"] for r in chimera_ids)):
        ids = [r for r in chimera_ids if truth[r]["chimera_class"] == cls]
        report(f"chimera_class={cls}", ids, arms)

    # cov_tier stratification of the hard case (clean reads fully spanning a
    # real IR locus) -- the direct test of the low-coverage mechanism.
    def locus_cov(rid):
        loci = truth[rid]["overlap_loci"]
        tiers = [cov_tier[l] for l in loci if l in cov_tier]
        return min(tiers) if tiers else None

    tier_buckets = [("cov_tier<0.3 (low)", lambda t: t is not None and t < 0.3),
                    ("0.3<=cov_tier<1.0 (mid)", lambda t: t is not None and 0.3 <= t < 1.0),
                    ("cov_tier>=1.0 (high)", lambda t: t is not None and t >= 1.0)]
    for label, pred_fn in tier_buckets:
        ids = [r for r in hard_ids if pred_fn(locus_cov(r))]
        if ids:
            report(f"hard case, {label}", ids, arms)

    print("\n".join(summary_lines))

    tsv_path = os.path.join(args.outdir, "classify_summary.tsv")
    cols = ["subset", "arm", "tp", "fp", "tn", "fn", "other", "precision", "recall", "f1", "fpr"]
    with open(tsv_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in summary_rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(f"\n[classify_eval] wrote {tsv_path}")


if __name__ == "__main__":
    main()
