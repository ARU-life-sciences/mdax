"""
test_pal.py — comparison of mdax calls against the Strobl et al. (2023) Perl baseline.

Reads (relative to the test_pal directory):
  mdax.calls.norm.tsv          — mdax output, normalised by test_pal.sh
  baseline.flagged.canon.ids   — baseline positive IDs (canonicalised)
  baseline.bp.tsv              — baseline breakpoints (id, bp)

Usage:
  python3 test_pal.py [--total-reads N] [--format {human,tsv}] [--tag LABEL]

  --total-reads N     Override the total read count (TN/specificity).
                      Inferred from out/subsample/ if omitted.
  --format tsv        Emit a single TSV data row (no header) for table assembly.
                      Column order: tag strict_P strict_R strict_F1 strict_J
                        generous_P generous_R generous_F1 generous_J
                        detected artefact low_ident real not_detected_in_base
  --tag LABEL         Label prepended to TSV row (useful in sweep scripts).
"""
import argparse, csv, glob, gzip, os, statistics, sys
from collections import Counter

# ── helpers ──────────────────────────────────────────────────────────────────

def div(a, b):
    return a / b if b else float("nan")


def fmt_pct(n, d):
    if not d:
        return "  n/a"
    return f"{100*n/d:5.1f}%"


def rule(char="─", width=60):
    return char * width


def infer_total_reads():
    """Try to count reads from the first FASTA/FASTQ found in out/subsample/."""
    patterns = [
        "out/subsample/*.fasta",
        "out/subsample/*.fa",
        "out/subsample/*.fastq.gz",
        "out/subsample/*.fq.gz",
    ]
    for pat in patterns:
        hits = glob.glob(pat)
        if hits:
            path = hits[0]
            n = 0
            opener = gzip.open if path.endswith(".gz") else open
            mode = "rt"
            with opener(path, mode) as f:
                for line in f:
                    if line.startswith(">") or line.startswith("@"):
                        n += 1
            return n
    return None


def confusion_metrics(tp, fp, fn, tn):
    prec = div(tp, tp + fp)
    rec  = div(tp, tp + fn)
    f1   = div(2 * tp, 2 * tp + fp + fn)
    jac  = div(tp, tp + fp + fn)
    spec = div(tn, tn + fp)
    return prec, rec, f1, jac, spec


def print_confusion(label, tp, fp, fn, tn, total):
    prec, rec, f1, jac, spec = confusion_metrics(tp, fp, fn, tn)
    print(f"\n{rule()}")
    print(f"  {label}")
    print(rule())
    print(f"  {'':20s}  {'baseline POS':>12}  {'baseline NEG':>12}")
    print(f"  {'mdax POS':20s}  {tp:>12,}  {fp:>12,}   = {tp+fp:,} called")
    print(f"  {'mdax NEG':20s}  {fn:>12,}  {tn:>12,}")
    print(f"  {'':20s}  {'─'*12}  {'─'*12}")
    print(f"  {'total':20s}  {tp+fn:>12,}  {fp+tn:>12,}   = {total:,}")
    print()
    print(f"  Precision  {prec:.3f}   (of mdax positives, fraction correct)")
    print(f"  Recall     {rec:.3f}   (of baseline positives, fraction found)")
    print(f"  F1         {f1:.3f}")
    print(f"  Jaccard    {jac:.3f}")
    print(f"  Specificity {spec:.3f}  (of baseline negatives, fraction left alone)")


# ── parse args ────────────────────────────────────────────────────────────────

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--total-reads", type=int, default=None,
                help="Override total read count (used for TN/specificity)")
ap.add_argument("--format", choices=["human", "tsv"], default="human",
                help="Output format: human-readable (default) or single TSV row")
ap.add_argument("--tag", default="",
                help="Label prepended to TSV row")
args = ap.parse_args()

# ── load data ─────────────────────────────────────────────────────────────────

mdax = {}          # read_id → {called, split, decision}
dec_counts = Counter()

with open("mdax.calls.norm.tsv") as f:
    for r in csv.reader(f, delimiter="\t"):
        mdax[r[0]] = {
            "called":   int(r[1]),
            "split":    int(r[2]) if r[2] else None,
            "decision": r[3] if len(r) > 3 else "",
        }
        dec_counts[mdax[r[0]]["decision"]] += 1

baseline_called = set()
with open("baseline.flagged.canon.ids") as f:
    for line in f:
        rid = line.strip()
        if rid:
            baseline_called.add(rid)

baseline_bp = {}
with open("baseline.bp.tsv") as f:
    for r in csv.reader(f, delimiter="\t"):
        baseline_bp[r[0]] = int(r[1])

# ── infer / validate total reads ─────────────────────────────────────────────

total_reads = args.total_reads
if total_reads is None:
    total_reads = infer_total_reads()
if total_reads is None:
    # fall back: mdax-detected + undetected approximation from mdax report
    # We can't know TN without total, so use the union of known IDs as a floor
    total_reads = len(set(mdax) | baseline_called)
    print(f"[warn] Could not find input FASTA/FASTQ; total_reads approximated as {total_reads:,}.",
          file=sys.stderr)

# ── summary ───────────────────────────────────────────────────────────────────

n_detected  = len(mdax)
n_artefact  = sum(1 for v in mdax.values() if v["decision"] == "artefact")
n_low_ident = sum(1 for v in mdax.values() if v["decision"] == "low_ident")
n_real      = sum(1 for v in mdax.values() if v["decision"] == "real")
n_other     = n_detected - n_artefact - n_low_ident - n_real
n_baseline  = len(baseline_called)

print(rule("═"))
print("  MDAX vs STROBL-2023 BASELINE  —  CALL AGREEMENT REPORT")
print(rule("═"))
print()
print(f"  Total reads in subsample          {total_reads:>10,}")
print(f"  Baseline flagged reads            {n_baseline:>10,}  ({fmt_pct(n_baseline, total_reads)} of total)")
print()
print(f"  mdax foldback detected            {n_detected:>10,}  ({fmt_pct(n_detected, total_reads)} of total)")
print(f"    → artefact (cut)                {n_artefact:>10,}  ({fmt_pct(n_artefact, total_reads)} of total)")
print(f"    → low_ident (no cut)            {n_low_ident:>10,}  ({fmt_pct(n_low_ident, total_reads)} of total)")
print(f"    → real (palindrome, no cut)     {n_real:>10,}  ({fmt_pct(n_real, total_reads)} of total)")
if n_other:
    print(f"    → other                         {n_other:>10,}")

# ── per-decision overlap with baseline ────────────────────────────────────────

print()
print(rule("─"))
print("  DECISION BUCKET BREAKDOWN  (overlap with baseline)")
print(rule("─"))
print(f"  {'decision':15s}  {'count':>7}  {'in baseline':>11}  {'not in base':>11}  "
      f"{'%pos in bkt':>11}  {'% of baseline':>13}")
print(f"  {'─'*15}  {'─'*7}  {'─'*11}  {'─'*11}  {'─'*11}  {'─'*13}")

decisions_ordered = ["artefact", "low_ident", "real"] + \
    sorted(k for k in dec_counts if k not in ("artefact", "low_ident", "real"))

for dec in decisions_ordered:
    if dec_counts[dec] == 0:
        continue
    ids_in_dec = {rid for rid, v in mdax.items() if v["decision"] == dec}
    in_base = len(ids_in_dec & baseline_called)
    not_in  = len(ids_in_dec) - in_base
    pct_of_bucket = fmt_pct(in_base, dec_counts[dec])  # fraction of this bucket that is baseline-pos
    pct_of_base   = fmt_pct(in_base, n_baseline)       # fraction of all baseline reads in this bucket
    print(f"  {dec:15s}  {dec_counts[dec]:>7,}  {in_base:>11,}  {not_in:>11,}  "
          f"{pct_of_bucket:>11}  {pct_of_base:>13}")

# reads with no foldback detected
undetected_ids = (set(mdax) | baseline_called) - set(mdax)
undetected_in_base = len(undetected_ids & baseline_called)
not_detected_count = total_reads - n_detected
not_detected_in_base = n_baseline - sum(
    len({rid for rid, v in mdax.items() if v["decision"] == dec} & baseline_called)
    for dec in decisions_ordered if dec_counts[dec] > 0
)
print(f"  {'(not detected)':15s}  {not_detected_count:>7,}  {not_detected_in_base:>11,}  "
      f"{not_detected_count - not_detected_in_base:>11,}  "
      f"{fmt_pct(not_detected_in_base, not_detected_count):>11}  "
      f"{fmt_pct(not_detected_in_base, n_baseline):>13}")

# ── confusion matrices ─────────────────────────────────────────────────────────

# Strict: only artefact = called
def compute_cm(called_decisions):
    tp = fp = fn = tn = 0
    called_ids = {rid for rid, v in mdax.items() if v["decision"] in called_decisions}
    # Iterate over all reads that appear in either set; reads absent from both
    # are pure TN and are added below using total_reads.
    known_ids = set(mdax) | baseline_called
    for rid in known_ids:
        m = 1 if rid in called_ids else 0
        b = 1 if rid in baseline_called else 0
        if m and b:     tp += 1
        elif m:         fp += 1
        elif b:         fn += 1
        else:           tn += 1
    # Reads not seen in either mdax output nor baseline are unambiguous TN.
    tn += total_reads - len(known_ids)
    return tp, fp, fn, tn

tp_s, fp_s, fn_s, tn_s = compute_cm({"artefact"})
tp_g, fp_g, fn_g, tn_g = compute_cm({"artefact", "low_ident"})

# ── TSV output (one data row, no header) ──────────────────────────────────────

if args.format == "tsv":
    ps, rs, f1s, js, _ = confusion_metrics(tp_s, fp_s, fn_s, tn_s)
    pg, rg, f1g, jg, _ = confusion_metrics(tp_g, fp_g, fn_g, tn_g)
    row = "\t".join([
        args.tag,
        f"{ps:.4f}", f"{rs:.4f}", f"{f1s:.4f}", f"{js:.4f}",
        f"{pg:.4f}", f"{rg:.4f}", f"{f1g:.4f}", f"{jg:.4f}",
        str(n_detected), str(n_artefact), str(n_low_ident), str(n_real),
        str(not_detected_in_base),
    ])
    print(row)
    sys.exit(0)

# ── human-readable output ─────────────────────────────────────────────────────

print_confusion("STRICT CALL  (decision == artefact only)",
                tp_s, fp_s, fn_s, tn_s, total_reads)

print_confusion("GENEROUS CALL  (artefact OR low_ident counted as positive)",
                tp_g, fp_g, fn_g, tn_g, total_reads)

# ── breakpoint deltas ─────────────────────────────────────────────────────────

abs_deltas = []
signed_deltas = []

for rid in baseline_called:
    if rid in mdax and mdax[rid]["called"] == 1:
        if mdax[rid]["split"] is not None and rid in baseline_bp:
            d = mdax[rid]["split"] - baseline_bp[rid]
            signed_deltas.append(d)
            abs_deltas.append(abs(d))

print()
print(rule())
print("  BREAKPOINT ACCURACY  (strict artefact calls vs baseline iter-2 bp)")
print(rule())

if abs_deltas:
    abs_deltas.sort()
    signed_deltas.sort()
    n = len(abs_deltas)

    # Bucket distribution
    buckets = [(0, 0), (1, 50), (51, 100), (101, 250), (251, 500), (501, 1000), (1001, 10**9)]
    bucket_counts = []
    for lo, hi in buckets:
        c = sum(1 for d in abs_deltas if lo <= d <= hi)
        bucket_counts.append(c)

    print(f"\n  n = {n}  (TP calls with both mdax and baseline breakpoints available)")
    print()
    print(f"  {'range (bp)':>12}  {'count':>6}  {'%':>6}  bar")
    for (lo, hi), c in zip(buckets, bucket_counts):
        bar = "█" * (c * 30 // n) if n else ""
        rng = f"≤{hi}" if lo == 0 else (f">{lo}" if hi >= 10**9 else f"{lo}–{hi}")
        print(f"  {rng:>12}  {c:>6}  {fmt_pct(c,n):>6}  {bar}")
    print()
    print(f"  abs min     {abs_deltas[0]:>6} bp")
    print(f"  abs median  {statistics.median(abs_deltas):>6.0f} bp")
    print(f"  abs mean    {statistics.mean(abs_deltas):>6.1f} bp")
    print(f"  abs max     {abs_deltas[-1]:>6} bp")
    print(f"  abs q90     {abs_deltas[int(0.90*(n-1))]:>6} bp")
    print(f"  abs q95     {abs_deltas[int(0.95*(n-1))]:>6} bp")
    print()
    print(f"  signed median  {statistics.median(signed_deltas):>+7.0f} bp  (positive = mdax downstream of baseline)")
    print(f"  signed mean    {statistics.mean(signed_deltas):>+7.1f} bp")
else:
    print("\n  No overlapping breakpoint calls.")

print()
print(rule("═"))
print()
