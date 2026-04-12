"""
report.py — summarise and compare mdax output across C. elegans MDA samples.

Usage (called by benchmark.sh, but also runnable standalone):
  python3 report.py \\
      --samples nrCaeEleg92 nrCaeEleg95 \\
      --labels  "nrCaeEleg92 (82.3% BUSCO)" "nrCaeEleg95 (64.5% BUSCO)" \\
      --reports out/nrCaeEleg92.report.tsv out/nrCaeEleg95.report.tsv \\
      [--max-reads 200000]

The report TSV has columns (0-indexed):
  0  read_id
  1  len
  2  event
  3  called
  4  coarse_split
  5  refined_split
  6  delta
  7  matches
  8  span_p1
  9  p2_span
  10 cross_frac
  11 coarse_score
  12 refined_score
  13 identity_est
  14 support_n
  15 support_rank_frac
  16 support_span
  17 decision  ← col index may vary; we look by header name
"""
import argparse, csv, os, statistics, sys
from collections import Counter, defaultdict

# ── helpers ───────────────────────────────────────────────────────────────────

def div(a, b):
    return a / b if b else float("nan")

def fmt_pct(n, d):
    if not d:
        return "  n/a"
    return f"{100*n/d:5.1f}%"

def rule(char="─", width=68):
    return char * width

def median_len(lengths):
    if not lengths:
        return "n/a"
    lengths = sorted(lengths)
    return f"{statistics.median(lengths):,.0f}"

def n50(lengths):
    if not lengths:
        return "n/a"
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    acc = 0
    for l in lengths:
        acc += l
        if acc >= half:
            return f"{l:,}"
    return "n/a"

# ── load one report TSV ────────────────────────────────────────────────────────

def load_report(path):
    """Return list of row dicts from an mdax report TSV."""
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows

# ── per-sample stats ──────────────────────────────────────────────────────────

def sample_stats(rows, total_reads):
    """Compute aggregate stats for one sample."""
    dec_counts = Counter(r["decision"] for r in rows)
    all_detected = len(rows)

    n_artefact  = dec_counts.get("artefact", 0)
    n_low_ident = dec_counts.get("low_ident", 0)
    n_real      = dec_counts.get("real", 0)

    # Read lengths — all reads in report (detected foldback)
    detected_lens = []
    for r in rows:
        try:
            detected_lens.append(int(r["len"]))
        except (ValueError, KeyError):
            pass

    # After-cut lengths: for artefact reads, the output length is refined_split
    cut_lens_before = []
    cut_lens_after  = []
    ident_vals = []
    support_vals = []

    for r in rows:
        if r["decision"] == "artefact":
            try:
                cut_lens_before.append(int(r["len"]))
                split = int(r["refined_split"]) if r.get("refined_split") else int(r["coarse_split"])
                cut_lens_after.append(split)
            except (ValueError, KeyError):
                pass
        if r["decision"] in ("artefact", "low_ident", "real"):
            try:
                ident_vals.append(float(r["identity_est"]))
            except (ValueError, KeyError):
                pass
        try:
            sn = int(r["support_n"])
            if sn > 0:
                support_vals.append(sn)
        except (ValueError, KeyError):
            pass

    return {
        "total":         total_reads,
        "detected":      all_detected,
        "artefact":      n_artefact,
        "low_ident":     n_low_ident,
        "real":          n_real,
        "dec_counts":    dec_counts,
        "detected_lens": detected_lens,
        "cut_before":    cut_lens_before,
        "cut_after":     cut_lens_after,
        "ident_vals":    ident_vals,
        "support_vals":  support_vals,
    }


def print_sample(label, stats):
    s = stats
    total   = s["total"]
    det     = s["detected"]
    art     = s["artefact"]
    low     = s["low_ident"]
    real    = s["real"]
    not_det = total - det

    print(f"\n  {label}")
    print(rule())
    print(f"  Reads processed               {total:>10,}")
    print(f"  Foldback detected             {det:>10,}  {fmt_pct(det, total)} of total")
    print(f"    → artefact  (cut)           {art:>10,}  {fmt_pct(art, total)} of total  |  {fmt_pct(art, det)} of detected")
    print(f"    → low_ident (below thresh)  {low:>10,}  {fmt_pct(low, total)} of total  |  {fmt_pct(low, det)} of detected")
    print(f"    → real      (palindrome)    {real:>10,}  {fmt_pct(real, total)} of total  |  {fmt_pct(real, det)} of detected")
    print(f"  Not detected                  {not_det:>10,}  {fmt_pct(not_det, total)} of total")
    print()

    if s["cut_before"]:
        bases_before = sum(s["cut_before"])
        bases_after  = sum(s["cut_after"])
        print(f"  Bases in artefact reads (before cut)  {bases_before:>14,} bp")
        print(f"  Bases retained after cut              {bases_after:>14,} bp  ({fmt_pct(bases_after, bases_before)} retained)")
        print(f"  Bases removed                         {bases_before - bases_after:>14,} bp")
        print(f"  Artefact read len  median={median_len(s['cut_before'])} bp,  N50={n50(s['cut_before'])}")
        print(f"  After-cut len      median={median_len(s['cut_after'])} bp")
        print()

    if s["ident_vals"]:
        ids = sorted(s["ident_vals"])
        n   = len(ids)
        print(f"  Identity of detected foldbacks  (n={n})")
        print(f"    median {statistics.median(ids):.3f}   mean {statistics.mean(ids):.3f}"
              f"   q10 {ids[int(0.10*(n-1))]:.3f}   q90 {ids[int(0.90*(n-1))]:.3f}")
        print()

        # identity histogram
        bins = [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01]
        labels = ["<0.50", "0.50-0.59", "0.60-0.69", "0.70-0.79",
                  "0.80-0.89", "0.90-0.94", "≥0.95"]
        cnts = [0] * len(labels)
        for v in ids:
            for i in range(len(bins) - 1):
                if bins[i] <= v < bins[i+1]:
                    cnts[i] += 1
                    break
        print(f"  Identity distribution:")
        for lbl, c in zip(labels, cnts):
            bar = "█" * (c * 25 // n) if n else ""
            print(f"    {lbl:10s}  {c:>6,}  {fmt_pct(c,n):>6}  {bar}")
        print()

    if s["support_vals"]:
        sv = s["support_vals"]
        print(f"  Support_n for confirmed junctions (n={len(sv)})")
        cnt = Counter(sv)
        for k in sorted(cnt)[:10]:
            print(f"    support={k:>3}: {cnt[k]:>5,} reads")
        if max(cnt) > 10:
            print(f"    …")


# ── comparison table ──────────────────────────────────────────────────────────

def print_comparison(labels, stats_list):
    print()
    print(rule("═"))
    print("  SAMPLE COMPARISON")
    print(rule("═"))
    print()

    metrics = [
        ("Total reads",        lambda s: f"{s['total']:>10,}"),
        ("Foldback detected",  lambda s: f"{s['detected']:>10,}  ({fmt_pct(s['detected'], s['total'])})"),
        ("  → artefact",       lambda s: f"{s['artefact']:>10,}  ({fmt_pct(s['artefact'], s['total'])})"),
        ("  → low_ident",      lambda s: f"{s['low_ident']:>10,}  ({fmt_pct(s['low_ident'], s['total'])})"),
        ("  → real",           lambda s: f"{s['real']:>10,}  ({fmt_pct(s['real'], s['total'])})"),
    ]

    col_w = 38
    header = f"  {'metric':28s}" + "".join(f"  {lbl[:col_w]:<{col_w}}" for lbl in labels)
    print(header)
    print("  " + "─" * (28 + len(labels) * (col_w + 2)))

    for (name, fn) in metrics:
        row = f"  {name:28s}"
        for s in stats_list:
            if s is not None:
                row += f"  {fn(s):<{col_w}}"
            else:
                row += f"  {'(no data)':<{col_w}}"
        print(row)

    # ratio row for artefact rate
    if len(stats_list) >= 2 and all(s is not None for s in stats_list):
        rates = [div(s["artefact"], s["total"]) for s in stats_list]
        det_rates = [div(s["detected"], s["total"]) for s in stats_list]
        low_rates  = [div(s["low_ident"], s["total"]) for s in stats_list]
        if rates[0] and rates[0] != float("nan"):
            ratio = div(rates[1], rates[0])
            if ratio == ratio:  # not nan
                print()
                print(f"  Artefact-rate ratio  [{labels[1]} / {labels[0]}] = {ratio:.2f}x")
                if ratio > 1.0:
                    print(f"  → {labels[1]} has {ratio:.2f}x MORE artefact reads, consistent with lower assembly quality.")
                else:
                    print(f"  → {labels[1]} has {ratio:.2f}x fewer artefact reads.")
                    print(f"     Note: if reads were read from the start of the file, PacBio CCS files are")
                    print(f"     sometimes length-sorted (shortest first) which can bias short-read subsets.")
                    print(f"     Consider MAX_READS=0 (full dataset) or a random subsample for a fairer comparison.")

    print()


# ── main ──────────────────────────────────────────────────────────────────────

def write_summary_tsv(path, labels, stats_list, max_reads):
    """Write a machine-readable per-sample summary TSV."""
    fields = [
        "sample", "max_reads", "total", "detected", "detected_frac",
        "artefact", "artefact_frac",
        "low_ident", "low_ident_frac",
        "real", "real_frac",
        "bases_before_bp", "bases_after_bp", "bases_removed_bp", "frac_retained",
        "identity_median", "identity_mean", "identity_q10", "identity_q90",
    ]
    rows = []
    for label, s in zip(labels, stats_list):
        if s is None:
            continue
        total = s["total"]
        det   = s["detected"]
        art   = s["artefact"]
        low   = s["low_ident"]
        real  = s["real"]

        bb = sum(s["cut_before"]) if s["cut_before"] else ""
        ba = sum(s["cut_after"])  if s["cut_after"]  else ""
        br = (bb - ba) if (bb != "" and ba != "") else ""
        fr = f"{ba/bb:.4f}" if bb else ""

        ids = sorted(s["ident_vals"])
        n   = len(ids)
        id_med  = f"{statistics.median(ids):.4f}" if ids else ""
        id_mean = f"{statistics.mean(ids):.4f}"   if ids else ""
        id_q10  = f"{ids[int(0.10*(n-1))]:.4f}"  if ids else ""
        id_q90  = f"{ids[int(0.90*(n-1))]:.4f}"  if ids else ""

        rows.append([
            label, str(max_reads), str(total), str(det), f"{div(det,total):.4f}",
            str(art), f"{div(art,total):.4f}",
            str(low), f"{div(low,total):.4f}",
            str(real), f"{div(real,total):.4f}",
            str(bb), str(ba), str(br), fr,
            id_med, id_mean, id_q10, id_q90,
        ])

    with open(path, "w") as f:
        f.write("\t".join(fields) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples",    nargs="+", required=True)
    ap.add_argument("--labels",     nargs="+", required=True)
    ap.add_argument("--reports",    nargs="+", required=True)
    ap.add_argument("--max-reads",  type=int,  default=0,
                    help="MAX_READS cap used (for display only; 0=all)")
    ap.add_argument("--out-tsv",    default=None,
                    help="Write machine-readable per-sample summary to this path")
    args = ap.parse_args()

    print()
    print(rule("═", 68))
    print("  MDAX BENCHMARK  —  C. elegans single-worm MDA (PacBio CCS)")
    print(rule("═", 68))
    cap = f"{args.max_reads:,}" if args.max_reads else "all"
    print(f"  Reads per sample: {cap}   |   samples: {len(args.samples)}")

    all_stats = []
    for sample, label, rpt_path in zip(args.samples, args.labels, args.reports):
        rows = load_report(rpt_path)
        if rows is None:
            print(f"\n  [SKIP] {rpt_path} not found — did benchmark.sh run successfully?")
            all_stats.append(None)
            continue

        # Total reads: if we capped at MAX_READS we know the total; otherwise
        # derive it from the fasta subset that sits next to the report.
        total = args.max_reads if args.max_reads else None
        if total is None:
            # Try the subset fasta
            subset_fa = rpt_path.replace(".report.tsv", ".subset.fasta").replace("out/", "out/")
            if os.path.exists(subset_fa):
                total = sum(1 for line in open(subset_fa) if line.startswith(">"))
            else:
                # Fall back: detected reads + undetected estimated as 0
                total = len(rows)
                print(f"  [warn] Could not find subset FASTA for {sample}; "
                      f"total_reads = detected reads only ({total:,})", file=sys.stderr)

        stats = sample_stats(rows, total)
        all_stats.append(stats)

        print_sample(f"{sample}  —  {label}", stats)

    print_comparison(args.labels, all_stats)

    if args.out_tsv:
        write_summary_tsv(args.out_tsv, args.labels, all_stats, args.max_reads)
        print(f"\n  Summary TSV written → {args.out_tsv}")


if __name__ == "__main__":
    main()
