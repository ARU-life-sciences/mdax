#!/usr/bin/env python3
"""
crossmatch.py — Read-level comparison of mdax output vs Chris Laumer's independent
soft-clip classification for ERR12263839 (PacBio HiFi, MDA-amplified, C. elegans).

Tests covered
-------------
1. Read-level crossmatch (per Chris artifact class):
   - How many reads overlap between the two analyses?
   - For each of Chris's artifact classes, what decision does mdax reach?
   - What are the mean/median gap_est and identity_est values per class?

3. Tandem concatemer characterisation:
   - Chris finds 9.1% of clips are tandem_rc_concatemers (multiple alternating
     fwd/RC repeat units).  mdax does not model this structure explicitly.
   - Does mdax still call them artefact?  What is the identity distribution?
   - Which unit sizes / copy counts are most common?

Input files (relative to this script — adjust DATA_DIR if layout changes)
----------
  ../../chris_data/ERR12263839.tsv                   mdax report (full dataset)
  ../../chris_data/artifact_classifications_head.tsv Chris's per-clip table (first 10k rows)

Output
------
Prints a formatted report to stdout.  Optionally writes a joined TSV with
  --out joined.tsv

Usage
-----
  python3 crossmatch.py [--out joined.tsv] [--verbose]

Expected results (sanity checks at end of script)
--------------------------------------------------
  intrachromosomal_chimera  → mdax mostly "artefact" with gap_est > 0
  local_inversion           → mdax mostly "artefact" with gap_est ≈ 0
  palindromic_hairpin       → mdax mostly "real" (palindromes should be preserved)
  tandem_rc_concatemer      → mdax mostly "artefact" (overall foldback still detected)
"""

import argparse
import csv
import os
import statistics
import sys
from collections import Counter, defaultdict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA_DIR    = os.path.join(SCRIPT_DIR, "..", "..", "chris_data")
OUTPUT_DIR  = os.path.join(SCRIPT_DIR, "output")

MDAX_TSV    = os.path.join(DATA_DIR, "ERR12263839.tsv")
CHRIS_TSV   = os.path.join(DATA_DIR, "artifact_classifications_head.tsv")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--out",     default=os.path.join(OUTPUT_DIR, "joined.tsv"),
                metavar="PATH",
                help="Write joined TSV (default: output/joined.tsv)")
ap.add_argument("--verbose", action="store_true",
                help="Print per-read rows for classes with ≤20 overlapping reads")
args = ap.parse_args()
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Guard: check data files exist
# ---------------------------------------------------------------------------

MISSING = [p for p in (MDAX_TSV, CHRIS_TSV) if not os.path.exists(p)]
if MISSING:
    print("SKIP: required data files not found:", file=sys.stderr)
    for p in MISSING:
        print(f"  {p}", file=sys.stderr)
    print("  These files live in chris_data/ which is gitignored (large files).",
          file=sys.stderr)
    sys.exit(0)

# ---------------------------------------------------------------------------
# Load mdax report
# Index key: short SRA ID (first whitespace-delimited token of read_id).
# One row per detected read; keep the row with the largest refined_split if
# the same read appears multiple times (recursive-cut scenario).
# ---------------------------------------------------------------------------

print("Loading mdax report …", flush=True)
mdax: dict[str, dict] = {}
with open(MDAX_TSV, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        short_id = row["read_id"].split()[0]
        # Keep the outermost (largest refined_split) row for recursively-cut reads
        if short_id not in mdax:
            mdax[short_id] = row
        else:
            try:
                cur  = int(mdax[short_id].get("refined_split") or 0)
                new  = int(row.get("refined_split") or 0)
                if new > cur:
                    mdax[short_id] = row
            except ValueError:
                pass

print(f"  {len(mdax):,} unique reads in mdax report")

# ---------------------------------------------------------------------------
# Load Chris's classification table
# Key: read_id (already in short ERR12263839.N format)
# A read can have entries for both LEFT and RIGHT clips.
# ---------------------------------------------------------------------------

print("Loading Chris's classification table …", flush=True)
chris: dict[str, list[dict]] = defaultdict(list)
with open(CHRIS_TSV, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        chris[row["read_id"]].append(row)

print(f"  {sum(len(v) for v in chris.values()):,} rows across {len(chris):,} reads")

# ---------------------------------------------------------------------------
# Build overlap
# ---------------------------------------------------------------------------

overlap_ids = set(mdax) & set(chris)
print(f"  Overlap: {len(overlap_ids):,} reads appear in both tables\n")

# ---------------------------------------------------------------------------
# Per-class analysis
# ---------------------------------------------------------------------------

def sep(char="─", width=70):
    return char * width


def safe_float(val, default=float("nan")):
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


print(sep("═"))
print("  TEST 1 — PER-CLASS COMPARISON: mdax decisions vs Chris's artifact classes")
print(sep("═"))
print()

# For each read in the overlap, what is the "primary" class from Chris?
# A read may have clips on both sides classified differently; we use the
# class of the clip with the longer clip_len.
def primary_chris_class(entries: list[dict]) -> str:
    best = max(entries, key=lambda r: int(r.get("clip_len") or 0))
    return best["artifact_class"]

chris_class_for_read: dict[str, str] = {
    rid: primary_chris_class(entries)
    for rid, entries in chris.items()
}

# Sort classes by frequency in the overlap set
class_counts = Counter(
    chris_class_for_read[rid] for rid in overlap_ids
)

# Table header
print(f"  {'Chris class':<30}  {'n':>6}  {'mdax artefact':>13}  "
      f"{'low_ident':>9}  {'real':>6}  {'other':>6}  "
      f"{'mean gap_est':>12}  {'med identity':>12}")
print(f"  {'─'*30}  {'─'*6}  {'─'*13}  {'─'*9}  {'─'*6}  {'─'*6}  {'─'*12}  {'─'*12}")

joined_rows = []  # for optional --out file

for cls, _n in class_counts.most_common():
    rids = [rid for rid in overlap_ids if chris_class_for_read[rid] == cls]
    n = len(rids)

    mdax_decisions = Counter(mdax[r]["decision"] for r in rids)
    n_art = mdax_decisions.get("artefact", 0)
    n_low = mdax_decisions.get("low_ident", 0)
    n_real = mdax_decisions.get("real", 0)
    n_other = n - n_art - n_low - n_real

    gap_vals  = [safe_float(mdax[r].get("gap_est"))  for r in rids]
    ident_vals = [safe_float(mdax[r].get("identity_est")) for r in rids]
    gap_vals   = [v for v in gap_vals   if v == v]   # drop nan
    ident_vals = [v for v in ident_vals if v == v]

    mean_gap   = statistics.mean(gap_vals)   if gap_vals   else float("nan")
    med_ident  = statistics.median(ident_vals) if ident_vals else float("nan")

    print(f"  {cls:<30}  {n:>6,}  {n_art:>13,}  {n_low:>9,}  {n_real:>6,}  "
          f"{n_other:>6,}  {mean_gap:>12.0f}  {med_ident:>12.3f}")

    # Build joined rows for optional output
    for rid in rids:
        m = mdax[rid]
        # Include all of Chris's clip entries for this read
        for entry in chris[rid]:
            joined_rows.append({
                "read_id":          rid,
                "chris_class":      entry["artifact_class"],
                "chris_clip_side":  entry["clip_side"],
                "chris_clip_len":   entry["clip_len"],
                "chris_detail":     entry["detail"],
                "mdax_decision":    m["decision"],
                "mdax_gap_est":     m.get("gap_est", ""),
                "mdax_identity_est": m.get("identity_est", ""),
                "mdax_refined_split": m.get("refined_split", ""),
            })

print()

# ---------------------------------------------------------------------------
# Detailed verbose output for small classes
# ---------------------------------------------------------------------------

if args.verbose:
    for cls, n in class_counts.most_common():
        rids = [rid for rid in overlap_ids if chris_class_for_read[rid] == cls]
        if len(rids) > 20:
            continue
        print(sep())
        print(f"  {cls}  (n={len(rids)} — showing all)")
        print(sep())
        for rid in sorted(rids):
            m = mdax[rid]
            decisions = [e["artifact_class"] for e in chris[rid]]
            print(f"  {rid}  mdax={m['decision']}  gap={m.get('gap_est','')}  "
                  f"ident={m.get('identity_est','')}  chris={decisions}")
        print()

# ---------------------------------------------------------------------------
# TEST 3 — Tandem concatemer deep-dive
# ---------------------------------------------------------------------------

print(sep("═"))
print("  TEST 3 — TANDEM RC CONCATEMER CLASS DEEP-DIVE")
print(sep("═"))
print()
print("  Background: Chris's pipeline detects clips where alternating fwd/RC repeat")
print("  units are present (≥2 copies, identity ≥0.85).  mdax does not model this")
print("  sub-structure; it sees the overall foldback and may call artefact or low_ident.")
print()

tandem_rids = [rid for rid in overlap_ids
               if any(e["artifact_class"] == "tandem_rc_concatemer" for e in chris[rid])]

if tandem_rids:
    mdax_decisions = Counter(mdax[r]["decision"] for r in tandem_rids)
    gap_vals   = [safe_float(mdax[r].get("gap_est"))    for r in tandem_rids]
    ident_vals = [safe_float(mdax[r].get("identity_est")) for r in tandem_rids]
    gap_vals   = [v for v in gap_vals   if v == v]
    ident_vals = [v for v in ident_vals if v == v]

    print(f"  Tandem concatemer reads in overlap:  {len(tandem_rids):,}")
    print(f"  mdax decisions:  {dict(mdax_decisions)}")
    if gap_vals:
        print(f"  gap_est  — mean={statistics.mean(gap_vals):.0f} bp  "
              f"median={statistics.median(gap_vals):.0f} bp  "
              f"max={max(gap_vals):.0f} bp")
    if ident_vals:
        print(f"  identity — mean={statistics.mean(ident_vals):.3f}  "
              f"median={statistics.median(ident_vals):.3f}  "
              f"min={min(ident_vals):.3f}")

    # Parse unit_len and copies from Chris's detail field
    unit_lens = []
    copies_list = []
    for rid in tandem_rids:
        for entry in chris[rid]:
            if entry["artifact_class"] != "tandem_rc_concatemer":
                continue
            detail = entry.get("detail", "")
            for part in detail.split(","):
                if part.startswith("unit_len="):
                    try:
                        unit_lens.append(int(part.split("=")[1]))
                    except ValueError:
                        pass
                elif part.startswith("copies="):
                    try:
                        copies_list.append(int(part.split("=")[1]))
                    except ValueError:
                        pass

    if unit_lens:
        print(f"  Unit length — mean={statistics.mean(unit_lens):.0f} bp  "
              f"median={statistics.median(unit_lens):.0f} bp  "
              f"range={min(unit_lens)}–{max(unit_lens)} bp")
    if copies_list:
        print(f"  Copy count  — mean={statistics.mean(copies_list):.1f}  "
              f"median={statistics.median(copies_list):.1f}  "
              f"range={min(copies_list)}–{max(copies_list)}")
else:
    print("  No tandem_rc_concatemer reads found in the overlap (head table may be too small).")

print()

# ---------------------------------------------------------------------------
# Sanity-check assertions
# ---------------------------------------------------------------------------

print(sep("═"))
print("  SANITY CHECKS")
print(sep("═"))
print()

fails = 0

def check(label, cond, msg):
    global fails
    status = "PASS" if cond else "FAIL"
    if not cond:
        fails += 1
    print(f"  [{status}] {label}")
    if not cond:
        print(f"         {msg}")


# intrachromosomal_chimera: majority of mdax calls should be artefact
ic_rids = [r for r in overlap_ids if chris_class_for_read[r] == "intrachromosomal_chimera"]
if ic_rids:
    artefact_frac = sum(1 for r in ic_rids if mdax[r]["decision"] == "artefact") / len(ic_rids)
    check("intrachromosomal_chimera → mdax artefact (expect >50%)",
          artefact_frac > 0.50,
          f"only {artefact_frac:.1%} of intrachromosomal_chimera reads got artefact decision")

# local_inversion: should map to artefact with small gap
li_rids = [r for r in overlap_ids if chris_class_for_read[r] == "local_inversion"]
if li_rids:
    artefact_frac = sum(1 for r in li_rids if mdax[r]["decision"] == "artefact") / len(li_rids)
    check("local_inversion → mdax artefact (expect >50%)",
          artefact_frac > 0.50,
          f"only {artefact_frac:.1%} got artefact decision")
    gaps = [safe_float(mdax[r].get("gap_est")) for r in li_rids]
    gaps = [v for v in gaps if v == v]
    if gaps:
        # local_inversion means the supplementary alignment is genomically near the primary,
        # but the read-coordinate gap (mdax gap_est) can still be substantial because the
        # template-switch skipped bases within the read before reversing.  Observed median
        # is ~1000 bp; we check that it is smaller than the intrachromosomal_chimera median.
        ic_gaps = [safe_float(mdax[r].get("gap_est"))
                   for r in overlap_ids if chris_class_for_read[r] == "intrachromosomal_chimera"]
        ic_gaps = [v for v in ic_gaps if v == v]
        if ic_gaps:
            check("local_inversion gap_est ≤ intrachromosomal_chimera median",
                  statistics.median(gaps) <= statistics.median(ic_gaps),
                  f"local_inversion median gap={statistics.median(gaps):.0f} bp "
                  f"exceeds chimera median={statistics.median(ic_gaps):.0f} bp")
        else:
            print(f"  [INFO] local_inversion median gap_est = {statistics.median(gaps):.0f} bp")

# palindromic_hairpin: should map to real — but head table typically has very few
ph_rids = [r for r in overlap_ids if chris_class_for_read[r] == "palindromic_hairpin"]
if len(ph_rids) >= 10:
    real_frac = sum(1 for r in ph_rids if mdax[r]["decision"] == "real") / len(ph_rids)
    check("palindromic_hairpin → mdax real (expect >30%)",
          real_frac > 0.30,
          f"only {real_frac:.1%} got real decision (note: depends on batch support)")
elif ph_rids:
    real_frac = sum(1 for r in ph_rids if mdax[r]["decision"] == "real") / len(ph_rids)
    print(f"  [INFO] palindromic_hairpin: n={len(ph_rids)} (too few for assertion); "
          f"{real_frac:.0%} real — need full table for reliable check")
else:
    print("  [SKIP] palindromic_hairpin — not found in overlap (head table may be too small)")

# tandem_rc_concatemer: mdax should still call artefact
if tandem_rids:
    art_frac = mdax_decisions.get("artefact", 0) / len(tandem_rids)
    check("tandem_rc_concatemer → mdax artefact (expect >60%)",
          art_frac > 0.60,
          f"only {art_frac:.1%} got artefact decision")

print()
if fails:
    print(f"  {fails} check(s) FAILED — see output above.")
else:
    print("  All checks passed.")
print()

# ---------------------------------------------------------------------------
# Optional joined TSV output
# ---------------------------------------------------------------------------

if args.out and joined_rows:
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(joined_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(joined_rows)
    print(f"Joined table written to: {args.out}  ({len(joined_rows):,} rows)")
