#!/usr/bin/env bash
# dotplots.sh — sample detected foldback reads by identity bin and plot self-dotplots.
#
# Usage:
#   ./dotplots.sh                          # both samples, 5 reads/bin
#   SAMPLE=nrCaeEleg92 ./dotplots.sh       # one sample only
#   N=10 ./dotplots.sh                     # more reads per bin
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$SCRIPT_DIR/out"
DOTPLOT_DIR="$SCRIPT_DIR/dotplots"
SAMPLES="${SAMPLE:-nrCaeEleg92 nrCaeEleg95}"
N="${N:-5}"
SEED="${SEED:-42}"

mkdir -p "$DOTPLOT_DIR"

for sample in $SAMPLES; do
    report="$OUT_DIR/${sample}.report.tsv"
    fasta="$OUT_DIR/${sample}.subset.fasta"

    [[ -f "$report" ]] || { echo "[SKIP] $report not found — run benchmark.sh first"; continue; }
    [[ -f "$fasta"  ]] || { echo "[SKIP] $fasta not found — run benchmark.sh first"; continue; }

    echo
    echo "══════════════════════════════════════════════"
    echo "  $sample"
    echo "══════════════════════════════════════════════"

    python3 - "$report" "$fasta" "$DOTPLOT_DIR/$sample" "$N" "$SEED" <<'PYEOF'
import csv, os, random, subprocess, sys
from collections import defaultdict

report_path, fasta_path, out_prefix, n_per_bin, seed = \
    sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

random.seed(seed)

BINS = [
    ("0.30-0.40", 0.30, 0.40),
    ("0.40-0.50", 0.40, 0.50),
    ("0.50-0.60", 0.50, 0.60),
    ("0.60-0.70", 0.60, 0.70),
    ("0.70-1.00", 0.70, 1.01),
]

# ── 1. read report, bucket by identity ───────────────────────────────────────
binned = defaultdict(list)
with open(report_path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        try:
            ident = float(row["identity_est"])
            rid   = row["read_id"].split()[0]
            dec   = row["decision"]
        except (ValueError, KeyError):
            continue
        for name, lo, hi in BINS:
            if lo <= ident < hi:
                binned[name].append((rid, ident, dec))
                break

sampled = {}
for name, lo, hi in BINS:
    pool = binned[name]
    s    = random.sample(pool, min(n_per_bin, len(pool)))
    sampled[name] = s
    print(f"  ident {name}: {len(s)}/{len(pool)} reads sampled"
          + (f"  (decisions: { {dec for _,_,dec in s} })" if s else ""))

all_ids = {rid for reads in sampled.values() for (rid, _, _) in reads}

# ── 2. extract sequences from FASTA ──────────────────────────────────────────
seqs = {}
current_id  = None
current_seq = []
with open(fasta_path) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if current_id and current_id in all_ids:
                seqs[current_id] = "".join(current_seq)
            current_id  = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id and current_id in all_ids:
        seqs[current_id] = "".join(current_seq)

print(f"  Sequences extracted: {len(seqs)}/{len(all_ids)}")

# ── 3. write per-bin FASTAs and run flexidot ──────────────────────────────────
# Use k=10 throughout — large enough to avoid noise on CCS reads, fast enough
# for ~10kb reads. Substitutions are too slow at this read length so kept at 0.
# The anti-diagonal foldback pattern is visible even at low identity because
# the repeated region spans thousands of bases.
bin_params = {
    "0.30-0.40": ("-k", "10", "-S", "0"),
    "0.40-0.50": ("-k", "10", "-S", "0"),
    "0.50-0.60": ("-k", "10", "-S", "0"),
    "0.60-0.70": ("-k", "10", "-S", "0"),
    "0.70-1.00": ("-k", "10", "-S", "0"),
}

for name, lo, hi in BINS:
    reads = sampled[name]
    if not reads:
        print(f"  [{name}] no reads — skipping")
        continue

    kw = bin_params[name]
    bin_dir = f"{out_prefix}_ident{name}"
    os.makedirs(bin_dir, exist_ok=True)

    print(f"  [{name}] flexidot {' '.join(kw)} on {len(reads)} reads …", flush=True)
    for i, (rid, ident, dec) in enumerate(reads):
        if rid not in seqs:
            continue
        seq = seqs[rid]
        safe_rid = rid.replace("/", "_").replace(" ", "_")
        fa_path = os.path.join(bin_dir, f"{i+1:02d}_{safe_rid}.fasta")
        out_name = f"{i+1:02d}_{safe_rid}"
        with open(fa_path, "w") as f:
            f.write(f">{safe_rid}  ident={ident:.3f}  {dec}  len={len(seq)}\n{seq}\n")
        cmd = [
            "flexidot",
            "-i", fa_path,
            "-o", out_name,
            "--outdir", bin_dir,
            "-m", "0",       # self dotplot
            "-f", "png",
            *kw,
            "-P", "15",
            "-E", "5",
            "-T", "40",
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"    WARN [{rid}]: {r.stderr.strip()[:200]}")
    print(f"    → {bin_dir}/")

PYEOF
done

echo
echo "Dotplots written to: $DOTPLOT_DIR"
