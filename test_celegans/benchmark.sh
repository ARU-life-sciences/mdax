#!/usr/bin/env bash
# benchmark.sh — run mdax on C. elegans MDA samples and compare artefact rates.
#
# Samples available (data/):
#   nrCaeEleg92  Repli-G DLB rep 2, 120 min  →  82.3% BUSCO, 1917 contigs (moderate artefact)
#   nrCaeEleg95  Repli-G DLB 120 min         →  64.5% BUSCO, 7944 contigs (high artefact)
#
# Usage:
#   ./benchmark.sh                  # process up to MAX_READS from each sample
#   MAX_READS=0 ./benchmark.sh      # process ALL reads (slow for nrCaeEleg95)
#   THREADS=16  ./benchmark.sh      # override thread count
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

MDAX_BIN="${MDAX_BIN:-$REPO_DIR/target/release/mdax}"
DATA_DIR="$SCRIPT_DIR/data"
OUT_DIR="$SCRIPT_DIR/out"

# Default: cap at 200k reads per sample for a quick but representative benchmark.
# Set MAX_READS=0 to process everything.
MAX_READS="${MAX_READS:-200000}"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR"

[[ -x "$MDAX_BIN" ]] || { echo "mdax binary not found at $MDAX_BIN — run cargo build --release first."; exit 1; }

# ── helper: optional subsetting ───────────────────────────────────────────────

# Extract up to MAX_READS FASTA entries from a gzip-compressed FASTA.
# If MAX_READS=0, just decompress in full.
subset_fasta() {
    local src="$1"
    local dst="$2"
    if [[ "$MAX_READS" -eq 0 ]]; then
        echo "  (no read cap — decompressing full file)"
        gunzip -c "$src" > "$dst"
    else
        echo "  (capping at $MAX_READS reads)"
        python3 - "$src" "$dst" "$MAX_READS" <<'PYEOF'
import sys, gzip
src, dst, n_max = sys.argv[1], sys.argv[2], int(sys.argv[3])
n = 0
with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
    for line in fi:
        if line.startswith(b">"):
            n += 1
            if n > n_max:
                break
        fo.write(line)
PYEOF
    fi
}

# ── process each sample ────────────────────────────────────────────────────────

declare -A SAMPLE_LABELS=(
    ["nrCaeEleg92"]="Repli-G DLB rep 2, 120 min  (assembly: 110 Mb, 1917 ctgs, 82.3% BUSCO)"
    ["nrCaeEleg95"]="Repli-G DLB 120 min          (assembly: 433 Mb, 7944 ctgs, 64.5% BUSCO)"
)

SAMPLES=("nrCaeEleg92" "nrCaeEleg95")

for sample in "${SAMPLES[@]}"; do
    src="$DATA_DIR/${sample}_ccs.trim.fasta.gz"
    if [[ ! -f "$src" ]]; then
        echo "[SKIP] $src not found"
        continue
    fi

    echo
    echo "════════════════════════════════════════════════════════════"
    echo "  Sample: $sample"
    echo "  ${SAMPLE_LABELS[$sample]}"
    echo "════════════════════════════════════════════════════════════"

    subset_fa="$OUT_DIR/${sample}.subset.fasta"
    report="$OUT_DIR/${sample}.report.tsv"
    out_fa="$OUT_DIR/${sample}.out.fasta"

    echo "  Preparing subset → $subset_fa"
    subset_fasta "$src" "$subset_fa"

    echo "  Running mdax (threads=$THREADS)…"
    "$MDAX_BIN" \
        --threads "$THREADS" \
        --cut-low-ident \
        --report  "$report" \
        --output  "$out_fa" \
        "$subset_fa" 2>&1 | grep -E '^\[|Foldback|real' | sed 's/^/    /'

    echo "  Report written → $report"
done

# ── generate combined report ───────────────────────────────────────────────────

echo
echo "════════════════════════════════════════════════════════════"
echo "  Generating comparison report…"
echo "════════════════════════════════════════════════════════════"

python3 "$SCRIPT_DIR/report.py" \
    --samples nrCaeEleg92 nrCaeEleg95 \
    --labels  "nrCaeEleg92 (82.3% BUSCO, 1917 ctgs)" \
              "nrCaeEleg95 (64.5% BUSCO, 7944 ctgs)" \
    --reports "$OUT_DIR/nrCaeEleg92.report.tsv" \
              "$OUT_DIR/nrCaeEleg95.report.tsv" \
    --max-reads "$MAX_READS" \
    --out-tsv "$OUT_DIR/summary.tsv"
