#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
BASELINE_DIR="out/baseline/SRR24201687_first20k"
INPUT_FQ="out/subsample/SRR24201687.first20k.fastq.gz"

MDAX_BIN="../target/release/mdax"

# Toggle fairness baseline:
#   FAIR=1 ./test_pal.sh
FAIR="${FAIR:-0}"

# Outputs produced by this script
MDAX_REPORT="mdax.report.tsv"
MDAX_OUT_FASTA="mdax.out.fa"
MDAX_CALLS="mdax_calls.tsv"

# Derived canonical baseline files
BASELINE_CANON_IDS="baseline.flagged.canon.ids"

# ---- sanity checks ----
for f in \
  "$BASELINE_DIR/baseline.flagged.norm.ids" \
  "$BASELINE_DIR/baseline.iter2.breakpoints.tsv" \
  "$INPUT_FQ" \
  "$MDAX_BIN"
do
  [[ -f "$f" ]] || { echo "Missing $f"; exit 1; }
done

echo "== Extracting baseline breakpoints (iteration 2) =="

# baseline.iter2.breakpoints.tsv: col1=id, col5=breakpoint
# Canonicalise baseline IDs here too (strip -dup...)
awk -F'\t' '
  NF>=5 {
    id=$1;
    sub(/-dup[0-9]+(\.[0-9]+)?$/, "", id);
    print id "\t" $5
  }
' "$BASELINE_DIR/baseline.iter2.breakpoints.tsv" \
  | LC_ALL=C sort -u > baseline.bp.tsv

echo "Baseline breakpoints: $(wc -l < baseline.bp.tsv)"

echo "== Canonicalising baseline flagged IDs =="

# baseline.flagged.norm.ids contains IDs, sometimes with -dup...
sed -E 's/-dup[0-9]+(\.[0-9]+)?$//' \
  "$BASELINE_DIR/baseline.flagged.norm.ids" \
  | LC_ALL=C sort -u > "$BASELINE_CANON_IDS"

echo "Baseline flagged (canonical): $(wc -l < "$BASELINE_CANON_IDS")"

# ---- run mdax ----
echo "== Running mdax (FAIR=${FAIR}) =="

# If FAIR=1, enable developer-only fairness baseline.
# Prefer the CLI flag (since you added it), but also set env for convenience.
if [[ "$FAIR" == "1" ]]; then
  export MDAX_FAIRNESS=1
  FAIR_FLAG="--fairness-baseline"
else
  unset MDAX_FAIRNESS || true
  FAIR_FLAG=""
fi

# Note:
# - mdax writes the TSV report to --report; do NOT rely on stdout parsing.
# - mdax writes corrected FASTA to --output.
"$MDAX_BIN" $FAIR_FLAG \
  --report "$MDAX_REPORT" \
  --output "$MDAX_OUT_FASTA" \
  "$INPUT_FQ" >/dev/null

echo "MDAX report: $MDAX_REPORT"
echo "MDAX output: $MDAX_OUT_FASTA"

# ---- normalise mdax calls ----
echo "== Normalising mdax calls =="

# mdax.report.tsv header:
# read_id len event called coarse_split refined_split ... identity_est support_n support_span decision
#
# We want: canonical_id, called(0/1), refined_split (blank if NA), decision
#
# Canonicalisation rules:
#  - keep only first whitespace token (needletail id can contain extra tokens)
#  - then strip any -dup... suffix (so it matches baseline canon)
awk -F'\t' '
  NR==1{
    for(i=1;i<=NF;i++){
      if ($i=="refined_split") rs_col=i;
      if ($i=="decision")      dec_col=i;
    }
    next
  }
  {
    # $1 may contain spaces; take first token
    split($1, a, " ");
    id=a[1];

    # strip perl-style dup suffix if present
    sub(/-dup[0-9]+(\.[0-9]+)?$/, "", id);

    rs  = rs_col  ? $rs_col  : "";
    dec = dec_col ? $dec_col : "";

    # Compare against baseline "flagged for chopping"
    called = (dec=="artefact") ? 1 : 0;

    print id "\t" called "\t" rs "\t" dec
  }
' "$MDAX_REPORT" | LC_ALL=C sort -u > mdax.calls.norm.tsv

echo "MDAX calls (rows in mdax.calls.norm.tsv): $(wc -l < mdax.calls.norm.tsv)"

# Keep a copy in the old name if you want
cp -f mdax.calls.norm.tsv "$MDAX_CALLS"

# ---- quick sanity check: do IDs overlap now? ----
echo "== Sanity check: ID overlap sample =="
comm -12 \
  <(cut -f1 mdax.calls.norm.tsv | LC_ALL=C sort -u) \
  <(LC_ALL=C sort -u "$BASELINE_CANON_IDS") \
  | head

echo "== Running comparison =="
python3 ./test_pal.py

echo
echo "== Disagreement lists written =="

# mdax_only = called by mdax but not by baseline
comm -23 \
  <(awk -F'\t' '$2==1{print $1}' mdax.calls.norm.tsv | LC_ALL=C sort -u) \
  <(LC_ALL=C sort -u "$BASELINE_CANON_IDS") \
  > mdax_only.ids

# baseline_only = called by baseline but not by mdax
comm -13 \
  <(awk -F'\t' '$2==1{print $1}' mdax.calls.norm.tsv | LC_ALL=C sort -u) \
  <(LC_ALL=C sort -u "$BASELINE_CANON_IDS") \
  > baseline_only.ids

echo "mdax_only.ids      $(wc -l < mdax_only.ids)"
echo "baseline_only.ids  $(wc -l < baseline_only.ids)"

echo
echo "Tip: run FAIR=1 ./test_pal.sh to benchmark under the locked baseline."

