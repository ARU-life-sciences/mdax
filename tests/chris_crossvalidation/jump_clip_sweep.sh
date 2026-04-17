#!/usr/bin/env bash
# =============================================================================
# jump_clip_sweep.sh — Test 2: max-jump-clip sensitivity analysis
#
# Background
# ----------
# mdax's delta-scan (refine_breakpoint_hamming / refine_breakpoint_banded_ed)
# tests offsets δ = 0 .. max_jump_clip/2 on each side of the estimated split
# to find the alignment that maximises arm identity.  If the true template-
# switch gap G exceeds max_jump_clip, the best δ is capped at max_jump_clip/2
# and the identity estimate is systematically too low.
#
# From the ERR12263839 gap_est histogram (chris_data/README.md):
#   - modal gap bin is 660–1320 bp
#   - ~61% of reads have gap_est > 660 bp
#   - the default --max-jump-clip is 1000, sitting inside the modal bin
#
# This script runs mdax on two inputs at four max-jump-clip values and shows
# how the artefact/low_ident split changes:
#
#   1. data/chris_artefacts.fasta   (57 ground-truth HiFi reads)
#      Focus: do the 8 "li_*" reads (high BLAST identity, currently low_ident
#      in mdax) flip to artefact when the clip window is wider?
#      NOTE: most li_* reads have BLAST gap_est ≈ 0, so this tests a different
#      hypothesis — see expected output comment below.
#
#   2. chris_data/ERR12263839.fastq.gz  (full ~1M-read HiFi dataset, optional)
#      Focus: population-level shift in artefact detection rate.
#      Skipped automatically if the file is not present.
#
# Usage
# -----
#   cd tests/chris_crossvalidation
#   ./jump_clip_sweep.sh
#
#   # Override binary path or number of threads:
#   MDAX_BIN=../../target/release/mdax THREADS=4 ./jump_clip_sweep.sh
#
# Expected output
# ---------------
# For the 57-read set most li_* reads are expected to STAY low_ident across all
# max-jump-clip values — their low identity is caused by the refine window being
# smaller than the arm (refine_arm=1200 < arm_len=5000–9000 bp), not by a large
# gap.  Confirmed if li_* decision is stable across columns.
#
# For the full dataset (if present), we expect the artefact count to increase
# when max_jump_clip grows from 1000 → 2000 (reads whose gap_est is in the
# 660–1320 bp modal bin are rescued), with diminishing returns above 3000 bp.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}/../.."

MDAX_BIN="${MDAX_BIN:-${REPO_ROOT}/target/release/mdax}"
THREADS="${THREADS:-4}"

SMALL_FASTA="${REPO_ROOT}/data/chris_artefacts.fasta"
LARGE_FASTQ="${REPO_ROOT}/chris_data/ERR12263839.fastq.gz"

JUMP_CLIPS=(500 1000 2000 4000)

# ── sanity checks ─────────────────────────────────────────────────────────────

if [[ ! -x "$MDAX_BIN" ]]; then
    echo "ERROR: mdax binary not found or not executable: $MDAX_BIN" >&2
    echo "  Build it with:  cargo build --release" >&2
    exit 1
fi

if [[ ! -f "$SMALL_FASTA" ]]; then
    echo "ERROR: 57-read test FASTA not found: $SMALL_FASTA" >&2
    exit 1
fi

TMPDIR_SWEEP="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_SWEEP"' EXIT

# ── helper: count decisions in a report TSV ───────────────────────────────────

count_decision() {
    local tsv="$1" decision="$2"
    awk -F'\t' -v d="$decision" 'NR>1 && $NF==d' "$tsv" | wc -l | tr -d ' '
}

# Extract decision for a specific read label (matches first token of read_id)
decision_for() {
    local tsv="$1" label="$2"
    awk -F'\t' -v l="$label" 'NR>1 && $1 ~ "^" l { print $NF; exit }' "$tsv"
}

print_header() {
    printf "\n  %-22s" "label"
    for jc in "${JUMP_CLIPS[@]}"; do
        printf "  %8s" "jc=$jc"
    done
    printf "\n  %s\n" "$(printf '─%.0s' {1..80})"
}

# =============================================================================
# Part A: 57-read ground-truth set
# =============================================================================

echo ""
echo "════════════════════════════════════════════════════════════════════════════"
echo "  PART A: 57-read ground-truth set (data/chris_artefacts.fasta)"
echo "  Config: --refine-mode hifi --min-identity 0.5 --mode balanced"
echo "  Varying: --max-jump-clip"
echo "════════════════════════════════════════════════════════════════════════════"

declare -A REPORT_A
for jc in "${JUMP_CLIPS[@]}"; do
    rpt="${TMPDIR_SWEEP}/small_jc${jc}.tsv"
    "$MDAX_BIN" \
        --threads "$THREADS" \
        --max-jump-clip "$jc" \
        --refine-mode hifi \
        --min-identity 0.5 \
        --mode balanced \
        --report "$rpt" \
        --output /dev/null \
        "$SMALL_FASTA" 2>/dev/null
    REPORT_A[$jc]="$rpt"
done

# Summary counts
echo ""
echo "  Overall decision counts:"
printf "  %-14s" "decision"
for jc in "${JUMP_CLIPS[@]}"; do printf "  %8s" "jc=$jc"; done
printf "\n  %s\n" "$(printf '─%.0s' {1..60})"

for dec in artefact low_ident real; do
    printf "  %-14s" "$dec"
    for jc in "${JUMP_CLIPS[@]}"; do
        n="$(count_decision "${REPORT_A[$jc]}" "$dec")"
        printf "  %8s" "$n"
    done
    printf "\n"
done

# li_* reads specifically
echo ""
echo "  li_* reads (high BLAST identity, expected low_ident in default config):"
echo "  Note: most li_* reads have BLAST gap_est ≈ 0; low_ident is not gap-driven."

LI_LABELS=(li_3279bp li_5619bp li_6296bp li_7152bp li_9340bp
           li_11367bp li_12397bp li_14193bp)

print_header
for label in "${LI_LABELS[@]}"; do
    printf "  %-22s" "$label"
    for jc in "${JUMP_CLIPS[@]}"; do
        d="$(decision_for "${REPORT_A[$jc]}" "$label")"
        printf "  %8s" "${d:-n/a}"
    done
    printf "\n"
done

# art_hi_* reads: should be stable artefact regardless
echo ""
echo "  art_hi_* reads (clean hairpins, should remain artefact throughout):"
ART_HI_LABELS=(art_hi_3820bp art_hi_5653bp art_hi_6156bp art_hi_9193bp
               art_hi_11012bp art_hi_13019bp)

print_header
for label in "${ART_HI_LABELS[@]}"; do
    printf "  %-22s" "$label"
    for jc in "${JUMP_CLIPS[@]}"; do
        d="$(decision_for "${REPORT_A[$jc]}" "$label")"
        printf "  %8s" "${d:-n/a}"
    done
    printf "\n"
done

# =============================================================================
# Part B: full dataset (optional)
# =============================================================================

echo ""
echo "════════════════════════════════════════════════════════════════════════════"
echo "  PART B: full ERR12263839 dataset (chris_data/ERR12263839.fastq.gz)"
echo "════════════════════════════════════════════════════════════════════════════"

if [[ ! -f "$LARGE_FASTQ" ]]; then
    echo ""
    echo "  SKIP: $LARGE_FASTQ not found."
    echo "  To enable, place the full FASTQ in chris_data/ and re-run."
    echo ""
    exit 0
fi

echo ""
echo "  Running on full dataset (this may take several minutes) …"
echo ""

declare -A REPORT_B
for jc in "${JUMP_CLIPS[@]}"; do
    echo "    max-jump-clip=$jc …"
    rpt="${TMPDIR_SWEEP}/large_jc${jc}.tsv"
    "$MDAX_BIN" \
        --threads "$THREADS" \
        --max-jump-clip "$jc" \
        --refine-mode hifi \
        --mode balanced \
        --report "$rpt" \
        --output /dev/null \
        "$LARGE_FASTQ" 2>/dev/null
    REPORT_B[$jc]="$rpt"
done

echo ""
echo "  Overall decision counts (full dataset):"
printf "  %-14s" "decision"
for jc in "${JUMP_CLIPS[@]}"; do printf "  %10s" "jc=$jc"; done
printf "\n  %s\n" "$(printf '─%.0s' {1..70})"

for dec in artefact low_ident real; do
    printf "  %-14s" "$dec"
    for jc in "${JUMP_CLIPS[@]}"; do
        n="$(count_decision "${REPORT_B[$jc]}" "$dec")"
        printf "  %10s" "$n"
    done
    printf "\n"
done

echo ""
echo "  Hypothesis: artefact count increases from jc=1000 → jc=2000 (reads in the"
echo "  660–1320 bp modal gap bin are rescued); diminishing returns above jc=3000."
echo ""
