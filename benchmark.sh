#!/usr/bin/env bash
# benchmark.sh — unified mdax benchmark runner.
#
# Runs two evaluations and prints a combined summary:
#
#   1. test_pal   — compare against the Strobl-2023 Perl baseline on
#                   20k SRR24201687 ONT reads with known artefact calls.
#
#   2. test_celegans — run on real single-worm C. elegans MDA (PacBio CCS)
#                   samples and compare artefact rates between:
#                     nrCaeEleg92 (Repli-G DLB rep 2, 82.3% BUSCO)
#                     nrCaeEleg95 (Repli-G DLB 120 min, 64.5% BUSCO)
#
# Usage:
#   ./benchmark.sh                     # default: cap C. elegans at 200k reads, 8 threads
#   MAX_READS=0 ./benchmark.sh         # process full C. elegans datasets (slow)
#   THREADS=16  ./benchmark.sh         # use more threads for C. elegans
#   SKIP_CELEGANS=1 ./benchmark.sh     # only run test_pal comparison
#   SKIP_PAL=1 ./benchmark.sh          # only run C. elegans benchmark
#
# Requirements:
#   cargo build --release       (builds ./target/release/mdax)
#   test_pal/out/ populated     (run test_pal/test_pal.sh first to get the baseline)
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDAX_BIN="$REPO_DIR/target/release/mdax"

MAX_READS="${MAX_READS:-200000}"
THREADS="${THREADS:-8}"
SKIP_PAL="${SKIP_PAL:-0}"
SKIP_CELEGANS="${SKIP_CELEGANS:-0}"

export MDAX_BIN MAX_READS THREADS

# ── build ─────────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════════"
echo "  Building mdax (release)…"
echo "════════════════════════════════════════════════════════════════════"
cd "$REPO_DIR"
RUSTFLAGS='-C target-feature=+aes' cargo build --release -q
echo "  Binary: $MDAX_BIN  ($(stat -f%z "$MDAX_BIN" 2>/dev/null || stat -c%s "$MDAX_BIN") bytes)"

# ── Part 1: test_pal ──────────────────────────────────────────────────────────
if [[ "$SKIP_PAL" != "1" ]]; then
    echo
    echo "════════════════════════════════════════════════════════════════════"
    echo "  PART 1: test_pal  (Strobl-2023 baseline comparison)"
    echo "════════════════════════════════════════════════════════════════════"
    PAL_DIR="$REPO_DIR/test_pal"

    if [[ ! -f "$PAL_DIR/mdax.calls.norm.tsv" ]]; then
        echo "  mdax.calls.norm.tsv not found — running test_pal.sh first…"
        cd "$PAL_DIR"
        bash test_pal.sh
        cd "$REPO_DIR"
    else
        echo "  Using existing mdax.calls.norm.tsv (re-run test_pal/test_pal.sh to refresh)"
    fi

    cd "$PAL_DIR"
    python3 test_pal.py
    cd "$REPO_DIR"
fi

# ── Part 2: C. elegans ────────────────────────────────────────────────────────
if [[ "$SKIP_CELEGANS" != "1" ]]; then
    echo
    echo "════════════════════════════════════════════════════════════════════"
    echo "  PART 2: test_celegans  (real single-worm MDA samples)"
    echo "════════════════════════════════════════════════════════════════════"
    cd "$REPO_DIR/test_celegans"
    bash benchmark.sh
    cd "$REPO_DIR"
fi

echo
echo "════════════════════════════════════════════════════════════════════"
echo "  Benchmark complete."
echo "════════════════════════════════════════════════════════════════════"
echo
