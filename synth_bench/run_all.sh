#!/usr/bin/env bash
# synth_bench/run_all.sh -- end-to-end synthetic-genome benchmark for mdax.
#
# Generates a reference genome with controlled, ground-truthed inverted-repeat
# loci (varying arm length / identity / local coverage tier), simulates a
# HiFi read set over it with MDA-style coverage bias plus injected foldback
# chimeras, then compares assembly contiguity across correction arms:
#
#   raw          -- no correction at all
#   mdax_*       -- mdax with a few parameter presets/overrides
#   refbased     -- reference-based dechimerization proxy (see refbased_clip.py)
#
# Usage:
#   ./run_all.sh [outdir] [genome_len] [coverage] [threads]
#
# Env overrides: MDAX_BIN, MINIMAP2, HIFIASM, GFATOOLS, SEED
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SB_DIR="$REPO_DIR/synth_bench"

OUTDIR="${1:-$SB_DIR/run1}"
GENOME_LEN="${2:-2000000}"
COVERAGE="${3:-30}"
THREADS="${4:-16}"

MDAX_BIN="${MDAX_BIN:-$REPO_DIR/target/release/mdax}"
MINIMAP2="${MINIMAP2:-/software/team301/minimap2/minimap2}"
HIFIASM="${HIFIASM:-/software/team301/hifiasm/hifiasm}"
GFATOOLS="${GFATOOLS:-/software/team301/gfatools/gfatools}"
SEED="${SEED:-1}"

mkdir -p "$OUTDIR"
echo "════════════════════════════════════════════════════════════════════"
echo "  synth_bench run: $OUTDIR"
echo "  genome=${GENOME_LEN}bp coverage=${COVERAGE}x threads=$THREADS seed=$SEED"
echo "════════════════════════════════════════════════════════════════════"

# ── 1) build mdax ────────────────────────────────────────────────────────────
if [[ ! -x "$MDAX_BIN" ]]; then
    echo "-- building mdax (release)"
    (cd "$REPO_DIR" && RUSTFLAGS='-C target-feature=+aes,+sse2,+avx2' cargo build --release -q)
fi

# ── 2) synthetic genome + reads ──────────────────────────────────────────────
GENOME_DIR="$OUTDIR/genome"
if [[ ! -f "$GENOME_DIR/genome.fa" ]]; then
    echo "-- generating synthetic genome"
    python3 "$SB_DIR/make_genome.py" --outdir "$GENOME_DIR" \
        --genome-len "$GENOME_LEN" --seed "$SEED"
fi

READS_DIR="$OUTDIR/reads"
if [[ ! -f "$READS_DIR/reads.fa" ]]; then
    echo "-- simulating reads"
    python3 "$SB_DIR/simulate_reads.py" \
        --genome "$GENOME_DIR/genome.fa" --bed "$GENOME_DIR/repeats_truth.bed" \
        --coverage-profile "$GENOME_DIR/coverage_profile.tsv" \
        --outdir "$READS_DIR" --target-coverage "$COVERAGE" --seed "$SEED"
fi
RAW_FA="$READS_DIR/reads.fa"

# ── 3) correction arms ────────────────────────────────────────────────────────
CORR_DIR="$OUTDIR/corrected"
mkdir -p "$CORR_DIR"

declare -A ARM_FASTA
ARM_FASTA[raw]="$RAW_FA"

run_mdax () {
    local name="$1"; shift
    local out_fa="$CORR_DIR/${name}.fa"
    local out_tsv="$CORR_DIR/${name}.report.tsv"
    if [[ ! -f "$out_fa" ]]; then
        echo "-- mdax arm: $name ($*)"
        "$MDAX_BIN" -o "$out_fa" -r "$out_tsv" -t "$THREADS" "$@" "$RAW_FA" \
            > "$CORR_DIR/${name}.log" 2>&1 || { echo "mdax arm $name FAILED, see $CORR_DIR/${name}.log"; return 1; }
    fi
    ARM_FASTA[$name]="$out_fa"
}

run_mdax mdax_balanced
run_mdax mdax_balanced_cutlow --cut-low-ident
run_mdax mdax_strict --mode strict
run_mdax mdax_permissive --mode permissive --cut-low-ident
run_mdax mdax_support1 --min-support 1 --cut-low-ident
run_mdax mdax_support8 --min-support 8 --cut-low-ident

# reference-based proxy
REFB_DIR="$OUTDIR/refbased"
if [[ ! -f "$REFB_DIR/reads.refbased.fa" ]]; then
    echo "-- reference-based clipping proxy"
    python3 "$SB_DIR/refbased_clip.py" --genome "$GENOME_DIR/genome.fa" --reads "$RAW_FA" \
        --outdir "$REFB_DIR" --minimap2 "$MINIMAP2" --threads "$THREADS"
fi
ARM_FASTA[refbased]="$REFB_DIR/reads.refbased.fa"

# ── 4) assemble each arm with hifiasm ────────────────────────────────────────
ASM_DIR="$OUTDIR/asm"
mkdir -p "$ASM_DIR"

declare -A ARM_ASSEMBLY
for name in "${!ARM_FASTA[@]}"; do
    fa="${ARM_FASTA[$name]}"
    prefix="$ASM_DIR/$name"
    p_ctg_gfa="${prefix}.bp.p_ctg.gfa"
    p_ctg_fa="${prefix}.p_ctg.fa"
    if [[ ! -f "$p_ctg_fa" ]]; then
        echo "-- hifiasm arm: $name"
        "$HIFIASM" -o "$prefix" -t "$THREADS" "$fa" \
            > "${prefix}.hifiasm.log" 2>&1
        "$GFATOOLS" gfa2fa "$p_ctg_gfa" > "$p_ctg_fa" 2>/dev/null
    fi
    ARM_ASSEMBLY[$name]="$p_ctg_fa"
done

# ── 5) evaluate ───────────────────────────────────────────────────────────────
EVAL_DIR="$OUTDIR/eval"
ASM_ARGS=()
for name in "${!ARM_ASSEMBLY[@]}"; do
    ASM_ARGS+=(--assembly "${name}=${ARM_ASSEMBLY[$name]}")
done

echo "-- evaluating assemblies"
python3 "$SB_DIR/evaluate.py" \
    --truth-genome "$GENOME_DIR/genome.fa" --repeats-bed "$GENOME_DIR/repeats_truth.bed" \
    "${ASM_ARGS[@]}" --outdir "$EVAL_DIR" --minimap2 "$MINIMAP2" --threads "$THREADS" \
    | tee "$OUTDIR/summary.txt"

echo
echo "Done. Summary: $OUTDIR/summary.txt"
echo "Per-locus repeat-breakage table: $EVAL_DIR/per_locus.tsv"
