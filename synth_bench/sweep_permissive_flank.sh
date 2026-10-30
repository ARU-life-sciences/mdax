#!/bin/bash
# Sweep permissive mode's sig_flank_bp on microlaimidae (the one confirmed
# signature-collision case) to find the smallest flank that resolves the
# blowup, holding sig_take=6 and sig_value_shift=8 (permissive's own
# defaults) fixed so we isolate flank_bp specifically.
set -uo pipefail
cd /software/team301/mdax
BIN=./target/release/mdax
READS=synth_bench/sweep_permissive/microlaimidae_sub200k.fastq
OUT=synth_bench/sweep_permissive
mkdir -p "$OUT"

for fb in 200 300 400 500 600 800; do
    log="$OUT/microlaimidae_fb${fb}.log"
    echo "=== sig_flank_bp=$fb ===" | tee -a "$OUT/summary.txt"
    RUSTFLAGS="-C target-cpu=native" "$BIN" "$READS" \
        --mode permissive --refine-mode ont --cut-low-ident --qc-only \
        --sig-flank-bp "$fb" \
        > "$log" 2>&1
    grep -E "SUPPORT UNIQUENESS|SUPPORT SANITY|QC SUMMARY|QC MODE ADVICE" "$log" | tee -a "$OUT/summary.txt"
done

# Also isolate value_shift at the *original* flank_bp=200 (permissive default),
# in case it's contributing more at this smaller window than it did at
# balanced's flank_bp=600.
for vs in 8 0; do
    log="$OUT/microlaimidae_fb200_vs${vs}.log"
    echo "=== sig_flank_bp=200 sig_value_shift=$vs ===" | tee -a "$OUT/summary.txt"
    RUSTFLAGS="-C target-cpu=native" "$BIN" "$READS" \
        --mode permissive --refine-mode ont --cut-low-ident --qc-only \
        --sig-value-shift "$vs" \
        > "$log" 2>&1
    grep -E "SUPPORT UNIQUENESS|SUPPORT SANITY|QC SUMMARY|QC MODE ADVICE" "$log" | tee -a "$OUT/summary.txt"
done

echo "DONE" | tee -a "$OUT/summary.txt"
