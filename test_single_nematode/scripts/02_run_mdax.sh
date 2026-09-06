#!/usr/bin/env bash
set -euo pipefail
FASTQ="$1"; OUT="$2"; CPU="${3:-16}"
MDAX=/software/team301/mdax/target/release/mdax
echo "[info] started: $(date -Iseconds)"
"${MDAX}" --mode permissive --refine-mode ont --cut-low-ident -t "${CPU}" \
    -o "${OUT}.fa" -r "${OUT}.report.tsv" "${FASTQ}"
echo "[info] finished: $(date -Iseconds)"
