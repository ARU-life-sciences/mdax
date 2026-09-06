#!/usr/bin/env bash
set -euo pipefail
INPUT="$1"; OUTDIR="$2"; CPU="${3:-16}"
source /nfs/users/nfs_m/mb39/miniconda3/etc/profile.d/conda.sh
conda activate flye_env
echo "[info] started: $(date -Iseconds)"
flye --nano-raw "${INPUT}" --meta --out-dir "${OUTDIR}" --threads "${CPU}"
echo "[info] finished: $(date -Iseconds)"
