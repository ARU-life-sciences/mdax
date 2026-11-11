#!/usr/bin/env bash
# Flye for raw PacBio subreads (not HiFi, not ONT) -- Lepidodermella squamata
# gastrotricha benchmark. --pacbio-raw matches the noisier, non-CCS-corrected
# subread error profile; mdax itself is run with --refine-mode ont on this
# same data (error-tolerance driven choice, not a platform claim).
set -euo pipefail
INPUT="$1"; OUTDIR="$2"; CPU="${3:-16}"
source /nfs/users/nfs_m/mb39/miniconda3/etc/profile.d/conda.sh
conda activate flye_env
echo "[info] started: $(date -Iseconds)"
flye --pacbio-raw "${INPUT}" --meta --out-dir "${OUTDIR}" --threads "${CPU}"
echo "[info] finished: $(date -Iseconds)"
