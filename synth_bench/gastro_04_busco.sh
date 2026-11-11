#!/usr/bin/env bash
# BUSCO for Lepidodermella squamata (Gastrotricha) -- no gastrotrich-specific
# odb10 lineage exists, so use metazoa_odb10 (the general animal lineage),
# unlike the nematode suite's nematoda_odb10.
set -euo pipefail
ASM="$1"; NAME="$2"; CPU="${3:-16}"
module load busco/5.4.6--pyhdfd78af_0
echo "[info] started: $(date -Iseconds)"
busco -m genome -l metazoa_odb10 -i "${ASM}" -o "${NAME}" -c "${CPU}" -f
echo "[info] finished: $(date -Iseconds)"
