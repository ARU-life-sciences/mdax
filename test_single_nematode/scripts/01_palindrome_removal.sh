#!/usr/bin/env bash
set -euo pipefail
FASTQ="$1"; ID="$2"; CPU="${3:-16}"
module load minimap2/2.28--he4a0461_1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="${SCRIPT_DIR}:${PATH}"
echo "[info] started: $(date -Iseconds)"
bash "${SCRIPT_DIR}/remove_palindrome.sh" "${FASTQ}" "${ID}" "${CPU}"
echo "[info] finished: $(date -Iseconds)"
ls -la "${ID}.palindrome_treated.fq.gz"
