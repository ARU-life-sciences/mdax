#!/usr/bin/env bash
# Download raw ONT reads for one species from ENA (fast, no auth) straight
# into the lustre reads/ dir, used as the upstream dependency for
# run_species_pipeline.sh when reads aren't already staged locally.
set -euo pipefail
URL="$1"   # https:// or ftp:// path to the fastq.gz on ENA
OUT="$2"   # destination path, e.g. .../reads/aphelenchoides.fastq.gz

echo "[info] started: $(date -Iseconds)"
echo "[info] fetching ${URL} -> ${OUT}"
curl --fail --retry 5 --retry-delay 10 -o "${OUT}.part" "${URL}"
mv "${OUT}.part" "${OUT}"

echo "[info] integrity check (gzip -t)"
gzip -t "${OUT}"

echo "[info] finished: $(date -Iseconds)"
ls -la "${OUT}"
