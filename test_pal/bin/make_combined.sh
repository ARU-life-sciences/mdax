#!/usr/bin/env bash
set -euo pipefail

in1="$1"   # e.g. data/SAMPLE_1.fastq.gz
in2="$2"   # e.g. data/SAMPLE_2.fastq.gz
out="$3"   # e.g. out/SAMPLE.combined.fastq.gz

mkdir -p "$(dirname "$out")"

# Concatenate gzip streams safely
cat "$in1" "$in2" > "$out"

echo "Wrote: $out"

