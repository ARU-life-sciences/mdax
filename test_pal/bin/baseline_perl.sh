#!/usr/bin/env bash
set -euo pipefail

fq="$1"        # e.g. out/SRR24201677.fastq.gz (relative to project root OR absolute)
id="$2"
cpu="${3:-8}"
minlen="${4:-1000}"

# Project root = directory containing this script, then ..
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$script_dir/.." && pwd)"

# Resolve input FASTQ to absolute path (works for relative or absolute)
if [[ "$fq" = /* ]]; then
  fq_abs="$fq"
else
  fq_abs="$root/$fq"
fi

# Ensure output dirs exist
outdir="$root/out/baseline/$id"
mkdir -p "$outdir"
cd "$outdir"

ln -sf "$fq_abs" input.fq.gz

echo "[baseline] root = $root"
echo "[baseline] input = $fq_abs"
echo "[baseline] work  = $outdir"

# ---------- 1st iteration ----------
~/Documents/software/minimap2/minimap2 -t "$cpu" -x ava-ont input.fq.gz input.fq.gz \
  | perl "$root/perl/pafIdentifyPalimdrom.pl" \
  1> "$id.palimProp.1stIte.list" \
  2> "$id.1stIte.log"

perl "$root/perl/fastq_partition.and.chop.palindrome.pl" \
  "$id.palimProp.1stIte.list" input.fq.gz "$minlen" \
  1> "$id.1stIte.chop.log" \
  2> "$id.1stIte.chop.err.log"

# ---------- 2nd iteration ----------
~/Documents/software/minimap2/minimap2 -t "$cpu" -x ava-ont \
  input.fq.gz.exclude.fq.gz input.fq.gz.exclude.fq.gz \
  | perl "$root/perl/pafIdentifyPalimdrom.pl" \
  1> "$id.palimProp.2ndIte.list" \
  2> "$id.2ndIte.log"

perl "$root/perl/fastq_partition.and.chop.palindrome.pl" \
  "$id.palimProp.2ndIte.list" input.fq.gz.exclude.fq.gz "$minlen" \
  1> "$id.2ndIte.chop.log" \
  2> "$id.2ndIte.chop.err.log"

# ---------- final treated reads ----------
cat \
  input.fq.gz.include.fq.gz \
  input.fq.gz.exclude.fq.gz.exclude.fq.gz \
  input.fq.gz.exclude.fq.gz.include.fq.gz \
  > "$id.palindrome_treated.fq.gz"

echo "[baseline] DONE -> $outdir/$id.palindrome_treated.fq.gz"

