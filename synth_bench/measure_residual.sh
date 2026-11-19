#!/usr/bin/env bash
# Run mdax on a full library, keep only the corrected versions of the reads in
# an existing reference-aligned truth sample (the full corrected FASTA goes
# through a FIFO and never lands on disk), re-align them, and score what was
# left behind and how well the cuts were placed (measure_residual.py).
#
# usage: measure_residual.sh NAME TRUTH_SAM READS PRESET MDAXBIN [mdax args...]
#   PRESET is the minimap2 preset used for TRUTH_SAM (map-hifi | map-ont)
# Needs minimap2 and seqkit on PATH (bq -l ...). Run from the output directory.
set -euo pipefail
NAME=$1; TRUTH=$2; READS=$3; PRESET=$4; MDAX=$5; shift 5
HERE=$(cd "$(dirname "$0")" && pwd)
REF=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_celegans_rerun/reference/GCF_000002985.6_WBcel235_genomic.fna
T=${THREADS:-16}

grep -v '^@' "$TRUTH" | cut -f1 | sort -u > "$NAME.ids"

rm -f "$NAME.fifo"; mkfifo "$NAME.fifo"
python3 "$HERE/filter_sampled_reads.py" "$NAME.ids" < "$NAME.fifo" > "$NAME.corr.fa" &
FILT=$!
"$MDAX" "$READS" "$@" -t "$T" -o "$NAME.fifo" -r "$NAME.report.full.tsv" 2> "$NAME.mdax.log"
wait $FILT
rm -f "$NAME.fifo"

# report rows for the sample only
awk -F'\t' 'NR==FNR{ids[$1]=1; next} FNR==1 || ($1 in ids)' "$NAME.ids" <(awk -F'\t' 'BEGIN{OFS="\t"}{split($1,a," "); $1=a[1]; print}' "$NAME.report.full.tsv") > "$NAME.report.tsv"
rm -f "$NAME.report.full.tsv"

case "$READS" in
    *.fastq|*.fq|*.fastq.gz|*.fq.gz) seqkit grep -f "$NAME.ids" "$READS" | seqkit fq2fa -w 0 > "$NAME.raw.fa" ;;
    *) seqkit grep -f "$NAME.ids" -w 0 "$READS" > "$NAME.raw.fa" ;;
esac

minimap2 -ax "$PRESET" --secondary=no -t "$T" "$REF" "$NAME.corr.fa" 2> /dev/null \
    | awk '$1~/^@/ || !and($2,256)' | cut -f1-6 > "$NAME.corr.sam"

python3 "$HERE/measure_residual.py" "$TRUTH" "$NAME.corr.sam" "$NAME.report.tsv" "$NAME.raw.fa" "$NAME.corr.fa" > "$NAME.result.txt"

if [ "${KEEP:-0}" != 1 ]; then
    rm -f "$NAME.corr.sam" "$NAME.corr.fa" "$NAME.raw.fa" "$NAME.ids"
fi
