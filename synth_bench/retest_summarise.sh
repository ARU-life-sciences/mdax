#!/usr/bin/env bash
# One line per assembly: contigs, total length, N50, largest, BUSCO string.
B=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_celegans_rerun; S=${1:-nrCaeEleg92}
row() { name=$1; fa=$2; bs=$3; [ -s "$fa" ] || return 0
  st=$(seqkit stats -a -T "$fa" | tail -1 | awk -F'\t' '{printf "%6d %7.1f %8.2f %7.2f", $4, $5/1e6, $13/1e6, $8/1e6}')
  bu=$(grep -h "C:" $bs 2>/dev/null | head -1 | tr -d '\t ')
  printf "%-20s %s  %s\n" "$name" "$st" "$bu"; }
printf "%-20s %6s %7s %8s %7s  %s\n" arm contigs Mb N50_Mb max_Mb BUSCO
for a in raw strict balanced permissive; do row "old_$a" $B/asm/$S/$a/${S}_$a.p_ctg.fa "$B/busco/${S}_$a/short_summary*.txt"; done
for d in $B/retest92/$S.*; do row "${d##*.}" $d/p_ctg.fa $d/busco_summary.txt; done
