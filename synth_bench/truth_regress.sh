#!/usr/bin/env bash
# Reference-based regression check for an mdax binary on the Sanger farm.
#
# Runs the binary on four full amplified libraries and scores a random read
# sample of each against a reference-derived truth (truth_join.py, copy in
# synth_bench/refbased_recall_eval.py):
#   nrCaeEleg92      C. elegans MDA HiFi, ~100x        (balanced)
#   nrCaeEleg94      C. elegans MDA HiFi, ~720x        (balanced)
#   celegans_wga     C. elegans MDA ONT, SRR24201716   (permissive, ONT refine)
#   aphelenchoides   A. besseyi MDA ONT, SRR24201688   (permissive, ONT refine;
#                    truth vs independent strain AORJ, GCA_024699835.1)
#
# usage:
#   truth_regress.sh submit MDAXBIN NAME   # one LSF job per library; prints job ids
#   truth_regress.sh check NAME            # table + exit 1 if any floor/ceiling fails
#
# Floors are set ~2 points under main at 2026-11-19 (frozen retest92/bin/mdax_v9) for
# the share of true foldbacks cut, and ~0.5 points over it for the share of linear
# reads cut. Raise them when a change lands that moves the numbers for good.
# Freeze a copy of the binary before submitting: the jobs read it when they start.
set -euo pipefail
U=/lustre/scratch122/tol/teams/blaxter/users/mb39
B=$U/mdax_celegans_rerun; D=$B/recall_check; O=$D/regress
SAMPLES="nrCaeEleg92 nrCaeEleg94 celegans_wga aphelenchoides"

# sample -> reads | mdax args | min foldback cut % | max linear cut %
spec() {
    case $1 in
        nrCaeEleg92)    echo "$B/reads/nrCaeEleg92_raw.fasta.gz|--mode balanced|77.0|2.8" ;;
        nrCaeEleg94)    echo "$B/reads/nrCaeEleg94_raw.fasta.gz|--mode balanced|72.0|2.3" ;;
        celegans_wga)   echo "$U/mdax_single_nematode_benchmark/reads/celegans_wga.fastq.gz|--mode permissive --refine-mode ont|75.5|2.7" ;;
        aphelenchoides) echo "$U/mdax_single_nematode_benchmark/reads/aphelenchoides.fastq.gz|--mode permissive --refine-mode ont|78.0|5.1" ;;
        *) echo "unknown sample $1" >&2; exit 2 ;;
    esac
}

case ${1:-} in
submit)
    M=$(readlink -f "$2"); NAME=$3; mkdir -p $O; cd $O
    for S in $SAMPLES; do
        IFS='|' read -r IN ARGS _ _ <<< "$(spec $S)"
        # shellcheck disable=SC2086
        bq -m 48G -n 16 -t 4h -J "tr_${NAME}_$S" -- bash -c \
            "$M $IN $ARGS -t 16 -o /dev/null -r $O/$NAME.$S.tsv 2> $O/$NAME.$S.log \
             && python3 $D/truth_join.py $D/$S.sam $O/$NAME.$S.tsv > $O/$NAME.$S.result.txt \
             && rm -f $O/$NAME.$S.tsv"
    done ;;
check)
    NAME=$2; fail=0
    printf "%-16s %10s %8s %10s %8s  %s\n" sample fb_cut floor lin_cut ceil status
    for S in $SAMPLES; do
        IFS='|' read -r _ _ FLOOR CEIL <<< "$(spec $S)"
        R=$O/$NAME.$S.result.txt
        if [ ! -s "$R" ]; then printf "%-16s %s\n" $S "MISSING ($R)"; fail=1; continue; fi
        fb=$(awk '$1=="foldback" && $2~/^n=/ {for(i=3;i<=NF;i++) if($i~/^artefact=/){sub(/artefact=/,"",$i); sub(/%/,"",$i); print $i}}' "$R")
        li=$(awk '$1=="linear"   && $2~/^n=/ {for(i=3;i<=NF;i++) if($i~/^artefact=/){sub(/artefact=/,"",$i); sub(/%/,"",$i); print $i}}' "$R")
        st=$(awk -v f="${fb:-0}" -v l="${li:-100}" -v F=$FLOOR -v C=$CEIL 'BEGIN{print (f>=F && l<=C) ? "ok" : "FAIL"}')
        [ "$st" = ok ] || fail=1
        printf "%-16s %9s%% %7s%% %9s%% %7s%%  %s\n" $S "${fb:-?}" $FLOOR "${li:-?}" $CEIL $st
    done
    exit $fail ;;
*)
    sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 2 ;;
esac
