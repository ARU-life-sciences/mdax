#!/usr/bin/env bash
# Resubmit only the theirs/mdax arms for the 11-species nematode suite,
# skipping raw -- raw's original flye_raw/busco_raw jobs are still alive
# (pending) and don't depend on the mdax binary, so they were deliberately
# left untouched when the mdax/theirs-dependent jobs were killed and are not
# resubmitted here to avoid duplicate/colliding output.
set -euo pipefail
ROOT="/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_single_nematode_benchmark"
SCRIPTS="${ROOT}/scripts"
LOGDIR="${ROOT}/logs"

SPECIES=(aphelenchoides celegans_wga enoplolaimus mesodorylaimus paralinhomoeus \
         ptycholaimellus rhynchonema sabatieria theristus trissonchulus_lat trissonchulus_wlg)

for SP in "${SPECIES[@]}"; do
    FASTQ="${ROOT}/reads/${SP}.fastq.gz"
    mkdir -p "${ROOT}/corrected/${SP}/theirs" "${ROOT}/corrected/${SP}/mdax" "${ROOT}/busco"

    JOBA=$(cd "${ROOT}/corrected/${SP}/theirs" && bsub \
        -J "${SP}_theirs" -q normal -n 16 \
        -R "select[mem>=64000] span[hosts=1] rusage[mem=64000]" -M 64000 \
        -o "${LOGDIR}/${SP}_theirs.out" -e "${LOGDIR}/${SP}_theirs.err" \
        "${SCRIPTS}/01_palindrome_removal.sh" "${FASTQ}" "${SP}" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

    JOBB=$(bsub \
        -J "${SP}_mdax" -q normal -n 16 \
        -R "select[mem>=32000] span[hosts=1] rusage[mem=32000]" -M 32000 \
        -o "${LOGDIR}/${SP}_mdax.out" -e "${LOGDIR}/${SP}_mdax.err" \
        "${SCRIPTS}/02_run_mdax.sh" "${FASTQ}" "${ROOT}/corrected/${SP}/mdax/${SP}_mdax" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

    JOBC1=$(bsub \
        -J "${SP}_flye_theirs" -w "done(${JOBA})" -q long -n 16 \
        -R "select[mem>=64000] span[hosts=1] rusage[mem=64000]" -M 64000 \
        -o "${LOGDIR}/${SP}_flye_theirs.out" -e "${LOGDIR}/${SP}_flye_theirs.err" \
        "${SCRIPTS}/03_flye.sh" "${ROOT}/corrected/${SP}/theirs/${SP}.palindrome_treated.fq.gz" "${ROOT}/asm/${SP}/theirs" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

    JOBC2=$(bsub \
        -J "${SP}_flye_mdax" -w "done(${JOBB})" -q long -n 16 \
        -R "select[mem>=64000] span[hosts=1] rusage[mem=64000]" -M 64000 \
        -o "${LOGDIR}/${SP}_flye_mdax.out" -e "${LOGDIR}/${SP}_flye_mdax.err" \
        "${SCRIPTS}/03_flye.sh" "${ROOT}/corrected/${SP}/mdax/${SP}_mdax.fa" "${ROOT}/asm/${SP}/mdax" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

    for pair in "theirs:${JOBC1}" "mdax:${JOBC2}"; do
        name="${pair%%:*}"
        dep="${pair##*:}"
        bsub \
            -J "${SP}_busco_${name}" -w "done(${dep})" -q long -n 16 \
            -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
            -o "${LOGDIR}/${SP}_busco_${name}.out" -e "${LOGDIR}/${SP}_busco_${name}.err" \
            bash -c "cd '${ROOT}/busco' && '${SCRIPTS}/04_busco.sh' '${ROOT}/asm/${SP}/${name}/assembly.fasta' '${SP}_${name}' 16" \
            > /dev/null
    done

    echo "[${SP}] theirs_correction=${JOBA} mdax_correction=${JOBB} flye_theirs=${JOBC1} flye_mdax=${JOBC2}"
done
