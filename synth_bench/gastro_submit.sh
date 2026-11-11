#!/usr/bin/env bash
# Submit the raw/theirs/mdax x flye x busco chain for one Lepidodermella
# squamata (Gastrotricha) sample -- two biosamples (SRR27483171/72) treated
# as separate individuals, matching the nrCaeEleg92-100 convention (multiple
# individual WGA preps of the same species, analysed separately) rather than
# pooling them into one assembly.
set -euo pipefail
SP="$1"           # sample slug, e.g. "SRR27483171"
FASTQ="$2"        # path to raw subreads fastq.gz
ROOT=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_gastrotricha_benchmark
SCRIPTS="${ROOT}/scripts"
LOGDIR="${ROOT}/logs"

mkdir -p "${ROOT}/corrected/${SP}/theirs" "${ROOT}/corrected/${SP}/mdax" \
         "${ROOT}/asm/${SP}/raw" "${ROOT}/asm/${SP}/theirs" "${ROOT}/asm/${SP}/mdax" \
         "${ROOT}/busco"

jid() { grep -oP '(?<=Job <)[0-9]+'; }

JOBA=$(cd "${ROOT}/corrected/${SP}/theirs" && bsub -J "${SP}_theirs" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${LOGDIR}/${SP}_theirs.out" -e "${LOGDIR}/${SP}_theirs.err" \
    "${SCRIPTS}/01_palindrome_removal.sh" "${FASTQ}" "${SP}" 16 | jid)

JOBB=$(bsub -J "${SP}_mdax" -q normal -n 16 \
    -R "select[mem>=32000] span[hosts=1] rusage[mem=32000]" -M 32000 \
    -o "${LOGDIR}/${SP}_mdax.out" -e "${LOGDIR}/${SP}_mdax.err" \
    "${SCRIPTS}/02_run_mdax.sh" "${FASTQ}" "${ROOT}/corrected/${SP}/mdax/${SP}_mdax" 16 | jid)

JOBC0=$(bsub -J "${SP}_flye_raw" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${LOGDIR}/${SP}_flye_raw.out" -e "${LOGDIR}/${SP}_flye_raw.err" \
    "${SCRIPTS}/03_flye.sh" "${FASTQ}" "${ROOT}/asm/${SP}/raw" 16 | jid)

JOBC1=$(bsub -J "${SP}_flye_theirs" -w "done(${JOBA})" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${LOGDIR}/${SP}_flye_theirs.out" -e "${LOGDIR}/${SP}_flye_theirs.err" \
    "${SCRIPTS}/03_flye.sh" "${ROOT}/corrected/${SP}/theirs/${SP}.palindrome_treated.fq.gz" "${ROOT}/asm/${SP}/theirs" 16 | jid)

JOBC2=$(bsub -J "${SP}_flye_mdax" -w "done(${JOBB})" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${LOGDIR}/${SP}_flye_mdax.out" -e "${LOGDIR}/${SP}_flye_mdax.err" \
    "${SCRIPTS}/03_flye.sh" "${ROOT}/corrected/${SP}/mdax/${SP}_mdax.fa" "${ROOT}/asm/${SP}/mdax" 16 | jid)

for pair in "raw:${JOBC0}" "theirs:${JOBC1}" "mdax:${JOBC2}"; do
    name="${pair%%:*}"; dep="${pair##*:}"
    bsub -J "${SP}_busco_${name}" -w "done(${dep})" -q long -n 16 \
        -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
        -o "${LOGDIR}/${SP}_busco_${name}.out" -e "${LOGDIR}/${SP}_busco_${name}.err" \
        bash -c "cd '${ROOT}/busco' && '${SCRIPTS}/04_busco.sh' '${ROOT}/asm/${SP}/${name}/assembly.fasta' '${SP}_${name}' 16" >/dev/null
done

echo "[${SP}] theirs=${JOBA} mdax=${JOBB} flye_raw=${JOBC0} flye_theirs=${JOBC1} flye_mdax=${JOBC2}"
