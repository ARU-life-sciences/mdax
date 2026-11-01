#!/usr/bin/env bash
# Submit the full 9-job LSF chain (correction x2 + Flye x3 + BUSCO x3) for one
# species' ONT read set. Mirrors the manually-run Epsilonema pilot, but as a
# reusable driver so scaling to more species doesn't mean hand-crafting bsub
# calls each time.
set -euo pipefail

SP="$1"           # species slug, e.g. "microlaimidae"
FASTQ="$2"        # path to raw ONT fastq.gz
DLJOB="${3:-}"    # optional: LSF job ID of an upstream download job FASTQ depends on
ROOT="/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_single_nematode_benchmark"
SCRIPTS="${ROOT}/scripts"
LOGDIR="${ROOT}/logs"

# When reads aren't staged locally yet, chain every job that reads FASTQ
# directly (both correction arms, plus raw Flye) onto the download job.
# Use ended(), not done(): done() is only satisfied on a *successful*
# download, so if it fails, a done() dependency would never be satisfied and
# these three jobs (plus everything chained behind them) would sit PEND
# forever with no visible error. ended() releases them either way, so a
# failed download still shows up the way it always did before this
# dependency existed: the job starts, can't read the missing/partial FASTQ,
# and fails loudly and immediately instead of hanging silently.
if [[ -n "${DLJOB}" ]]; then
    DLDEP="-w ended(${DLJOB})"
else
    DLDEP=""
fi

mkdir -p "${ROOT}/corrected/${SP}/theirs" "${ROOT}/corrected/${SP}/mdax" \
         "${ROOT}/asm/${SP}/raw" "${ROOT}/asm/${SP}/theirs" "${ROOT}/asm/${SP}/mdax" \
         "${ROOT}/busco"

# --- correction stage ---
JOBA=$(cd "${ROOT}/corrected/${SP}/theirs" && bsub \
    -J "${SP}_theirs" ${DLDEP} -q normal -n 16 \
    -R "select[mem>=64000] span[hosts=1] rusage[mem=64000]" -M 64000 \
    -o "${LOGDIR}/${SP}_theirs.out" -e "${LOGDIR}/${SP}_theirs.err" \
    "${SCRIPTS}/01_palindrome_removal.sh" "${FASTQ}" "${SP}" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

JOBB=$(bsub \
    -J "${SP}_mdax" ${DLDEP} -q normal -n 16 \
    -R "select[mem>=32000] span[hosts=1] rusage[mem=32000]" -M 32000 \
    -o "${LOGDIR}/${SP}_mdax.out" -e "${LOGDIR}/${SP}_mdax.err" \
    "${SCRIPTS}/02_run_mdax.sh" "${FASTQ}" "${ROOT}/corrected/${SP}/mdax/${SP}_mdax" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

# --- Flye stage (raw depends only on the optional download job, ${DLDEP};
# theirs/mdax additionally depend on their own correction job) ---
# NOTE: `long` queue (48h limit), not `normal` (12h) -- Flye on ~6-15GB of ONT
# WGA reads routinely needs more than 12h even at 16 threads (confirmed: every
# -n 8 Flye run in this project's normal-queue attempts hit TERM_RUNLIMIT
# without finishing). Don't move this back to `normal` without re-verifying
# runtime headroom.
JOBC0=$(bsub \
    -J "${SP}_flye_raw" ${DLDEP} -q long -n 16 \
    -R "select[mem>=64000] span[hosts=1] rusage[mem=64000]" -M 64000 \
    -o "${LOGDIR}/${SP}_flye_raw.out" -e "${LOGDIR}/${SP}_flye_raw.err" \
    "${SCRIPTS}/03_flye.sh" "${FASTQ}" "${ROOT}/asm/${SP}/raw" 16 2>&1 | grep -oP '(?<=Job <)[0-9]+')

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

# --- BUSCO stage, each depending on its Flye job ---
for pair in "raw:${JOBC0}" "theirs:${JOBC1}" "mdax:${JOBC2}"; do
    name="${pair%%:*}"
    dep="${pair##*:}"
    bsub \
        -J "${SP}_busco_${name}" -w "done(${dep})" -q long -n 16 \
        -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
        -o "${LOGDIR}/${SP}_busco_${name}.out" -e "${LOGDIR}/${SP}_busco_${name}.err" \
        bash -c "cd '${ROOT}/busco' && '${SCRIPTS}/04_busco.sh' '${ROOT}/asm/${SP}/${name}/assembly.fasta' '${SP}_${name}' 16" \
        > /dev/null
done

echo "[${SP}] submitted: theirs_correction=${JOBA} mdax_correction=${JOBB} flye_raw=${JOBC0} flye_theirs=${JOBC1} flye_mdax=${JOBC2}"
