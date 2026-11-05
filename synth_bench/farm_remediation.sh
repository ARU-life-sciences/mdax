#!/usr/bin/env bash
# Remediate the farm batch after the 2026-11-05 audit: every failure was a
# resource-sizing issue (memory/runtime limits too small), not an mdax
# correctness problem, except celegans_wga whose raw read download was
# corrupted (now re-fetched and verified) and needs both correction steps
# (theirs + mdax) rerun from scratch since their outputs were built on the
# bad data.
set -euo pipefail

NEMROOT=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_single_nematode_benchmark
CELROOT=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_celegans_rerun
NEMSCRIPTS="${NEMROOT}/scripts"
CELSCRIPTS="${CELROOT}/scripts"
NEMLOG="${NEMROOT}/logs"
CELLOG="${CELROOT}/logs"

jid() { grep -oP '(?<=Job <)[0-9]+'; }

##############################################################################
# Part A: celegans_wga -- raw read corruption fixed; rerun both correction
# arms from scratch, then their downstream flye+busco. "raw" arm never
# actually ran (was stuck PEND the whole time), so it just needs a fresh
# submission against the now-fixed fastq.
##############################################################################
SP=celegans_wga
FASTQ="${NEMROOT}/reads/${SP}.fastq.gz"

JOB_FLYE_RAW=$(bsub -J "${SP}_flye_raw" -q long -n 16 \
    -R "select[mem>=64000] span[hosts=1] rusage[mem=64000]" -M 64000 \
    -o "${NEMLOG}/${SP}_flye_raw.out" -e "${NEMLOG}/${SP}_flye_raw.err" \
    "${NEMSCRIPTS}/03_flye.sh" "${FASTQ}" "${NEMROOT}/asm/${SP}/raw" 16 | jid)
bsub -J "${SP}_busco_raw" -w "done(${JOB_FLYE_RAW})" -q long -n 16 \
    -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
    -o "${NEMLOG}/${SP}_busco_raw.out" -e "${NEMLOG}/${SP}_busco_raw.err" \
    bash -c "cd '${NEMROOT}/busco' && '${NEMSCRIPTS}/04_busco.sh' '${NEMROOT}/asm/${SP}/raw/assembly.fasta' '${SP}_raw' 16" >/dev/null

JOB_THEIRS=$(cd "${NEMROOT}/corrected/${SP}/theirs" && bsub -J "${SP}_theirs" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${NEMLOG}/${SP}_theirs.out" -e "${NEMLOG}/${SP}_theirs.err" \
    "${NEMSCRIPTS}/01_palindrome_removal.sh" "${FASTQ}" "${SP}" 16 | jid)
JOB_FLYE_THEIRS=$(bsub -J "${SP}_flye_theirs" -w "done(${JOB_THEIRS})" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${NEMLOG}/${SP}_flye_theirs.out" -e "${NEMLOG}/${SP}_flye_theirs.err" \
    "${NEMSCRIPTS}/03_flye.sh" "${NEMROOT}/corrected/${SP}/theirs/${SP}.palindrome_treated.fq.gz" "${NEMROOT}/asm/${SP}/theirs" 16 | jid)
bsub -J "${SP}_busco_theirs" -w "done(${JOB_FLYE_THEIRS})" -q long -n 16 \
    -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
    -o "${NEMLOG}/${SP}_busco_theirs.out" -e "${NEMLOG}/${SP}_busco_theirs.err" \
    bash -c "cd '${NEMROOT}/busco' && '${NEMSCRIPTS}/04_busco.sh' '${NEMROOT}/asm/${SP}/theirs/assembly.fasta' '${SP}_theirs' 16" >/dev/null

JOB_MDAX=$(bsub -J "${SP}_mdax" -q normal -n 16 \
    -R "select[mem>=32000] span[hosts=1] rusage[mem=32000]" -M 32000 \
    -o "${NEMLOG}/${SP}_mdax.out" -e "${NEMLOG}/${SP}_mdax.err" \
    "${NEMSCRIPTS}/02_run_mdax.sh" "${FASTQ}" "${NEMROOT}/corrected/${SP}/mdax/${SP}_mdax" 16 | jid)
JOB_FLYE_MDAX=$(bsub -J "${SP}_flye_mdax" -w "done(${JOB_MDAX})" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${NEMLOG}/${SP}_flye_mdax.out" -e "${NEMLOG}/${SP}_flye_mdax.err" \
    "${NEMSCRIPTS}/03_flye.sh" "${NEMROOT}/corrected/${SP}/mdax/${SP}_mdax.fa" "${NEMROOT}/asm/${SP}/mdax" 16 | jid)
bsub -J "${SP}_busco_mdax" -w "done(${JOB_FLYE_MDAX})" -q long -n 16 \
    -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
    -o "${NEMLOG}/${SP}_busco_mdax.out" -e "${NEMLOG}/${SP}_busco_mdax.err" \
    bash -c "cd '${NEMROOT}/busco' && '${NEMSCRIPTS}/04_busco.sh' '${NEMROOT}/asm/${SP}/mdax/assembly.fasta' '${SP}_mdax' 16" >/dev/null

echo "[${SP}] flye_raw=${JOB_FLYE_RAW} theirs=${JOB_THEIRS} flye_theirs=${JOB_FLYE_THEIRS} mdax=${JOB_MDAX} flye_mdax=${JOB_FLYE_MDAX}"

##############################################################################
# Part B: nrCaeEleg93/94/95 hifiasm hit TERM_MEMLIMIT at 131GB (medium tier)
# across all 3 modes each. Corrected .fa inputs are fine (correction itself
# succeeded) -- just resubmit hifiasm with much more memory, then busco.
##############################################################################
for SAMPLE in nrCaeEleg93 nrCaeEleg94 nrCaeEleg95; do
    for MODE in strict balanced permissive; do
        JOB_HIFI=$(bsub -J "${SAMPLE}_hifiasm_${MODE}_retry" -q basement -n 32 \
            -R "select[mem>=250000] span[hosts=1] rusage[mem=250000]" -M 250000 \
            -o "${CELLOG}/${SAMPLE}_hifiasm_${MODE}_retry.out" -e "${CELLOG}/${SAMPLE}_hifiasm_${MODE}_retry.err" \
            "${CELSCRIPTS}/02_hifiasm.sh" "${CELROOT}/corrected/${SAMPLE}/${SAMPLE}_${MODE}.fa" "${CELROOT}/asm/${SAMPLE}/${MODE}/${SAMPLE}_${MODE}" 32 | jid)
        bsub -J "${SAMPLE}_busco_${MODE}_retry" -w "done(${JOB_HIFI})" -q long -n 16 \
            -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
            -o "${CELLOG}/${SAMPLE}_busco_${MODE}_retry.out" -e "${CELLOG}/${SAMPLE}_busco_${MODE}_retry.err" \
            bash -c "cd '${CELROOT}/busco' && '${CELSCRIPTS}/03_busco.sh' '${CELROOT}/asm/${SAMPLE}/${MODE}/${SAMPLE}_${MODE}.p_ctg.fa' '${SAMPLE}_${MODE}' 16" >/dev/null
        echo "[${SAMPLE}/${MODE}] hifiasm=${JOB_HIFI}"
    done
done

##############################################################################
# Part C: nrCaeEleg33 permissive hifiasm aborted on a malloc_consolidate
# heap-corruption crash -- not a resource limit. Retry as-is first.
##############################################################################
JOB_HIFI33=$(bsub -J "nrCaeEleg33_hifiasm_permissive_retry" -q basement -n 32 \
    -R "select[mem>=196000] span[hosts=1] rusage[mem=196000]" -M 196000 \
    -o "${CELLOG}/nrCaeEleg33_hifiasm_permissive_retry.out" -e "${CELLOG}/nrCaeEleg33_hifiasm_permissive_retry.err" \
    "${CELSCRIPTS}/02_hifiasm.sh" "${CELROOT}/corrected/nrCaeEleg33/nrCaeEleg33_permissive.fa" "${CELROOT}/asm/nrCaeEleg33/permissive/nrCaeEleg33_permissive" 32 | jid)
bsub -J "nrCaeEleg33_busco_permissive_retry" -w "done(${JOB_HIFI33})" -q long -n 16 \
    -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
    -o "${CELLOG}/nrCaeEleg33_busco_permissive_retry.out" -e "${CELLOG}/nrCaeEleg33_busco_permissive_retry.err" \
    bash -c "cd '${CELROOT}/busco' && '${CELSCRIPTS}/03_busco.sh' '${CELROOT}/asm/nrCaeEleg33/permissive/nrCaeEleg33_permissive.p_ctg.fa' 'nrCaeEleg33_permissive' 16" >/dev/null
echo "[nrCaeEleg33/permissive retry] hifiasm=${JOB_HIFI33}"

##############################################################################
# Part D: 6 nematode-suite "theirs" corrections hit TERM_MEMLIMIT (64GB) or
# TERM_RUNLIMIT (normal queue's 12h cap) -- resubmit with more memory and on
# `long` (48h) instead of `normal`. Chain flye_theirs + busco_theirs after.
##############################################################################
for SP in mesodorylaimus ptycholaimellus sabatieria theristus trissonchulus_lat trissonchulus_wlg; do
    FASTQ="${NEMROOT}/reads/${SP}.fastq.gz"
    JOB_THEIRS=$(cd "${NEMROOT}/corrected/${SP}/theirs" && bsub -J "${SP}_theirs_retry" -q long -n 16 \
        -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
        -o "${NEMLOG}/${SP}_theirs_retry.out" -e "${NEMLOG}/${SP}_theirs_retry.err" \
        "${NEMSCRIPTS}/01_palindrome_removal.sh" "${FASTQ}" "${SP}" 16 | jid)
    JOB_FLYE=$(bsub -J "${SP}_flye_theirs_retry" -w "done(${JOB_THEIRS})" -q long -n 16 \
        -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
        -o "${NEMLOG}/${SP}_flye_theirs_retry.out" -e "${NEMLOG}/${SP}_flye_theirs_retry.err" \
        "${NEMSCRIPTS}/03_flye.sh" "${NEMROOT}/corrected/${SP}/theirs/${SP}.palindrome_treated.fq.gz" "${NEMROOT}/asm/${SP}/theirs" 16 | jid)
    bsub -J "${SP}_busco_theirs_retry" -w "done(${JOB_FLYE})" -q long -n 16 \
        -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
        -o "${NEMLOG}/${SP}_busco_theirs_retry.out" -e "${NEMLOG}/${SP}_busco_theirs_retry.err" \
        bash -c "cd '${NEMROOT}/busco' && '${NEMSCRIPTS}/04_busco.sh' '${NEMROOT}/asm/${SP}/theirs/assembly.fasta' '${SP}_theirs' 16" >/dev/null
    echo "[${SP} theirs retry] theirs=${JOB_THEIRS} flye=${JOB_FLYE}"
done

##############################################################################
# Part E: trissonchulus_lat flye_mdax hit TERM_MEMLIMIT (64GB). mdax
# correction itself succeeded; just resubmit flye with more memory.
##############################################################################
SP=trissonchulus_lat
JOB_FLYE=$(bsub -J "${SP}_flye_mdax_retry" -q long -n 16 \
    -R "select[mem>=128000] span[hosts=1] rusage[mem=128000]" -M 128000 \
    -o "${NEMLOG}/${SP}_flye_mdax_retry.out" -e "${NEMLOG}/${SP}_flye_mdax_retry.err" \
    "${NEMSCRIPTS}/03_flye.sh" "${NEMROOT}/corrected/${SP}/mdax/${SP}_mdax.fa" "${NEMROOT}/asm/${SP}/mdax" 16 | jid)
bsub -J "${SP}_busco_mdax_retry" -w "done(${JOB_FLYE})" -q long -n 16 \
    -R "select[mem>=16000] span[hosts=1] rusage[mem=16000]" -M 16000 \
    -o "${NEMLOG}/${SP}_busco_mdax_retry.out" -e "${NEMLOG}/${SP}_busco_mdax_retry.err" \
    bash -c "cd '${NEMROOT}/busco' && '${NEMSCRIPTS}/04_busco.sh' '${NEMROOT}/asm/${SP}/mdax/assembly.fasta' '${SP}_mdax' 16" >/dev/null
echo "[${SP} flye_mdax retry] flye=${JOB_FLYE}"
