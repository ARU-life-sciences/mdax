#!/usr/bin/env bash
# Resubmit only the mdax-dependent arms (strict/balanced/permissive) for every
# C. elegans rerun sample, skipping raw -- raw's original jobs are still
# alive (running or pending) and don't depend on the mdax binary, so they
# were deliberately left untouched when the mdax-dependent jobs were killed
# and are not resubmitted here to avoid duplicate/colliding output.
set -euo pipefail
ROOT="/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_celegans_rerun"
SCRIPTS="${ROOT}/scripts"
LOGDIR="${ROOT}/logs"
MANIFEST="${ROOT}/manifest.tsv"

while IFS=$'\t' read -r SAMPLE RAW TIER; do
    case "$TIER" in
        small)  MDAX_MEM=32000;  MDAX_CPU=16; HIFI_MEM=65000;  HIFI_CPU=32 ;;
        medium) MDAX_MEM=65000;  MDAX_CPU=16; HIFI_MEM=131000; HIFI_CPU=32 ;;
        large)  MDAX_MEM=98000;  MDAX_CPU=16; HIFI_MEM=196000; HIFI_CPU=32 ;;
        *) echo "unknown tier $TIER for $SAMPLE" >&2; exit 1 ;;
    esac
    BUSCO_MEM=24000; BUSCO_CPU=16

    mkdir -p "${ROOT}/corrected/${SAMPLE}"
    mkdir -p "${ROOT}/asm/${SAMPLE}/strict" "${ROOT}/asm/${SAMPLE}/balanced" "${ROOT}/asm/${SAMPLE}/permissive"
    mkdir -p "${ROOT}/busco"

    declare -A MDAXJOB
    for MODE in strict balanced permissive; do
        OUT="${ROOT}/corrected/${SAMPLE}/${SAMPLE}_${MODE}"
        JID=$(bsub \
            -J "${SAMPLE}_mdax_${MODE}" -q normal -n "${MDAX_CPU}" \
            -R "select[mem>=${MDAX_MEM}] span[hosts=1] rusage[mem=${MDAX_MEM}]" -M "${MDAX_MEM}" \
            -o "${LOGDIR}/${SAMPLE}_mdax_${MODE}.out" -e "${LOGDIR}/${SAMPLE}_mdax_${MODE}.err" \
            "${SCRIPTS}/01_mdax_correct.sh" "${RAW}" "${OUT}" "${MODE}" "${MDAX_CPU}" \
            2>&1 | grep -oP '(?<=Job <)[0-9]+')
        MDAXJOB[$MODE]=$JID
    done

    declare -A HIFIJOB
    for MODE in strict balanced permissive; do
        JID=$(bsub \
            -J "${SAMPLE}_hifiasm_${MODE}" -w "done(${MDAXJOB[$MODE]})" -q basement -n "${HIFI_CPU}" \
            -R "select[mem>=${HIFI_MEM}] span[hosts=1] rusage[mem=${HIFI_MEM}]" -M "${HIFI_MEM}" \
            -o "${LOGDIR}/${SAMPLE}_hifiasm_${MODE}.out" -e "${LOGDIR}/${SAMPLE}_hifiasm_${MODE}.err" \
            "${SCRIPTS}/02_hifiasm.sh" "${ROOT}/corrected/${SAMPLE}/${SAMPLE}_${MODE}.fa" "${ROOT}/asm/${SAMPLE}/${MODE}/${SAMPLE}_${MODE}" "${HIFI_CPU}" \
            2>&1 | grep -oP '(?<=Job <)[0-9]+')
        HIFIJOB[$MODE]=$JID
    done

    declare -A BUSCOJOB
    for ARM in strict balanced permissive; do
        JID=$(bsub \
            -J "${SAMPLE}_busco_${ARM}" -w "done(${HIFIJOB[$ARM]})" -q long -n "${BUSCO_CPU}" \
            -R "select[mem>=${BUSCO_MEM}] span[hosts=1] rusage[mem=${BUSCO_MEM}]" -M "${BUSCO_MEM}" \
            -o "${LOGDIR}/${SAMPLE}_busco_${ARM}.out" -e "${LOGDIR}/${SAMPLE}_busco_${ARM}.err" \
            bash -c "cd '${ROOT}/busco' && '${SCRIPTS}/03_busco.sh' '${ROOT}/asm/${SAMPLE}/${ARM}/${SAMPLE}_${ARM}.p_ctg.fa' '${SAMPLE}_${ARM}' ${BUSCO_CPU}" \
            2>&1 | grep -oP '(?<=Job <)[0-9]+')
        BUSCOJOB[$ARM]=$JID
    done

    echo "[${SAMPLE}] mdax: strict=${MDAXJOB[strict]} balanced=${MDAXJOB[balanced]} permissive=${MDAXJOB[permissive]} | hifiasm: strict=${HIFIJOB[strict]} balanced=${HIFIJOB[balanced]} permissive=${HIFIJOB[permissive]} | busco: strict=${BUSCOJOB[strict]} balanced=${BUSCOJOB[balanced]} permissive=${BUSCOJOB[permissive]}"
done < "${MANIFEST}"
