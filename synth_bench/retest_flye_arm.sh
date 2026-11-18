#!/usr/bin/env bash
# mdax -> Flye -> BUSCO -> tidy, for noisy-read benchmarks.
# usage: flye_arm.sh NAME READS FLYE_READ_FLAG BUSCO_LINEAGE MDAXBIN [mdax args...]
set -euo pipefail
source /etc/profile.d/modules.sh
NAME=$1; READS=$2; FLAG=$3; LIN=$4; MDAX=$5; shift 5
W=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_retest_v4/$NAME; mkdir -p $W; cd $W
echo "[cmd] $MDAX $READS $*" > cmd.txt
$MDAX "$READS" "$@" -t 16 -o reads.fa -r report.tsv 2> mdax.log
gzip -f report.tsv
source /nfs/users/nfs_m/mb39/miniconda3/etc/profile.d/conda.sh; conda activate flye_env
flye $FLAG reads.fa --meta --out-dir flye --threads 16 > flye.log 2>&1
cp flye/assembly.fasta assembly.fasta; rm -rf flye reads.fa
conda deactivate
module load busco/5.4.6--pyhdfd78af_0
busco -m genome -l $LIN -i assembly.fasta -o busco -c 16 -f > busco.log 2>&1
cp busco/short_summary*.txt busco_summary.txt; rm -rf busco busco_downloads
