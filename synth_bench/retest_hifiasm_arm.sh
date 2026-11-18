#!/usr/bin/env bash
# One arm end to end: mdax -> hifiasm -> BUSCO -> tidy. usage: arm.sh SAMPLE NAME MDAXBIN [mdax args...]
set -euo pipefail
source /etc/profile.d/modules.sh
B=/lustre/scratch122/tol/teams/blaxter/users/mb39/mdax_celegans_rerun
S=$1; NAME=$2; MDAX=$3; shift 3
IN=$(ls $B/reads/${S}_raw.fasta* | head -1)
W=$B/retest92/$S.$NAME; mkdir -p $W; cd $W
echo "[cmd] $MDAX $IN $* --cut-low-ident" > cmd.txt
$MDAX "$IN" "$@" --cut-low-ident -t 32 -o reads.fa -r report.tsv 2> mdax.log
/software/team301/hifiasm/hifiasm -o asm -t 32 reads.fa 2> hifiasm.log
/software/team301/gfatools/gfatools gfa2fa asm.bp.p_ctg.gfa > p_ctg.fa
rm -f reads.fa asm.*
module load busco/5.4.6--pyhdfd78af_0
busco -m genome -l nematoda_odb10 -i p_ctg.fa -o busco -c 32 -f > busco.log 2>&1
cp busco/short_summary*.txt busco_summary.txt; rm -rf busco busco_downloads
gzip -f report.tsv
