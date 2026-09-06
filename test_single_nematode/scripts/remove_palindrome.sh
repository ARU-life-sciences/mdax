#!/bin/bash

SAMPLE=$1 ### fastq_file
ID=$2 ### sample_ID 
CPU=$3

minimap2 -t $CPU -x ava-ont $SAMPLE $SAMPLE | cat |pafIdentifyPalimdrom.pl 1> $ID.palimProp.1stIte.list 2> $ID.1stIte.log
fastq_partition.and.chop.palindrome.pl $ID.palimProp.1stIte.list $SAMPLE 1000
minimap2 -t $CPU -x ava-ont $SAMPLE.exclude.fq.gz $SAMPLE.exclude.fq.gz | cat | pafIdentifyPalimdrom.pl 1> $ID.palimProp.2ndIte.list 2> $ID.2ndIte.log
fastq_partition.and.chop.palindrome.pl $ID.palimProp.2ndIte.list $SAMPLE.exclude.fq.gz 1000
cat $SAMPLE.include.fq.gz $SAMPLE.exclude.fq.gz.exclude.fq.gz $SAMPLE.exclude.fq.gz.include.fq.gz > $ID.palindrome_treated.fq.gz