#!/bin/bash

# a bash script that makes an input for QuASAR analysis

set -o errexit
set -o nounset

BAMF=$1

GENOMEF=~/genomes/hg19/hg19all.fa

SNPF="1KG_SNPs_filt.bed"
SNPID=${SNPF%.bed}

MPILEUPF=${BAMF%.bam}.mpileup.txt.gz
MPILEUPF=`basename $MPILEUPF`
MPILEUPBEDF=${MPILEUPF%.txt.gz}.bed.gz

samtools mpileup -f $GENOMEF -l $SNPF $BAMF |gzip -c >$MPILEUPF

less ${MPILEUPF} | awk -v OFS='\t' '{ if ($4>0 && $5 !~ /[^\^][<>]/ && $5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $5 !~ /-[0-9]+[ACGTNacgtn]+/ && $5 !~ /[^\^]\*/) print $1,$2-1,$2,$3,$4,$5,$6}'| sortBed -i stdin | intersectBed -a stdin -b $SNPF -wo|cut -f 1-7,11-14 | gzip -c >${MPILEUPBEDF}

R --vanilla --args ${MPILEUPBEDF} < convertPileupToQuasar2.R
