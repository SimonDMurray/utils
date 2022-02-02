#!/bin/bash

set -euo pipefail

# This script removes all contigs from a bam file
# If nomenclature is chr1, chr10, chr11 etc uncomment first samtools line
# If nomenclature is 1, 2,3 outcomment second samtools line
# Uncomment grep line always

SAMPLE=$1
FILE=$2
INFILE=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1280/results
OUTFILE=/lustre/scratch117/cellgen/cellgeni/TIC-10x/tic-1303/data

# SAMPLE is a sanger sample id i.e. FCA_GND10375779
# FILE is the bam file i.e. possorted_bam.bam
# Example:
# if a data directory cantains samples HCA1, HCA2, HCA3 then
# for i in *; do ../actions/remove_contigs.sh $i possorted_bam.bam
# this would run script 3 times
# ../actions/remove_contigs.sh HCA1 possorted_bam.bam
# ../actions/remove_contigs.sh HCA2 possorted_bam.bam
# ../actions/remove_contigs.sh HCA3 possorted_bam.bam

cd $INFILE
mkdir -p $OUTFILE/$SAMPLE
cd $SAMPLE
#samtools index $FILE
chrs=$(samtools view -H $FILE | grep chr | cut -f 2 | sed 's/[^:]*://')
samtools view -h -b $FILE $chrs > $OUTFILE/$SAMPLE/tmp 
samtools reheader -c 'grep -v SN:[A-Z]' $OUTFILE/$SAMPLE/tmp > $OUTFILE/$SAMPLE/removed.bam
rm $OUTFILE/$SAMPLE/tmp
cd ../
