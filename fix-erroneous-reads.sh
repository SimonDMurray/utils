#!/bin/bash

set -euo pipefail

## A FEW NOTES:
## 1) Install bbmap tools with:
# `wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && tar -xvf bbmap.tar.gz && rm bbmap.tar.gz`
## 2) Ensure fastqs are gzipped and in cellranger format i.e. SAMPLEID_S1_L00*_R[1|2]_001.fastq.gz 

FQDIR=${1:?"Please provide path to directory containing fastqs"}
SAMPLEID=${2:?"Please provide the sampleID for the broken fastqs"}
MIN=${3:?"Please provide the minimum length a read needs to be to be kept in the fixed fastq"}
MAX=${4:?"Please provide the maximum length a read needs to be to be kept in the fixed fastq"}

cd $FQDIR

for i in "${SAMPLEID}"*R1*fastq.gz; do
  zcat $i | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= $MIN && length(seq) <= $MAX) {print header, seq, qheader, qseq}}' > filtered.fastq
  gzip filtered.fastq
  mv filtered.fastq.gz "filtered_${i}"
done

for i in "filtered_${SAMPLEID}"*R1*.fastq.gz; do
  name1=`echo $i | cut -f 2- -d _`
  name2=`echo $name1 | sed 's/R1/R2/g'`
  /path/to/bbmap/repair.sh in1=$i in2=$name2 out1="repair_filtered_${name1}" out2="repair_filtered_${name2}" outs="${SAMPLEID}_singletons.fastq" repair 1>bbduk.log 2>bbduk.err
done

rm "${SAMPLEID}"*R1*fastq.gz
rm filtered_${SAMPLEID}"*R1*.fastq.gz
