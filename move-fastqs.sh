#!/bin/bash

set -euo pipefail

data_dir=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-XXXX/data
fastq_dir=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-XXXX/fastqs

cd $data_dir
for sample in *; do 
  if [[ $sample == "library_info" ]]; then
    continue
  else
    cd $sample
    for fq in *fastq*; do 
      file="${sample}_${fq}"
      mv $data_dir/$sample/$fq $fastq_dir/$file 
    done
    cd ..
  fi
done
