#!/bin/bash

set -euo pipefail

BAM=$1
OUTPUT="${2:-bam-to-fastq-output}"

mkdir -p logs

bsub -n 1 -R"span[hosts=1]" -M 32000 -R"select[mem>32000] rusage[mem=32000]" -G cellgeni -q long -o logs/ooo.$BAM.%J.txt -e logs/eee.$BAM.%J.txt \
  /software/cellgeni/bamtofastq $BAM $OUTPUT
