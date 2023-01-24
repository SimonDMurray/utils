#!/bin/bash

set -euo pipefail

#This script enables dynamic job array submission based on the number of samples submitted
#Assumes each line of sample file is a different sample
#Depending on inputs to script a second script may need to be submitted to utilise $LSB_JOBINDEX variable
#(i.e. if the input is only a single sample from the file and the line of the sample is equal to number in $LSB_JOBINDEX 
#then a second script is needed as $LSB_JOBINDEX won't be initialised until after the job is submitted

SCRIPT="${1:?/path/to/work/dir}"
OUTDIR="${2:?/path/to/out/dir"}
SAMPLEFILE="${3:?/path/to/sample/file"}

GROUP="cellgeni"
CPUS=16
MEM=128000
QUEUE="long"

#### DO NOT CHANGE BELOW THIS LINE ####

FILELENGTH=`wc -l $SAMPLEFILE | cut -f 1 -d " "`

bsub \
  -G ${GROUP} \
  -J "job[1-${FILELENGTH}]" \
  -n ${CPUS} \
  -R "span[hosts=1] select[mem>${MEM}] rusage[mem=${MEM}]" \
  -M ${MEM} \
  -o ${OUTDIR}/%J.%I.bsub.log \
  -e ${OUTDIR}/%J.%I.bsub.err \
  -q ${QUEUE} \
  $SCRIPT $OUTDIR $SAMPLEFILE
