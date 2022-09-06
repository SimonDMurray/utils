#!/bin/bash

##Generic submission script that illustrates how to submit, taking a script to be ran and a sample to run the script on as well as all additional arguments

set -euo pipefail

#User inputs
SCRIPT=$1
SAMPLE=$2
#All arguments from third onwards are put into the variable ARG
ARGS=${@:3}

#BSUB ARGUMENTS

CPUS=16
MEM=4000
QUE="long"
GROUP="cellgeni"
IMAGE="/path/to/singularity/image.sif"


###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

bsub -G $GROUP -n $CPU -R"span[hosts=1] select[mem>${MEM}] rusage[mem=${MEM}]" -M $MEM -o $WDIR/$SAMPLE-%J.bsub.log -e $WDIR/$SAMPLE-%J.bsub.err -q $QUE \
  singularity exec -B /lustre,/nfs $IMAGE $SCRIPT $SAMPLE $ARG
