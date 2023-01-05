#!/bin/bash

set -euo pipefail

# Different versions of cellranger label chromosomes differently
# Some label chr1, chr2, chr3 etc
# Some label 1, 2, 3 etc
# This script lets you alter a vcf file between the two

input_vcf=""
output_vcf="out.vcf"
remove_chr=false

while getopts :i:o:ch opt
do
    case "$opt" in
    i)
      input_vcf=$OPTARG
      ;;
    o)
      output_vcf=$OPTARG
      ;;
    c)
      remove_chr=true
      ;;
    h)
      cat <<EOU
-i path to input vcf
-o path to output location (default CWD)
-c choose to remove chr instead of add it
-h displays this message!
EOU
      exit
            ;;
    :) echo "Flag $OPTARG needs argument"
        exit 1;;
    ?) echo "Flag $OPTARG unknown"              
        exit 1;;
   esac
done

if [[  $input_vcf == "" ]]; then
   echo "Need input vcf! (see -h)"
   false
fi

grep "#" $input_vcf > $output_vcf

if [[  $remove_chr == true ]]; then
  grep -v "#" $input_vcf | sed 's/^chr//' >> $output_vcf
else
  grep -v "#" $input_vcf | sed 's/^/chr/' >> $output_vcf
fi
