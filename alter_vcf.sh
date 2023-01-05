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

if [[  $remove_chr == true ]]; then
  #remove any occurrence of chr in the header
  grep "##" $input_vcf | sed 's/chr//' > $output_vcf
  #add field identifier line
  grep "#CHROM" $input_vcf >> $output_vcf
  #remove chr from the beginning of each sample line
  grep -v "#" $input_vcf | sed 's/^chr//' >> $output_vcf
else
  #add chr to the beginning of any chromosomes starting with a number or X,Y,Z in the header 
  grep "##" $input_vcf | sed -e '/ID=[0-9]/ s/ID=/ID=chr/ ; /ID=[0-9]/! s///' | sed -e '/ID=[M,X,Y]/ s/ID=/ID=chr/ ; /ID=[M,X,Y]/! s///' > $output_vcf
  #add field identifier line
  grep "#CHROM" $input_vcf >> $output_vcf
  #add chr to any sample line starting with a number or X,Y,Z
  grep -v "#" $input_vcf | egrep "^[0-9]|^[M,X,Y]" | sed 's/^/chr/' >> $output_vcf
  #add any sample lines not starting with a number or X,Y,Z (contigs)
  grep -v "#" $input_vcf | egrep -v "^[0-9]|^[M,X,Y]" >> $output_vcf
fi
