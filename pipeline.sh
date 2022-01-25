#!/bin/bash

set -euo pipefail

input_file=""
output_file="out.bam"
input_vcf=""
input_fasta=""
input_sam=true
input_bam=false
input_cram=false
remove_chr=true
remove_tmp=true

while getopts :i:o:v:f:bcrkh opt
do
    case "$opt" in
    i)
      input_file=$OPTARG
      ;;
    o)
      output_file=$OPTARG
      ;;
    v)
      input_vcf=$OPTARG
      ;;
    f)
      input_fasta=$OPTARG
      ;;
    b)
      input_sam=false
      input_bam=true
      input_cram=false
      ;; 
    c)
      input_sam=false
      input_bam=false
      input_cram=true
      ;;
    r)
      remove_chr=false
      ;;
    k)
      remove_tmp=false
      ;;
    h)
      cat <<EOU
-i path to input file
-o path to output location (default CWD)
-v path to input vcf for chrom comparison
-f path to reference fasta file (needed with CRAM)
-b input file type is bam
-c input file type is cram
-r add "chr" to file (default is to remove)
-k does not remove temporary bam files (default is to remove)
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

if [[  $input_file == "" ]]; then
   echo "Need input file! (see -h)"
   false
fi

if [[  $input_vcf == "" ]]; then
   echo "Need input vcf! (see -h)"
   false
fi

if $input_cram && [[ $input_fasta == "" ]]; then
   echo "Need input fasta! (see -h)"
   false
fi

#convert BAM to SAM
if $input_bam; then
  echo "input is BAM file"
  samtools view -h -o tmp.sam $input_file
  input_file=tmp.sam
fi

#conver CRAM to SAM
if $input_cram; then
  echo "input is CRAM file (default)"
  samtools view -h -T $input_fasta -o tmp.sam $input_file
  input_file=tmp.sam
fi

if $input_sam; then
  echo "input is SAM file (default)"
fi

# get chromosome nomenclature of vcf and input file
grep -v '##' $input_vcf | cut -f 1 | uniq > vcf_chroms.txt
samtools view $input_file | cut -f 3 | uniq > input_chroms.txt

# remove any contigs present in input file that are not named chr
if $remove_chr; then
  grep "@" $input_file > filtered.sam 
    cat input_chroms.txt | while read i; do
    if [[ $i == chr* ]]; then
      grep $i $input_file >> filtered.sam
    fi;
  done
else
  cat $input_file > filtered.sam
fi

if $remove_chr; then
  sed s/chr/""/g filtered.sam > altered.sam
else
  awk '$3="chr"$3' filtered.sam > altered.sam
samtools sort -o sorted.sam altered.sam

if $remove_tmp; then
  if !$input_sam; then
    rm tmp.sam
  fi
  rm filtered.sam
  rm altered.sam
fi
