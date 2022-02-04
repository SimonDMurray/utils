#!/bin/bash

set -euo pipefail

input_file=""
output_file="out.bam"
input_vcf=""
input_fasta=""
input_sam=false
input_bam=true
input_cram=false
remove_chr=true
remove_tmp=true

while getopts :i:o:v:f:scrkh opt
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
    s)
      input_sam=true
      input_bam=false
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
-s input file type is sam
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

#convert SAM to BAM
if $input_sam; then
  echo "input is SAM file"
  samtools view -b -o tmp.bam $input_file
  input_file=tmp.bam
fi

#conver CRAM to SAM
if $input_cram; then
  echo "input is CRAM file (default)"
  samtools view -b -T $input_fasta -o tmp.bam $input_file
  input_file=tmp.bam
fi

if $input_bam; then
  echo "input is BAM file (default)"
fi

# get chromosome nomenclature of vcf and input file
grep -v '#' $input_vcf | cut -f 1 | uniq > vcf_chroms.txt
samtools view -H $input_file | grep @SQ | cut -f 2 | cut -c 4- > input_chroms.txt

# reorder header for sorting later and so chromosome order goes 1, 2, 3 not chr1, chr10, chr11
if $remove_chr; then
  samtools view -H $input_file | grep @HD > filtered.sam
  samtools view -H $input_file | grep @SQ | grep "chr*" > sq.tmp
  > reorder.tmp
  > store.tmp
  cat sq.tmp | while read i; do
  if [[ $i == *chr[1,2][0-9]* ]]; then
    echo $i >> store.tmp
    elif [[ $i == *chr[A-Z]* ]]; then
      echo $i >> store.tmp
    else
      echo $i >> reorder.tmp
    fi
  done
  cat store.tmp >> reorder.tmp
  sed -e 's/  */\t/g' reorder.tmp > tab.tmp
  cat tab.tmp >> filtered.sam
  rm sq.tmp
  rm reorder.tmp
  rm store.tmp
  rm tab.tmp
  samtools view -H $input_file | grep @PG >> filtered.sam
  samtools view -H $input_file | grep @CO >> filtered.sam
fi

# remove any contigs present in input file that are not named chr
if $remove_chr; then
  samtools index $input_file
  chroms=$(grep "chr" input_chroms.txt)
  samtools view $input_file $chroms >> filtered.sam 
else
  cat $input_file > filtered.sam
fi

# removing or adding chr nomenclature to chromosomes
if $remove_chr; then
  sed s/chr/""/g filtered.sam > altered.sam
else
  samtools view -H filtered.sam | grep "@HD" > altered.sam
  samtools view -H filtered.sam | grep "@SQ" | sed -r -e 's/^.{7}/&chr/' >> altered.sam
  samtools view -H filtered.sam | grep "@PG" >> altered.sam
  samtools view -H filtered.sam | grep "@CO" >> altered.sam
  samtools view filtered.sam | awk '$3="chr"$3' | sed -e 's/  */\t/g' >> altered.sam
fi

# sorting SAM file based on new chromosome order
if $remove_chr; then
  samtools sort -o sorted.sam altered.sam
fi

# convertin SAM to BAM for output
if $remove_chr; then
  samtools view -b -o out.bam sorted.sam
else
  samtools view -b -o out.bam altered.sam
fi

# removing temporary files
if $remove_tmp; then
  if [[ $input_bam == false ]]; then
    rm tmp.bam
  fi
  rm filtered.sam
  rm altered.sam
  if $remove_chr; then
    rm $input_file.bai
    rm sorted.sam
  fi
fi
