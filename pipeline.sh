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
grep -v '#' $input_vcf | cut -f 1 | uniq > vcf_chroms.txt
samtools view $input_file | cut -f 3 | uniq > input_chroms.txt

# reorder header for sorting later and so chromosome order goes 1, 2, 3 not chr1, chr10, chr11
if $remove_chr; then
	grep "@HD" $input_file > filtered.sam
	grep "@SQ" $input_file | grep "chr*" > sq.tmp
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
	rm store.tmp
	rm sq.tmp
	sed -e 's/  */\t/g' reorder.tmp > tab.tmp
	cat tab.tmp >> filtered.sam
	rm reorder.tmp
	rm tab.tmp
	grep "@PG" $input_file >> filtered.sam
	grep "@CO" $input_file >> filtered.sam
fi

# remove any contigs present in input file that are not named chr
if $remove_chr; then
  cat input_chroms.txt | while read i; do
    if [[ $i == chr* ]]; then
      grep -v "@" $input_file | grep $i >> filtered.sam
    fi;
  done
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

# sorting BAM file based on new chromosome order
if $remove_chr; then
  samtools sort -o sorted.sam altered.sam
fi

# removing temporary files
if $remove_tmp; then
  if [[ $input_sam == false ]]; then
    rm tmp.sam
  fi
  rm filtered.sam
  if $remove_chr; then
    rm altered.sam
  fi
fi
