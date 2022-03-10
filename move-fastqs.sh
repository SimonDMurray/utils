#!/bin/bash

data_dir=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1313/data5
fastq_dir=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1313/fastqs5

cd $data_dir
for i in *;
  do cd $i
  for j in *fastq*;
    do k=$i'_'$j
    mv $data_dir/$i/$j $fastq_dir/$k ;
  done;
  cd ../
done
