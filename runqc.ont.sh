#!/bin/bash

file1=$1
trim1=$1
fastp1=$1
fastphtml1=$1
stats1=$1
file1+=".fastq"
trim1+="_fastp_trim.fastq"
fastp1+="_fastp.fastq"
fastphtml1+="_fastp.html"
stats1+=".trim.txt"
echo $file1 $file2 ;

if test -f $file1; then
  #  fastqc $file1 &
    echo "fastqc $file1 "

 #   ~/FASTQC/fastp --overrepresentation_analysis  --html $fastphtml1 -i $file1 -o $fastp1
    echo " ~/FASTQC/fastp --overrepresentation_analysis  --html $fastphtml1 -i $file1 -o $fastp1"
   /mnt/lustre/RDS-live/downing/.cargo/bin/nanoq -i $fastp1 -q 10 -l 100 -S 5 -E 5 -O u -vvv -r $stats1 -o $trim1

    fastqc $trim1 
    echo "fastqc $trim1"
fi

echo $file1  $trim1  $fastp1 $fastphtml1
