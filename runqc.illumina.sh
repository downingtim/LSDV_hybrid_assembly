#!/bin/bash

file1=$1
file2=$1
trim1=$1
trim2=$1
fastp1=$1
fastp2=$1
fastphtml1=$1
fastphtml2=$1
stats1=$1
stats2=$1
file1+="_1.fastq.gz"
file2+="_2.fastq.gz"
trim1+="_1_fastp_trim.fastq"
trim2+="_2_fastp_trim.fastq"
fastp1+="_1_fastp.fastq"
fastp2+="_2_fastp.fastq"
fastphtml1+="_1_fastp.html"
fastphtml2+="_2_fastp.html"
stats1+="_1.trim.txt"
stats2+="_2.trim.txt"
echo $file1 $file2 ;

if test -f $file1; then
    fastqc $file1 &
    fastqc $file2 &
    echo "fastqc $file1 "

    ~/FASTQC/fastp --overrepresentation_analysis  --html $fastphtml1 -i $file1 -o $fastp1
    echo " ~/FASTQC/fastp --overrepresentation_analysis  --html $fastphtml1 -i $file1 -o $fastp1"
    ~/FASTQC/fastp --overrepresentation_analysis  --html $fastphtml2 -i $file2 -o $fastp2 
    
    ~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 70 -i $fastp1 -o $trim1 -v  $stats1
    echo "~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 70 -i $fastp1 -o $trim1 -v  $stats1"
    ~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 70 -i $fastp2 -o $trim2 -v  $stats2

    fastqc $trim1 &
    echo "fastqc $trim1"
    fastqc $trim2
fi

echo $file1 $file2 $trim1 $trim2 $fastp1 $fastp2 $fastphtml1 $fastphtml2
