#!/bin/bash

file1="merge"
file2=$file1
file1+=".fastp.chop.fastq.gz"
sam=$file1
sam+=".sam"
fasta="pypolca_c.fa"

samtools faidx $fasta
minimap2 -ax map-ont $fasta $file1  > $sam
echo "minimap2 -ax lr:hq $fasta $file1  > $sam"
bam=$file1
bam+=".bam"
samtools view -bS $sam > $bam 
echo "samtools view -bS $sam > $bam "
sorted=$bam
sorted+=".sorted"
~/samtools-1.9/samtools sort -n $bam -o $sorted
echo "~/samtools-1.9/samtools sort -n $bam -o $sorted"
sorted2=$sorted
sorted2+=".2"
rmdup=$file1
rmdup+=".rmdup.bam"
~/samtools-1.9/samtools fixmate -m $sorted $sorted2
echo "~/samtools-1.9/samtools fixmate -m $sorted $sorted2"
~/samtools-1.9/samtools sort $sorted2 -o $sorted
echo "~/samtools-1.9/samtools sort $sorted2 -o $sorted"
~/samtools-1.9/samtools markdup -r -s $sorted $rmdup
echo "~/samtools-1.9/samtools markdup -r -s $sorted $rmdup"
~/samtools-1.9/samtools index $rmdup
echo "~/samtools-1.9/samtools index $rmdup"
rmdup2=$file1
rmdup2+=".flagstat"
cov=$file1
cov+=".coverage"
~/samtools-1.9/samtools flagstat $rmdup >$rmdup2
echo "~/samtools-1.9/samtools flagstat $rmdup >$rmdup2"
/usr/bin/samtools coverage $rmdup > $cov
echo "/usr/bin/samtools coverage $rmdup > $cov"
rm -rf SAM_FILES/$sam BAM_FILES/$bam B*/$sorted B*/$sorted2 
samtools depth -a $rmdup  > cov.ont.txt
