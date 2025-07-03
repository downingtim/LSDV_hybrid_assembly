#!/bin/bash

# Check if the file stem argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <file_stem>"
  file_stem="Oman_2009_S1_trim_fastp"
else
  file_stem="$1"
fi

module load samtools/samtools1.15
module load katuali
module load medaka/medaka1.7.1

## # # # # #  Pilon iteration 1

# Index the first polished assembly
echo "Indexing the Medaka-polished assembly for Pilon iteration 1..."
minimap2 -d Flye_ONT_ASM/Flye_ONT_ASM.assembly1.mmi Flye_ONT_ASM/medaka1/consensus.fasta
ls -lt Flye_ONT_ASM/Flye_ONT_ASM.assembly1.mmi

# Map Illumina reads to the Medaka-polished assembly
echo "Mapping Illumina reads to the Medaka-polished assembly..."
minimap2 -ax sr Flye_ONT_ASM/medaka1/consensus.fasta  KRAKEN_VALID_FILES/${file_stem}_1.fastq  KRAKEN_VALID_FILES/${file_stem}_2.fastq > SAM_FILES/${file_stem}.sam 2> ERROR_FILES/${file_stem}.errors.txt  
ls -lt SAM_FILES/${file_stem}.sam

# Convert SAM to BAM, sort, fix mates, mark duplicates, index, flagstat, coverage, and mpileup in one pipeline
echo "Processing BAM files: sorting, fixing mates, marking duplicates, indexing, and generating stats..."
samtools view -bS SAM_FILES/${file_stem}.sam | \
samtools sort -n - | \
samtools fixmate -m - - | \
samtools sort - | \
samtools markdup -r -s - BAM_FILES/${file_stem}.sort.rmdup.bam && \
samtools index BAM_FILES/${file_stem}.sort.rmdup.bam && \
samtools flagstat BAM_FILES/${file_stem}.sort.rmdup.bam > BAM_FILES/${file_stem}.flagstat && \
samtools coverage BAM_FILES/${file_stem}.sort.rmdup.bam > COVERAGE/${file_stem}.coverage.txt 2> ERROR_FILES/${file_stem}.coverage.txt
ls -lt COVERAGE/${file_stem}.coverage.txt BAM_FILES/${file_stem}.sort.rmdup.bam 

# Variant calling with bcftools
echo "Running bcftools mpileup and variant calling..."
bcftools mpileup -A -d 2500 -Ob -f Flye_ONT_ASM/assembly.fasta BAM_FILES/${file_stem}.sort.rmdup.bam | \
bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/${file_stem}.vcf 
bgzip BCF_VCF_FILES/${file_stem}.vcf && tabix -f -p vcf BCF_VCF_FILES/${file_stem}.vcf.gz
ls -lt BCF_VCF_FILES/${file_stem}.vcf.gz

# Normalize the VCF
echo "Normalizing the VCF..."
bcftools norm -c w -f Flye_ONT_ASM/assembly.fasta -m-both -Oz -o BCF_VCF_FILES/${file_stem}.norm.vcf.gz BCF_VCF_FILES/${file_stem}.vcf.gz && tabix -f -p vcf BCF_VCF_FILES/${file_stem}.norm.vcf.gz
ls -lt  BCF_VCF_FILES/${file_stem}.norm.vcf.gz

# Run Pilon polishing
echo "Running Pilon polishing for iteration 1..."
rm -rf Flye_ONT_ASM/${file_stem}.polished_assembly*
java -Xmx32G -jar pilon-1.24.jar --genome Flye_ONT_ASM/assembly.fasta --bam BAM_FILES/${file_stem}.sort.rmdup.bam --output Flye_ONT_ASM/${file_stem}.polished_assembly --fix all --changes --threads 28 --tracks --vcf
ls -lt Flye_ONT_ASM/${file_stem}.polished_assembly*  BAM_FILES/${file_stem}.sort.rmdup.bam

## # # # # #  Pilon iteration 2

# Index the polished assembly
echo "Indexing the polished assembly for Pilon iteration 2..."
minimap2 -d Flye_ONT_ASM/Flye_ONT_ASM.${file_stem}.polished.mmi Flye_ONT_ASM/${file_stem}.polished_assembly.fasta
ls -lt Flye_ONT_ASM/Flye_ONT_ASM.${file_stem}.polished.mmi

# Map reads to the polished assembly
echo "Mapping Illumina reads to the polished assembly (iteration 2)..."
minimap2 -ax sr Flye_ONT_ASM/${file_stem}.polished_assembly.fasta KRAKEN_VALID_FILES/${file_stem}_1.fastq KRAKEN_VALID_FILES/${file_stem}_2.fastq > SAM_FILES/${file_stem}.polished.sam 2> ERROR_FILES/${file_stem}.polished.errors.txt
ls -lt SAM_FILES/${file_stem}.polished.sam

# Convert SAM to BAM, sort, fix mates, mark duplicates, index, flagstat, coverage, and mpileup in one pipeline
echo "Processing BAM files for Pilon iteration 2..."
samtools view -bS SAM_FILES/${file_stem}.polished.sam | \
samtools sort -n - | \
samtools fixmate -m - - | \
samtools sort - | \
samtools markdup -r -s - BAM_FILES/${file_stem}.polished.sort.rmdup.bam && \
samtools index BAM_FILES/${file_stem}.polished.sort.rmdup.bam && \
samtools flagstat BAM_FILES/${file_stem}.polished.sort.rmdup.bam > BAM_FILES/${file_stem}.polished.flagstat && \
samtools coverage BAM_FILES/${file_stem}.polished.sort.rmdup.bam > COVERAGE/${file_stem}.polished.coverage.txt 2> ERROR_FILES/${file_stem}.polished.coverage.txt
ls -lt COVERAGE/${file_stem}.polished.coverage.txt BAM_FILES/${file_stem}.polished.sort.rmdup.bam 

# bcftools mpileup and variant calling
echo "Running bcftools mpileup and variant calling (iteration 2)..."
bcftools mpileup -A -d 2500 -Ob -f Flye_ONT_ASM/polished_assembly.fasta BAM_FILES/${file_stem}.polished.sort.rmdup.bam | \
bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/${file_stem}.polished.vcf 
bgzip BCF_VCF_FILES/${file_stem}.polished.vcf && tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished.vcf.gz
ls -lt BCF_VCF_FILES/${file_stem}.polished.vcf.gz

# Normalize the polished VCF
echo "Normalizing the polished VCF..."
bcftools norm -c w -f Flye_ONT_ASM/polished_assembly.fasta -m-both -Oz -o BCF_VCF_FILES/${file_stem}.polished.norm.vcf.gz BCF_VCF_FILES/${file_stem}.polished.vcf.gz && tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished.norm.vcf.gz
ls -lt BCF_VCF_FILES/${file_stem}.polished.norm.vcf.gz

# Run Pilon polishing iteration 2
echo "Running Pilon polishing iteration 2..."
rm -rf Flye_ONT_ASM/${file_stem}.polished_assembly2*
java -Xmx32G -jar pilon-1.24.jar --genome Flye_ONT_ASM/${file_stem}.polished_assembly.fasta --bam BAM_FILES/${file_stem}.polished.sort.rmdup.bam --output Flye_ONT_ASM/${file_stem}.polished_assembly2 --fix all --changes --threads 28 --tracks --vcf

# Final outputs
echo "Pilon iteration 2 complete. Final output files:"
ls -lt Flye_ONT_ASM/${file_stem}.polished_assembly2*

## # # # # #  Pilon iteration 3

# Index the polished2 assembly
minimap2 -d Flye_ONT_ASM/Flye_ONT_ASM.${file_stem}.polished2.mmi Flye_ONT_ASM/${file_stem}.polished_assembly2.fasta

# Map reads to the polished2 assembly
minimap2 -ax sr Flye_ONT_ASM/${file_stem}.polished_assembly2.fasta KRAKEN_VALID_FILES/${file_stem}__1.fastq KRAKEN_VALID_FILES/${file_stem}__2.fastq > SAM_FILES/${file_stem}.polished2.sam 2> ERROR_FILES/${file_stem}.polished2.errors.txt

# Convert to BAM
samtools view -bS SAM_FILES/${file_stem}.polished2.sam > BAM_FILES/${file_stem}.polished2.bam

# Sort, fix mates, sort again, mark duplicates, index, get flagstat stats, get coverage, and run mpileup
~/samtools-1.9/samtools sort -n BAM_FILES/${file_stem}.polished2.bam -o BAM_FILES/${file_stem}.polished2.sort.bam
~/samtools-1.9/samtools fixmate -m BAM_FILES/${file_stem}.polished2.sort.bam BAM_FILES/${file_stem}.polished2.sort.bam.2
~/samtools-1.9/samtools sort BAM_FILES/${file_stem}.polished2.sort.bam.2 -o BAM_FILES/${file_stem}.polished2.sort.bam
~/samtools-1.9/samtools markdup -r -s BAM_FILES/${file_stem}.polished2.sort.bam BAM_FILES/${file_stem}.polished2.sort.rmdup.bam
~/samtools-1.9/samtools index BAM_FILES/${file_stem}.polished2.sort.rmdup.bam
~/samtools-1.9/samtools flagstat BAM_FILES/${file_stem}.polished2.sort.rmdup.bam > BAM_FILES/${file_stem}.polished2.flagstat
/usr/bin/samtools coverage BAM_FILES/${file_stem}.polished2.sort.rmdup.bam > COVERAGE/${file_stem}.polished2.coverage.txt 2> ERROR_FILES/${file_stem}.polished2.coverage.txt

# bcftools mpileup and variant calling
bcftools mpileup -A -d 2500 -Ob -f Flye_ONT_ASM/polished_assembly2.fasta BAM_FILES/${file_stem}.polished2.sort.rmdup.bam | \
bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/${file_stem}.polished2.vcf

# Compress and index the polished2 VCF
bgzip BCF_VCF_FILES/${file_stem}.polished2.vcf
tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished2.vcf.gz

# Normalize the polished2 VCF
bcftools norm -c w -f Flye_ONT_ASM/polished_assembly2.fasta -m-both -Oz -o BCF_VCF_FILES/${file_stem}.polished2.norm.vcf.gz BCF_VCF_FILES/${file_stem}.polished2.vcf.gz
tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished2.norm.vcf.gz

rm -rf Flye_ONT_ASM/${file_stem}.polished_assembly3*
java -Xmx32G -jar pilon-1.24.jar  --genome Flye_ONT_ASM/${file_stem}.polished_assembly.fasta  --bam BAM_FILES/${file_stem}.polished2.sort.rmdup.bam  --output Flye_ONT_ASM/${file_stem}.polished_assembly3  --fix all --changes --threads 28  --tracks --vcf
ls -lt Flye_ONT_ASM/${file_stem}.polished_assembly3*

# medaka 2

# Index the assembly
minimap2 -d Flye_ONT_ASM/Flye_ONT_ASM.polished_assembly3.mmi Flye_ONT_ASM/${file_stem}.polished_assembly3.fasta
ls -lt  Flye_ONT_ASM/Flye_ONT_ASM.polished_assembly3.mmi 

# Run Medaka polishing
rm -rf  Flye_ONT_ASM/medaka2/
/mnt/lustre/RDS-ephemeral/downing/LSDV/HYBRID_ASSEMBLY/Flye_ONT_ASM/medaka/bin/medaka_consensus -i KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz  -d Flye_ONT_ASM/${file_stem}.polished_assembly3.fasta  -o Flye_ONT_ASM/medaka2/ -m r941_min_sup_g507
ls -lt  Flye_ONT_ASM/medaka2/

