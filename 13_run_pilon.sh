#!/bin/bash

# Check if the file stem argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <optional_file_prefix>"
    echo "If no prefix is provided, sample names will be used directly"
fi

prefix="$1"

# Load required modules
module load samtools/samtools1.15
module load katuali
module load medaka/medaka1.7.1

# Define sample array
# samples=(449_05k_S1 449_05p_S4 449_07k_S2 449_07p_S5 449_11k_S3 449_11p_S6 Mongolia1 Mongolia2)
samples=( <add here>)

for sample in "${samples[@]}"; do
    echo "===== Processing sample: $sample ====="
    
    # Define directories with sample name
    SAMPLE_DIR="SPAdes_ILLUMINA_FASTQS/${sample}"
    
    # Define file stem (sample name or prefix+sample)
    if [ -n "$prefix" ]; then
        file_stem="${prefix}_${sample}"
    else
        file_stem="${sample}"
    fi
    
    echo "Using file stem: $file_stem"
    
    ##########  Pilon iteration 1 ##########
    echo "===== Starting Pilon iteration 1 ====="
    
    # Index the assembly
    echo "Indexing the assembly..."
    minimap2 -d ${SAMPLE_DIR}/assembly1.mmi ${SAMPLE_DIR}/scaffolds.fasta
    
    # Map Illumina reads to the assembly
    echo "Mapping Illumina reads to the assembly..."
    minimap2 -ax sr ${SAMPLE_DIR}/scaffolds.fasta \
        ILLUMINA_FASTQS/${file_stem}_1.fastq \
        ILLUMINA_FASTQS/${file_stem}_2.fastq \
        > SAM_FILES/${file_stem}.sam 2> ERROR_FILES/${file_stem}.errors.txt
    
    ls -lt SAM_FILES/${file_stem}.sam
    
    # Process BAM files
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
    bcftools mpileup -A -Ob -f ${SAMPLE_DIR}/scaffolds.fasta BAM_FILES/${file_stem}.sort.rmdup.bam | \
        bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/${file_stem}.vcf
    
    bgzip -f BCF_VCF_FILES/${file_stem}.vcf && tabix -f -p vcf BCF_VCF_FILES/${file_stem}.vcf.gz
    ls -lt BCF_VCF_FILES/${file_stem}.vcf.gz
    
    # Normalize the VCF
    echo "Normalizing the VCF..."
    bcftools norm -c w -f ${SAMPLE_DIR}/scaffolds.fasta -m-both -Oz \
        -o BCF_VCF_FILES/${file_stem}.norm.vcf.gz BCF_VCF_FILES/${file_stem}.vcf.gz && \
        tabix -f -p vcf BCF_VCF_FILES/${file_stem}.norm.vcf.gz
    
    ls -lt BCF_VCF_FILES/${file_stem}.norm.vcf.gz
    
    # Run Pilon polishing
    echo "Running Pilon polishing for iteration 1..."
    rm -rf ${SAMPLE_DIR}/${file_stem}.polished*
    
    java -Xmx32G -jar pilon-1.24.jar  --genome ${SAMPLE_DIR}/scaffolds.fasta \
        --bam BAM_FILES/${file_stem}.sort.rmdup.bam --output ${SAMPLE_DIR}/${file_stem}.polished_scaffolds \
        --fix all --changes --threads 28 --tracks --vcf

    ##########  Pilon iteration 2 ##########
    echo "===== Starting Pilon iteration 2 ====="
    
    # Index the polished assembly
    echo "Indexing the polished assembly for Pilon iteration 2..."
    minimap2 -d ${SAMPLE_DIR}/${file_stem}.polished.mmi ${SAMPLE_DIR}/${file_stem}.polished_scaffolds.fasta
    
    ls -lt ${SAMPLE_DIR}/${file_stem}.polished.mmi
    
    # Map reads to the polished assembly
    echo "Mapping Illumina reads to the polished assembly (iteration 2)..."
    minimap2 -ax sr ${SAMPLE_DIR}/${file_stem}.polished_scaffolds.fasta \
        ILLUMINA_FASTQS/${file_stem}_1.fastq \
        ILLUMINA_FASTQS/${file_stem}_2.fastq \
        > SAM_FILES/${file_stem}.polished.sam 2> ERROR_FILES/${file_stem}.polished.errors.txt
    
    ls -lt SAM_FILES/${file_stem}.polished.sam
    
    # Process BAM files for iteration 2
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
    
    # Variant calling for iteration 2
    echo "Running bcftools mpileup and variant calling (iteration 2)..."
    bcftools mpileup -A -d 10000 -Ob -f ${SAMPLE_DIR}/${file_stem}.polished_scaffolds.fasta \
        BAM_FILES/${file_stem}.polished.sort.rmdup.bam | \
        bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/${file_stem}.polished.vcf
    
    bgzip -f BCF_VCF_FILES/${file_stem}.polished.vcf && tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished.vcf.gz
    ls -lt BCF_VCF_FILES/${file_stem}.polished.vcf.gz
    
    # Normalize the polished VCF
    echo "Normalizing the polished VCF..."
    bcftools norm -c w -f ${SAMPLE_DIR}/${file_stem}.polished_scaffolds.fasta -m-both -Oz \
        -o BCF_VCF_FILES/${file_stem}.polished.norm.vcf.gz BCF_VCF_FILES/${file_stem}.polished.vcf.gz && \
        tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished.norm.vcf.gz
    
    ls -lt BCF_VCF_FILES/${file_stem}.polished.norm.vcf.gz
    
    # Run Pilon polishing iteration 2
    echo "Running Pilon polishing iteration 2..."
    rm -rf ${SAMPLE_DIR}/${file_stem}.polished_scaffold2*
    
    java -Xmx32G -jar pilon-1.24.jar \
        --genome ${SAMPLE_DIR}/${file_stem}.polished_scaffolds.fasta \
        --bam BAM_FILES/${file_stem}.polished.sort.rmdup.bam \
        --output ${SAMPLE_DIR}/${file_stem}.polished_scaffold2 \
        --fix all --changes --threads 28 --tracks --vcf
    
    # Final outputs
    echo "Pilon iteration 2 complete. Final output files:"
    ls -lt ${SAMPLE_DIR}/${file_stem}.polished_scaffold2*
    
    
    ##########  Pilon iteration 3 ##########
    echo "===== Starting Pilon iteration 3 ====="
    
    # Index the polished2 assembly
    echo "Indexing the polished2 assembly..."
    minimap2 -d ${SAMPLE_DIR}/${file_stem}.polished2.mmi ${SAMPLE_DIR}/${file_stem}.polished_scaffold2.fasta
    
    # Map reads to the polished2 assembly
    echo "Mapping Illumina reads to the polished2 assembly..."
    minimap2 -ax sr ${SAMPLE_DIR}/${file_stem}.polished_scaffold2.fasta \
        ILLUMINA_FASTQS/${file_stem}_1.fastq      ILLUMINA_FASTQS/${file_stem}_2.fastq \
        > SAM_FILES/${file_stem}.polished2.sam 2> ERROR_FILES/${file_stem}.polished2.errors.txt
    
    # Process BAM files for iteration 3
    echo "Processing BAM files for Pilon iteration 3..."
    samtools view -bS SAM_FILES/${file_stem}.polished2.sam > BAM_FILES/${file_stem}.polished2.bam
    samtools sort -n BAM_FILES/${file_stem}.polished2.bam -o BAM_FILES/${file_stem}.polished2.sort.bam
    samtools fixmate -m BAM_FILES/${file_stem}.polished2.sort.bam BAM_FILES/${file_stem}.polished2.sort.bam.2
    samtools sort BAM_FILES/${file_stem}.polished2.sort.bam.2 -o BAM_FILES/${file_stem}.polished2.sort.bam
    samtools markdup -r -s BAM_FILES/${file_stem}.polished2.sort.bam BAM_FILES/${file_stem}.polished2.sort.rmdup.bam
    samtools index BAM_FILES/${file_stem}.polished2.sort.rmdup.bam
    samtools flagstat BAM_FILES/${file_stem}.polished2.sort.rmdup.bam > BAM_FILES/${file_stem}.polished2.flagstat
    samtools coverage BAM_FILES/${file_stem}.polished2.sort.rmdup.bam > COVERAGE/${file_stem}.polished2.coverage.txt 2> ERROR_FILES/${file_stem}.polished2.coverage.txt
    
    # Variant calling for iteration 3
    echo "Running bcftools mpileup and variant calling (iteration 3)..."
    bcftools mpileup -A -Ob -f ${SAMPLE_DIR}/${file_stem}.polished_scaffold2.fasta \
        BAM_FILES/${file_stem}.polished2.sort.rmdup.bam | \
        bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/${file_stem}.polished2.vcf
    
    # Compress and index the polished2 VCF
    bgzip -f BCF_VCF_FILES/${file_stem}.polished2.vcf
    tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished2.vcf.gz
    
    # Normalize the polished2 VCF
    echo "Normalizing the polished2 VCF..."
    bcftools norm -c w -f ${SAMPLE_DIR}/${file_stem}.polished_scaffold2.fasta -m-both -Oz \
        -o BCF_VCF_FILES/${file_stem}.polished2.norm.vcf.gz BCF_VCF_FILES/${file_stem}.polished2.vcf.gz
    tabix -f -p vcf BCF_VCF_FILES/${file_stem}.polished2.norm.vcf.gz
    
    # Run Pilon polishing iteration 3
    echo "Running Pilon polishing iteration 3..."
    rm -rf ${SAMPLE_DIR}/${file_stem}.polished_scaffold3*
    
    java -Xmx32G -jar pilon-1.24.jar --genome ${SAMPLE_DIR}/${file_stem}.polished_scaffold2.fasta \
        --bam BAM_FILES/${file_stem}.polished2.sort.rmdup.bam  --output ${SAMPLE_DIR}/${file_stem}.polished_scaffold3 \
        --fix all --changes --threads 28 --tracks --vcf
    
    echo "Pilon iteration 3 complete. Final output files:"
    ls -lt ${SAMPLE_DIR}/${file_stem}.polished_scaffold3*
    
    echo "===== All processing completed for sample: $sample ====="
    echo ""
done

echo "Script completed successfully!"
