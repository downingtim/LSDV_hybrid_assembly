#!/bin/bash

# Directory containing the FASTP files
fastp_files_dir="ONT_FASTQS"

# Base command components
base_command="kraken2"
db_path="/home/share/"
threads=3
output_dir_base="/mnt/lustre/RDS-ephemeral/downing/ASFV_ASSEMBLY/KRAKEN_FILES2"
output_dir_base2="/mnt/lustre/RDS-ephemeral/downing/ASFV_ASSEMBLY/KRAKEN_REPORT2"

# Get a list of sample names by listing the files in the directory and extracting the base names
sample_names=$(ls $fastp_files_dir | grep -E '_fastp_trim.fastq$' | sed -E 's/_fastp_trim.fastq$//' | sort | uniq)

# Iterate over sample names and construct the command for each
for sample_name in $sample_names; do
    # Construct output directories
    output_dir="${output_dir_base}/${sample_name}"
    output_dir2="${output_dir_base2}/${sample_name}"

    # Construct the full command
    command="${base_command} --db ${db_path} --threads ${threads} --report ${output_dir2} --use-names --output ${output_dir} ${fastp_files_dir}/${sample_name}_fastp_trim.fastq "

    # Print the command for debugging purposes
    echo "${command}"
done
