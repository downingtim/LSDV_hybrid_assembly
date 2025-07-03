#!/bin/bash

fastp_files_dir="ONT_FASTQS"

# Base command components  merged_25_905_3_1_fastp_trim.fastq
base_command="python /mnt/lustre/RDS-live/downing/KrakenTools/extract_kraken_reads.py -k "
output_dir_base="/mnt/lustre/RDS-ephemeral/downing/ASFV_ASSEMBLY/KRAKEN_FILES2"
output_dir_base2="/mnt/lustre/RDS-ephemeral/downing/ASFV_ASSEMBLY/KRAKEN_REPORT2"
sample_names=$(ls $fastp_files_dir | grep -E '_fastp_trim.fastq$' | sed -E 's/_fastp_trim.fastq$//' | sort | uniq)

for sample_name in $sample_names; do
    output_dir="${output_dir_base}/${sample_name}"
    output_dir2="${output_dir_base2}/${sample_name}"
    out_file=${sample_name}

    taxon=2732526 # ASFV
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_fastp_trim.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents "
    echo "${command}"

    taxon=9903 # cow
    command="${base_command} $output_dir_base/${sample_name}  -s1 $fastp_files_dir/${sample_name}_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_cow.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents "
    echo "${command}"

    taxon=9822 # pig
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_pig.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents "
    echo "${command}"

    taxon=9935 # sheep
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_sheep.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name}  --include-children --include-parents "
    echo "${command}"

    taxon=9605 # human
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_human.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name}   --include-children --include-parents "
    echo "${command}"
done
