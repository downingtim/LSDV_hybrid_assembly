#!/bin/bash

fastp_files_dir="ILLUMINA_FASTQS"

# Base command components  merged_25_905_3_1_fastp_trim.fastq
base_command="python /mnt/lustre/RDS-live/downing/KrakenTools/extract_kraken_reads.py -k "
output_dir_base="KRAKEN_FILES"
output_dir_base2="KRAKEN_REPORT"
sample_names=$(ls $fastp_files_dir | grep -E '_[12]_fastp_trim.fastq$' | sed -E 's/_([12])_fastp_trim.fastq$//' | sort | uniq)

for sample_name in $sample_names; do
    output_dir="${output_dir_base}/${sample_name}"
    output_dir2="${output_dir_base2}/${sample_name}"
    out_file=${sample_name}

    taxon=10497 
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_1_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_1_fastp_trim.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents "
    echo "${command}"
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_2_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_2_fastp_trim.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name}  --include-children --include-parents "
    echo "${command}"

    taxon=9822 # pig
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_1_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_1_pig.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents "
#     echo "${command}"
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_2_fastp_trim.fastq -o KRAKEN_VALID_FILES/${out_file}_2_pig.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name}  --include-children --include-parents "
#     echo "${command}"

done
