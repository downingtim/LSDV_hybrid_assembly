________________________________
###########  [2] SPADES

SPAdes v3.13.1 was run in careful mode using a range of input: the ONT reads, the Illumina reads, and the Flye assembly. 

# long only
spades -t 32 --careful  -s KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz -o SPAdes_ONT_2

# illumina only
spades -t 22 --careful  -1 KRAKEN_VALID_FILES/Oman.repair.1.fastq.gz   -2 KRAKEN_VALID_FILES/Oman.repair.2.fastq.gz     -o SPAdes_Illumina_1  

# hybrid
spades -t 22 --careful  -1 KRAKEN_VALID_FILES/Oman.repair.1.fastq.gz   -2 KRAKEN_VALID_FILES/Oman.repair.2.fastq.gz  --nanopore KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz -o SPAdes_Hybrid_1  

# with Flye assembly hybrid
spades -t 22 --careful  -1 KRAKEN_VALID_FILES/Oman.repair.1.fastq  -2 KRAKEN_VALID_FILES/Oman.repair.2.fastq --nanopore KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz --trusted-contigs Flye_ONT_ASM/assembly.fasta -o SPAdes_Hybrid_Flye &> o 

# with Flye assembly, Illumina
spades -t 22 --careful  -1 KRAKEN_VALID_FILES/Oman.repair.1.fastq  -2 KRAKEN_VALID_FILES/Oman.repair.2.fastq  --trusted-contigs Flye_ONT_ASM/assembly.fasta -o SPAdes_Illumina_Flye
