########### [3] Flye v2.9.1-b1784
# Flye v2.9.1-b1784 using a subset of the longest reads for initial disjointig assembly, 10 polishing iterations, an expected genome size of 151 Kb, and high-quality ONT reads.
# examine 'Flye_ONT_noAlt' -> ONT: 182 Kb, 973-fold depth - no asm
# examine 'Flye_ONT_ASM2' -> ONT: 151 Kb 1,107-fold depth - with asm
# depth 20 or 30 fold works less well
# no-alt-contigs  - no alternative contigs
# asm-coverage X  - reduced coverage for initial disjointig assembly 

/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  --no-alt-contigs --asm-coverage 30 -g 0.151m -i 10 --nano-hq  KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz   --out-dir Flye_ONT_ASM30 

/mnt/lustre/RDS-live/downing/Flye/bin/flye --threads 32 --no-alt-contigs --asm-coverage 50 -g 0.151m -i 10 --nano-hq KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz --out-dir Flye_ONT_ASM50

# 4 Kb threshold
# cp -r KRAKEN_VALID_FILES/merge_trim_fastp.fastq Flye_ONT_ASM/MAP_OMAN_OMAN/merge.fastp.chop.fastq 
seqkit seq --min-len 4000 KRAKEN_VALID_FILES/merge_trim_fastp.fastq  -o merge.4kb.fastq &
/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  --no-alt-contigs --asm-coverage 50  -g 0.151m -i 10 --nano-hq  merge.4kb.fastq  --out-dir Flye_ONT_ASM_4KB &> 4.o & 

# 2 Kb threshold
seqkit seq --min-len 2000 KRAKEN_VALID_FILES/merge_trim_fastp.fastq  -o merge.2kb.fastq &
/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  --no-alt-contigs  --asm-coverage 50 -g 0.151m -i 10 --nano-hq  merge.2kb.fastq  --out-dir Flye_ONT_ASM_2KB &> 2.o & 

# 1 Kb threshold
seqkit seq --min-len 1000 KRAKEN_VALID_FILES/merge_trim_fastp.fastq  -o merge.1kb.fastq &
/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  --no-alt-contigs --asm-coverage 50 -g 0.151m -i 10 --nano-hq  merge.1kb.fastq  --out-dir Flye_ONT_ASM_1KB &> 1.o & 

# 500 bp threshold
seqkit seq --min-len 500 KRAKEN_VALID_FILES/merge_trim_fastp.fastq  -o merge.500bp.fastq &
/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  --no-alt-contigs -g 0.151m -i 10 --nano-hq  merge.500bp.fastq  --out-dir Flye_ONT_ASM_500BP &> 500bp.o

/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  --no-alt-contigs  -g 0.151m -i 10 --nano-hq  KRAKEN_VALID_FILES/merge_trim_fastp.fastq.gz   --out-dir Flye_ONT_ASM

# Flye with alternative contigs
/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32  -g 0.151m -i 10 --nano-hq  KRAKEN_VALID_FILES/merge_trim_fastp.fastq  --out-dir Flye_ONT_ASM_With_Alt

--keep-haplotypes
/mnt/lustre/RDS-live/downing/Flye/bin/flye   --threads 32 --keep-haplotypes  -g 0.151m -i 10 --nano-hq  KRAKEN_VALID_FILES/merge_trim_fastp.fastq  --out-dir Flye_ONT_ASM_With_Alt2
