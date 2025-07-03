########### medaka v2.0.1 -> before pilon
# /ephemeral/downing/LSDV/HYBRID_ASSEMBLY/Flye_ONT_ASM/medaka/bin/medaka_consensus

module load medaka

/ephemeral/downing/LSDV/HYBRID_ASSEMBLY/Flye_ONT_ASM/medaka/bin/medaka_consensus -i ../KRAKEN_VALID_FILES/merge_trim_fastp.fastq -d polished_assembly2.fasta -o medaka_1/ -t 24  -m r1041_e82_400bps_sup_v5.0.0 

cd Flye_ONT_ASM_2KB
/ephemeral/downing/LSDV/HYBRID_ASSEMBLY/Flye_ONT_ASM/medaka/bin/medaka_consensus -i ../merge.2kb.fastq -d assembly.fasta -o medaka/ -t 32 -m r1041_e82_400bps_sup_v5.0.0 

#  -m r941_min_sup_g507: This specifies the model used for consensus. Choose a model based on your sequencing platform (e.g., r941_min_sup_g507 is a model for the ONT R9.4.1 chemistry).
# -m r10_min_sup_g507 : - is a model for the ONT R10.3 chemistry.
nucmer medaka_1/consensus.fasta polished_assembly.fasta
show-diff -q out.delta > medaka_changes.txt
cat medaka_1/consensus.fasta polished_assembly.fasta > medaka_1/test.fa
mafft --thread 50 --auto  test.fa  >  test.aln
