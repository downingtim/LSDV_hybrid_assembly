################### Assembly quality - Quast  

# conda install -c bioconda quast 
quast Flye_ONT_ASM/assembly.fasta  -o Flye_ONT_ASM/QUAST 
~/quast/quast.py Flye_ONT_ASM2/assembly.fasta  -o Flye_ONT_ASM2/QUAST 
~/quast/quast.py Flye_ONT_ASM_With_Alt/assembly.fasta  -o Flye_ONT_ASM_With_Alt/QUAST 
~/quast/quast.py Flye_ONT_ASM_With_Alt2/assembly.fasta  -o Flye_ONT_ASM_With_Alt2/QUAST 
~/quast/quast.py Flye_ONT_ASM20/assembly.fasta  -o Flye_ONT_ASM20/QUAST 
~/quast/quast.py Flye_ONT_ASM30/assembly.fasta  -o Flye_ONT_ASM30/QUAST 
~/quast/quast.py Flye_ONT_ASM50/assembly.fasta  -o Flye_ONT_ASM50/QUAST 
~/quast/quast.py  Flye_ONT_ASM_4KB/assembly.fasta  -o Flye_ONT_ASM_4KB/QUAST 
~/quast/quast.py  Flye_ONT_ASM_4KB/assembly.fasta  -o Flye_ONT_ASM_4KB/QUAST 
~/quast/quast.py  Flye_ONT_ASM_2KB/assembly.fasta  -o Flye_ONT_ASM_2KB/QUAST & 
~/quast/quast.py  Flye_ONT_ASM_1_KB/assembly.fasta  -o Flye_ONT_ASM_1_KB/QUAST & 

seqtk seq -r in.fa > out.fa

quast SPAdes_ONT_2/scaffolds.fasta  -o SPAdes_ONT_2/QUAST 
quast SPAdes_Illumina_1/scaffolds.fasta  -o SPAdes_Illumina_1/QUAST 
quast SPAdes_Hybrid_1/scaffolds.fasta  -o SPAdes_Hybrid_1/QUAST 
quast SPAdes_Illumina_Flye/scaffolds.fasta  -o SPAdes_Illumina_Flye/QUAST 
quast SPAdes_Hybrid_Flye/scaffolds.fasta  -o SPAdes_Hybrid_Flye/QUAST 

grep "Total length" */QUAST/report.txt 
grep "auN" */QUAST/report.txt 

# Flye with ONT worked best 150,935 bp vs SPAdes with ONT, 136 Kb
# Flye ASM2 (10-fold coverage) has coverage of 1,107 vs 1,100 -> better
