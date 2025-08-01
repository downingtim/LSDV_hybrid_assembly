Polypolish v0.6.0 was used to correct errors in repeats in the Pilon output. To run this, BWA v0.7.17 was used to map the Illumina to all possible matches. These reads were filter by insert size and then used for polishing.

bwa index pilon_assembly1.fasta
bwa mem -t 34 -a pilon_assembly1.fasta ../KRAKEN_VALID_FILES/Oman_2009_S1_trim_fastp_1.fastq > a1.sam
bwa mem -t 34 -a pilon_assembly1.fasta ../KRAKEN_VALID_FILES/Oman_2009_S1_trim_fastp_2.fastq  > a2.sam
~/Polypolish/target/release/polypolish filter --in1 a1.sam --in2 a2.sam --out1 f1.sam --out2 f2.sam
~/Polypolish/target/release/polypolish polish pilon.fasta f1.sam f2.sam > poly.fasta
 rm *.amb *.ann *.bwt *.pac *.sa *.sam

__________
Pypolca v0.3.1 in careful mode was used to correct using the Illumina reads. 

Pypolca v0.3.1 

#  [1] repair.sh if needed 
repair.sh  in1=../KRAKEN_VALID_FILES/Oman_2009_S1__1_trim_fastp.fastq.gz  in2=../KRAKEN_VALID_FILES/Oman_2009_S1__2_trim_fastp.fastq.gz  out1=../KRAKEN_VALID_FILES/Oman.repair.1.fastq out2=../KRAKEN_VALID_FILES/Oman.repair.2.fastq outs=../KRAKEN_VALID_FILES/Oman.repair.0.gz ain=t
# ain = allow identical names = TRUE

#  conda create -n pypolca_env pypolca    # creates conda environment with pypolca 
conda activate pypolca_env    # activates conda environment
pypolca run -t 24 -a pilon.fasta -1 ../KRAKEN_VALID_FILES/Oman.repair.1.fastq -2 ../KRAKEN_VALID_FILES/Oman.repair.2.fastq   -o PYPOLCA_OUT --careful -f
# runs pypolca with --careful
