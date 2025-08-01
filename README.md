Construction & analysis of the LSDV Oman 2009 hybrid genome assembly.

We used a combination of long ONT reads and short Illumina reads to create a new reference genome sequence for the Oman 2009 sample. The optimised genome produced de novo was a single contig spanning 151,091 bp with a GC content of 25.9%. This came from Flye v2.9.1-b1784 (Kolmogorov et al 2019) output followed by polishing with Medaka v2.0.1 (https://github.com/nanoporetech/medaka) that made 15 changes, Pilon v1.24 (Walker et al 2014) that fixed two SNPs, and correction with Polypolish v0.6.0 (Bouras et al 2024) that fixed nine indel and two substitution errors. In total, the Illumina data only affected 28 bp (0.019%) of the 151,091 bp genome. Alternative approaches produced diffuse assemblies with large numbers of contigs and inaccurate total assembly length.

--Technical details--

Sequencing data:

The Illumina FASTQ files are in FASTQ_ILLUMINA
The ONT FASTQ file links and merged file are in FASTQ_ONT

QC:

The script "1_runqc.illumina.sh" was used to run Fastqc, Fastp, and fastq_quality_trimmer on Illumina data.
The script "1_runqc.ont.sh" was used to run Fastqc, Fastp, and nanoq on ONT data.
 
Viral read extraction:

The scripts "3_kraken.illumina.sh" and "4_kraken.ont.sh" were used to get the Kraken2 taxonomic assignment and report files.
The scripts "5_kraken.illumina.extract.sh" and "6_kraken.ont.extract.sh" were used to extract the viral reads from the above Kraken2 analyses.

Read counts:

The script "7_get_read_stats.py" was used to measure read number and rates.

Assembly & optimisation:

Flye v2.9.1-b1784 using a subset of the longest reads for initial disjointig assembly, 10 polishing iterations, an expected genome size of 151 Kb, and high-quality ONT reads. See "8_flye.sh".
SPAdes v3.13.1 was run in careful mode using a range of input: the ONT reads, the Illumina reads, and the Flye assembly. See "9_spades.sh".
Medaka was run with comamnds in 10_medaka.sh.
Polypolish and pypolca were run with commands in 11_polypolish_pypolca.sh.
Quast was used to assess assembly quality with 12_quast.sh.
Pilon: The script "13_run_pilon.sh" implemented the Pilon iterations.
GC cotent plot: The script "14_gc.R" was used to visualise the GC content in R (Figure S1).

Read mapping to the reference & visualisation:

The scripts "15_run.sh" and "16_run_ont.sh" were used to map Illumin and ONT reads (respectively) to the Oman FASTA reference. Variations on these were used for mapping to the alternative FASTA versions.
The script "17_cov.R" was used to visualise read depth and coverage levels.
The script "18_plot_cov_combined.R" was used to make Figure 1.

Notes:
These scripts are GAIT-assisted.