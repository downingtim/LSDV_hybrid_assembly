Construction & analysis of the LSDV Oman 2009 hybrid genome assembly.

We used a combination of long ONT reads and short Illumina reads to create a new reference genome sequence for the Oman 2009 sample. The optimised genome produced de novo was a single contig spanning 151,091 bp with a GC content of 25.9%. This came from Flye v2.9.1-b1784 (Kolmogorov et al 2019) output followed by polishing with Medaka v2.0.1 (https://github.com/nanoporetech/medaka) that made 15 changes, Pilon v1.24 (Walker et al 2014) that fixed two SNPs, and correction with Polypolish v0.6.0 (Bouras et al 2024) that fixed nine indel and two substitution errors. In total, the Illumina data only affected 28 bp (0.019%) of the 151,091 bp genome. Alternative approaches produced diffuse assemblies with large numbers of contigs and inaccurate total assembly length.

Technical details:

The sequencing data:
The Illumina FASTQ files are in FASTQ_ILLUMINA
The ONT FASTQ file links and merged file are in FASTQ_ONT

QC:
The script "1_runqc.sh" was used to run Fastqc, Fastp, the Fastq_quality_trimmer, nanoq and multiqc.

Viral read extraction:
The scripts "kraken.illumina.sh" and "kraken.ont.sh" were used to get the Kraken2 taxonomic assignment and report files.
The scripts "kraken.illumina.extract.2.sh" and "kraken.ont.extract.sh" were used to extract the viral reads from the above Kraken2 analyses.

Read counts:
The script "2_get_read_stats.py" was used to measure read number and rates.

Flye v2.9.1-b1784 using a subset of the longest reads for initial disjointig assembly, 10 polishing iterations, an expected genome size of 151 Kb, and high-quality ONT reads.
SPAdes v3.13.1 was run in careful mode using a range of input: the ONT reads, the Illumina reads, and the Flye assembly. 

Pilon:
The script "run_pilon.sh" implemented the Pilon iterations.

GC cotent plot:
The script "6_gc.R" was used to visualise the GC content in R.

Read mapping to the reference:
The scripts "10_run.sh" and "11_run_ont.sh" were used to map Illumin and ONT reads (respectively) to the Oman FASTA reference. Variations on these were used for mapping to the alternative FASTA versions.
The script "12_cov3.R" was used to visualise read depth and coverage levels.
The script "13_plot_cov_combined.R" was used to make Figure 1.

Notes:
These scripts are GAIT-assisted.