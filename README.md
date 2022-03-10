# utils
Repo filled with little useful bioinformatics scripts

alter_vcf.sh - alters vcf to add or remove chr naming convention used by CellRanger so that BAM file and VCF match (much quicker to edit VCF file)

pipeline.sh - manipulates SAM, BAM and CRAM files so that the chr naming convention used by CellRanger is added or removed, reorders the header as chr naming convention orders head chr1, chr10, chr11 etc, resorts BAM file to match new header and removes contigs. Also generates files listing chromosomes and contigs present in BAM and VCF file to check naming convention matches.

remove_contigs.sh - removes extra contigs added to BAM files by processes such as Starsolo 

move-fastqs.sh - combines fastq label with cram label and moves to a new directory so that starsolo can execute on uniquely labelled fastqs
