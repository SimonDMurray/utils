# utils
Repo filled with little useful bioinformatics scripts

`alter_vcf.sh` - alters vcf to add or remove chr naming convention used by CellRanger so that BAM file and VCF match (much quicker to edit VCF file)

`pipeline.sh` - manipulates SAM, BAM and CRAM files so that the chr naming convention used by CellRanger is added or removed, reorders the header as chr naming convention orders head chr1, chr10, chr11 etc, resorts BAM file to match new header and removes contigs. Also generates files listing chromosomes and contigs present in BAM and VCF file to check naming convention matches.

`remove_contigs.sh` - removes extra contigs added to BAM files by processes such as Starsolo 

`move-fastqs.sh` - combines fastq label with cram label and moves to a new directory so that starsolo can execute on uniquely labelled fastqs

`add_assignee.sh` - adds a new file to each ticket in our teams ticket structure with the information on who the jira ticket is assigned to

`wget_parallel.sh` - parallelises downloading files using wget and xargs

`bz2_to_gz.sh` - parallelises the conversion of bz2 to gz

`generic-bsub.sh` - illustrates the concept of how to bsub properly

`bsub_bam_to_fastq.sh` - script that converts bams to fastqs on Sanger FARM

`generic-job-array-bsub.sh` - illustrates dynamic bsub job array submission utilising inputted sample file
