#!/bin/bash

######################
# Preliminary steps: #
######################
# change working directory
cd /Users/nickwinters/Desktop/DS_Projects/R/Bioinformatics/BulkRNAseqWF

# get the genome indices for step 3
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# get gtf for step 4
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz

######################
# STEP 1: Run fastqc #
######################

# Obtain quality control metrics for FASTQ file
fastqc data/demo.fastq -o data/

###########################
# STEP 2: Run trimmomatic #
###########################

# Trim the last 10 sequences
trimmomatic SE -threads 4 data/demo.fastq data/demo_trimmed.fastq TRAILING:10 -phred33
echo "Trimmomatic finished running!"

# STEP 2a: Run fastqc again
fastqc data/demo_trimmed.fastq -o data/

######################
# STEP 3: Run HISAT2 #
######################

# Map reads using HISAT2
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U data/demo_trimmed.fastq | samtools sort -o HISAT2/demo_trimmed.bam

##############################################
# STEP 4: Run featureCounts - Quantification #
##############################################

# Obtain a counts matrix
featureCounts -S 2 -a HISAT2/grch38/Homo_sapiens.GRCh38.106.gtf -o quants/demo_featurecounts.txt HISAT2/demo_trimmed.bam

#####################
# STEP 4: View data #
#####################

cat quants/demo_featurecounts.txt.summary

cat quants/demo_featurecounts.txt

cat quants/demo_featurecounts.txt | cut -f1,7 | less
