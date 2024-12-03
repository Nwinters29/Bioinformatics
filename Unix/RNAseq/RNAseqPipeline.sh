#!/bin/bash

SECONDS=0

######################
# Preliminary steps: #
######################
# change working directory
cd /Users/nickwinters/Desktop/DS_Projects/R/Bioinformatics/BulkRNAseqWF
# get the genome indices for step 3
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# get gtf for step 4
# wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz

######################
# STEP 1: Run fastqc #
######################
# Arguments: -o was used to specify where the qc report would be stored on the computer
fastqc data/demo.fastq -o data/

# Observation: 
#      0 were flagged as poor quality
#      No adapter sequences were present
# note: 
#      Recommended to trim reads with scores less than 20 and any represeneted adapter sequences.
#      If trimming is necessary it is recommended to keep at least 80% of sequence length
#      Quality of reads tend to drop off towards the end of the seqence due to signal decay during the end of the sequencing run.

##################################################################
# STEP 2: Run trimmomatic to trim the last reads in the seqeunce #
##################################################################
trimmomatic SE -threads 4 data/demo.fastq data/demo_trimmed.fastq TRAILING:10 -phred33
echo "Trimmomatic finished running!"

# STEP 2a: Run fastqc again
fastqc data/demo_trimmed.fastq -o data/

## Observation: 
##      Median base quality has improved

######################
# STEP 3: Run HISAT2 #
######################
## note: was RNAseq obtained using a stranded protocol?
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U data/demo_trimmed.fastq | samtools sort -o HISAT2/demo_trimmed.bam
echo "HISAT2 finished running!"

##############################################
# STEP 4: Run featureCounts - Quantification #
##############################################
featureCounts -S 2 -a HISAT2/grch38/Homo_sapiens.GRCh38.106.gtf -o quants/demo_featurecounts.txt HISAT2/demo_trimmed.bam
echo "featureCounts finished running!"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

############################
# OPTIONAL STEP: View data #
############################

#cat quants/demo_featurecounts.txt.summary

#cat quants/demo_featurecounts.txt

#cat quants/demo_featurecounts.txt | cut -f1,7 | less
