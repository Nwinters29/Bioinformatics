#!/bin/bash

SECONDS=0

# Assign variables
bwa_ref="/home/nwinters29/tiny-test-data/genomes/Hsapiens/hg19/bwa/hg19.fa"
sam_ref="/home/nwinters29/tiny-test-data/genomes/Hsapiens/hg19/seq/hg19.fa"

sample_name=$1
echo $sample_name

# align to reference genome
bwa mem $bwa_ref ${sample_name}_1.fq.gz.bak ${sample_name}_2.fq.gz.bak > $sample_name.sam
echo "alignment finished! SAM file created."

# sort and index
samtools sort ${sample_name}.sam > ${sample_name}_s.bam

# call variants
bcftools mpileup -Ou -f $sam_ref ${sample_name}_s.bam | bcftools call -vmO v -o $sample_name.vcf
echo "variants called! VCF file created."

duration=$SECONDS
echo "$(($duration/60)) minutes and $(($duration%60)) seconds elapsed."
