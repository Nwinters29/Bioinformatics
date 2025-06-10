#!/usr/bin/env python3

# Import necessary packages
import sys
import argparse
import subprocess

def run(command):
    """ Function used to run the subprocess module for each created command in the main function """
    process = subprocess.run(command, capture_output = True, shell = True)
    return process

def parse_args():
    """ Argument parser function that is used to extract the base level name of the input FASTQ files """
    parser = argparse.ArgumentParser(description="Call variants and return statistics")
    parser.add_argument('-i', '--input_sample_name', type = str, required=True, help='Base sample name')
    return parser.parse_args()
    
def main():
    """ Function that runs bash commands for aligning reads, sorting the alignment, and calling variants """

    # Assign variables
    args = parse_args()
    sample = args.input_sample_name
    bwa_ref = "~/tiny-test-data/genomes/Hsapiens/hg19/bwa/hg19.fa"
    sam_ref = "~/tiny-test-data/genomes/Hsapiens/hg19/seq/hg19.fa"

    # 1: Align FASTQ files to the referencde genome using BWA
    align = f"bwa mem {bwa_ref} {sample}_1.fq.gz.bak {sample}_2.fq.gz.bak > {sample}.sam"
    aligned = run(align)

    # 2: Sort aligned reads by name
    collate = f"samtools collate -o {sample}-c.bam {sample}.sam"
    collated = run(collate)

    # 3: Fill in mate coordinates
    fixmate = f"samtools fixmate -m {sample}-c.bam {sample}-f.bam"
    fixmated = run(fixmate)

    # 4: Sort reads by coordinate
    sort = f"samtools sort {sample}-f.bam.sam > {sample}-s.bam"
    sorted = run(sort)

    # 5: Remove duplicates
    markdup = f"samtools markdup -r {sample}-s.bam {sample}-nodup.bam"
    marked = run(markdup)

    # 6: Call variants
    call = f"bcftools mpileup -Ou -f {sam_ref} {sample}-nodup.bam | bcftools call -vmO v -o {sample}.vcf"
    called = run(call)

    # 7: return Statistics
    stats = f"bcftools stats {sample}.vcf"
    statistics = run(stats)

    # Output
    print(f"Pipeline completed. Output files: {sample}.sam, {sample}-c.bam, {sample}-f.bam, {sample}-s.bam, {sample}-nodup.bam, {sample}.vcf")
    return(statistics.stdout.decode())

if __name__ == '__main__':
    sys.exit(main())