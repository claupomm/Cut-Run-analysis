#!/bin/bash

#SBATCH -c 8               # number of core to be used
#SBATCH -t 0-06:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=10000        # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}


# get paths, tools
homer=/path/to/homer_v4.8/bin/

date
echo "Prepare tags for homer, only common chromosomes, without scaffolds..."
samtools view -b -@ 8 ../../bam/$sample.s.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > ../../bam/$sample.chr.bam
echo "Done."

date
echo "Make tag directory required for findPeak..."
$homer/makeTagDirectory ../../homer/$sample ../../bam/$sample.chr.bam -tbp 1
echo "Done."
date
