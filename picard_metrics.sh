#!/bin/bash

#SBATCH -c 3                # number of core to be used
#SBATCH -t 0-06:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=20G         # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}

picard=/path/to/picard-2.9.2/picard.jar
cd ../../bam


date
echo "Estimate duplicate reads via picard..."
java -XX:ParallelGCThreads=2 -Xmx20g -jar $picard MarkDuplicates I=$sample.s.bam O=$sample.rmdup.bam METRICS_FILE=$sample.metrics ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE
echo "Done."
date

rm $sample.rmdup.bam
