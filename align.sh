#!/bin/bash

#SBATCH -c 14               # number of core to be used
#SBATCH -t 0-06:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=100000        # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}

# get paths, tools
star=/path/to/STAR-2.5.3a/bin/Linux_x86_64/STAR
genome=/path/to/gencode/26/star/g_all # => ensembl_v88+89

date
echo "Alignment..."
$star --outFileNamePrefix ../../star/$sample. --genomeDir $genome --runThreadN 14 --readFilesIn ../../trim/$sample.R*.gz --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outSAMmode Full --outSAMmapqUnique 60 --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat > ../../bam/$sample.s.bam
echo "Done."

samtools index ../../bam/$sample.s.bam
rm -r ../../star/$sample._STARtmp

