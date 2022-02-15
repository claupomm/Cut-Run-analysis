#!/bin/bash

#SBATCH -c 4               # number of core to be used
#SBATCH -t 0-01:00          # estimated run-time in D-HH:MM
#SBATCH -p ultrashort            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=1G        # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}

# Setting all samples to one read number (lowest sample size); take out scaffold mapping before; subsampling via samtools; peak calling with macs2 each with corresponding IgG control
# samtools view -q 60 => unique mappings only
echo "Filter alignments on chr1-22XYM and uniquely mapped reads only."
samtools view -@ 4 -q 60 -b ../../bam/$sample.s.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > ../../bam/$sample.chr.bam
samtools index ../../bam/$sample.chr.bam
echo "DONE."

# count filtered mappings per sample
sample_size=$(samtools view -@ 4 -h ../../bam/$sample.chr.bam | wc -l)
echo "Unqiue mapping read counts for $sample: $sample_size"

# calculate scaling factor
factor=$(echo "scale=2; 11602310 / $sample_size" | bc)
echo "scaling factor: $factor"

# downsize to smallest sample size
echo "Downsize to smallest sample size 11602310 (G6_MU_2_H3K4me3)..."
samtools view -@ 4 -b -s 1$factor ../../bam/$sample.chr.bam > ../../bam/$sample.chr.12M.bam
samtools index ../../bam/$sample.chr.12M.bam
echo "DONE."

