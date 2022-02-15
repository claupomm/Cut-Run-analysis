# Cut&Run-analysis
Instructions for analysis of Cut&Run-sequencing data. Cut&Run-seq data were generated via CUT&RUN (Cleavage Under Target and Release Using Nuclease), for which 40-200 bp fragment length are expected. Sequencing was done on a NextSeq500 in-house. This Cut&Run workflow includes demultiplexing, trimming, mapping, and peak calling. More descriptions and interpretations have been published (https://www.nature.com/articles/s41598-020-66224-1#Sec9) and stored at ArrayExpress/ENA (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8203). 


# 1. trimming
# 2. mapping
# 3. keep all duplicates, but count those
# 4. access to bam files


## Prepare the project folder and fastq files
Set the Project folder:
```
DIR=/path/to/Project_folder
mkdir $DIR
cd $DIR
# work on sample sheet: Sample_Name later as folder name in Project folder "raw"
```

## Demultiplex sequences
Prepare the NextSeq_SampleSheet.csv to your needs. Creating the fastq files of the NextSeq run in the following will allow for one mismatch in the barcode and be submitted to job management system SLURM.

```
# set variables
DIR=/path/to/Project_folder
nextseq_folder=date_xxx_xxx_xxx
samples=NextSeq_Samplesheet.csv
p=ultrashort # p=short <6h, p=mid <2d, p=long <4d
c=20 # number of core to be used

# demultiplex with 1 mismatch in barcode
mkdir /path/to/NextSeq/demux/$nextseq_folder/
cd /path/to/NextSeq/demux/$nextseq_folder/
# submit job on slurm
sbatch -c $c -p $p -J bcl2fq -o bcl2fq.out -e bcl2fq.err <<EOF
#!/bin/sh
bcl2fastq --fastq-compression-level 6 -r 20 -d 20 -p 20 -w 20 --barcode-mismatches 1 --runfolder-dir /path/to/NextSeq/raw/$nextseq_folder/ --output-dir /path/to/NextSeq/demux/$nextseq_folder/ --no-lane-splitting --sample-sheet $DIR/$samples
EOF
```

## Create the subfolders with samples
Save or link fastq files in the corresponding sample folder.
```
mkdir $DIR $DIR/raw
cd $DIR/raw
mkdir D96_HG3_IgG D96_HG3_CTCF D96_HG3_H3K4me3 D96_HG3_H3K27ac
# take project name from samplesheet
raw=/path/to/NextSeq/demux/$nextseq_folder/D96 
for i in $(find $DIR/raw/* -maxdepth 1 -type d); do
cd $i
ln -s $raw/*/${i#${DIR}/raw/}*fastq.gz .
done

# Get sample names (string)
SAMPLES=$(find $DIR/raw/* -maxdepth 1 -type d)
cd $DIR
mkdir trim fastqc star plots bam macs2 homer
```

## Start pipeline on SLURM, iterate each sample
```
for SAMPLE in $SAMPLES
do
cd $SAMPLE
sample=${PWD##*/}
# trim sequences
# quality control via fastQC, plots, 50min
RES=$(sbatch -J $sample.1 -o $sample.1.out -e $sample.1.err ../../trim.sh)
# alignment without joining sequences, 30min, include 
RES2=$(sbatch --dependency=afterok:${RES##* } -J $sample.2 -o $sample.2.out -e $sample.2.err ../../align.sh)
# count duplicates => bam/$sample.metrics
RES3=$(sbatch --dependency=afterok:${RES2##* } -J $sample.3 -o $sample.3.out -e $sample.3.err ../../picard_metrics.sh)
done
```

## Alignment statistics
```
cd $DIR/star
grep "mapped reads %" *final* # 85-89
grep "mapped reads n" *final* # 30-42M
grep "Number of input reads" *final* # 35-47M
grep "% of reads mapped to multiple loci" *.Log.final.out # 7-8%

# Graphical info on sequencing + mapping
cd $DIR
R --file="stats_seq_rna.R"
```

## Peak calling via macs2
Cut&Run-seq peak calling was done via macs2 for D06_HG3_CTCF + D96_HG3_H3K4me3 + D96_HG3_H3K27ac whereas D96_HG3_IgG served as control.
```
cd $DIR
for sample in D06_HG3_CTCF D96_HG3_H3K4me3 D96_HG3_H3K27ac; do
sample=D96_HG3_H3K4me3
cd $DIR/raw/$sample
sbatch -c 1 -p ultrashort -J $sample.4 -o $sample.4.out -e $sample.4.err <<EOF
#!/bin/sh
macs2 callpeak -t ../../bam/$sample.s.bam -n ../../macs2/$sample --gsize hs -c ../../bam/D96_HG3_IgG.s.bam
EOF
done
```

## Peak calling via HOMER
Peak calling preparation
```
for SAMPLE in $SAMPLES
do
cd $SAMPLE
sample=${PWD##*/}
sbatch -J $sample.5 -o $sample.5.out -e $sample.5.err ../../homer_preparation.sh
done
```
Actual peak calling
```
homer=~/Programme/homer_v4.8/bin
cd $DIR/homer
for sample in D06_HG3_CTCF D96_HG3_H3K4me3 D96_HG3_H3K27ac; do
sbatch -c 1 -p ultrashort -J $sample.h -o $sample.h.out -e $sample.h.err <<EOF
#!/bin/sh
$homer/findPeaks $sample/ -style histone -o auto -i D96_HG3_IgG/ -size 500 -minDist 1000
EOF
done
```
Playing around with parameters...
```
for sample in D06_HG3_CTCF D96_HG3_H3K4me3 D96_HG3_H3K27ac; do
cd $DIR/raw/$sample
sbatch -c 1 -p ultrashort -J $sample.3m -o $sample.3m.out -e $sample.3m.err <<EOF
#!/bin/sh
macs2 callpeak -t ../../bam/$sample.s.bam -n ../../macs2/$sample.IgG --gsize hs -c ../../bam/D96_HG3_IgG.s.bam --nomodel
EOF
done
```

Take out scaffold mapping before; subsampling via samtools; peak calling with macs2 each with corresponding IgG control; count filtered mappings per sample and use unique mappings only; calculate scaling factor; setting all samples to one read number (lowest sample size).

```
for sample in D96_HG3_IgG D06_HG3_CTCF D96_HG3_H3K4me3 D96_HG3_H3K27ac; do
cd $DIR/raw/$sample
RES=$(sbatch -c 4 -p ultrashort -J $sample.3.2 -o $sample.3.2.out -e $sample.3.2.err ../../count_sample_size_scaling.sh)
sbatch -c 1 -p ultrashort --dependency=afterok:${RES##* } -J $sample.8.6 -o $sample.8.6.out -e $sample.8.6.err <<EOF
#!/bin/sh
macs2 callpeak --qvalue 0.1 --keep-dup all --broad -t ../../bam/$sample.chr.12M.bam -n ../../macs2/$sample.12M.q01.dupAll.broad. --gsize hs -c ../../bam/D96_HG3_IgG.chr.12M.bam 
EOF
done
```


## Upload to arrayexpress/ENA
```
nextseq_folder=date_xxx_xxx_xxx
cd /path/to/NextSeq/demux/$nextseq_folder/D96
ln -s */D96_HG3_CTCF_S3_R1_001.fastq.gz . # all fastq files in one folder
ln -s */D96_HG3_H3K4me3_S4_R1_001.fastq.gz .
ln -s */D96_HG3_H3K27ac_S5_R1_001.fastq.gz .
ln -s */D96_HG3_IgG_S1_R1_001.fastq.gz .
ftp ftp-private.ebi.ac.uk
Name: aexpress # => ftp field in ftp upload, ftp/aspera upload
Password: aexpress1
cd xxxxxxxxxxx/ # as given in ftp/aspera upload
prompt # interactive mode off
mput *gz
exit
# calculate md5sum
md5sum *gz
8d50a41f204243ad8a272bde164dd044  D96_HG3_CTCF_S3_R1_001.fastq.gz
ad272792e07b2304758f71d5d4d53b7b  D96_HG3_H3K27ac_S5_R1_001.fastq.gz
7ae683327946c3c987e92f132552cd8f  D96_HG3_H3K4me3_S4_R1_001.fastq.gz
4c841e868e305e90405c3320bded7c3c  D96_HG3_IgG_S1_R1_001.fastq.gz
rm *gz
```




