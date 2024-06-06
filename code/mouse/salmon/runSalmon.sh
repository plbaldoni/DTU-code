#!/bin/bash
#SBATCH --mem=66g
#SBATCH --cpus-per-task=11
#SBATCH --time=12:00:00

module load salmon/1.10.2

salmonindex=../../../output/mouse/salmon-index/transcripts_index/

mkdir -p ../../../output/mouse/salmon/

for filename in ../../../data/mouse/fastq/*_R1.fastq.gz; do

outname=${filename%_R1.fastq.gz}
outname=${outname##*/}

salmon quant \
-i $salmonindex \
-l A \
-1 $filename \
-2 ${filename/_R1.fastq.gz/_R2.fastq.gz} \
-p 10 \
--numGibbsSamples 100 \
--dumpEq \
-o ../../../output/mouse/salmon/$outname

done
