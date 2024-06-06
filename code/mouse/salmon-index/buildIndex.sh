#!/bin/bash

#SBATCH --mem=200g
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# Creating decoy-aware transcriptome as descripted in the first bullet point of
# https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode

anndir=../../../data/annotation/mm39/
destdir=../../../output/mouse/salmon-index/
salmon=~/lab_smyth/baldoni.p/software/salmon-1.10.0_linux_x86_64/bin/salmon
salmondecoy=~/lab_smyth/baldoni.p/software/SalmonTools/scripts/generateDecoyTranscriptome.sh
mashmap=~/lab_smyth/baldoni.p/software/MashMap/mashmap
bedtools=/stornext/System/data/apps/bedtools/bedtools-2.26.0/bin/bedtools

mkdir -p $destdir

cp ${anndir}GRCm39.genome.fa.gz \
  ${anndir}gencode.vM35.annotation.gtf.gz \
  ${anndir}gencode.vM35.transcripts.fa.gz ./

gunzip GRCm39.genome.fa.gz
gunzip gencode.vM35.annotation.gtf.gz
gunzip gencode.vM35.transcripts.fa.gz

$salmondecoy \
-j 18 \
-g GRCm39.genome.fa \
-t gencode.vM35.transcripts.fa \
-a gencode.vM35.annotation.gtf \
-m $mashmap \
-b $bedtools \
-o $destdir

# Building index

$salmon index \
-t ${destdir}gentrome.fa \
-i ${destdir}transcripts_index \
-d ${destdir}decoys.txt \
-k 31 \
-p 18

# Removing extra files

rm -rf GRCm39.genome.fa gencode.vM35.transcripts.fa gencode.vM35.annotation.gtf ${destdir}gentrome.fa ${destdir}decoys.txt
