#!/bin/bash

#SBATCH --mem=32g
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH -p long

prj=/vast/projects/diffSplice/DTU-code

module load nextflow

nextflow pull plbaldoni/nextflow-rnaseq

nextflow run plbaldoni/nextflow-rnaseq -with-report report.html -with-trace -with-timeline timeline.html -resume -r main \
  --align --subjunc \
  --reads "$prj/data/mouse/fastq/*{R1,R2}.fastq.gz" \
  --outdir "$prj/output/mouse/subread/" \
  --subreadAnnoType "GTF" \
  --subreadAnno "$prj/data/annotation/mm39/gencode.vM35.annotation.gtf.gz" \
  --subreadIndex "$prj/output/mouse/subread-index/GRCm39.genome" \
  --subreadGenome "$prj/data/annotation/mm39/GRCm39.genome.fa.gz" \
  --gsize 2410055689

