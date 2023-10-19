#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-10-19 16:21:18
# @DESCRIPTION:

# Number of input parameters
param=$#

#! 1. Align fastq to hte MT rCRS
cellranger count \
  --id=R1571846 \
  --fastqs=/home/liuc9/scratch/tmp/scmochatest/R1571846 \
  --sample=R1571846 \
  --transcriptome=/mnt/isilon/xing_lab/liuc9/refdata/rCRS/rCRS_cellranger \
  --nosecondary \
  --disable-ui \
  --localcores 10

#! 2. Create MD tags
samtools calmd -e -Q --output-fmt BAM input_MT.bam >input_MT_MD.bam
