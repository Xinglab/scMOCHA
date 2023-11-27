#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-10-19 16:21:18
# @DESCRIPTION:

# Number of input parameters
param=$#

# #! 1. Align fastq to hte MT rCRS
cellranger count \
  --id=Pei-1 \
  --fastqs=/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-1 \
  --sample=Pei-1 \
  --transcriptome=/mnt/isilon/xing_lab/liuc9/refdata/rCRS/rCRS_cellranger \
  --nosecondary \
  --disable-ui \
  --localcores 10

# #! 2. Create MD tags
# samtools calmd -e -Q --output-fmt BAM input_MT.bam >input_MT_MD.bam

#! mapping use STAR
cd /home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-2
~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
  --genomeDir /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/rCRS_starindex2.7.9a \
  --readFilesIn Pei-2_S2_L001_R2_001.fastq.gz Pei-2_S2_L001_R1_001.fastq.gz \
  --soloType Droplet \
  --readFilesCommand zcat \
  --soloCBwhitelist /scr1/users/liuc9/tools/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt \
  --outFileNamePrefix cj_mapping_ \
  --runThreadN 10 \
  --soloUMIlen 12 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR
