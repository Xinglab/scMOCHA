#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-03-19 13:53:11
# @DESCRIPTION:

# Number of input parameters
param=$#

output_id="Flu2"
fastqs="/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2"
sample_id="Flu2"
transcriptome="/home/liuc9/data/refdata/cellranger/refdata-gex-GRCh38-2020-A"
cpu=20

cellranger count \
  --id=${output_id} \
  --fastqs=${fastqs} \
  --sample=${sample_id} \
  --transcriptome=${transcriptome} \
  --nosecondary \
  --disable-ui \
  --localcores ${cpu}