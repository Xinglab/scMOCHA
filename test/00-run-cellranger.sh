#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2022-11-29 05:59:14
# @DESCRIPTION:

# Number of input parameters
param=$#

/home/liuc9/tools/cellranger-7.0.1/bin/cellranger count \
  --id=pbmc_10k_v3_a \
  --fastqs=/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/pbmc_10k_v3_fastqs \
  --sample=pbmc_10k_v3 \
  --transcriptome=/mnt/isilon/xing_lab/liuc9/refdata/cellranger/refdata-gex-GRCh38-2020-A \
  --nosecondary \
  --disable-ui \
  --localcores 20
