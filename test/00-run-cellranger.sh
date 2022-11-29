#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2022-11-29 05:59:14
# @DESCRIPTION:

# Number of input parameters
param=$#

/home/liuc9/tools/cellranger-7.0.1/bin/cellranger count \
  --id=pbmc_10k_v3 \
  --transcriptome=/home/liuc9/data/refdata/cellranger/refdata-gex-GRCh38-2020-A \
  --fastqs=/home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/pbmc_10k_v3_fastqs \
  --sample=pbmc_10k_v3 \
  --nosecondary
