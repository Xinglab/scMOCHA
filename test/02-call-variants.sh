#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2022-11-29 11:34:42
# @DESCRIPTION:

# Number of input parameters
param=$#

/home/liuc9/tools/cellsnp-lite-1.2.2/bin/bin/cellsnp-lite \
  -s /home/liuc9/github/scRNAseq-MitoVariant/data/pbmc_10k_v3_possorted_genome_bam.bam \
  -b /home/liuc9/github/scRNAseq-MitoVariant/data/pbmc_10k_v3_filtered_feature_bc_matrix/barcodes.tsv.gz \
  -O /home/liuc9/github/scRNAseq-MitoVariant/data/test/cellsnp \
  --chrom chrMT \
  --UMItag Auto \
  --minMAF 0.1 \
  --minCOUNT 20 \
  --genotype \
  --gzip \
  -p 20