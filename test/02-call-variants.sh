#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2022-11-29 11:34:42
# @DESCRIPTION:

# Number of input parameters
param=$#

# /home/liuc9/tools/cellsnp-lite-1.2.2/bin/bin/cellsnp-lite \
#   -s /home/liuc9/github/scMOCHA/data/PBMC_10k_v3_10x/pbmc_10k_v3_possorted_genome_bam.bam \
#   -b /home/liuc9/github/scMOCHA/data/PBMC_10k_v3_10x/pbmc_10k_v3_filtered_feature_bc_matrix/barcodes.tsv.gz \
#   -O /home/liuc9/github/scMOCHA/data/PBMC_10k_v3_10x/test/cellsnp \
#   --chrom chrMT \
#   --UMItag Auto \
#   --minMAF 0.1 \
#   --minCOUNT 20 \
#   --genotype \
#   --gzip \
#   -p 20

/home/liuc9/tools/cellsnp-lite-1.2.2/bin/bin/cellsnp-lite \
  -s /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/MT.bam \
  -b /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/barcodes.tsv \
  -O /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/cellsnp \
  --chrom MT \
  --UMItag Auto \
  --minMAF 0.05 \
  --minCOUNT 10 \
  --genotype \
  --gzip \
  -p 20