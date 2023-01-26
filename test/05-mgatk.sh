#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-01-23 16:05:43
# @DESCRIPTION:

# Number of input parameters
param=$#

# ~/tools/anaconda3/bin/mgatk tenx -i /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/pbmc_10k_v3_possorted_genome_bam.bam \
#   -n CRR_test1 -o CRR_test1_mgatk -c 12 -ub UB \
#   -bt CB -b /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/pbmc_10k_v3_filtered_feature_bc_matrix/barcodes.tsv.gz \
#   --mito-genome /home/liuc9/scratch/mitochondrial/rCRS/NC_012920.rCRS.fasta

# ~/tools/anaconda3/bin/mgatk tenx -i /scr1/users/liuc9/tools/mgatk/tests/barcode/test_barcode.bam \
#   -n CRR_test1 -o CRR_test1_mgatk -c 12 -ub UB \
#   -bt CB -b /scr1/users/liuc9/tools/mgatk/tests/barcode/test_barcodes.txt


mgatk tenx -i /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/MT.bam \
  -o /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk \
  -n mgatk \
  -c 12 -ub UB  -bt CB \
  -b /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/barcodes.tsv \
  --mito-genome /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/fasta/rCRS.fasta