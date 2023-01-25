#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2022-11-29 13:14:32
# @DESCRIPTION:

# Number of input parameters
param=$#
INPUT_DIR=/home/liuc9/github/scRNAseq-MitoVariant/data/PBMC_10k_v3_10x/test/cellsnp
OUT_DIR=/home/liuc9/github/scRNAseq-MitoVariant/data/PBMC_10k_v3_10x/test/mquad
~/tools/anaconda3/bin/mquad -c $INPUT_DIR -o $OUT_DIR -p 20 --minDP 5