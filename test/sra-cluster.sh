#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-03-19 17:49:23
# @DESCRIPTION:

# Number of input parameters
param=$#


h5file="/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/filtered_feature_bc_matrix.h5"

Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/cellcluster_10x.R ${h5file}

gunzip -c ${gzipped_barcodes} > ${barcodes}

mgatk tenx -i ${possorted_genome_bam} \
  -n sc \
  -c ${cpu} -ub UB  -bt CB \
  -b ${barcodes} \
  --mito-genome ${rCRS}

tar czf mgatk_single_cell_level.tar.gz "final"