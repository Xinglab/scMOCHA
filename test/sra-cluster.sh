#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-03-19 17:49:23
# @DESCRIPTION:

# Number of input parameters
param=$#

# Flu2
# cell_cluster_annotation
h5file="/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/filtered_feature_bc_matrix.h5"

Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/cellcluster_10x.R ${h5file}

# call_variant_on_single_cell_level

possorted_genome_bam=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/possorted_genome_bam.bam

gzipped_barcodes=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
barcodes=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/filtered_feature_bc_matrix/barcodes.tsv

cpu=20
rCRS=/home/liuc9/github/scRNAseq-MitoVariant/fasta/rCRS.chrM.fasta

barcode_cluster=/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/barcode_cluster.tsv
barcode_bulk=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu2/Flu2/outs/barcode_bulk.tsv

gunzip -c ${gzipped_barcodes} > ${barcodes}
mgatk tenx -i ${possorted_genome_bam} \
  -n sc \
  -c ${cpu} -ub UB  -bt CB \
  -b ${barcodes} \
  --mito-genome ${rCRS}

tar czf mgatk_single_cell_level.tar.gz "mgatk_out/final"


# call_variant_on_cell_cluster_level
samtools view -hb ${possorted_genome_bam} chrM > MT.bam
samtools index MT.bam

sinto addtags \
  -b MT.bam \
  -f ${barcode_cluster} \
  -o MT_cluster.bam \
  -p ${cpu}

samtools index MT_cluster.bam

mgatk bcall -i MT_cluster.bam \
    -o mgatk_cluster \
    -n mgatk_cluster \
    -c ${cpu} -bt CJ \
    --mito-genome ${rCRS} \
    --keep-temp-files

python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling.py mgatk_cluster/final/ mgatk_cluster 16569 10 chrM

sinto addtags \
  -b MT.bam \
  -f ${barcode_bulk} \
  -o MT_bulk.bam \
  -p ${cpu}

samtools index MT_bulk.bam

mgatk bcall -i MT_bulk.bam \
  -o mgatk_bulk \
  -n mgatk_bulk \
  -c ${cpu} -bt CJ \
  --mito-genome ${rCRS} \
  --keep-temp-files

python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling.py mgatk_bulk/final/ mgatk_bulk 16569 10 chrM


# flu5
# cell_cluster_annotation
h5file="/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs/filtered_feature_bc_matrix.h5"

Rscript /home/liuc9/github/scRNAseq-MitoVariant/bin/cellcluster_10x.R ${h5file}

# call_variant_on_single_cell_level


possorted_genome_bam=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs/possorted_genome_bam.bam

gzipped_barcodes=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
barcodes=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs/filtered_feature_bc_matrix/barcodes.tsv

cpu=20
rCRS=/home/liuc9/github/scRNAseq-MitoVariant/fasta/rCRS.chrM.fasta

barcode_cluster=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs/barcode_cluster.tsv
barcode_bulk=/scr1/users/liuc9/mitochondrial/testdata/1_old_donor_pbmc/flu5/Flu5/outs/barcode_bulk.tsv

gunzip -c ${gzipped_barcodes} > ${barcodes}
mgatk tenx -i ${possorted_genome_bam} \
  -n sc \
  -c ${cpu} -ub UB  -bt CB \
  -b ${barcodes} \
  --mito-genome ${rCRS}

tar czf mgatk_single_cell_level.tar.gz "mgatk_out/final"


# call_variant_on_cell_cluster_level
samtools view -hb ${possorted_genome_bam} chrM > MT.bam
samtools index MT.bam

sinto addtags \
  -b MT.bam \
  -f ${barcode_cluster} \
  -o MT_cluster.bam \
  -p ${cpu}
samtools index MT_cluster.bam

mgatk bcall -i MT_cluster.bam \
    -o mgatk_cluster \
    -n mgatk_cluster \
    -c ${cpu} -bt CJ \
    --mito-genome ${rCRS} \
    --keep-temp-files

python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling.py mgatk_cluster/final/ mgatk_cluster 16569 10 chrM

sinto addtags \
  -b MT.bam \
  -f ${barcode_bulk} \
  -o MT_bulk.bam \
  -p ${cpu}

samtools index MT_bulk.bam

mgatk bcall -i MT_bulk.bam \
  -o mgatk_bulk \
  -n mgatk_bulk \
  -c ${cpu} -bt CJ \
  --mito-genome ${rCRS} \
  --keep-temp-files

python /home/liuc9/github/scRNAseq-MitoVariant/bin/variant_calling.py mgatk_bulk/final/ mgatk_bulk 16569 10 chrM
