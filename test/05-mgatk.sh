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

# for single cell variant calling
mgatk tenx -i /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/MT.bam \
  -o /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk \
  -n mgatk \
  -c 12 -ub UB  -bt CB \
  -b /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/barcodes.tsv \
  --mito-genome /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/fasta/rCRS.fasta



# for cell cluster variant calling

sinto addtags \
  -b /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT.bam \
  -f /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/barcode_cluster.tsv \
  -o /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT_cluster.bam \
  -p 40

samtools index /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT_cluster.bam

mgatk bcall -i /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT_cluster.bam \
  -o /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_cluster \
  -n mgatk_cluster \
  -c 40 -bt CJ \
  --mito-genome /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/fasta/rCRS.fasta \
  --keep-temp-files

python /home/liuc9/github/scMOCHA/test/variant_calling.py /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_cluster/final/ mgatk_cluster 16569 10 MT

# for bulk cell clustering

sinto addtags \
  -b /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT.bam \
  -f /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/barcode_bulk.tsv \
  -o /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT_bulk.bam \
  -p 40

samtools index /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT_bulk.bam


mgatk bcall -i /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/MT_bulk.bam \
  -o /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_bulk \
  -n mgatk_bulk \
  -c 40 -bt CJ \
  --mito-genome /scr1/users/liuc9/mitochondrial/PBMC_10k_v3_10x/MTbam/fasta/rCRS.fasta \
  --keep-temp-files

python /home/liuc9/github/scMOCHA/test/variant_calling.py /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/mgatkmtbam_cluster/mgatk_bulk/final/ mgatk_bulk 16569 10 MT


# test
python /home/liuc9/github/scMOCHA/test/variant_calling.py /home/liuc9/scratch/mitochondrial/PBMC_10k_v3_10x/MTbam/mgatk/final mgatk_bulk 16569 10 MT
