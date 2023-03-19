#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-03-19 13:47:48
# @DESCRIPTION:

# Number of input parameters
param=$#


sra_dir="/home/liuc9/scratch/mitochondrial/testdata/1_old_donor_pbmc"


fasterq-dump --mem 50G --threads 20 --split-3 --gzip ${sra_dir}/flu2/SRR11680210.sra

fasterq-dump --mem 50G --threads 20 --split-3 --gzip ${sra_dir}/flu5/SRR11680214.sra


# fasterq-dump -S --include-technical

mv SRR11680210_1.fastq SRR11680210_S1_L001_R1_001.fastq

mv SRR11680210_2.fastq SRR11680210_S1_L001_R2_001.fastq


mv SRR11680214_1.fastq SRR11680214_S1_L001_R1_001.fastq

mv SRR11680214_2.fastq SRR11680214_S1_L001_R2_001.fastq