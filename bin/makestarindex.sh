#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-11-26 19:38:37
# @DESCRIPTION:

# Number of input parameters
param=$#

~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 80 \
  --runMode genomeGenerate \
  --genomeDir /home/liuc9/data/refdata/mitoscape/rCRS_starindex2.7.9a \
  --genomeFastaFiles /mnt/isilon/xing_lab/liuc9/refdata/rCRS/NC_012920.rCRS.fasta \
  --sjdbGTFfile /mnt/isilon/xing_lab/liuc9/refdata/rCRS/Homo_sapiens.GRCh38.107.MT.gtf \
  --sjdbOverhang 100

seqkit grep -vrp "MT" genome.fa >genome.no_MT.fa

~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 80 \
  --runMode genomeGenerate \
  --genomeDir /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/genome_no_MT_starindex2.7.9a \
  --genomeFastaFiles /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/genome.no_MT.fa \
  --sjdbGTFfile /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/Homo_sapiens.GRCh38.107.new.no_MT.gtf \
  --sjdbOverhang 100
