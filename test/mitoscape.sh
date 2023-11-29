#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-10-19 16:21:18
# @DESCRIPTION:

# Number of input parameters
param=$#

# STAR
# https://www.biostars.org/p/9476222/
# https://divingintogeneticsandgenomics.com/post/understand-10x-scrnaseq-and-scatac-fastqs/
# https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
# https://help.geneiousbiologics.com/hc/en-us/articles/4781289585300-Understanding-Single-Cell-technologies-Barcodes-and-UMIs

# #! 1. Align fastq to hte MT rCRS
cellranger count \
  --id=Pei-1 \
  --fastqs=/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-1 \
  --sample=Pei-1 \
  --transcriptome=/mnt/isilon/xing_lab/liuc9/refdata/rCRS/rCRS_cellranger \
  --nosecondary \
  --disable-ui \
  --localcores 10

samtools calmd -e --output-fmt BAM --threads 20 input_MT.bam /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/rCRS_cellranger/fasta/genome.fa >input_MT_MD.bam

/scr1/users/liuc9/tools/bamtofastq input_MT_MD.bam

# cellranger count \
#   --id=Pei-1 \
#   --fastqs=/home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-1 \
#   --sample=Pei-1 \
#   --transcriptome=/mnt/isilon/xing_lab/liuc9/refdata/mitoscape/genome_no_MT_cellranger \
#   --nosecondary \
#   --disable-ui \
#   --localcores 20

# #! 2. Create MD tags
# samtools calmd -e -Q --output-fmt BAM input_MT.bam >input_MT_MD.bam

#! mapping use STAR
cd /home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-2
~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
  --genomeDir /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/rCRS_starindex2.7.9a \
  --readFilesIn Pei-2_S2_L001_R2_001.fastq.gz Pei-2_S2_L001_R1_001.fastq.gz --readFilesCommand zcat \
  --outFileNamePrefix cj_mapping_ \
  --outSAMtype BAM SortedByCoordinate \
  --soloType Droplet \
  --soloCBwhitelist /scr1/users/liuc9/tools/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt \
  --runThreadN 10 \
  --soloUMIlen 12 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR

cd /home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-3

~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
  --genomeDir /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/genome_no_MT_starindex2.7.9a \
  --readFilesIn Pei-3_S3_L001_R2_001.fastq.gz Pei-3_S3_L001_R1_001.fastq.gz \
  --soloType Droplet \
  --readFilesCommand zcat \
  --soloCBwhitelist /scr1/users/liuc9/tools/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt \
  --outFileNamePrefix cj_mapping_no_MT_ \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 10 \
  --soloUMIlen 12 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR

# Concating gz file into one
cd /home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-4

cat Pei-4_S4_L001_R2_001.fastq.gz Pei-4_S4_L002_R2_001.fastq.gz >Pei-4_S1_L001_R2_001.fastq.gz
cat Pei-4_S4_L001_R1_001.fastq.gz Pei-4_S4_L002_R1_001.fastq.gz >Pei-4_S1_L001_R1_001.fastq.gz
cat Pei-4_S4_L001_I1_001.fastq.gz Pei-4_S4_L002_I1_001.fastq.gz >Pei-4_S1_L001_I1_001.fastq.gz

~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
  --genomeDir /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/genome_no_MT_starindex2.7.9a \
  --readFilesIn Pei-4_S4_L001_R2_001.fastq.gz Pei-4_S4_L001_R1_001.fastq.gz \
  --soloType Droplet \
  --readFilesCommand zcat \
  --soloCBwhitelist /scr1/users/liuc9/tools/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt \
  --outFileNamePrefix cj_mapping_no_MT_ \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 10 \
  --soloUMIlen 12 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR

cd /home/liuc9/github/scMOCHA/05-Liming/cellline/torun/Pei-5
cat *R1* >Pei-5_S1_L001_R1_001.fastq.gz
cat *R2* >Pei-5_S1_L001_R2_001.fastq.gz
cat *I1* >Pei-5_S1_L001_I1_001.fastq.gz

~/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
  --genomeDir /mnt/isilon/xing_lab/liuc9/refdata/mitoscape/genome_no_MT_starindex2.7.9a \
  --readFilesIn Pei-5_S1_L001_R2_001.fastq.gz Pei-5_S1_L001_R1_001.fastq.gz \
  --soloType Droplet \
  --readFilesCommand zcat \
  --soloCBwhitelist /scr1/users/liuc9/tools/cellranger-7.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt \
  --outFileNamePrefix cj_mapping_no_MT_ \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 10 \
  --soloUMIlen 12 \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR
