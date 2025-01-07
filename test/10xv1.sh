#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-04-16 14:50:29
# @DESCRIPTION:

# Number of input parameters
param=$#

# https://kb.10xgenomics.com/hc/en-us/articles/360043386291-How-to-format-v1-chemistry-datasets-to-work-with-current-Cell-Ranger-versions

# rename
# I2 -> I1 8bp sample index
# RA -> R1 98bp transcript
# I1 -> R2 14bp barcode
# R3 -> 10bp UMI

datadir=/mnt/isilon/u01_project/liuc9/From_PT/10X_V1/pbmc3k_fastqs
targetdir=/mnt/isilon/u01_project/liuc9/From_PT/10X_V1/rename
for filename in $(ls ${datadir}/*.fastq.gz); do
  # echo ${filename}
  basefilename=$(basename ${filename})
  # echo ${basefilename}
  readtype=$(echo ${basefilename} | awk -F '[_-]' '{print $2}')
  case ${readtype} in
  I2)
    newreadtype=I1
    ;;
  RA)
    newreadtype=R1
    ;;
  I1)
    newreadtype=R2
    ;;
  esac
  si_name=$(echo ${basefilename} | awk -F '[_-]' '{print $4}')
  lane=$(echo ${basefilename} | awk -F '[_-]' '{print $6}')
  newname=si_${si_name}_S1_L${lane}_${newreadtype}_001.fastq.gz
  # ln -s ${filename} ${targetdir}/${newname}
  if [ ${newreadtype} = "R1" ]; then
    # echo ${basefilename} ${newname}
    # echo ${newreadtype} ${basefilename} R1
    # split RA into R1 and R3
    r1filename=si_${si_name}_S1_L${lane}_R1_001.fastq.gz
    r3filename=si_${si_name}_S1_L${lane}_R3_001.fastq.gz
    zcat "$filename" | awk 'NR % 8 == 1 || NR % 8 == 2 || NR % 8 == 3 || NR % 8 == 4 {print $0}' | gzip >"${targetdir}/${r1filename}" &
    zcat "$filename" | awk 'NR % 8 == 5 || NR % 8 == 6 || NR % 8 == 7 || NR % 8 == 0 {print $0}' | gzip >"${targetdir}/${r3filename}" &
  else
    # echo ${newreadtype} ${basefilename} ${newname}
    ln -s ${filename} ${targetdir}/${newname}
  fi
done

output_id=pbmc_v1data_si_ACGCGGAA
fastqs=/mnt/isilon/u01_project/liuc9/From_PT/10X_V1/rename
sample_id=si_ACGCGGAA
transcriptome=/home/liuc9/data/refdata/mgatk_index/Human
cpu=10

cellranger count \
  --id=${output_id} \
  --fastqs=${fastqs} \
  --sample=${sample_id} \
  --transcriptome=${transcriptome} \
  --nosecondary \
  --disable-ui \
  --localcores ${cpu} &
