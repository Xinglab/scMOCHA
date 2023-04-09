#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-02-23 16:54:40
# @DESCRIPTION:

# Number of input parameters
# param=$#
# cd /home/liuc9/tmp/mito/flu2
module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scRNAseq-MitoVariant/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scRNAseq-MitoVariant/scmtah.wdl \
    -i /home/liuc9/github/scRNAseq-MitoVariant/cellranger-a.json &

# java -jar /home/liuc9/tools/cromwell-78.jar \
#     run /home/liuc9/github/scRNAseq-MitoVariant/cellranger.wdl \
#     -i /home/liuc9/github/scRNAseq-MitoVariant/cellranger.json

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scRNAseq-MitoVariant/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scRNAseq-MitoVariant/scmtah.wdl \
    -i /home/liuc9/github/scRNAseq-MitoVariant/cellranger-flu5-a.json &

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scRNAseq-MitoVariant/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scRNAseq-MitoVariant/scmtah.wdl \
    -i /home/liuc9/github/scRNAseq-MitoVariant/cellranger-mm068.json &

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scRNAseq-MitoVariant/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scRNAseq-MitoVariant/scmtah.wdl \
    -i /home/liuc9/github/scRNAseq-MitoVariant/cellranger-mm188.json &
