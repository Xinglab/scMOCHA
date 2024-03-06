#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-02-23 16:54:40
# @DESCRIPTION:

# Number of input parameters
# param=$#
# cd /home/liuc9/tmp/mito/flu2
module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scMOCHA/scmtah.wdl \
    -i /home/liuc9/github/scMOCHA/cellranger-a.json &

# java -jar /home/liuc9/tools/cromwell-78.jar \
#     run /home/liuc9/github/scMOCHA/cellranger.wdl \
#     -i /home/liuc9/github/scMOCHA/cellranger.json

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scMOCHA/scmtah.wdl \
    -i /home/liuc9/github/scMOCHA/cellranger-flu5-a.json &

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scMOCHA/scmtah.wdl \
    -i /home/liuc9/github/scMOCHA/cellranger-mm068.json &

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/ref.conf \
    -jar /home/liuc9/tools/cromwell-78.jar \
    run /home/liuc9/github/scMOCHA/scmtah.wdl \
    -i /home/liuc9/github/scMOCHA/cellranger-mm188.json &
