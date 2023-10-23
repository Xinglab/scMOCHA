#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-08-10 16:46:23.96468

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/singularity.slurm.conf \
  -jar /home/liuc9/tools/cromwell-78.jar \
  run /home/liuc9/github/scMOCHA/scMOCHA-image.wdl \
  -i /scr1/users/liuc9/mitochondrial/realdata/03-ADKP/torun/R9443041/R9443041.json 1>/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/torun/R9443041/R9443041.log 2>/scr1/users/liuc9/mitochondrial/realdata/03-ADKP/torun/R9443041/R9443041.err &
