#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-08 16:50:08

module load Java/15.0.1
nohup java -Dconfig.file=/home/liuc9/github/scMOCHA/config/ref.conf \
-jar /home/liuc9/tools/cromwell-78.jar \
run /home/liuc9/github/scMOCHA/scMOCHA.wdl \
-i /scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/SRR11680207/SRR11680207.json 1>/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/SRR11680207/SRR11680207.log 2>/scr1/users/liuc9/mitochondrial/realdata/01-Sci_Immunol_32651212/SRR11680207/SRR11680207.err &
