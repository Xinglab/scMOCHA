#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-20 15:14:09
# @DESCRIPTION:

# Number of input parameters
param=$#
# Input parameters

targetpath="/home/liuc9/github/scMOCHA/03-ADKP/output"

for tarball in $(ls $targetpath/*.tar.gz); do
  cmd="tar -xzf $tarball -C $targetpath &"
  echo $cmd
  eval $cmd
done
