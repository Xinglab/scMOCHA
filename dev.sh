#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-20 15:14:09
# @DESCRIPTION:

# Number of input parameters
param=$#
# Input parameters
conda env create --file=environment.yaml
conda activate scmocha
