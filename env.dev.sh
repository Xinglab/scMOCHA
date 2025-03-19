#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2023-05-20 15:14:09
# @DESCRIPTION:

# mamba env create --name scmocha --file scmocha.env.yaml
# Number of input parameters
param=$#
# Input parameters
mamba create -n scmocha python -y
mamba activate scmocha
mamba env update -n scmocha --file env.dev.yaml

# mamba env remove -n scmocha
