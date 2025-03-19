#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-19 15:44:32
# @DESCRIPTION:
# @VERSION: v0.0.1

param=$#
# Input parameters
if [ "$param" -lt 1 ]; then
  echo "Usage: $0 <environment_name>"
  echo "This script creates a conda environment from env.prod.yaml file"
  exit 1
fi

env_name=$1

echo "Creating conda environment '$env_name' from env.prod.yaml..."
mamba env create -n "$env_name" -f env.prod.yaml

if [ $? -eq 0 ]; then
  echo "Environment '$env_name' successfully created."
  echo "Activate with: conda activate $env_name"
else
  echo "Failed to create environment '$env_name'."
  exit 1
fi
