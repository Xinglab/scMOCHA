#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-19 15:24:18
# @DESCRIPTION:
# @VERSION: v0.0.1

# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Navigate to the script's directory
cd "$SCRIPT_DIR" || exit 1

# Now we're in the script's directory (/home/liuc9/github/scMOCHA/)
echo "Current directory: $(pwd)"

# Export the conda environment to a YAML file
echo "Exporting 'scmocha' environment to env.prod.yaml..."
mamba env export -n scmocha --no-builds | grep -v "prefix:" >env.prod.yaml

# Confirm the export was successful
if [ $? -eq 0 ]; then
  echo "Environment successfully exported to: $SCRIPT_DIR/env.prod.yaml"
else
  echo "Failed to export environment"
  exit 1
fi
