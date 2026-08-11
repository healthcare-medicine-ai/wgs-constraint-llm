#!/bin/bash
#SBATCH --job-name=manifest
#SBATCH --partition=mrivas
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
export PYTHONUNBUFFERED=1
/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python pipelines/00_input_manifest.py
