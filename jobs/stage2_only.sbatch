#!/bin/bash
#SBATCH --partition=mrivas
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
export PYTHONUNBUFFERED=1
/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python \
  pipelines/02_run_metaregression.py --input-suffix "$1" --label "$2"
