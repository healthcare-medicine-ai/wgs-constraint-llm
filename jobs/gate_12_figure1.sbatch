#!/bin/bash
#SBATCH --job-name=gate12b
#SBATCH --partition=mrivas
#SBATCH --time=02:00:00
#SBATCH --mem=96G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/12_figure1_scn1a.py --suffix _REGEN
echo
$PY genentech_fix/compare_one_figure.py "Figure 1: observed vs predicted mutations for SCN1A" _REGEN
