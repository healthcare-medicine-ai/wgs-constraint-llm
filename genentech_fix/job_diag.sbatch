#!/bin/bash
#SBATCH --job-name=gerp-diag
#SBATCH --partition=mrivas
#SBATCH --time=00:45:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/diag_%j.out
#SBATCH --error=logs/diag_%j.err

set -euo pipefail

WORKDIR=/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/genentech_fix
PYTHON=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python

cd "$WORKDIR"
export PYTHONUNBUFFERED=1

echo "host      : $(hostname)"
echo "job id    : ${SLURM_JOB_ID:-none}"
echo "started   : $(date -Is)"
echo

"$PYTHON" diag_duplication.py

echo
echo "finished  : $(date -Is)"
