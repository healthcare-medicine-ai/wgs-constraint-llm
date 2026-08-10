#!/bin/bash
#SBATCH --job-name=verify-liftover
#SBATCH --partition=mrivas
#SBATCH --time=01:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/verifylift_%j.out
#SBATCH --error=logs/verifylift_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/genentech_fix
export PYTHONUNBUFFERED=1
echo "host: $(hostname)"; echo "started: $(date -Is)"; echo
/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python verify_liftover.py \
  --chroms chr1,chr7,chr12,chr20 --per-chrom 150
echo; echo "finished: $(date -Is)"
