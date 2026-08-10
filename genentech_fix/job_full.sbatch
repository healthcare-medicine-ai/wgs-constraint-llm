#!/bin/bash
#SBATCH --job-name=gerp-full
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/full_%x_%j.out
#SBATCH --error=logs/full_%x_%j.err

# Usage:
#   sbatch --job-name=repro job_full.sbatch control   # reproduces submitted numbers
#   sbatch --job-name=fixed job_full.sbatch fixed     # both corrections applied
#
# The control run exists to prove the port is faithful. If it does not
# reproduce results/epilepsy_unified_model_pvalues.tsv from baseline_R2,
# nothing downstream should be believed.

set -euo pipefail

MODE="${1:?must pass 'control' or 'fixed'}"

WORKDIR=/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/genentech_fix
PYTHON=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python

cd "$WORKDIR"
export PYTHONUNBUFFERED=1

case "$MODE" in
  control) BUILD_ARGS="--no-fix --out-suffix _REPRO"; SUFFIX="_REPRO"; LABEL="REPRO (no fix)" ;;
  fixed)   BUILD_ARGS="--out-suffix _FIXED";          SUFFIX="_FIXED"; LABEL="FIXED" ;;
  *) echo "unknown mode: $MODE" >&2; exit 2 ;;
esac

echo "host    : $(hostname)"
echo "mode    : $MODE"
echo "started : $(date -Is)"
echo

"$PYTHON" 01_build_regression_input.py $BUILD_ARGS
echo
"$PYTHON" 02_run_metaregression.py --input-suffix "$SUFFIX" --label "$LABEL"

echo
echo "finished: $(date -Is)"
