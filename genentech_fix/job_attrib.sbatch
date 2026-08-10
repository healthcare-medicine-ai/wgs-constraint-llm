#!/bin/bash
#SBATCH --job-name=gerp-attrib
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/attrib_%x_%j.out
#SBATCH --error=logs/attrib_%x_%j.err

# Attribution runs: isolate each bug so the response letter can state which
# correction moved which gene, rather than inferring it from n_variants.
#
#   sbatch --job-name=gerponly job_attrib.sbatch gerponly
#   sbatch --job-name=amonly   job_attrib.sbatch amonly

set -euo pipefail

MODE="${1:?must pass 'gerponly' or 'amonly'}"

WORKDIR=/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/genentech_fix
PYTHON=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python

cd "$WORKDIR"
export PYTHONUNBUFFERED=1

case "$MODE" in
  gerponly) ARGS="--gerp-fix on  --am-fix off --out-suffix _GERPONLY"; SUF="_GERPONLY" ;;
  amonly)   ARGS="--gerp-fix off --am-fix on  --out-suffix _AMONLY";   SUF="_AMONLY" ;;
  *) echo "unknown mode: $MODE" >&2; exit 2 ;;
esac

echo "host    : $(hostname)"
echo "mode    : $MODE"
echo "started : $(date -Is)"
echo

"$PYTHON" 01_build_regression_input.py $ARGS
echo
"$PYTHON" 02_run_metaregression.py --input-suffix "$SUF" --label "$MODE"

echo
echo "finished: $(date -Is)"
