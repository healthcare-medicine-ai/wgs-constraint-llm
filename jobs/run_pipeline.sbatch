#!/bin/bash
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Run the epilepsy meta-regression end to end.
#
#   sbatch --job-name=v2repro jobs/run_pipeline.sbatch control
#   sbatch --job-name=v2fixed jobs/run_pipeline.sbatch fixed
#   sbatch --job-name=v2gerp  jobs/run_pipeline.sbatch gerponly
#   sbatch --job-name=v2am    jobs/run_pipeline.sbatch amonly
#
# `control` reproduces the Revision 2 submission and must match
# genentech_fix/baseline_R2 exactly. If it does not, nothing else is
# interpretable -- treat that as a hard failure, not a discrepancy.
set -euo pipefail
MODE="${1:?pass control|fixed|gerponly|amonly}"
REPO=/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PYTHON=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
cd "$REPO"
export PYTHONUNBUFFERED=1
case "$MODE" in
  control)  ARGS="--no-fix --out-suffix _REPRO";                    SUF=_REPRO ;;
  fixed)    ARGS="--out-suffix _FIXED";                             SUF=_FIXED ;;
  gerponly) ARGS="--gerp-fix on --am-fix off --out-suffix _GERPONLY"; SUF=_GERPONLY ;;
  amonly)   ARGS="--gerp-fix off --am-fix on --out-suffix _AMONLY";   SUF=_AMONLY ;;
  *) echo "unknown mode $MODE" >&2; exit 2 ;;
esac
echo "host: $(hostname)  mode: $MODE  started: $(date -Is)"; echo
$PYTHON pipelines/01_build_regression_input.py $ARGS
echo
$PYTHON pipelines/02_run_metaregression.py --input-suffix "$SUF" --label "$MODE"
echo; echo "finished: $(date -Is)"
