#!/bin/bash
#SBATCH --job-name=gerp-smoke
#SBATCH --partition=mrivas
#SBATCH --time=01:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/smoke_%j.out
#SBATCH --error=logs/smoke_%j.err

set -euo pipefail

WORKDIR=/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/genentech_fix
PYTHON=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python

cd "$WORKDIR"
export PYTHONUNBUFFERED=1

echo "host    : $(hostname)"
echo "started : $(date -Is)"
echo

# Smoke test on two small chromosomes only. Proves the port runs end to end
# and that the +/-1 shift actually changes GERP values, before committing to
# the full genome-wide run.
echo "########## SMOKE: no-fix (reproduces submitted behaviour) ##########"
"$PYTHON" 01_build_regression_input.py --no-fix --out-suffix _SMOKE_NOFIX --chroms chr21,chr22

echo
echo "########## SMOKE: fixed ##########"
"$PYTHON" 01_build_regression_input.py --out-suffix _SMOKE_FIXED --chroms chr21,chr22

echo
echo "########## compare GERP columns ##########"
"$PYTHON" - <<'PY'
import pandas as pd
R = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
     "osthoag/wgs-constraint-llm/results/")
key = ["chr", "pos", "ref", "alt", "gene_id", "group"]

a = pd.read_csv(R + "constraint_gerp_am_epi25_variants_SMOKE_NOFIX.tsv.gz", sep="\t")
b = pd.read_csv(R + "constraint_gerp_am_epi25_variants_SMOKE_FIXED.tsv.gz", sep="\t")

print(f"no-fix rows : {len(a):,}")
print(f"fixed  rows : {len(b):,}")
print(f"row delta   : {len(b) - len(a):,}  (expected negative: AM dedupe)")

m = a[key + ["GERP_RS"]].drop_duplicates(key).merge(
    b[key + ["GERP_RS"]].drop_duplicates(key), on=key,
    suffixes=("_nofix", "_fixed"))
diff = (m["GERP_RS_nofix"] != m["GERP_RS_fixed"]).mean()
corr = m["GERP_RS_nofix"].corr(m["GERP_RS_fixed"])
print(f"\ncommon variants compared : {len(m):,}")
print(f"GERP values that changed : {100 * diff:.2f}%")
print(f"correlation old vs new   : {corr:.4f}")
print("\n(A high correlation with most values changed is the expected")
print(" signature of a 1 bp shift in an autocorrelated track.)")
PY

echo
echo "finished: $(date -Is)"
