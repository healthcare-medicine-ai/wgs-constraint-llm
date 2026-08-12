#!/bin/bash
#SBATCH --job-name=gate01f
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# gate_01_refactor gates the --no-fix reproduce path of stage 1. Nothing gated
# the CORRECTED path, and that is exactly where a stale artifact survived the
# whole overhaul: constraint_gerp_am_epi25_variants_FIXED.tsv.gz was written
# 2026-08-10 02:35, twenty-one hours before commit a15fc7a reverted epi25.py
# from .isin() to the notebook's chained !=. Every corrected epilepsy number
# published up to 2026-08-11 was fitted from it. The perturbation was 3e-15 and
# changed no significance call, but nothing was checking.
#
# This rebuilds the corrected input and demands it match the committed one.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/01_build_regression_input.py --out-suffix _FIXEDCHK || exit 1
echo
A=$(zcat results/constraint_gerp_am_epi25_variants_FIXED.tsv.gz    | md5sum | cut -d' ' -f1)
B=$(zcat results/constraint_gerp_am_epi25_variants_FIXEDCHK.tsv.gz | md5sum | cut -d' ' -f1)
echo "  committed (_FIXED)   : $A"
echo "  rebuilt   (_FIXEDCHK): $B"
if [ "$A" = "$B" ]; then echo; echo "GATE 01-FIXED: PASS"; exit 0; fi
echo; echo "GATE 01-FIXED: FAIL -- the committed corrected input is not what the code produces"
$PY - <<'PY'
import numpy as np, pandas as pd
R = "results/"
a = pd.read_csv(R + "constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t")
b = pd.read_csv(R + "constraint_gerp_am_epi25_variants_FIXEDCHK.tsv.gz", sep="\t")
print(f"  rows {len(a):,} vs {len(b):,}")
for c in a.columns:
    if not (np.issubdtype(a[c].dtype, np.number) and np.issubdtype(b[c].dtype, np.number)):
        continue
    same = ((a[c] == b[c]) | (a[c].isna() & b[c].isna())).mean()
    if same < 1.0:
        x, y = a[c].to_numpy(float), b[c].to_numpy(float)
        f = np.isfinite(x) & np.isfinite(y) & (x != 0)
        print(f"    {c:<18} bit-identical {100*same:6.2f}%  max rel "
              f"{(np.abs(x[f]-y[f])/np.abs(x[f])).max():.3e}")
PY
exit 1
