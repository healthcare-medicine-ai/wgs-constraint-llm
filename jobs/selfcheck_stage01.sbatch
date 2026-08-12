#!/bin/bash
#SBATCH --job-name=selfchk
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Control for the VIAFIG comparison.
#
# Rebuilding through the published constraint+GERP table left effect_size only
# 81.8% bit-identical to results/constraint_gerp_am_epi25_variants_FIXED.tsv.gz.
# effect_size is computed from allele counts before GERP is touched, so the
# route cannot explain it. Either stage 01 is not deterministic run to run, or
# the committed _FIXED file predates some change.
#
# This reruns stage 01 by the ORIGINAL bigWig route, changing nothing, and
# compares. If effect_size differs here too, the route is exonerated and the
# nondeterminism is the real finding -- the same one gate 11 could not explain
# for schizophrenia.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/01_build_regression_input.py --out-suffix _SELFCHK || exit 1
echo
echo "=== bigWig route rerun vs the committed _FIXED file ==="
$PY - <<'PY'
import numpy as np
import pandas as pd
R = "results/"
a = pd.read_csv(R + "constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t")
b = pd.read_csv(R + "constraint_gerp_am_epi25_variants_SELFCHK.tsv.gz", sep="\t")
print(f"  rows {len(a):,} vs {len(b):,}")
for c in ("prob_0", "GERP_RS", "effect_size", "var_effect_size", "am_pathogenicity"):
    x, y = a[c], b[c]
    same = ((x == y) | (x.isna() & y.isna())).mean()
    xv, yv = x.to_numpy(float), y.to_numpy(float)
    f = np.isfinite(xv) & np.isfinite(yv) & (xv != 0)
    rel = np.abs(xv[f] - yv[f]) / np.abs(xv[f])
    print(f"  {c:<18} bit-identical {100*same:6.2f}%   max rel {rel.max():.3e}")
PY
