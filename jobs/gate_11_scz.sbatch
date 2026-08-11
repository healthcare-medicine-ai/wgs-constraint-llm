#!/bin/bash
#SBATCH --job-name=gate11
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# First execution of the extracted schizophrenia pipeline, in reproduce mode.
# Must match the published SCZ outputs from October 2025.
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/11_schizophrenia.py --no-fix
echo
echo "=== gate: compare against the published SCZ results ==="
$PY - <<'PY'
import numpy as np, pandas as pd
R = "results/"
a = pd.read_csv(R + "schizophrenia_variants_unified_gerp_model_pvalues.tsv", sep="\t")
b = pd.read_csv(R + "schizophrenia_variants_unified_gerp_model_pvalues_REPRO.tsv", sep="\t")
print(f"  published rows : {len(a):,}")
print(f"  regenerated    : {len(b):,}")
m = a.merge(b, on=["gene_id", "group"], suffixes=("_a", "_b"), how="outer", indicator=True)
print(f"  overlap        : {(m['_merge']=='both').sum():,}   "
      f"only published {(m['_merge']=='left_only').sum():,}   "
      f"only new {(m['_merge']=='right_only').sum():,}")
both = m[m["_merge"] == "both"]
ok = True
for c in ["n_variants", "p_constraint", "p_gerp", "p_pathogenicity",
          "p_pLoF", "p_missense", "p_unified"]:
    x = pd.to_numeric(both[c + "_a"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(both[c + "_b"], errors="coerce").to_numpy(float)
    f = np.isfinite(x) & np.isfinite(y)
    same = np.isclose(x[f], y[f], rtol=1e-9, atol=0).mean() if f.sum() else 0
    print(f"  {c:<18} identical {100*same:6.2f}%")
    ok &= same > 0.9999
print("\nGATE 11:", "PASS" if ok and (m['_merge']=='both').all() else "REVIEW")
PY
