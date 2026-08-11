#!/bin/bash
#SBATCH --job-name=gate01def
#SBATCH --partition=mrivas
#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
# A failing gate exits nonzero. The first version of these scripts printed
# the verdict and exited 0, so a --dependency=afterok job downstream launched
# against a failing gate. A gate that cannot stop anything is not a gate.
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/02_run_metaregression.py --input-suffix _REGATE --label regate
echo
echo "=== do the p-values still match the SUBMITTED results exactly? ==="
$PY - <<'PY'
import sys
import numpy as np, pandas as pd
sub = pd.read_csv("genentech_fix/baseline_R2/results/epilepsy_unified_model_pvalues.tsv", sep="\t")
new = pd.read_csv("results/epilepsy_unified_model_pvalues_REGATE.tsv", sep="\t")
m = sub.merge(new, on=["gene_id","group"], suffixes=("_a","_b"))
print(f"  gene-groups compared: {len(m):,}")
ok = True
for c in ["n_variants","p_constraint","p_gerp","p_pathogenicity","p_pLoF","p_missense","p_unified"]:
    x = pd.to_numeric(m[c+"_a"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(m[c+"_b"], errors="coerce").to_numpy(float)
    f = np.isfinite(x) & np.isfinite(y)
    exact = (x[f] == y[f]).mean()
    close = np.isclose(x[f], y[f], rtol=1e-12).mean()
    print(f"  {c:<18} bit-identical {100*exact:6.2f}%   within 1e-12 {100*close:6.2f}%")
    ok &= close > 0.99999
print("\nVERDICT:", "p-values unaffected by the ULP difference" if ok else "REVIEW")
sys.exit(0 if ok else 1)
PY
