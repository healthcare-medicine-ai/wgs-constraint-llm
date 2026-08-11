#!/bin/bash
#SBATCH --job-name=scz_fix
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# The corrected schizophrenia analysis: GERP indexed 1-based, AlphaMissense
# collapsed per variant, and the full ten-consequence filter that the submitted
# analysis actually used. Run only after gate_11 confirms the reproduce path
# still matches the September 2025 cache -- otherwise the corrected numbers
# conflate the corrections with a broken input, which is exactly what happened
# on the first attempt.
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/11_schizophrenia.py
echo
echo "=== corrected vs published schizophrenia results ==="
$PY - <<'PY'
import pandas as pd
R = "results/"
thr = 3.4e-7
pub = pd.read_csv(R + "table_A2_schizophrenia_full_p_value_comparison_n25.csv") \
    if __import__("os").path.exists(R + "table_A2_schizophrenia_full_p_value_comparison_n25.csv") else None
fix = pd.read_csv(R + "table_A2_schizophrenia_full_p_value_comparison_n25_FIXED.csv")
print(f"  corrected gene-groups tested : {len(fix):,}")
sig = fix[fix["Unified Model p-value"] < thr]
print(f"  corrected significant        : {len(sig)}")
for _, r in sig.sort_values("Unified Model p-value").iterrows():
    print(f"    {str(r['Gene Name']):<12} unified={r['Unified Model p-value']:.3e}"
          f"  SCHEMA={r['SCHEMA p-value']:.3e}")
if pub is not None:
    ps = pub[pub["Unified Model p-value"] < thr]
    print(f"\n  published significant        : {len(ps)}")
    for _, r in ps.sort_values("Unified Model p-value").iterrows():
        print(f"    {str(r['Gene Name']):<12} unified={r['Unified Model p-value']:.3e}")
PY
