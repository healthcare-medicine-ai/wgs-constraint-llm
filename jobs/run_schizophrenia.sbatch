#!/bin/bash
#SBATCH --job-name=scz
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
#   sbatch --job-name=scz-fixed jobs/run_schizophrenia.sbatch
#   sbatch --job-name=scz-repro jobs/run_schizophrenia.sbatch --no-fix
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
export PYTHONUNBUFFERED=1
/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python \
  pipelines/11_schizophrenia.py "$@"
echo
echo "=== corrected vs published schizophrenia results ==="
/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python - <<'PY'
import numpy as np, pandas as pd
from pathlib import Path
R = Path("results")
pub = pd.read_csv(R/"schizophrenia_variants_unified_gerp_model_pvalues.tsv", sep="\t")
new_path = R/"schizophrenia_variants_unified_gerp_model_pvalues_FIXED.tsv"
if not new_path.exists():
    raise SystemExit("no corrected output to compare")
new = pd.read_csv(new_path, sep="\t")
THR = 3.4e-7
def sig(df):
    d = df[df["n_variants"] >= 25]
    g = d.groupby("gene_id", as_index=False).agg({"gene_name":"first","p_unified":"min"})
    return g[g["p_unified"] < THR]
a, b = sig(pub), sig(new)
print(f"  published gene-groups : {len(pub):,}")
print(f"  corrected gene-groups : {len(new):,}")
print(f"  significant, published: {len(a)}  -> {sorted(a['gene_name'].astype(str))}")
print(f"  significant, corrected: {len(b)}  -> {sorted(b['gene_name'].astype(str))}")
lost = set(a["gene_id"]) - set(b["gene_id"]); gained = set(b["gene_id"]) - set(a["gene_id"])
nm = dict(zip(a["gene_id"], a["gene_name"])); nm.update(dict(zip(b["gene_id"], b["gene_name"])))
print(f"  lost   : {[nm[g] for g in lost]}")
print(f"  gained : {[nm[g] for g in gained]}")
PY
