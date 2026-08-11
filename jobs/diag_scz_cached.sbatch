#!/bin/bash
#SBATCH --job-name=sczcached
#SBATCH --partition=mrivas
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY - <<'PY'
import sys, numpy as np, pandas as pd
sys.path.insert(0, "src")
from wgs_constraint import get_config, fit_per_gene, prepare_regression_input
cfg = get_config()
cached = cfg.result("constraint_gerp_am_scz_variants.tsv.gz")   # 2025-09-15
print(f"using cached input: {cached}", flush=True)
df = pd.read_csv(cached, sep="\t")
print(f"  rows: {len(df):,}", flush=True)
out = fit_per_gene(prepare_regression_input(df, clip_epsilon=cfg.param("clip_epsilon")))
print(f"  models fitted: {len(out):,}", flush=True)
pub = pd.read_csv(cfg.result("schizophrenia_variants_unified_gerp_model_pvalues.tsv"), sep="\t")
print(f"  published rows: {len(pub):,}")
m = pub.merge(out, on=["gene_id","group"], suffixes=("_a","_b"))
print(f"  overlap: {len(m):,}")
for c in ["n_variants","p_constraint","p_gerp","p_unified"]:
    x = pd.to_numeric(m[c+"_a"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(m[c+"_b"], errors="coerce").to_numpy(float)
    f = np.isfinite(x) & np.isfinite(y)
    print(f"  {c:<14} bit-identical {100*(x[f]==y[f]).mean():6.2f}%")
PY
