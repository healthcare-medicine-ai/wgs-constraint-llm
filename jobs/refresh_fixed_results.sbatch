#!/bin/bash
#SBATCH --job-name=refixed
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# The committed constraint_gerp_am_epi25_variants_FIXED.tsv.gz was written
# 2026-08-10 02:35, before commit a15fc7a reverted epi25.py from .isin() back to
# the notebook's chained != at 23:22 the same day. So it carries the ULP
# difference that revert removed, and every corrected epilepsy number published
# so far derives from it.
#
# _SELFCHK is the same pipeline under current code and is bit-identical to an
# independent rerun. This refits from it and checks whether any gene changes
# significance. Expected: none. Verified, not assumed.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/02_run_metaregression.py --input-suffix _SELFCHK --label selfchk || exit 1
echo
echo "=== corrected results: stale _FIXED input vs current-code _SELFCHK ==="
$PY - <<'PY'
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from wgs_constraint import get_config
cfg = get_config(); thr = cfg.param("exome_wide_threshold")
R = "results/"
a = pd.read_csv(R + "epilepsy_unified_model_pvalues_FIXED.tsv", sep="\t")
b = pd.read_csv(R + "epilepsy_unified_model_pvalues_SELFCHK.tsv", sep="\t")
m = a.merge(b, on=["gene_id", "group"], suffixes=("_a", "_b"), how="outer", indicator=True)
print(f"  models {len(a):,} vs {len(b):,}   overlap {(m['_merge']=='both').sum():,}")
both = m[m["_merge"] == "both"]
x = pd.to_numeric(both["p_unified_a"], errors="coerce").to_numpy(float)
y = pd.to_numeric(both["p_unified_b"], errors="coerce").to_numpy(float)
n = pd.to_numeric(both["n_variants_a"], errors="coerce").to_numpy(float)
f = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
d = np.abs(np.log10(x[f]) - np.log10(y[f]))
print(f"  n_variants bit-identical: {100*(pd.to_numeric(both['n_variants_a'],errors='coerce').to_numpy(float)==pd.to_numeric(both['n_variants_b'],errors='coerce').to_numpy(float)).mean():.2f}%")
print(f"  p_unified  bit-identical: {100*(x[f]==y[f]).mean():.2f}%   median |dlog| {np.median(d):.2e}   max {d.max():.2e}")
floor = cfg.param("min_variants")
for label, mask in (("all models", f), ("n>=%d" % floor, f & (n >= floor))):
    sa, sb = x[mask] < thr, y[mask] < thr
    print(f"  {label:<12} significant {int(sa.sum())} vs {int(sb.sum())}   "
          f"DISAGREEMENTS {int((sa!=sb).sum())}")
PY
