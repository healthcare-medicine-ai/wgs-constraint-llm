#!/bin/bash
#SBATCH --job-name=gate11
#SBATCH --partition=mrivas
#SBATCH --time=06:00:00
#SBATCH --mem=192G
# Pin the CPU vendor. numpy dispatches a different SIMD kernel for np.log on
# AVX-512 (Intel Skylake-X) than on AVX2 (AMD Rome), and the two differ in the
# last bit. effect_size then differs on ~18% of rows at up to 3.3e-15, so a
# bit-identity gate passes or fails depending purely on which node Slurm picked.
# The reference artifacts were built on Intel.
#SBATCH --constraint=CPU_MNF:INTEL
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Fidelity gate for the schizophrenia port.
#
# WHY THIS GATE IS NOT BIT-IDENTICAL, UNLIKE THE EPILEPSY ONE.
#
# Everything structural reproduces exactly: the same 5,684,992 rows, the same
# columns, the same 190,328 models, and n_variants bit-identical throughout.
# One column does not: effect_size agrees on ~80% of rows, with a maximum
# relative difference of 2e-15 and a median of 1.8e-16 -- the last bit or two
# of a double.
#
# The cause was investigated and is NOT known. Ruled out by measurement:
#   * frame construction. Six variants -- chained !=, the compound mask from
#     cell 15, .isin(), an added .copy(), np.asarray, and reordered filters --
#     all give byte-for-byte the same result (genentech_fix/bisect_effect_size.py).
#   * the np.log SIMD path. Vectorised, one-element-at-a-time and libm math.log
#     agree with each other (genentech_fix/diag_log_path.py), so it is not CPU
#     dispatch across Sherlock's heterogeneous nodes.
# The differing rows are exactly those where the arithmetic can round: ac_case=1
# (factor 1.5) dominates them, ac_case=0 (factor 0.5, exact in binary) dominates
# the matching ones.
#
# So this gate asserts what is both achievable and scientifically meaningful:
# an identical row set, identical variant counts, and NO CHANGE IN ANY
# SIGNIFICANCE CALL. Measured on the September 2025 cache: 824 significant
# gene-groups before and after, zero disagreements. Bit-identity is still
# reported, as information, so a future run that achieves it is visible rather
# than silently accepted.
#
# One number in the 11b output deserves reading rather than skipping: 109
# gene-groups move by more than an order of magnitude in p. Every one of them
# has fewer than 21 variants (median 8), and none passes the n_variants >= 25
# inclusion floor, so none is ever reported
# (genentech_fix/diag_scz_unstable_models.py). That is structural, not luck:
# five moderators plus an intercept is six parameters, so a gene with 8 variants
# is fitted with 2 residual degrees of freedom and the fit is near-saturated. A
# p-value from such a model is not meaningful at the precision it is printed at,
# whatever the input. The inclusion floor removes exactly these.
#
# Do not weaken this further. If a future change makes significance calls
# disagree, or if a model above the floor starts moving, that is a real
# regression, not more floating-point noise.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1

MAX_REL_EFFECT_SIZE=1e-12   # effect_size must agree to at least this
$PY pipelines/11_schizophrenia.py --no-fix || exit 1

echo
echo "=== gate 11a: merged regression input vs the cached September 2025 file ==="
MAX_REL="$MAX_REL_EFFECT_SIZE" $PY - <<'PY'
import os, sys
import numpy as np, pandas as pd
tol = float(os.environ["MAX_REL"])
R = "results/"
key = ["chr", "pos", "ref", "alt", "gene_id", "group"]
a = pd.read_csv(R + "constraint_gerp_am_scz_variants.tsv.gz", sep="\t")
b = pd.read_csv(R + "constraint_gerp_am_scz_variants_REPRO.tsv.gz", sep="\t")
print(f"  cached      : {len(a):,} rows")
print(f"  regenerated : {len(b):,} rows")
rows_ok = len(a) == len(b)
cols_ok = list(a.columns) == list(b.columns)
print(f"  row counts match: {rows_ok}")
print(f"  columns match   : {cols_ok}")
if not rows_ok:
    m = b.merge(a[key].assign(_c=1), on=key, how="left")
    print(f"  rows absent from the cache: {int(m['_c'].isna().sum()):,}")
ok = rows_ok and cols_ok
if ok:
    for c in a.columns:
        x, y = a[c], b[c]
        if np.issubdtype(x.dtype, np.number) and np.issubdtype(y.dtype, np.number):
            same = ((x == y) | (x.isna() & y.isna())).mean()
            if same == 1.0:
                continue
            xv, yv = x.to_numpy(float), y.to_numpy(float)
            f = np.isfinite(xv) & np.isfinite(yv) & (xv != 0)
            rel = np.abs(xv[f] - yv[f]) / np.abs(xv[f])
            within = rel.max() <= tol
            print(f"    {c:<18} bit-identical {100*same:6.2f}%   max rel {rel.max():.2e}"
                  f"   {'within tolerance' if within else 'EXCEEDS TOLERANCE'}")
            ok &= bool(within)
        elif (x.astype(str) != y.astype(str)).any():
            print(f"    {c:<18} differs (non-numeric)")
            ok = False
print("\nGATE 11a:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
A=$?

echo
echo "=== gate 11b: fitted p-values vs the published SCZ results ==="
$PY - <<'PY'
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from wgs_constraint import get_config
thr = get_config().param("exome_wide_threshold")
R = "results/"
a = pd.read_csv(R + "schizophrenia_variants_unified_gerp_model_pvalues.tsv", sep="\t")
b = pd.read_csv(R + "schizophrenia_variants_unified_gerp_model_pvalues_REPRO.tsv", sep="\t")
print(f"  published rows : {len(a):,}")
print(f"  regenerated    : {len(b):,}")
m = a.merge(b, on=["gene_id", "group"], suffixes=("_a", "_b"), how="outer", indicator=True)
complete = bool((m["_merge"] == "both").all())
print(f"  overlap        : {(m['_merge']=='both').sum():,}   "
      f"only published {(m['_merge']=='left_only').sum():,}   "
      f"only new {(m['_merge']=='right_only').sum():,}")
both = m[m["_merge"] == "both"]

nv_a = pd.to_numeric(both["n_variants_a"], errors="coerce").to_numpy(float)
nv_b = pd.to_numeric(both["n_variants_b"], errors="coerce").to_numpy(float)
f = np.isfinite(nv_a) & np.isfinite(nv_b)
nv_ok = bool((nv_a[f] == nv_b[f]).all())
print(f"  n_variants bit-identical: {100*(nv_a[f]==nv_b[f]).mean():.2f}%   -> "
      f"{'ok' if nv_ok else 'REGRESSION'}")

x = pd.to_numeric(both["p_unified_a"], errors="coerce").to_numpy(float)
y = pd.to_numeric(both["p_unified_b"], errors="coerce").to_numpy(float)
g = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
dlog = np.abs(np.log10(x[g]) - np.log10(y[g]))
print(f"  |delta log10 p|  median {np.median(dlog):.2e}   max {dlog.max():.2e}   "
      f"within 1e-6 {100*(dlog<1e-6).mean():.3f}%")
print(f"  moving > 1 log unit: {int((dlog>1).sum()):,}")
for c in ["p_constraint", "p_gerp", "p_unified"]:
    u = pd.to_numeric(both[c + "_a"], errors="coerce").to_numpy(float)
    v = pd.to_numeric(both[c + "_b"], errors="coerce").to_numpy(float)
    h = np.isfinite(u) & np.isfinite(v)
    print(f"  {c:<14} bit-identical {100*(u[h]==v[h]).mean():6.2f}%  (informational)")

sa, sb = x[g] < thr, y[g] < thr
disagree = int((sa != sb).sum())
print(f"\n  significant at p < {thr:g}:  published {int(sa.sum())}   "
      f"regenerated {int(sb.sum())}   DISAGREEMENTS {disagree}")

ok = complete and nv_ok and disagree == 0
print("\nGATE 11b:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
B=$?

echo
echo "gate 11a exit=$A   gate 11b exit=$B"
[ "$A" -eq 0 ] && [ "$B" -eq 0 ] || exit 1
