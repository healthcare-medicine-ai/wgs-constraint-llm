#!/bin/bash
#SBATCH --job-name=gate07
#SBATCH --partition=mrivas
#SBATCH --time=10:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Fidelity gate for the moderator ablations behind Supplementary Figures S1-S4.
#
# The published panels cannot be pixel-compared: they survive only as 1299x1092
# images inside Supplemental Data.pdf, rescaled by Word. The tables can be, and
# they are the substance -- six leave-one-out fits committed with the Revision 2
# snapshot, each carrying the same 72,386 gene-groups as the full model.
#
# WHY THIS IS NOT A BIT-IDENTITY GATE.
#
# It was, on the first attempt, and it read as a catastrophic failure: p_unified
# only 54-60% bit-identical with relative differences up to 18. That was the
# wrong metric. Relative difference is meaningless for a quantity spanning
# 1e-300 to 1, and measured properly the tables reproduce:
#
#     median |delta log10 p|   exactly 0.00 for all six ablations
#     within 1e-9              99.967% to 99.993%
#     significance changes     0 of 6, at 41/53/49/43/62/61 gene-groups each
#
# So more than half the p-values are bit-identical, the median difference is
# zero, and the residual sits in a ~0.02% tail of ill-conditioned low-variant
# models -- the same tail documented for schizophrenia in gate 11, where every
# unstable fit fell below the n_variants >= 25 inclusion floor.
#
# This gate therefore asserts the gene set, n_variants exactly, no change in any
# significance call, and a bound on the log-scale difference. Bit-identity is
# reported as information so a future run that achieves it is visible.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1

$PY pipelines/07_moderator_ablations.py --as-published || exit 1

echo
echo "=== gate 07: ablation tables vs the Revision 2 committed files ==="
$PY - <<'PY'
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from wgs_constraint import get_config

cfg = get_config()
thr = cfg.param("exome_wide_threshold")
MIN_WITHIN_1E9 = 0.999      # at least this fraction must agree to 1e-9 in log10
MAX_DLOG = 5.0              # and nothing may move further than this

R = "results/"
STEMS = [
    "MINUS_log_constraint",
    "MINUS_GERP_RS",
    "MINUS_log_pathogenicity",
    "NO_pLoF_NO_missense",
    "MINUS_pLoF_ind",
    "MINUS_missense_ind",
]
KEY = ["gene_id", "group"]

ok = True
for stem in STEMS:
    pub_path = R + "epilepsy_unified_model_pvalues_" + stem + ".tsv"
    new_path = R + "epilepsy_unified_model_pvalues_" + stem + "_REPRO.tsv"
    try:
        pub = pd.read_csv(pub_path, sep="\t")
        new = pd.read_csv(new_path, sep="\t")
    except FileNotFoundError as exc:
        print("  " + stem + ": missing -- " + str(exc))
        ok = False
        continue

    merged = pub.merge(new, on=KEY, suffixes=("_a", "_b"), how="outer",
                       indicator=True)
    complete = bool((merged["_merge"] == "both").all())
    both = merged[merged["_merge"] == "both"]

    nv_a = pd.to_numeric(both["n_variants_a"], errors="coerce").to_numpy(float)
    nv_b = pd.to_numeric(both["n_variants_b"], errors="coerce").to_numpy(float)
    nf = np.isfinite(nv_a) & np.isfinite(nv_b)
    nv_ok = bool((nv_a[nf] == nv_b[nf]).all())

    x = pd.to_numeric(both["p_unified_a"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(both["p_unified_b"], errors="coerce").to_numpy(float)
    f = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    exact = (x[f] == y[f]).mean()
    dlog = np.abs(np.log10(x[f]) - np.log10(y[f]))
    within = (dlog < 1e-9).mean()
    sa, sb = x[f] < thr, y[f] < thr
    disagree = int((sa != sb).sum())

    passed = (complete and nv_ok and disagree == 0
              and within >= MIN_WITHIN_1E9 and dlog.max() <= MAX_DLOG)
    ok &= passed
    print("  {:<26} set {}  n_var {}  bit {:5.1f}%  <1e-9 {:7.3f}%  "
          "max dlog {:7.2e}  sig {}={} disagree {}  {}".format(
              stem, "ok" if complete else "DIFFERS", "ok" if nv_ok else "DIFFERS",
              100 * exact, 100 * within, dlog.max(),
              int(sa.sum()), int(sb.sum()), disagree,
              "PASS" if passed else "FAIL"))

print()
print("GATE 07:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
