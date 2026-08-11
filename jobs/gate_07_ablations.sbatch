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
# This runs stage 07 from the pre-correction input and demands bit-identity
# against all six.
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

    line = "  {:<26} published {:>7,}  new {:>7,}".format(
        stem, len(pub), len(new))
    if not complete:
        line += "  GENE SET DIFFERS"
        ok = False
    for col in ("n_variants", "p_unified"):
        x = pd.to_numeric(both[col + "_a"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(both[col + "_b"], errors="coerce").to_numpy(float)
        f = np.isfinite(x) & np.isfinite(y)
        exact = (x[f] == y[f]).mean() if f.sum() else 0.0
        line += "  {} {:6.2f}%".format(col, 100 * exact)
        if exact < 1.0:
            nz = f & (x != 0)
            rel = np.abs(x[nz] - y[nz]) / np.abs(x[nz]) if nz.sum() else np.array([0.0])
            line += "(maxrel {:.1e})".format(rel.max())
            ok = False
    print(line)

print()
print("GATE 07:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
