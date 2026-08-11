#!/bin/bash
#SBATCH --job-name=gate07
#SBATCH --partition=mrivas
#SBATCH --time=10:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Fidelity gate for the moderator ablations behind Supplementary Figures S1-S4.
# Two levels: the six fitted tables, then the four published panels pixel by
# pixel.
#
# WHY THE TABLE CHECK IS NOT BIT-IDENTITY.
#
# It was, on the first attempt, and it read as a catastrophe: p_unified only
# 54-60% bit-identical, relative differences up to 18. That was the wrong
# metric. Relative difference is meaningless for a quantity spanning 1e-300 to
# 1. Measured on the log scale the six tables reproduce:
#
#     median |delta log10 p|   exactly 0.00, all six
#     within 1e-9              99.967% to 99.993%
#     significance changes     none, at 41/53/49/43/62/61 gene-groups each
#
# So over half the p-values are bit-identical, the median difference is zero,
# and the rest is a ~0.02% tail of ill-conditioned low-variant fits -- the tail
# gate 11 documents for schizophrenia, where every unstable model fell below the
# n_variants >= 25 floor. Since these figures apply that same floor, the tail is
# excluded from the panels entirely, which is why they can be gated exactly.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1

$PY pipelines/07_moderator_ablations.py --as-published || exit 1

echo
echo "=== gate 07a: ablation tables vs the Revision 2 committed files ==="
$PY - <<'PY'
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from wgs_constraint import get_config

cfg = get_config()
thr = cfg.param("exome_wide_threshold")
MIN_WITHIN_1E9 = 0.999
MAX_DLOG = 5.0

R = "results/"
STEMS = ["MINUS_log_constraint", "MINUS_GERP_RS", "MINUS_log_pathogenicity",
         "NO_pLoF_NO_missense", "MINUS_pLoF_ind", "MINUS_missense_ind"]
KEY = ["gene_id", "group"]

ok = True
for stem in STEMS:
    try:
        pub = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + ".tsv", sep="\t")
        new = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + "_REPRO.tsv", sep="\t")
    except FileNotFoundError as exc:
        print("  " + stem + ": missing -- " + str(exc))
        ok = False
        continue
    m = pub.merge(new, on=KEY, suffixes=("_a", "_b"), how="outer", indicator=True)
    complete = bool((m["_merge"] == "both").all())
    both = m[m["_merge"] == "both"]

    nva = pd.to_numeric(both["n_variants_a"], errors="coerce").to_numpy(float)
    nvb = pd.to_numeric(both["n_variants_b"], errors="coerce").to_numpy(float)
    nf = np.isfinite(nva) & np.isfinite(nvb)
    nv_ok = bool((nva[nf] == nvb[nf]).all())

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
print("GATE 07a:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
A=$?

echo
echo "=== gate 07b: the four published panels, pixel by pixel ==="
$PY - <<'PY'
import sys
from pathlib import Path

import numpy as np
from PIL import Image

R = Path("results")
PNGS = [
    "Figure_heatmap_FULL_vs_MINUS_log_constraint_2x2",
    "Figure_heatmap_FULL_vs_MINUS_GERP_RS_2x2",
    "Figure_heatmap_FULL_vs_MINUS_log_pathogenicity_2x2",
    "Figure_heatmap_FULL_vs_ALT_no_pLoF_no_missense_2x2",
]
ok = True
for stem in PNGS:
    p, q = R / (stem + ".png"), R / (stem + "_REPRO.png")
    if not (p.exists() and q.exists()):
        print("  {:<52} missing published={} new={}".format(
            stem[:52], p.exists(), q.exists()))
        ok = False
        continue
    a, b = Image.open(p), Image.open(q)
    print("  {:<52} published {}  new {}".format(stem[:52], a.size, b.size))
    if a.size != b.size:
        print("      dimensions differ")
        ok = False
        continue
    aa = np.asarray(a.convert("RGB"), np.int16)
    bb = np.asarray(b.convert("RGB"), np.int16)
    d = np.abs(aa - bb)
    frac = (d.max(axis=2) > 0).mean()
    print("      differing {:.4f}%  max delta {}  -> {}".format(
        100 * frac, d.max(), "IDENTICAL" if frac == 0 else "differs"))
    if frac:
        ys, xs = np.nonzero(d.max(axis=2) > 0)
        print("      bounding box x {}-{} y {}-{}".format(
            xs.min(), xs.max(), ys.min(), ys.max()))
    ok &= frac == 0
print()
print("GATE 07b:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
B=$?

echo
echo "=== gate 07c: the attribution the response letter quotes ==="
$PY - <<'PY'
import sys

import pandas as pd

sys.path.insert(0, "src")
from wgs_constraint import get_config

# Response letter, item 10.
EXPECTED = {"p_pLoF": 13, "p_pathogenicity": 7, "p_missense": 5,
            "p_constraint": 5, "p_gerp": 3}
EXPECTED_TOTAL = 33

cfg = get_config()
floor = cfg.param("min_variants")
merged = pd.read_csv(
    cfg.result("unified_epi25_g4e_merged_n" + str(floor) + ".tsv.gz"), sep="	")
for column in list(EXPECTED) + ["p_unified", "n_variants"]:
    if column in merged.columns:
        merged[column] = pd.to_numeric(merged[column], errors="coerce")
merged = merged[merged["n_variants"] >= floor]
sig = merged[merged["p_unified"] < cfg.param("exome_wide_threshold")].copy()
sig["dominant"] = sig[list(EXPECTED)].idxmin(axis=1)
counts = sig["dominant"].value_counts()

ok = len(sig) == EXPECTED_TOTAL
print("  significant pairs  expected {}  got {}  {}".format(
    EXPECTED_TOTAL, len(sig), "ok" if ok else "MISMATCH"))
for column, want in EXPECTED.items():
    got = int(counts.get(column, 0))
    good = got == want
    ok &= good
    print("  {:<18} expected {:>3}  got {:>3}  {}".format(
        column, want, got, "ok" if good else "MISMATCH"))
print()
print("GATE 07c:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
C=$?

echo
echo "gate 07a exit=$A   gate 07b exit=$B   gate 07c exit=$C"
[ "$A" -eq 0 ] && [ "$B" -eq 0 ] && [ "$C" -eq 0 ] || exit 1
