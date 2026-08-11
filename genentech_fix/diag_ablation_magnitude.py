"""Is the ablation p-value difference the benign ULP effect, or something real?

Relative difference is the wrong metric for a p-value spanning 1e-300 to 1.
Measure agreement in -log10(p), and whether any gene-group changes significance.
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from wgs_constraint import get_config

cfg = get_config()
thr = cfg.param("exome_wide_threshold")
floor = cfg.param("min_variants")
R = "results/"
STEMS = ["MINUS_log_constraint", "MINUS_GERP_RS", "MINUS_log_pathogenicity",
         "NO_pLoF_NO_missense", "MINUS_pLoF_ind", "MINUS_missense_ind"]
KEY = ["gene_id", "group"]

print(f"{'ablation':<26} {'med dlog':>10} {'max dlog':>10} {'<1e-9':>8} "
      f"{'<1e-6':>8} {'sig a':>6} {'sig b':>6} {'disag':>6} {'disag n>=25':>12}")
for stem in STEMS:
    a = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + ".tsv", sep="\t")
    b = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + "_REPRO.tsv", sep="\t")
    m = a.merge(b, on=KEY, suffixes=("_a", "_b"))
    x = pd.to_numeric(m["p_unified_a"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(m["p_unified_b"], errors="coerce").to_numpy(float)
    n = pd.to_numeric(m["n_variants_a"], errors="coerce").to_numpy(float)
    f = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    d = np.abs(np.log10(x[f]) - np.log10(y[f]))
    sa, sb = x[f] < thr, y[f] < thr
    dis = sa != sb
    dis_above = dis & (n[f] >= floor)
    print(f"{stem:<26} {np.median(d):10.2e} {d.max():10.2e} "
          f"{100*(d<1e-9).mean():7.3f}% {100*(d<1e-6).mean():7.3f}% "
          f"{int(sa.sum()):6} {int(sb.sum()):6} {int(dis.sum()):6} "
          f"{int(dis_above.sum()):12}")

print()
print("n_variants of the gene-groups that disagree, if any:")
for stem in STEMS:
    a = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + ".tsv", sep="\t")
    b = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + "_REPRO.tsv", sep="\t")
    m = a.merge(b, on=KEY, suffixes=("_a", "_b"))
    x = pd.to_numeric(m["p_unified_a"], errors="coerce")
    y = pd.to_numeric(m["p_unified_b"], errors="coerce")
    d = m[(x < thr) != (y < thr)]
    if len(d):
        print(f"  {stem}: {len(d)} disagree, n_variants "
              f"min {d['n_variants_a'].min()} max {d['n_variants_a'].max()}")
    else:
        print(f"  {stem}: none")
