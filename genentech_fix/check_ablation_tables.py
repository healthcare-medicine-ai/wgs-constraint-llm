"""Compare the regenerated ablation tables against the Revision 2 committed ones."""
import sys

import numpy as np
import pandas as pd

R = "results/"
STEMS = ["MINUS_log_constraint", "MINUS_GERP_RS", "MINUS_log_pathogenicity",
         "NO_pLoF_NO_missense", "MINUS_pLoF_ind", "MINUS_missense_ind"]
KEY = ["gene_id", "group"]

ok = True
for stem in STEMS:
    pub = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + ".tsv", sep="\t")
    new = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + "_REPRO.tsv", sep="\t")
    m = pub.merge(new, on=KEY, suffixes=("_a", "_b"), how="outer", indicator=True)
    complete = bool((m["_merge"] == "both").all())
    both = m[m["_merge"] == "both"]
    line = "  {:<26} pub {:>7,} new {:>7,}".format(stem, len(pub), len(new))
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
            rel = np.abs(x[nz] - y[nz]) / np.abs(x[nz])
            line += "(maxrel {:.1e})".format(rel.max())
            ok = False
    print(line, flush=True)
print()
print("GATE 07 (tables):", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
