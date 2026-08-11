"""Does the response letter's dominant-moderator attribution reproduce?

Letter (item 10): "Across the 33 gene-group pairs reaching exome-wide
significance under the unified model, the dominant moderator is the predicted
loss-of-function term for 13, predicted pathogenicity for 7, the missense
annotation for 5, HMM-derived constraint for 5, and GERP RS for 3."
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
SINGLES = {
    "MINUS_pLoF_ind": "pLoF",
    "MINUS_log_pathogenicity": "pathogenicity",
    "MINUS_missense_ind": "missense",
    "MINUS_log_constraint": "HMM constraint",
    "MINUS_GERP_RS": "GERP RS",
}

full = pd.read_csv(R + "epilepsy_unified_model_pvalues.tsv", sep="\t")
full["p_unified"] = pd.to_numeric(full["p_unified"], errors="coerce")
full["n_variants"] = pd.to_numeric(full["n_variants"], errors="coerce")

for label, sel in (("no floor", full),
                   (f"n_variants >= {floor}", full[full["n_variants"] >= floor])):
    sig = sel[sel["p_unified"] < thr].copy()
    for stem in SINGLES:
        abl = pd.read_csv(R + "epilepsy_unified_model_pvalues_" + stem + ".tsv",
                          sep="\t")
        abl["p_unified"] = pd.to_numeric(abl["p_unified"], errors="coerce")
        j = sig.merge(abl[["gene_id", "group", "p_unified"]],
                      on=["gene_id", "group"], suffixes=("", "_abl"))
        sig[stem] = (np.log10(j["p_unified_abl"].to_numpy(dtype=float))
                     - np.log10(j["p_unified"].to_numpy(dtype=float)))
    sig["dominant"] = sig[list(SINGLES)].idxmax(axis=1)
    counts = sig["dominant"].value_counts()
    print(f"{label}: {len(sig)} exome-wide significant gene-groups")
    for stem, name in SINGLES.items():
        print(f"    {name:<16} {int(counts.get(stem, 0)):>3}")
    print()

print("letter says: 33 pairs -- pLoF 13, pathogenicity 7, missense 5, "
      "HMM constraint 5, GERP RS 3")
