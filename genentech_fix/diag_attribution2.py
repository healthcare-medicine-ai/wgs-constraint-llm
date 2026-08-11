"""Second reading of "dominant moderator": the most significant coefficient."""
import sys

import pandas as pd

sys.path.insert(0, "src")
from wgs_constraint import get_config

cfg = get_config()
thr = cfg.param("exome_wide_threshold")
floor = cfg.param("min_variants")
R = "results/"
COEF = {
    "p_pLoF": "pLoF",
    "p_pathogenicity": "pathogenicity",
    "p_missense": "missense",
    "p_constraint": "HMM constraint",
    "p_gerp": "GERP RS",
}

sources = [
    ("raw p-value table",
     R + "epilepsy_unified_model_pvalues.tsv", "\t"),
    ("merged epi25+g4e table",
     R + f"unified_epi25_g4e_merged_n{floor}.tsv.gz", "\t"),
]
for label, path, sep in sources:
    try:
        df = pd.read_csv(path, sep=sep)
    except FileNotFoundError:
        print(f"{label}: not present ({path})\n")
        continue
    for c in list(COEF) + ["p_unified", "n_variants"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "n_variants" in df.columns:
        df = df[df["n_variants"] >= floor]
    sig = df[df["p_unified"] < thr].copy()
    have = [c for c in COEF if c in sig.columns]
    if not have:
        print(f"{label}: no coefficient columns\n")
        continue
    sig["dominant"] = sig[have].idxmin(axis=1)
    counts = sig["dominant"].value_counts()
    print(f"{label}: {len(sig)} significant gene-groups (n>={floor})")
    for c in COEF:
        if c in have:
            print(f"    {COEF[c]:<16} {int(counts.get(c, 0)):>3}")
    print()

print("letter: 33 pairs -- pLoF 13, pathogenicity 7, missense 5, "
      "HMM constraint 5, GERP RS 3")
