"""Which schizophrenia models are numerically unstable, and do they matter?

Gate 11b passes: no gene-group changes significance. But 109 move by more than
an order of magnitude in p from input differences of order 1e-15. A p-value that
mobile is not meaningful to the precision it is reported at, so the question is
whether any of them reach the table the manuscript reports.
"""
import sys
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from wgs_constraint import get_config

cfg = get_config()
R = cfg.results_dir
thr = cfg.param("exome_wide_threshold")
floor = cfg.param("min_variants")

a = pd.read_csv(R / "schizophrenia_variants_unified_gerp_model_pvalues.tsv", sep="\t")
b = pd.read_csv(R / "schizophrenia_variants_unified_gerp_model_pvalues_REPRO.tsv", sep="\t")
m = a.merge(b, on=["gene_id", "group"], suffixes=("_a", "_b"))
for c in ("p_unified_a", "p_unified_b", "n_variants_a"):
    m[c] = pd.to_numeric(m[c], errors="coerce")
m = m[(m["p_unified_a"] > 0) & (m["p_unified_b"] > 0)].copy()
m["dlog"] = np.abs(np.log10(m["p_unified_a"]) - np.log10(m["p_unified_b"]))

mov = m[m["dlog"] > 1]
print(f"gene-groups compared      : {len(m):,}")
print(f"moving > 1 log unit       : {len(mov):,}  ({100*len(mov)/len(m):.4f}%)")
print()
print("n_variants of the movers vs everyone else:")
print(f"  movers   median {mov['n_variants_a'].median():.0f}   "
      f"min {mov['n_variants_a'].min():.0f}   max {mov['n_variants_a'].max():.0f}")
print(f"  others   median {m.loc[m['dlog']<=1,'n_variants_a'].median():.0f}")
print()
print("where do the movers sit in p?")
for lo, hi, lbl in [(0, thr, f"p < {thr:g} (significant)"),
                    (thr, 1e-4, "3.4e-7 to 1e-4 (suggestive)"),
                    (1e-4, 1e-2, "1e-4 to 1e-2"),
                    (1e-2, 1.01, "p > 1e-2")]:
    n = int(((mov["p_unified_a"] >= lo) & (mov["p_unified_a"] < hi)).sum())
    print(f"  {lbl:<30} {n:>5}")
print()
above = mov[mov["n_variants_a"] >= floor]
print(f"movers passing the n_variants >= {floor} floor: {len(above):,}")
if len(above):
    schema = pd.read_csv(cfg.data_dir / "SCHEMA_gene_results.tsv.bgz", sep="\t",
                         compression="gzip")
    rep = above.merge(schema, on=["gene_id", "group"])
    print(f"  of those, reaching the reported table (SCHEMA merge): {len(rep):,}")
    if len(rep):
        cols = ["gene_id", "group", "n_variants_a", "p_unified_a", "p_unified_b", "dlog"]
        print(rep.sort_values("dlog", ascending=False)[cols].head(15).to_string(index=False))
else:
    print("  -> every unstable model is below the inclusion floor and never reported")
