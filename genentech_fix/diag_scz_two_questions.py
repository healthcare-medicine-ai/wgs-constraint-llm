"""Separate two questions the failing gate conflates.

A. Is the REGRESSION faithful? Run fit_per_gene on the CACHED input -- the exact
   bytes the notebook fed its own regression -- and compare to the published
   p-values. Any difference here is mine, independent of the input.

B. Does the ULP-level effect_size difference MATTER? Relative difference is the
   wrong metric for a p-value spanning 1e-300 to 1. What matters is agreement in
   -log10(p) and whether any gene crosses the significance threshold.

C. What characterises the 19.6% of rows whose effect_size differs?
"""
import sys
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from wgs_constraint import get_config, fit_per_gene, prepare_regression_input

cfg = get_config()
eps = cfg.param("clip_epsilon")
thr = cfg.param("exome_wide_threshold")
R = cfg.results_dir

pub = pd.read_csv(R / "schizophrenia_variants_unified_gerp_model_pvalues.tsv", sep="\t")

print("=" * 70)
print("A. REGRESSION FIDELITY -- my fit on the notebook's own cached input")
print("=" * 70, flush=True)
cached = pd.read_csv(R / "constraint_gerp_am_scz_variants.tsv.gz", sep="\t")
fit_cached = fit_per_gene(prepare_regression_input(cached, clip_epsilon=eps))
m = pub.merge(fit_cached, on=["gene_id", "group"], suffixes=("_a", "_b"))
print(f"  gene-groups compared: {len(m):,}")
for c in ["n_variants", "p_constraint", "p_gerp", "p_unified"]:
    x = pd.to_numeric(m[c + "_a"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(m[c + "_b"], errors="coerce").to_numpy(float)
    f = np.isfinite(x) & np.isfinite(y)
    print(f"  {c:<14} bit-identical {100*(x[f]==y[f]).mean():6.2f}%")

print()
print("=" * 70)
print("B. DOES THE INPUT DIFFERENCE MATTER? published vs REPRO")
print("=" * 70, flush=True)
rep = pd.read_csv(R / "schizophrenia_variants_unified_gerp_model_pvalues_REPRO.tsv", sep="\t")
m = pub.merge(rep, on=["gene_id", "group"], suffixes=("_a", "_b"))
x = pd.to_numeric(m["p_unified_a"], errors="coerce").to_numpy(float)
y = pd.to_numeric(m["p_unified_b"], errors="coerce").to_numpy(float)
f = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
dlog = np.abs(np.log10(x[f]) - np.log10(y[f]))
print(f"  gene-groups: {f.sum():,}")
print(f"  |delta log10(p)|   max {dlog.max():.3e}   median {np.median(dlog):.3e}")
print(f"    within 1e-6 : {100*(dlog < 1e-6).mean():.4f}%")
print(f"    within 1e-3 : {100*(dlog < 1e-3).mean():.4f}%")
sa, sb = x[f] < thr, y[f] < thr
print(f"  significant published {sa.sum()}   REPRO {sb.sum()}   disagreements {int((sa!=sb).sum())}")
if (sa != sb).any():
    d = m[f].iloc[np.where(sa != sb)[0]]
    for _, r in d.iterrows():
        print(f"    {r.get('gene_name_a', r['gene_id']):<14} pub={r['p_unified_a']:.3e} rep={r['p_unified_b']:.3e}")

print()
print("=" * 70)
print("C. WHICH ROWS HAVE A DIFFERENT effect_size?")
print("=" * 70, flush=True)
key = ["chr", "pos", "ref", "alt", "gene_id", "group"]
new = pd.read_csv(R / "constraint_gerp_am_scz_variants_REPRO.tsv.gz", sep="\t")
j = cached[key + ["ac_case", "an_case", "ac_ctrl", "an_ctrl", "effect_size"]].merge(
    new[key + ["effect_size"]], on=key, suffixes=("_a", "_b"))
j["same"] = j["effect_size_a"] == j["effect_size_b"]
print(f"  rows: {len(j):,}   differing: {int((~j['same']).sum()):,} "
      f"({100*(~j['same']).mean():.2f}%)")
for col in ["ac_case", "ac_ctrl", "an_case", "an_ctrl"]:
    print(f"  {col:<9} differing-rows mean {j.loc[~j['same'], col].mean():12.2f}   "
          f"same-rows mean {j.loc[j['same'], col].mean():12.2f}")
print("\n  ac_case value counts among differing vs same rows:")
print(pd.DataFrame({
    "differing": j.loc[~j["same"], "ac_case"].value_counts().head(8),
    "same": j.loc[j["same"], "ac_case"].value_counts().head(8)}).fillna(0).astype(int).to_string())
