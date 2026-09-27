"""Are pLoF_ind and missense_ind collinear with the intercept?

If every variant in a gene-group is either pLoF or missense, then
pLoF_ind + missense_ind = 1 = the intercept column, so their individual
coefficients (and p-values) are not identified. Removing one of them then
leaves the fit, and its F-test, unchanged.
"""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from wgs_constraint import get_config
cfg = get_config(); R = cfg.results_dir
floor = cfg.param("min_variants")
df = pd.read_csv(R / "constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t",
                 usecols=["gene_name", "group", "pLoF_ind", "missense_ind"])
cls = np.select([df.pLoF_ind == 1, df.missense_ind == 1], ["pLoF", "missense"], "neither")
print("consequence classes across the whole table:")
print(pd.Series(cls).value_counts().to_string())
print(f"  rows flagged both: {int(((df.pLoF_ind == 1) & (df.missense_ind == 1)).sum())}")
df["neither"] = ((df.pLoF_ind == 0) & (df.missense_ind == 0)).astype(int)
g = df.groupby(["gene_name", "group"]).agg(n=("neither", "size"), pLoF=("pLoF_ind", "sum"),
                                           mis=("missense_ind", "sum"), neither=("neither", "sum"))
g = g[g.n >= floor]
print(f"\ngene-groups with n >= {floor}: {len(g):,}")
print(f"  no 'neither' variant (pLoF, missense and intercept collinear): "
      f"{int((g.neither == 0).sum()):,} ({100 * (g.neither == 0).mean():.1f}%)")
print(f"  all one class (a moderator is constant): "
      f"{int(((g.pLoF == 0) | (g.mis == 0)).sum()):,}")
