"""Is the n_variants >= 25 floor applied to true variant counts after the fix?"""
import sys
import pandas as pd
sys.path.insert(0, "src")
from wgs_constraint import get_config, prepare_regression_input
cfg = get_config(); R = cfg.results_dir
floor, thr = cfg.param("min_variants"), cfg.param("exome_wide_threshold")
K = ["gene_id", "gene_name", "group"]

# 1. corrected n_variants vs distinct variants actually entering each regression
d = pd.read_csv(R / "constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t")
d = prepare_regression_input(d, clip_epsilon=cfg.param("clip_epsilon"))
rows = d.groupby(K).size().rename("rows")
distinct = d.drop_duplicates(K + ["chr", "pos", "ref", "alt"]).groupby(K).size().rename("distinct")
fx = pd.read_csv(R / "epilepsy_unified_model_pvalues_FIXED.tsv", sep="\t")
m = fx.set_index(K)[["n_variants"]].join(rows).join(distinct)
print("1. corrected n_variants vs distinct variants per gene-group")
print(f"   gene-groups                      : {len(m):,}")
print(f"   n_variants == rows               : {(m.n_variants == m.rows).mean():.4%}")
print(f"   n_variants == distinct variants  : {(m.n_variants == m.distinct).mean():.4%}")
del d

# 2. R2 inflation at the floor
r2 = pd.read_csv(R / "epilepsy_unified_model_pvalues.tsv", sep="\t")
j = r2.merge(fx, on=K, suffixes=("_r2", "_fix"), how="outer")
crossed = j[(j.n_variants_r2 >= floor) & (j.n_variants_fix < floor)]
ratio = (j.n_variants_r2 / j.n_variants_fix)
print(f"\n2. R2 inflation")
print(f"   median R2/true n ratio           : {ratio.median():.2f}   max {ratio.max():.1f}")
print(f"   passed the floor in R2 only because of inflation: {len(crossed):,} gene-groups")
print(f"   ... of which exome-wide significant in R2: "
      f"{int((pd.to_numeric(crossed.p_unified_r2, errors='coerce') < thr).sum())}")

# 3. why each lost R2 significant pair was lost
lost = j[(pd.to_numeric(j.p_unified_r2, errors="coerce") < thr) & (j.n_variants_r2 >= floor)
         & ~((pd.to_numeric(j.p_unified_fix, errors="coerce") < thr) & (j.n_variants_fix >= floor))]
print(f"\n3. R2-significant pairs lost after correction: {len(lost)}")
for _, r in lost.sort_values("gene_name").iterrows():
    why = "below floor" if r.n_variants_fix < floor else "weaker signal"
    print(f"   {r.gene_name:<9} {r.group:<5} n {int(r.n_variants_r2):>5} -> {int(r.n_variants_fix):>4}"
          f"   p {float(r.p_unified_r2):.1e} -> {float(r.p_unified_fix):.1e}   [{why}]")
