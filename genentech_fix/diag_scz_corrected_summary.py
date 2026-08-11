"""Published vs corrected schizophrenia significance, and where PCDHA4 sits."""
import sys
sys.path.insert(0, "src")
import pandas as pd
from wgs_constraint import get_config
cfg = get_config(); R = cfg.results_dir
thr = cfg.param("exome_wide_threshold")
pub = pd.read_csv(R / "table_A2_schizophrenia_full_p_value_comparison_n25.csv")
fix = pd.read_csv(R / "table_A2_schizophrenia_full_p_value_comparison_n25_FIXED.csv")
U = "Unified Model p-value"
print(f"published rows {len(pub):,}   corrected rows {len(fix):,}")
print()
for nm, d in (("PUBLISHED", pub), ("CORRECTED", fix)):
    s = d[d[U] < thr].sort_values(U)
    print(f"{nm}: {len(s)} significant at p < {thr:g}")
    for _, r in s.iterrows():
        print(f"    {str(r['Gene Name']):<12} {r[U]:.3e}   n_variants={r['# of Variants']}")
print()
m = pub.merge(fix, on=["Gene ID", "Group"], suffixes=("_pub", "_fix"))
lost = m[(m[U + "_pub"] < thr) & (m[U + "_fix"] >= thr)]
gained = m[(m[U + "_pub"] >= thr) & (m[U + "_fix"] < thr)]
print(f"lost significance : {len(lost)}")
for _, r in lost.iterrows():
    print(f"    {str(r['Gene Name_pub']):<12} {r[U+'_pub']:.3e} -> {r[U+'_fix']:.3e}"
          f"   n_variants {r['# of Variants_pub']} -> {r['# of Variants_fix']}")
print(f"gained significance: {len(gained)}")
for _, r in gained.iterrows():
    print(f"    {str(r['Gene Name_pub']):<12} {r[U+'_pub']:.3e} -> {r[U+'_fix']:.3e}")
print()
for g in ("PCDHA4", "PCDHGA5", "XPO7"):
    sub = m[m["Gene Name_pub"] == g]
    if sub.empty:
        print(f"  {g:<9} not present in the n>=25 comparison table")
    for _, r in sub.iterrows():
        print(f"  {g:<9} {r['Group']:<6} pub {r[U+'_pub']:.3e} -> fix {r[U+'_fix']:.3e}"
              f"   n {r['# of Variants_pub']} -> {r['# of Variants_fix']}")
