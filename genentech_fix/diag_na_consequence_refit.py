"""Refit every gene-group that contains missing-consequence Epi25 variants, with
those variants excluded (as the Methods intend), and compare p_unified."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config, prepare_regression_input, fit_per_gene, GROUP_KEY
cfg = get_config()
floor, thr = cfg.param("min_variants"), cfg.param("exome_wide_threshold")
df = prepare_regression_input(pd.read_csv(cfg.result("constraint_gerp_am_epi25_variants_FIXED.tsv.gz"), sep="\t"),
                              clip_epsilon=cfg.param("clip_epsilon"))
na = (df.pLoF_ind == 0) & (df.missense_ind == 0)
affected = df.loc[na, GROUP_KEY].drop_duplicates()
sub = df.merge(affected, on=GROUP_KEY)
print(f"affected gene-groups: {len(affected):,}; rows before {len(sub):,}, missing-consequence rows removed {int(((sub.pLoF_ind==0)&(sub.missense_ind==0)).sum()):,}")
sub = sub[~((sub.pLoF_ind == 0) & (sub.missense_ind == 0))]
refit = fit_per_gene(sub, progress=False)[GROUP_KEY + ["n_variants", "p_unified"]].rename(
    columns={"n_variants": "n_new", "p_unified": "p_new"})
full = pd.read_csv(cfg.result("epilepsy_unified_model_pvalues_FIXED.tsv"), sep="\t")[GROUP_KEY + ["n_variants", "p_unified"]]
c = full.merge(refit, on=GROUP_KEY)
c["p_unified"] = c.p_unified.to_numpy(dtype=float); c["p_new"] = c.p_new.to_numpy(dtype=float)
c = c[c.n_variants >= floor]
c["dlog10"] = np.log10(c.p_new) - np.log10(c.p_unified)
print(f"compared (n>={floor} before): {len(c):,}; now below the floor: {int((c.n_new < floor).sum())}")
print("median |change in -log10 p|:", round(float(c.dlog10.abs().median()), 3))
for lab, t in (("exome-wide", thr), ("suggestive", 1e-4)):
    was, now = c.p_unified < t, (c.p_new < t) & (c.n_new >= floor)
    print(f"\n{lab}: significant before {int(was.sum())}, after {int(now.sum())}; lost {int((was & ~now).sum())}, gained {int((now & ~was).sum())}")
    x = c[was | now].sort_values("p_unified")
    print(x[["gene_name", "group", "n_variants", "n_new", "p_unified", "p_new"]].to_string(index=False, float_format=lambda v: f"{v:.2e}"))
