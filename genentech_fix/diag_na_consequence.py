"""Epi25 variants with a missing consequence: do they reach the meta-regression?

load_epi25_variants excludes consequence == "NA", but pandas reads the literal
"NA" as NaN, and NaN != "NA" is True, so those rows may pass the filter and
enter the model with both variant-class indicators set to 0.
"""
import pandas as pd
raw = pd.read_csv("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/epi25_variant_results.tsv.gz", sep="\t",
                  usecols=["variant_id", "group", "consequence"])
print("raw rows with consequence read as missing (NaN):", int(raw.consequence.isna().sum()),
      "| literal 'NA' strings:", int((raw.consequence == "NA").sum()))
inp = pd.read_csv("results/constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t",
                  usecols=["chr", "pos", "ref", "alt", "gene_id", "gene_name", "group", "pLoF_ind", "missense_ind"],
                  dtype={"chr": str})
nei = inp[(inp.pLoF_ind == 0) & (inp.missense_ind == 0)].copy()
print("regression-input rows with neither indicator:", len(nei))
raw[["chr", "pos", "ref", "alt"]] = raw.variant_id.str.split(":", expand=True)
raw["pos"] = raw.pos.astype(int)
nei["pos"] = nei.pos.astype(int)
m = nei.merge(raw[["chr", "pos", "ref", "alt", "group", "consequence"]], on=["chr", "pos", "ref", "alt", "group"], how="left")
print("their raw Epi25 consequence:", m.consequence.fillna("<missing>").value_counts().to_dict())
key = ["gene_id", "gene_name", "group"]
cnt = nei.groupby(key).size().rename("n_neither").reset_index()
tot = inp.groupby(key).size().rename("n_total").reset_index()
cnt = cnt.merge(tot, on=key)
full = pd.read_csv("results/epilepsy_unified_model_pvalues_FIXED.tsv", sep="\t")
full = full[full.n_variants >= 25][key + ["p_unified"]]
c = cnt.merge(full, on=key)
print("gene-groups (n>=25) containing such variants:", len(c))
for lab, t in (("exome-wide", 3.4e-7), ("suggestive", 1e-4)):
    x = c[c.p_unified < t].sort_values("p_unified")
    print(f"  {lab}: {len(x)} gene-groups")
    if len(x):
        print(x.head(40).to_string(index=False))
