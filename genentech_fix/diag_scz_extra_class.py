import pandas as pd
R = "results/"
key = ["chr", "pos", "ref", "alt", "gene_id", "group"]
cols = key + ["pLoF_ind", "missense_ind", "am_pathogenicity"]
old = pd.read_csv(R + "constraint_gerp_am_scz_variants.tsv.gz", sep="\t", usecols=key)
new = pd.read_csv(R + "constraint_gerp_am_scz_variants_REPRO.tsv.gz", sep="\t", usecols=cols)
m = new.merge(old.assign(_c=1), on=key, how="left")
extra = m[m["_c"].isna()]

def show(name, d):
    ct = d.groupby(["pLoF_ind", "missense_ind"]).size()
    tot = len(d)
    print(f"  {name} (n={tot:,})")
    for (p, ms), n in ct.items():
        label = {(1,0): "pLoF", (0,1): "missense", (0,0): "NEITHER", (1,1): "both"}[(p, ms)]
        print(f"    pLoF={p} missense={ms}  {label:<9} {n:>10,}  ({100*n/tot:5.2f}%)")

print("indicator breakdown")
show("extra rows only", extra)
show("whole regenerated file", new)
print()
# ref/alt length -> SNV vs indel
for nm, d in (("extra", extra), ("all", new)):
    snv = ((d["ref"].str.len() == 1) & (d["alt"].str.len() == 1)).mean() * 100
    print(f"  {nm:<6} SNVs: {snv:5.2f}%")
