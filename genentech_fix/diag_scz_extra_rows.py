"""Characterise the rows present when the SCZ input is regenerated but absent
from the September 2025 cache."""
import pandas as pd
R = "results/"
key = ["chr", "pos", "ref", "alt", "gene_id", "group"]
old = pd.read_csv(R + "constraint_gerp_am_scz_variants.tsv.gz", sep="\t")
new = pd.read_csv(R + "constraint_gerp_am_scz_variants_REPRO.tsv.gz", sep="\t")
print(f"cached      : {len(old):,} rows, columns {list(old.columns)}")
print(f"regenerated : {len(new):,} rows")
print()

merged = new.merge(old[key].assign(_in_cache=1), on=key, how="left")
extra = merged[merged["_in_cache"].isna()]
print(f"rows only in the regenerated file: {len(extra):,}")
print()
print("null counts among those extra rows (vs the file as a whole):")
for c in new.columns:
    e = extra[c].isna().mean() * 100
    a = new[c].isna().mean() * 100
    flag = "   <-- distinctive" if abs(e - a) > 20 else ""
    print(f"  {c:<20} extra {e:6.2f}%   overall {a:6.2f}%{flag}")
print()
# Are these duplicate (chr,pos,ref,alt,gene,group) keys, i.e. extra fan-out?
dup_new = new.duplicated(subset=key).sum()
dup_old = old.duplicated(subset=key).sum()
print(f"duplicate keys, cached      : {dup_old:,}")
print(f"duplicate keys, regenerated : {dup_new:,}")
print()
print("distinct keys:")
print(f"  cached      : {old[key].drop_duplicates().shape[0]:,}")
print(f"  regenerated : {new[key].drop_duplicates().shape[0]:,}")
