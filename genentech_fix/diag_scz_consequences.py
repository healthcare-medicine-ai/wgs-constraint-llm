"""Recover which SCHEMA consequences the September 2025 cached merge kept."""
import sys

import pandas as pd
sys.path.insert(0, "src")
from wgs_constraint import get_config

cfg = get_config()
key = ["chr", "pos", "ref", "alt", "gene_id", "group"]

src = pd.read_csv(cfg.data_dir / "SCHEMA_variant_results_hg38.tsv.gz", sep="\t",
                  compression="gzip", usecols=key + ["consequence"])
src = src.drop_duplicates(subset=key)
print(f"SCHEMA distinct keys: {len(src):,}", flush=True)

old = pd.read_csv(cfg.result("constraint_gerp_am_scz_variants.tsv.gz"),
                  sep="\t", usecols=key).drop_duplicates()
new = pd.read_csv(cfg.result("constraint_gerp_am_scz_variants_REPRO.tsv.gz"),
                  sep="\t", usecols=key).drop_duplicates()

both = new.merge(old.assign(_c=1), on=key, how="left").merge(src, on=key, how="left")
kept = both[both["_c"].notna()]
dropped = both[both["_c"].isna()]
print(f"distinct keys regenerated : {len(both):,}")
print(f"  present in cache        : {len(kept):,}")
print(f"  absent from cache       : {len(dropped):,}\n")

tab = pd.DataFrame({
    "in_cache": kept["consequence"].value_counts(),
    "not_in_cache": dropped["consequence"].value_counts(),
}).fillna(0).astype(int)
tab["pct_absent"] = (100 * tab["not_in_cache"]
                     / (tab["in_cache"] + tab["not_in_cache"])).round(2)
pd.set_option("display.width", 130)
print(tab.sort_values("not_in_cache", ascending=False).to_string())
print()
fully = tab[(tab["in_cache"] == 0) & (tab["not_in_cache"] > 0)].index.tolist()
print("consequences ENTIRELY absent from the cache (filtered in Sept 2025):")
for c in fully:
    print(f"    {c!r},")
