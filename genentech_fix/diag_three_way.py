"""Do the two fresh runs agree with each other, or is each run different?"""
import numpy as np
import pandas as pd

R = "results/"
files = {
    "FIXED (committed)": R + "constraint_gerp_am_epi25_variants_FIXED.tsv.gz",
    "SELFCHK (rerun, bigWig)": R + "constraint_gerp_am_epi25_variants_SELFCHK.tsv.gz",
    "VIAFIG (rerun, published tbl)": R + "constraint_gerp_am_epi25_variants_VIAFIG.tsv.gz",
}
data = {k: pd.read_csv(v, sep="\t", usecols=["effect_size"])["effect_size"]
              .to_numpy(dtype=float) for k, v in files.items()}
for k, v in data.items():
    print(f"  {k:<32} {len(v):,} rows")
print()
keys = list(data)
for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
        a, b = data[keys[i]], data[keys[j]]
        same = (a == b).mean()
        f = np.isfinite(a) & np.isfinite(b) & (a != 0)
        rel = np.abs(a[f] - b[f]) / np.abs(a[f])
        print(f"  {keys[i]:<32} vs {keys[j]:<32} "
              f"bit-identical {100*same:8.4f}%   max rel {rel.max():.3e}")
