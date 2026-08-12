"""Where does the published-table route differ from the bigWig route?"""
import numpy as np
import pandas as pd

R = "results/"
a = pd.read_csv(R + "constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t")
b = pd.read_csv(R + "constraint_gerp_am_epi25_variants_VIAFIG.tsv.gz", sep="\t")
print(f"rows {len(a):,} vs {len(b):,}")
print(f"columns identical order: {list(a.columns) == list(b.columns)}")
print()
for c in a.columns:
    x, y = a[c], b[c]
    if np.issubdtype(x.dtype, np.number) and np.issubdtype(y.dtype, np.number):
        same = ((x == y) | (x.isna() & y.isna())).mean()
        if same == 1.0:
            print(f"  {c:<20} identical")
            continue
        xv, yv = x.to_numpy(float), y.to_numpy(float)
        f = np.isfinite(xv) & np.isfinite(yv) & (xv != 0)
        rel = np.abs(xv[f] - yv[f]) / np.abs(xv[f])
        print(f"  {c:<20} bit-identical {100*same:6.2f}%   max rel {rel.max():.3e}"
              f"   median rel {np.median(rel[rel > 0]) if (rel > 0).any() else 0:.3e}")
    else:
        same = (x.astype(str) == y.astype(str)).mean()
        print(f"  {c:<20} {'identical' if same == 1 else f'differs {100*same:.2f}% same'}")

print()
print("row order identical on the key?",
      (a[["chr", "pos", "ref", "alt", "gene_id", "group"]].astype(str).agg("|".join, axis=1)
       == b[["chr", "pos", "ref", "alt", "gene_id", "group"]].astype(str).agg("|".join, axis=1)
       ).all())
