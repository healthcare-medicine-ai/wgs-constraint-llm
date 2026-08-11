"""Report where two renders of the same figure differ."""
import sys
import numpy as np
from pathlib import Path
from PIL import Image
base, suffix = sys.argv[1], sys.argv[2]
R = Path("results")
a = np.asarray(Image.open(R/f"{base}.png").convert("RGB"), np.int16)
b = np.asarray(Image.open(R/f"{base}{suffix}.png").convert("RGB"), np.int16)
d = (np.abs(a-b).max(axis=2) > 0)
ys, xs = np.nonzero(d)
print(f"  differing pixels : {d.sum():,} of {d.size:,} ({100*d.mean():.4f}%)")
print(f"  bounding box     : x {xs.min()}-{xs.max()}, y {ys.min()}-{ys.max()}")
print(f"  image size       : {a.shape[1]} x {a.shape[0]}")
print(f"  region           : {100*xs.min()/a.shape[1]:.0f}-{100*xs.max()/a.shape[1]:.0f}% across, "
      f"{100*ys.min()/a.shape[0]:.0f}-{100*ys.max()/a.shape[0]:.0f}% down")
# Is it one contiguous blob (a text label) or scattered (data)?
print(f"  distinct rows    : {len(np.unique(ys))}   distinct cols: {len(np.unique(xs))}")
