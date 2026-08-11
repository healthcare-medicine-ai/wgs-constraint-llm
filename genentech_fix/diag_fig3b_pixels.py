import numpy as np
from pathlib import Path
from PIL import Image
R = Path("results")
base = "Figure 3b: HMM vs MTR per gene"
a = np.asarray(Image.open(R / f"{base}.png").convert("RGB"), np.int16)
b = np.asarray(Image.open(R / f"{base}_ASPUBLISHED.png").convert("RGB"), np.int16)
d = (np.abs(a - b).max(axis=2) > 0)
ys, xs = np.nonzero(d)
print(f"image {a.shape[1]}x{a.shape[0]}   differing pixels: {d.sum():,} ({100*d.mean():.4f}%)")
print(f"  bounding box : x {xs.min()}-{xs.max()}   y {ys.min()}-{ys.max()}")
print(f"  box size     : {xs.max()-xs.min()+1} x {ys.max()-ys.min()+1}")
verdict = ("localised -- consistent with a single text annotation"
           if (ys.max() - ys.min()) < 60 else "SCATTERED -- not just a text label")
print(f"  verdict      : {verdict}")
