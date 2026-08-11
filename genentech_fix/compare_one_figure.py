"""Pixel-compare one regenerated figure against its published counterpart.

Exits nonzero if the figures are missing, differently sized, or differ in any
pixel, so callers using --dependency=afterok are actually gated by the result.
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image

base, suffix = sys.argv[1], sys.argv[2]
R = Path("results")
p, q = R / f"{base}.png", R / f"{base}{suffix}.png"
if not (p.exists() and q.exists()):
    print(f"  missing: published={p.exists()} regenerated={q.exists()}")
    sys.exit(1)
a, b = Image.open(p), Image.open(q)
print(f"  published {a.size}   regenerated {b.size}")
if a.size != b.size:
    print("  dimensions differ -> cannot pixel-compare")
    sys.exit(1)
aa = np.asarray(a.convert("RGB"), np.int16)
bb = np.asarray(b.convert("RGB"), np.int16)
d = np.abs(aa - bb)
frac = (d.max(axis=2) > 0).mean()
print(f"  differing {100*frac:.4f}%   max channel delta {d.max()}")
print(f"  -> {'IDENTICAL' if frac == 0 else 'differs'}")
sys.exit(0 if frac == 0 else 1)
