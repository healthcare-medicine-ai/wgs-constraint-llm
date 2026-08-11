"""Pixel-compare a regenerated figure against its published counterpart.

Byte comparison is useless for PNGs (metadata differs per write), so this
compares decoded pixels. Used to gate figure-producing stages whose original
save call was never recorded in the notebook.
"""
import sys
import numpy as np
from pathlib import Path
from PIL import Image

R = Path("results")
PAIRS = [
    ("Figure 4a: percent difference distribution of GERP RS vs HMM constraint", "4a"),
    ("Figure 4b: chi-square statistics of GERP RS vs HMM constraint", "4b"),
]
suffix = sys.argv[1] if len(sys.argv) > 1 else "_ASPUBLISHED"
ok = True
for base, label in PAIRS:
    p, q = R / f"{base}.png", R / f"{base}{suffix}.png"
    if not (p.exists() and q.exists()):
        print(f"  {label}: missing ({p.exists()=} {q.exists()=})"); ok = False; continue
    a, b = Image.open(p), Image.open(q)
    print(f"  {label}: published {a.size}   regenerated {b.size}")
    if a.size != b.size:
        print("      dimensions differ -> cannot pixel-compare"); ok = False; continue
    aa = np.asarray(a.convert("RGB"), dtype=np.int16)
    bb = np.asarray(b.convert("RGB"), dtype=np.int16)
    d = np.abs(aa - bb)
    frac = (d.max(axis=2) > 0).mean()
    print(f"      pixels differing {100*frac:.4f}%   max channel delta {d.max()}")
    verdict = "IDENTICAL" if frac == 0 else ("near-identical" if frac < 0.01 else "differs")
    print(f"      -> {verdict}")
    ok &= frac < 0.01
print("\nGATE:", "pass" if ok else "not established")
