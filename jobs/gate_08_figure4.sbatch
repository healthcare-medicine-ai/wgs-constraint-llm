#!/bin/bash
#SBATCH --job-name=gate08
#SBATCH --partition=mrivas
#SBATCH --time=04:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Figures 4a and 4b had no gate job. They were only ever checked by
# genentech_fix/compare_figures.py, which passes anything under a 1% pixel
# tolerance and exits 0 regardless -- so "Figures 4a/4b reproduce" was never an
# enforced claim. This demands pixel-exactness and exits nonzero otherwise.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1

$PY pipelines/08_hmm_gerp_joint_distribution.py --as-published || exit 1

echo
echo "=== pixel comparison against the published panels ==="
$PY - <<'PY'
import sys
from pathlib import Path

import numpy as np
from PIL import Image

R = Path("results")
BASES = [
    "Figure 4a: percent difference distribution of GERP RS vs HMM constraint",
    "Figure 4b: chi-square statistics of GERP RS vs HMM constraint",
]
ok = True
for base in BASES:
    p, q = R / f"{base}.png", R / f"{base}_ASPUBLISHED.png"
    if not (p.exists() and q.exists()):
        print(f"  {base[:38]}: missing published={p.exists()} new={q.exists()}")
        ok = False
        continue
    a, b = Image.open(p), Image.open(q)
    print(f"  {base[:38]}: published {a.size}  new {b.size}")
    if a.size != b.size:
        print("      dimensions differ")
        ok = False
        continue
    aa = np.asarray(a.convert("RGB"), np.int16)
    bb = np.asarray(b.convert("RGB"), np.int16)
    d = np.abs(aa - bb)
    frac = (d.max(axis=2) > 0).mean()
    print(f"      differing {100*frac:.4f}%  max delta {d.max()}  "
          f"-> {'IDENTICAL' if frac == 0 else 'differs'}")
    if frac:
        ys, xs = np.nonzero(d.max(axis=2) > 0)
        print(f"      bounding box x {xs.min()}-{xs.max()} y {ys.min()}-{ys.max()}")
    ok &= frac == 0
print()
print("GATE 08:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
