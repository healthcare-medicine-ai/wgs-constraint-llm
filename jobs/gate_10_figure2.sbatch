#!/bin/bash
#SBATCH --job-name=gate10b
#SBATCH --partition=mrivas
#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/10_rgc_aou_joint.py --suffix _REGEN
echo
echo "=== pixel comparison against the published panels ==="
$PY - <<'PY'
import numpy as np
from pathlib import Path
from PIL import Image
R = Path("results")
for base in ["Figure 2a: observed joint distribution of RGC vs AoU WES constraint predictions",
             "Figure 2b: chi-square statistics of RGC vs AoU WES constraint predictions"]:
    p, q = R/f"{base}.png", R/f"{base}_REGEN.png"
    if not (p.exists() and q.exists()):
        print(f"  {base[:30]}: missing"); continue
    a, b = Image.open(p), Image.open(q)
    print(f"  {base[:30]}: published {a.size}  new {b.size}")
    if a.size != b.size: print("      dimensions differ"); continue
    aa=np.asarray(a.convert("RGB"),np.int16); bb=np.asarray(b.convert("RGB"),np.int16)
    d=np.abs(aa-bb); frac=(d.max(axis=2)>0).mean()
    print(f"      differing {100*frac:.4f}%  max delta {d.max()}  -> {'IDENTICAL' if frac==0 else 'differs'}")
PY
