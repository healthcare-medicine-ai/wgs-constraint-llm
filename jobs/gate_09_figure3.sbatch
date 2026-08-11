#!/bin/bash
#SBATCH --job-name=gate09b
#SBATCH --partition=mrivas
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
# A failing gate exits nonzero. The first version of these scripts printed
# the verdict and exited 0, so a --dependency=afterok job downstream launched
# against a failing gate. A gate that cannot stop anything is not a gate.
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/09_gene_constraint_figures.py --as-published
echo
echo "=== pixel comparison against the published panels ==="
$PY - <<'PY'
import sys
import numpy as np
from pathlib import Path
from PIL import Image
R = Path("results")
ok = True
for base in ["Figure 3a: HMM vs gnomAD missense z-score per gene",
             "Figure 3b: HMM vs MTR per gene"]:
    p, q = R/f"{base}.png", R/f"{base}_ASPUBLISHED.png"
    if not (p.exists() and q.exists()):
        print(f"  {base[:34]}: missing")
        ok = False
        continue
    a, b = Image.open(p), Image.open(q)
    print(f"  {base[:34]}: published {a.size}  new {b.size}")
    if a.size != b.size:
        print("      dimensions differ")
        ok = False
        continue
    aa = np.asarray(a.convert("RGB"), np.int16); bb = np.asarray(b.convert("RGB"), np.int16)
    d = np.abs(aa-bb); frac = (d.max(axis=2) > 0).mean()
    print(f"      differing {100*frac:.4f}%  max delta {d.max()}  "
          f"-> {'IDENTICAL' if frac==0 else 'differs'}")
    ok &= frac == 0
print("
GATE 09:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
