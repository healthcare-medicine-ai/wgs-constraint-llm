#!/usr/bin/env python
"""Patch the GERP merge out of the notebooks and into the shared module.

The same ~30 lines of bigWig lookup were pasted into three notebooks, each
carrying the identical coordinate defect. This replaces every copy with a call
to `wgs_constraint.annotate_with_gerp`, so the logic exists once.

Also clears cell outputs. They were produced by the defective code, so leaving
them would ship wrong numbers rendered as authoritative output, and the
embedded PNGs are most of the repository's size.

Run with --dry-run first; it reports what it would change without writing.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

MARKER = "gerp_vals = gerp_vec[pos_idx]"
ASSIGN = re.compile(r"^\s*(\w+)\s*=\s*pd\.concat\(merged_chunks", re.M)
PREDVAR = re.compile(r"^\s*(\w+)\s*=\s*pd\.read_csv\(\s*(\w+)\s*,\s*sep", re.M)

TEMPLATE = '''\
# GERP RS annotation.
#
# This cell previously inlined a bigWig lookup that indexed a 0-based value
# vector with 1-based HMM positions, shifting every score one base downstream.
# The lookup now lives in src/wgs_constraint/gerp.py, so the three notebooks
# that need it cannot drift apart again; the coordinate convention is an
# explicit argument there rather than an assumption.
#
# Loading the predictions is kept here because later cells use `pred`.
import sys
sys.path.insert(0, "../src")

import numpy as np
import pandas as pd

from wgs_constraint import annotate_with_gerp

# --- load HMM predictions (once) ---
pred = pd.read_csv(hmm_predictions_path, sep="\\t", dtype={{"chr": "string"}})
if "position" in pred.columns and "pos" not in pred.columns:
    pred = pred.rename(columns={{"position": "pos"}})
pred["pos"] = pred["pos"].astype(np.int64)
pred["chr"] = pred["chr"].astype("string")
pred_cols = pred.columns.tolist()

{target} = annotate_with_gerp(
    pred,
    gerp_file_path,
    position_base=1,      # HMM positions are 1-based
)
{target}["chr"] = {target}["chr"].astype("category")
{target}
'''


def patch(path: Path, dry_run: bool) -> bool:
    nb = json.loads(path.read_text(encoding="utf-8"))
    changed = False
    n_outputs_cleared = 0
    patched_cells = []

    for i, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))

        if MARKER in source:
            m = ASSIGN.search(source)
            if not m:
                print(f"  !! cell {i}: found the defect but no "
                      f"'X = pd.concat(merged_chunks' assignment; skipping")
                continue
            target = m.group(1)

            # The replacement re-creates `pred` itself, so the original cell
            # must have loaded it from `hmm_predictions_path` under that name.
            # Bail out rather than silently produce a notebook that cannot run.
            if "pred = pd.read_csv(hmm_predictions_path" not in source:
                print(f"  !! cell {i}: unexpected prediction loading; "
                      f"refusing to patch")
                continue

            cell["source"] = TEMPLATE.format(target=target).splitlines(
                keepends=True)
            cell["outputs"] = []
            cell["execution_count"] = None
            patched_cells.append((i, target))
            changed = True

        if cell.get("outputs"):
            n_outputs_cleared += len(cell["outputs"])
            cell["outputs"] = []
            cell["execution_count"] = None
            changed = True

    size_before = path.stat().st_size
    if changed and not dry_run:
        path.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n",
                        encoding="utf-8")

    print(f"\n{path.name}")
    if patched_cells:
        for i, target in patched_cells:
            print(f"  cell {i:>3}: GERP merge -> "
                  f"annotate_with_gerp(pred, ...) -> {target}")
    else:
        print("  no GERP merge cell found")
    print(f"  outputs cleared: {n_outputs_cleared}")
    if not dry_run and changed:
        print(f"  size: {size_before/1024:.0f} KB -> "
              f"{path.stat().st_size/1024:.0f} KB")
    return changed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("notebooks", nargs="+")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    any_changed = False
    for name in args.notebooks:
        p = Path(name)
        if not p.exists():
            print(f"missing: {p}", file=sys.stderr)
            continue
        any_changed |= patch(p, args.dry_run)

    print("\n" + ("DRY RUN -- nothing written" if args.dry_run else "written"))
    return 0 if any_changed else 1


if __name__ == "__main__":
    sys.exit(main())
