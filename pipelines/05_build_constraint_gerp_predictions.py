#!/usr/bin/env python
"""Stage 5 -- rebuild the genome-wide constraint + GERP prediction table.

`HMM_rgc_ALL_RS_merged_predictions.tsv.gz` is consumed by
`notebooks/Constraint Measures Comparison.ipynb`, which produces the HMM-versus-
GERP joint distribution panels -- published Figure 4, becoming Figure 3C/3D in
Revision 3.

That table was built with the GERP coordinate defect baked in, so the figure
derived from it is affected as well. This regenerates it through the corrected
shared annotation. It is separate from the epilepsy pipeline, which builds its
own constraint table in memory and never reads this file.

  python pipelines/05_build_constraint_gerp_predictions.py
  python pipelines/05_build_constraint_gerp_predictions.py --no-fix   # as published
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import annotate_with_gerp, get_config  # noqa: E402

OUTPUT_STEM = "HMM_rgc_ALL_RS_merged_predictions"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fix", action="store_true",
                    help="reproduce the published (shifted) table")
    ap.add_argument("--out-suffix", default=None,
                    help="default: '' when fixed, '_ASPUBLISHED' when --no-fix")
    ap.add_argument("--chroms", default=None)
    args = ap.parse_args()

    fix = not args.no_fix
    suffix = args.out_suffix if args.out_suffix is not None else (
        "" if fix else "_ASPUBLISHED")

    cfg = get_config()
    print(cfg.describe(), flush=True)
    print(f"GERP off-by-one fix : {fix}\n", flush=True)

    print("[1/3] reading HMM predictions ...", flush=True)
    pred = pd.read_csv(cfg.derived("hmm_predictions"), sep="\t",
                       dtype={"chr": "string"})
    if "position" in pred.columns and "pos" not in pred.columns:
        pred = pred.rename(columns={"position": "pos"})
    pred["pos"] = pred["pos"].astype("int64")
    pred["chr"] = pred["chr"].astype("string")
    print(f"      {len(pred):,} positions", flush=True)

    print("[2/3] annotating with GERP ...", flush=True)
    merged = annotate_with_gerp(
        pred, cfg.data("gerp_bigwig"),
        position_base=1 if fix else 0,
        only_chromosomes=args.chroms.split(",") if args.chroms else None)
    print(f"      {len(merged):,} positions with GERP "
          f"(out-of-bounds dropped: "
          f"{merged.attrs.get('gerp_out_of_bounds_dropped', 0):,})", flush=True)

    print("[3/3] writing ...", flush=True)
    out = cfg.result(f"{OUTPUT_STEM}{suffix}.tsv.gz")
    merged.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"      wrote {out} ({os.path.getsize(out) / 1e6:.1f} MB)", flush=True)

    print("\nGERP summary (sanity check against the published table):")
    print(merged["GERP_RS"].describe().to_string())


if __name__ == "__main__":
    main()
