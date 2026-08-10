#!/usr/bin/env python
"""Stage 1 -- build the variant-level input for the Epi25 meta-regression.

Equivalent to `Epilepsy Analysis.ipynb` cells 4-16, with two corrections
relative to the Revision 2 submission. Both are toggleable so the submitted
behaviour can be reproduced exactly, which is what makes any difference
attributable to the corrections rather than to this rewrite.

  FIX 1  GERP coordinate off-by-one   (see src/wgs_constraint/gerp.py)
  FIX 2  AlphaMissense transcript fan-out (see src/wgs_constraint/alphamissense.py)

Usage
-----
  python pipelines/01_build_regression_input.py                    # both fixes
  python pipelines/01_build_regression_input.py --no-fix           # as submitted
  python pipelines/01_build_regression_input.py --gerp-fix on --am-fix off
  python pipelines/01_build_regression_input.py --chroms chr21,chr22   # smoke test
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wgs_constraint import (  # noqa: E402
    annotate_with_gerp, get_config, load_collapsed,
)
from wgs_constraint.alphamissense import KEY_COLUMNS, SCORE_COLUMN  # noqa: E402
from wgs_constraint.epi25 import (  # noqa: E402
    add_effect_sizes, load_epi25_variants, load_gene_annotation,
)

VARIANT_COLUMNS = [
    "chr", "pos", "ref", "alt", "gene_id", "gene_name", "group",
    "ac_case", "an_case", "ac_ctrl", "an_ctrl",
    "pLoF_ind", "missense_ind", "effect_size", "var_effect_size",
]


def log(msg):
    print(msg, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-suffix", default="_FIXED")
    ap.add_argument("--no-fix", action="store_true",
                    help="reproduce the submitted behaviour: both fixes off")
    ap.add_argument("--gerp-fix", choices=["on", "off"], default=None)
    ap.add_argument("--am-fix", choices=["on", "off"], default=None)
    ap.add_argument("--chroms", default=None,
                    help="comma-separated chromosome subset, for smoke tests")
    args = ap.parse_args()

    base = not args.no_fix
    fix_gerp = base if args.gerp_fix is None else (args.gerp_fix == "on")
    fix_am = base if args.am_fix is None else (args.am_fix == "on")

    cfg = get_config()
    log(cfg.describe())
    log(f"GERP off-by-one fix : {fix_gerp}")
    log(f"AlphaMissense collapse : {fix_am}")
    log("")

    # -- 1. annotation and variants -------------------------------------
    log("[1/5] GENCODE annotation ...")
    genes = load_gene_annotation(cfg.data("gene_annotation"))
    log(f"      protein-coding CDS rows: {len(genes):,}")

    log("[2/5] Epi25 variants ...")
    variants = load_epi25_variants(
        cfg.data("epi25_variants"), genes,
        excluded_chromosomes=cfg.param("excluded_chromosomes"),
        excluded_consequences=cfg.param("excluded_consequences"))
    variants = add_effect_sizes(
        variants, max_allele_count=cfg.param("max_allele_count"))
    log(f"      qualifying variant rows: {len(variants):,}")
    variants = variants[VARIANT_COLUMNS]

    # -- 2. GERP ---------------------------------------------------------
    log("[3/5] GERP annotation of HMM predictions ...")
    predictions = pd.read_csv(cfg.derived("hmm_predictions"), sep="\t",
                              dtype={"chr": "string"})
    if "position" in predictions.columns and "pos" not in predictions.columns:
        predictions = predictions.rename(columns={"position": "pos"})
    predictions["pos"] = predictions["pos"].astype("int64")
    predictions["chr"] = predictions["chr"].astype("string")

    constraint = annotate_with_gerp(
        predictions,
        cfg.data("gerp_bigwig"),
        # This argument IS the defect. position_base=0 reproduces the submitted
        # behaviour by reading 1-based positions as though they were 0-based.
        position_base=1 if fix_gerp else 0,
        only_chromosomes=args.chroms.split(",") if args.chroms else None,
    )
    log(f"      positions with GERP: {len(constraint):,} "
        f"(out-of-bounds dropped: "
        f"{constraint.attrs.get('gerp_out_of_bounds_dropped', 0):,})")

    # -- 3. AlphaMissense -------------------------------------------------
    log("[4/5] AlphaMissense ...")
    wanted = pd.MultiIndex.from_arrays(
        [variants[c] for c in KEY_COLUMNS]).unique()
    log(f"      distinct Epi25 variants: {len(wanted):,}")
    if fix_am:
        alphamissense = load_collapsed(
            cfg.data("alphamissense"), restrict_to=wanted,
            how=cfg.param("alphamissense_collapse"))
    else:
        alphamissense = _load_uncollapsed(cfg.data("alphamissense"), wanted)

    # -- 4. merge and write -----------------------------------------------
    log("[5/5] merging and writing ...")
    merged = pd.merge(
        constraint[["chr", "pos", "prob_0", "GERP_RS"]],
        pd.merge(variants, alphamissense, on=KEY_COLUMNS, how="left"),
        on=["chr", "pos"], how="inner")

    distinct = merged[["chr", "pos", "ref", "alt", "gene_id",
                       "group"]].drop_duplicates().shape[0]
    log(f"      rows: {len(merged):,}  distinct (variant,gene,group): "
        f"{distinct:,}  inflation: {len(merged) / max(distinct, 1):.4f}")

    out = cfg.result(f"constraint_gerp_am_epi25_variants{args.out_suffix}.tsv.gz")
    merged.to_csv(out, index=False, compression="gzip", sep="\t")
    log(f"      wrote {out} ({os.path.getsize(out) / 1e6:.1f} MB)")


def _load_uncollapsed(path, restrict_to):
    """Submitted behaviour: no de-duplication, so isoforms fan out on merge."""
    from wgs_constraint.alphamissense import _DEFAULT_CHUNK, _HEADER_ROW, _RENAME
    kept = []
    reader = pd.read_csv(path, sep="\t", header=_HEADER_ROW,
                         usecols=list(_RENAME) + [SCORE_COLUMN],
                         chunksize=_DEFAULT_CHUNK)
    for chunk in reader:
        chunk = chunk.rename(columns=_RENAME)
        index = pd.MultiIndex.from_arrays([chunk[c] for c in KEY_COLUMNS])
        chunk = chunk[index.isin(restrict_to)]
        if len(chunk):
            kept.append(chunk)
    frame = pd.concat(kept, ignore_index=True)
    print(f"      AlphaMissense: {len(frame):,} rows, NOT collapsed "
          f"(reproducing submitted behaviour)", flush=True)
    return frame


if __name__ == "__main__":
    main()
