#!/usr/bin/env python
"""Build the variant-level input for the Epi25 meta-regression.

Steps
  1. Annotate the genome-wide HMM constraint predictions with GERP RS from the
     hg38 bigWig (see src/wgs_constraint/gerp.py for the coordinate convention).
  2. Load Epi25 variant results, keep autosomal coding variants with
     ac_case + ac_ctrl <= 5, and compute Haldane-Anscombe effect sizes.
  3. Attach AlphaMissense pathogenicity, one score per variant (maximum over
     transcripts).
  4. Join everything on (chr, pos).

Inputs (config.yaml)
  data/    gencode.v44.basic.annotation.gtf.gz, epi25_variant_results.tsv.gz,
           AlphaMissense_hg38.tsv.gz, All_hg38_RS.bw
  results/ HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz

Outputs (results/)
  HMM_rgc_ALL_RS_merged_predictions.tsv.gz   constraint + GERP, every position
  epilepsy_regression_input.tsv.gz           one row per (variant, gene, group)

Usage
  python pipelines/02_regression_input.py
  python pipelines/02_regression_input.py --constraint-gerp HMM_rgc_ALL_RS_merged_predictions.tsv.gz

--constraint-gerp skips step 1 and reads a constraint + GERP table written by
an earlier run of this stage. That avoids reading the 16.4 GB bigWig again, and
gives an identical regression input. Stages 06 and 09 read that table from
results/, so keep it there.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wgs_constraint import annotate_with_gerp, get_config, load_collapsed  # noqa: E402
from wgs_constraint.alphamissense import KEY_COLUMNS  # noqa: E402
from wgs_constraint.epi25 import (  # noqa: E402
    add_effect_sizes, load_epi25_variants, load_gene_annotation,
)

REGRESSION_INPUT = "epilepsy_regression_input.tsv.gz"

VARIANT_COLUMNS = [
    "chr", "pos", "ref", "alt", "gene_id", "gene_name", "group",
    "ac_case", "an_case", "ac_ctrl", "an_ctrl",
    "pLoF_ind", "missense_ind", "effect_size", "var_effect_size",
]


def log(msg):
    print(msg, flush=True)


def constraint_with_gerp(cfg):
    """Annotate the HMM predictions with GERP and write the genome-wide table."""
    predictions = pd.read_csv(cfg.derived("hmm_predictions"), sep="\t",
                              dtype={"chr": "string"})
    predictions["pos"] = predictions["pos"].astype("int64")
    log(f"      HMM predictions: {len(predictions):,} positions")

    table = annotate_with_gerp(predictions, cfg.data("gerp_bigwig"))
    log(f"      with GERP: {len(table):,} (out of range: "
        f"{table.attrs['gerp_out_of_bounds_dropped']:,})")

    out = cfg.derived("constraint_gerp")
    table.to_csv(out, sep="\t", index=False, compression="gzip")
    log(f"      wrote {out}")
    return table


def read_constraint_with_gerp(path):
    # round_trip parsing recovers exactly the values the table was written
    # from; pandas' default parser can be off by one unit in the last place.
    table = pd.read_csv(path, sep="\t",
                        usecols=["chr", "pos", "prob_0", "GERP_RS"],
                        dtype={"chr": "string"}, float_precision="round_trip")
    table["pos"] = table["pos"].astype("int64")
    log(f"      read {path}: {len(table):,} positions")
    return table


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--constraint-gerp", type=Path, default=None,
                    help="existing constraint + GERP table (chr, pos, prob_0, "
                         "GERP_RS) to use instead of annotating from the bigWig")
    args = ap.parse_args()

    cfg = get_config()
    log(cfg.describe())

    log("[1/4] constraint + GERP ...")
    if args.constraint_gerp:
        constraint = read_constraint_with_gerp(args.constraint_gerp)
    else:
        constraint = constraint_with_gerp(cfg)

    log("[2/4] Epi25 variants ...")
    genes = load_gene_annotation(cfg.data("gene_annotation"))
    variants = load_epi25_variants(
        cfg.data("epi25_variants"), genes,
        excluded_chromosomes=cfg.param("excluded_chromosomes"),
        excluded_consequences=cfg.param("excluded_consequences"))
    variants = add_effect_sizes(
        variants, max_allele_count=cfg.param("max_allele_count"))
    variants = variants[VARIANT_COLUMNS]
    log(f"      qualifying variant rows: {len(variants):,}")

    log("[3/4] AlphaMissense ...")
    wanted = pd.MultiIndex.from_arrays(
        [variants[c] for c in KEY_COLUMNS]).unique()
    alphamissense = load_collapsed(
        cfg.data("alphamissense"), restrict_to=wanted,
        how=cfg.param("alphamissense_collapse"))

    log("[4/4] merging ...")
    merged = pd.merge(
        constraint[["chr", "pos", "prob_0", "GERP_RS"]],
        pd.merge(variants, alphamissense, on=KEY_COLUMNS, how="left"),
        on=["chr", "pos"], how="inner")
    out = cfg.result(REGRESSION_INPUT)
    merged.to_csv(out, index=False, compression="gzip", sep="\t")
    log(f"      wrote {out} ({len(merged):,} rows)")


if __name__ == "__main__":
    main()
