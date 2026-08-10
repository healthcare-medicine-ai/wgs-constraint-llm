#!/usr/bin/env python
"""
Step 1: build the variant-level regression input for the Epi25 unified
meta-regression.

Faithful port of `Epilepsy Analysis.ipynb` cells 4-16, with exactly two
corrections relative to the submitted version:

  FIX 1 -- GERP coordinate off-by-one.
      The notebook indexed a 0-based bigWig value vector with 1-based HMM
      positions:
          gerp_vec = bw.values(chrom, 0, chr_len)   # element i == 1-based pos i+1
          gerp_vals = gerp_vec[pos_idx]             # returns pos+1
      HMM positions are 1-based: `get_sequence()` in HMM Mutation Predictions
      indexes its masks directly with 1-based GTF/VCF coordinates. The correct
      index is therefore pos-1. The bounds guard is corrected in step with it;
      the old guard (pos >= 0) combined with a -1 index would wrap pos==0 to
      the last base of the chromosome, and (pos < chr_len) silently dropped
      the final base.

  FIX 2 -- AlphaMissense transcript fan-out.
      AlphaMissense_hg38.tsv.gz is the isoform-level release: 71,697,556 rows
      for 71,236,463 variants. Merging it un-deduplicated multiplied 0.77% of
      the regression rows, which inflated n_variants AND gave those variants
      extra weight in the WLS. Collapsed to one row per variant by MAX
      am_pathogenicity (decision: 2026-08-08; duplicates are almost entirely
      alternative transcripts of the same protein -- only 1,854 of 454,089
      span more than one uniprot id).

Everything else -- filters, thresholds, imputation, transformations -- is
unchanged from the submitted pipeline.

Output: constraint_gerp_am_epi25_variants.tsv.gz (path set by --out-suffix)
"""

import argparse
import os

import numpy as np
import pandas as pd
import pyBigWig
from tqdm import tqdm

DATA = "/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/"
RESULTS = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
           "osthoag/wgs-constraint-llm/results/")

GENE_ANNOTATION = DATA + "gencode.v44.basic.annotation.gtf.gz"
EPI25_VARIANTS = DATA + "epi25_variant_results.tsv.gz"
GERP_BIGWIG = DATA + "All_hg38_RS.bw"
ALPHAMISSENSE = DATA + "AlphaMissense_hg38.tsv.gz"
HMM_PREDICTIONS = RESULTS + "HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz"

AM_CHUNK = 5_000_000


def log(msg):
    print(msg, flush=True)


def chr_key(c):
    s = str(c).replace("chr", "")
    if s.isdigit():
        return (0, int(s))
    return (1, {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(s, 1000))


def load_genes():
    log("[1/6] reading GENCODE v44 ...")
    gene_df = pd.read_csv(
        GENE_ANNOTATION, sep="\t", comment="#", header=None,
        names=["chr", "source", "feature", "start", "end", "score",
               "strand", "frame", "attribute"],
        dtype={"start": int, "end": int},
    )
    gene_df["gene_id"] = gene_df["attribute"].str.extract(r'gene_id "(.*?)"')
    gene_df["gene_type"] = gene_df["attribute"].str.extract(r'gene_type "(.*?)"')
    gene_df["gene_name"] = gene_df["attribute"].str.extract(r'gene_name "(.*?)"')
    gene_df = gene_df.drop("attribute", axis=1)
    gene_df = gene_df[(gene_df["gene_type"] == "protein_coding")
                      & (gene_df["feature"] == "CDS")]
    gene_df["std_gene_id"] = gene_df["gene_id"].str.split(".").str[0]
    log(f"      protein-coding CDS rows: {len(gene_df):,}")
    return gene_df


def load_variants(gene_df):
    log("[2/6] reading Epi25 variant results ...")
    df = pd.read_csv(EPI25_VARIANTS, sep="\t")
    df[["chr", "pos", "ref", "alt"]] = df["variant_id"].str.split(":", expand=True)
    df["pos"] = df["pos"].astype(int)

    df = df[(df["chr"] != "chrX") & (df["chr"] != "chrY") & (df["chr"] != "chrMT")]
    df = df[(df["consequence"] != "synonymous")
            & (df["consequence"] != "non_coding")
            & (df["consequence"] != "NA")]

    df = pd.merge(
        df,
        gene_df[["std_gene_id", "gene_name"]].drop_duplicates(),
        left_on="gene_id", right_on="std_gene_id", how="left",
    ).drop("std_gene_id", axis=1)

    # AC <= 5 rare-variant cut and both-arm coverage requirement (cell 15).
    df = df[(df["ac_ctrl"] + df["ac_case"] <= 5)
            & (df["an_case"] > 0) & (df["an_ctrl"] > 0)]

    alt_case = df["ac_case"]
    alt_ctrl = df["ac_ctrl"]
    ref_case = df["an_case"] - df["ac_case"]
    ref_ctrl = df["an_ctrl"] - df["ac_ctrl"]

    # Haldane-Anscombe corrected log odds ratio and its Woolf variance.
    df["effect_size"] = np.log(
        ((0.5 + alt_case) * (0.5 + ref_ctrl))
        / ((0.5 + ref_case) * (0.5 + alt_ctrl))
    )
    df["var_effect_size"] = (1 / (0.5 + ref_case) + 1 / (0.5 + ref_ctrl)
                             + 1 / (0.5 + alt_case) + 1 / (0.5 + alt_ctrl))

    df["pLoF_ind"] = (df["consequence"] == "pLoF").astype("int32")
    df["missense_ind"] = ((df["consequence"] == "damaging_missense")
                          | (df["consequence"] == "other_missense")).astype("int32")

    log(f"      qualifying variant rows: {len(df):,}")
    return df


def merge_gerp(apply_fix, only_chroms=None):
    log(f"[3/6] annotating HMM predictions with GERP (fix={apply_fix}) ...")
    pred = pd.read_csv(HMM_PREDICTIONS, sep="\t", dtype={"chr": "string"})
    if "position" in pred.columns and "pos" not in pred.columns:
        pred = pred.rename(columns={"position": "pos"})
    pred["pos"] = pred["pos"].astype(np.int64)
    pred["chr"] = pred["chr"].astype("string")

    if only_chroms:
        pred = pred[pred["chr"].isin(only_chroms)]
        log(f"      SMOKE TEST: restricted to {sorted(only_chroms)} "
            f"-> {len(pred):,} positions")

    bw = pyBigWig.open(GERP_BIGWIG)
    bw_chroms = bw.chroms()
    chroms = sorted(pred["chr"].unique().tolist(), key=chr_key)

    merged_chunks = []
    n_dropped_bounds = 0

    for chrom in tqdm(chroms, desc="GERP merge"):
        if chrom not in bw_chroms:
            continue
        chr_len = int(bw_chroms[chrom])

        sub = pred.loc[pred["chr"] == chrom].copy()
        if sub.empty:
            continue

        if apply_fix:
            # 1-based positions, valid range [1, chr_len]
            mask = (sub["pos"] >= 1) & (sub["pos"] <= chr_len)
        else:
            mask = (sub["pos"] >= 0) & (sub["pos"] < chr_len)
        if not mask.all():
            n_dropped_bounds += int((~mask).sum())
            sub = sub.loc[mask]
        if sub.empty:
            continue

        gerp_vec = np.array(bw.values(chrom, 0, chr_len, numpy=True),
                            dtype=np.float32)
        pos_idx = sub["pos"].to_numpy(dtype=np.int64)

        # FIX 1: bigWig vector is 0-based; HMM positions are 1-based.
        gerp_vals = gerp_vec[pos_idx - 1] if apply_fix else gerp_vec[pos_idx]

        sub["GERP_RS"] = gerp_vals
        sub = sub.dropna(subset=["GERP_RS"])
        merged_chunks.append(sub)

        del gerp_vec

    bw.close()

    merged = pd.concat(merged_chunks, axis=0, ignore_index=True)
    log(f"      positions with GERP: {len(merged):,} "
        f"(dropped out-of-bounds: {n_dropped_bounds:,})")
    return merged


def load_alphamissense(variant_keys, apply_fix):
    """Load AlphaMissense restricted to variants we actually need.

    Restricting to the Epi25 key set before collapsing is exactly equivalent
    to the notebook's left-merge (non-matching AM rows never contributed) and
    turns a 71M-row groupby into a ~3M-row one.
    """
    log(f"[4/6] reading AlphaMissense (collapse={apply_fix}) ...")
    kept = []
    n_seen = 0

    reader = pd.read_csv(
        ALPHAMISSENSE, sep="\t", header=3,
        usecols=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity"],
        chunksize=AM_CHUNK,
    )
    for chunk in reader:
        n_seen += len(chunk)
        chunk = chunk.rename(columns={"#CHROM": "chr", "POS": "pos",
                                      "REF": "ref", "ALT": "alt"})
        idx = pd.MultiIndex.from_arrays(
            [chunk["chr"], chunk["pos"], chunk["ref"], chunk["alt"]])
        kept.append(chunk[idx.isin(variant_keys)])

    am = pd.concat(kept, ignore_index=True)
    log(f"      scanned {n_seen:,} rows; {len(am):,} match Epi25 variants")

    if apply_fix:
        before = len(am)
        # FIX 2: one row per variant, MAX pathogenicity across transcripts.
        am = (am.groupby(["chr", "pos", "ref", "alt"], sort=False,
                         as_index=False)["am_pathogenicity"].max())
        log(f"      collapsed {before:,} -> {len(am):,} rows "
            f"({before - len(am):,} transcript duplicates removed)")
    return am


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-suffix", default="_FIXED",
                    help="suffix for the output filename")
    ap.add_argument("--no-fix", action="store_true",
                    help="reproduce the submitted (buggy) behaviour: both fixes off")
    ap.add_argument("--gerp-fix", choices=["on", "off"], default=None,
                    help="override the GERP off-by-one fix independently")
    ap.add_argument("--am-fix", choices=["on", "off"], default=None,
                    help="override the AlphaMissense collapse independently")
    ap.add_argument("--chroms", default=None,
                    help="comma-separated chromosome subset, for smoke tests")
    args = ap.parse_args()

    # Default: both fixes on. --no-fix turns both off. The two --*-fix flags
    # override either one, so a single bug can be isolated for attribution.
    base = not args.no_fix
    fix_gerp = base if args.gerp_fix is None else (args.gerp_fix == "on")
    fix_am = base if args.am_fix is None else (args.am_fix == "on")
    log(f"configuration: GERP off-by-one fix = {fix_gerp}, "
        f"AlphaMissense collapse = {fix_am}")

    only_chroms = set(args.chroms.split(",")) if args.chroms else None

    gene_df = load_genes()
    variants_df = load_variants(gene_df)
    merged_constraint_df = merge_gerp(fix_gerp, only_chroms)

    variants_subset = variants_df[[
        "chr", "pos", "ref", "alt", "gene_id", "gene_name", "group",
        "ac_case", "an_case", "ac_ctrl", "an_ctrl",
        "pLoF_ind", "missense_ind", "effect_size", "var_effect_size"]]

    variant_keys = pd.MultiIndex.from_arrays([
        variants_subset["chr"], variants_subset["pos"],
        variants_subset["ref"], variants_subset["alt"]]).unique()
    log(f"      distinct Epi25 variants: {len(variant_keys):,}")

    am_subset = load_alphamissense(variant_keys, fix_am)

    log("[5/6] merging constraint + variants + pathogenicity ...")
    constraint_subset = merged_constraint_df[["chr", "pos", "prob_0", "GERP_RS"]]
    out = pd.merge(
        constraint_subset,
        pd.merge(variants_subset, am_subset,
                 on=["chr", "pos", "ref", "alt"], how="left"),
        on=["chr", "pos"], how="inner",
    )

    n_rows = len(out)
    n_distinct = out[["chr", "pos", "ref", "alt", "gene_id",
                      "group"]].drop_duplicates().shape[0]
    log(f"      rows: {n_rows:,}  distinct (variant,gene,group): {n_distinct:,}"
        f"  inflation: {n_rows / max(n_distinct, 1):.4f}")

    out_path = (RESULTS + "constraint_gerp_am_epi25_variants"
                + args.out_suffix + ".tsv.gz")
    log(f"[6/6] writing {out_path} ...")
    out.to_csv(out_path, index=False, compression="gzip", sep="\t")
    log(f"      done ({os.path.getsize(out_path) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
