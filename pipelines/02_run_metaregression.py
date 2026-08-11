#!/usr/bin/env python
"""Stage 2 -- unified meta-regression, comparison table, and reported counts.

Equivalent to `Epilepsy Analysis.ipynb` cells 18-24 plus the Table S1 export
from cell 22. No methodological changes: same clipping, mean imputation, log
transforms, WLS specification, inclusion floor and thresholds. Numbers move
only because of the corrections applied upstream in stage 1.

Emits the three counts the manuscript actually reports, so a corrected run is
directly comparable to the printed claims:

    "N significant gene-group associations (M unique genes) that the Epi25
     analysis does not, while failing to identify only one, STX1B"

For the submitted pipeline those are 28, 24 and STX1B respectively.
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import (  # noqa: E402
    GROUP_KEY, MODERATORS, fit_per_gene, get_config, prepare_regression_input,
)

DISPLAY_NAMES = {
    "gene_id": "Gene ID", "gene_name": "Gene Name", "group": "Group",
    "n_variants": "# of Variants", "p_constraint": "Constraint p-value",
    "p_gerp": "GERP p-value", "p_pathogenicity": "Pathogenicity p-value",
    "p_pLoF": "pLoF p-value", "p_missense": "Missense p-value",
    "p_unified": "Unified Model p-value",
}


def log(msg):
    print(msg, flush=True)


def prepare(path, cfg):
    log(f"[1/5] loading {Path(path).name} ...")
    df = prepare_regression_input(
        pd.read_csv(path, sep="\t"), clip_epsilon=cfg.param("clip_epsilon"))
    log(f"      regression rows: {len(df):,}")
    return df


def fit(input_df):
    log("[2/5] fitting per-gene WLS models ...")
    out = fit_per_gene(input_df)
    log(f"      models fitted: {len(out):,}")
    return out


def join_reference_sets(unified, cfg):
    log("[3/5] joining Genes4Epilepsy and Epi25 ...")
    min_variants = cfg.param("min_variants")

    g4e = pd.read_csv(cfg.data("genes4epilepsy"), sep="\t", dtype=str)
    candidates = [c for c in g4e.columns
                  if any(k in c.lower() for k in ("gene", "symbol", "hgnc"))]
    preferred = [c for c in candidates if c.lower() in (
        "gene", "gene_name", "gene symbol", "genesymbol", "hgnc",
        "hgnc_symbol", "symbol")]
    column = preferred[0] if preferred else candidates[0]
    g4e_set = set(g4e[column].dropna().astype(str).str.strip().str.upper())
    log(f"      G4E symbols: {len(g4e_set):,} (column {column!r})")

    df = unified[unified["n_variants"] >= min_variants].copy()
    log(f"      gene-groups with n_variants >= {min_variants}: {len(df):,}")
    df["p_unified"] = pd.to_numeric(df["p_unified"], errors="coerce")
    df["G4E"] = df["gene_name"].astype(str).str.strip().str.upper().isin(g4e_set)

    merged = []
    for group in df["group"].dropna().unique():
        files = sorted(glob.glob(
            str(cfg.data_dir / f"{group}_results_2023_12_30_14_36_*.csv")))
        if not files:
            log(f"      WARNING: no Epi25 file for group {group}")
            continue
        epi25 = pd.read_csv(files[0], sep=",")
        dm = next((c for c in epi25.columns
                   if "Damaging Missense" in c and "p" in c), None)
        ptv = next((c for c in epi25.columns if "PTV" in c and "p" in c), None)
        for c in (dm, ptv):
            epi25[c] = pd.to_numeric(epi25[c], errors="coerce")
        part = pd.merge(df[df["group"] == group], epi25,
                        left_on="gene_id", right_on="Gene", how="inner").copy()
        part["epi25_min_p"] = part[[dm, ptv]].min(axis=1, skipna=True)
        part = part.rename(columns={dm: "Epi25 DM p-value",
                                    ptv: "Epi25 PTV p-value"})
        merged.append(part)

    out = pd.concat(merged, ignore_index=True)
    out.replace([np.inf, -np.inf], np.nan, inplace=True)
    log(f"      merged rows: {len(out):,}")
    return out


def report(merged, cfg, label):
    thr = cfg.param("exome_wide_threshold")
    log("[4/5] significance summary")

    usig = merged["p_unified"] < thr
    esig = merged["epi25_min_p"] < thr
    only_unified = merged[usig & ~esig]
    only_epi25 = merged[esig & ~usig]

    print()
    print("=" * 72)
    print(f"REPORTED COUNTS [{label}]   threshold p < {thr:g}")
    print("=" * 72)
    print(f"unified-significant gene-group pairs        : {int(usig.sum())}")
    print(f"  ...unique genes                           : "
          f"{merged[usig]['gene_name'].nunique()}")
    print()
    print("The manuscript's Table 1 claim, recomputed:")
    print(f"  significant pairs Epi25 does not identify : {len(only_unified)}"
          f"   (submitted: 28)")
    print(f"  ...unique genes                           : "
          f"{only_unified['gene_name'].nunique()}   (submitted: 24)")
    names = sorted(only_epi25["gene_name"].dropna().unique())
    print(f"  identified by Epi25 but not unified       : {len(only_epi25)}"
          f" -> {names}   (submitted: STX1B)")
    print(f"  significant in both                       : "
          f"{int((usig & esig).sum())}")
    print("=" * 72)

    # Column name must not start with an underscore: itertuples() renames
    # those positionally, which is how an earlier version crashed here.
    genes = (merged.assign(g4e=merged["G4E"].astype(bool))
             .groupby("gene_id", as_index=False)
             .agg({"gene_name": "first", "g4e": "max",
                   "p_unified": "min", "epi25_min_p": "min"}))
    genes["meta_sig"] = genes["p_unified"] < thr
    sig = genes[genes["meta_sig"]].sort_values("p_unified")
    print(f"\nGene-level significant ({len(sig)}):")
    for r in sig.itertuples():
        print(f"  {r.gene_name:<16}{r.gene_id:<20}{r.p_unified:.3e}"
              f"  G4E={bool(r.g4e)}")
    print()
    return genes


def export_tables(merged, cfg, suffix):
    """Table S1 (suggestive set) and the full comparison table -- cell 22."""
    log("[5/5] exporting comparison tables ...")
    min_variants = cfg.param("min_variants")
    order = ["gene_id", "gene_name", "group", "n_variants",
             "p_constraint", "p_gerp", "p_pathogenicity", "p_pLoF",
             "p_missense", "p_unified"]
    table = merged[[c for c in order if c in merged.columns]
                   + ["Epi25 DM p-value", "Epi25 PTV p-value", "G4E"]].copy()
    table = table.rename(columns=DISPLAY_NAMES)
    table["Epi25 min p-value"] = table[
        ["Epi25 DM p-value", "Epi25 PTV p-value"]].min(axis=1, skipna=True)

    full = cfg.result(
        f"table_A1_epilepsy_full_p_value_comparison_n{min_variants}{suffix}.csv")
    table.to_csv(full, index=False)
    log(f"      wrote {full.name} ({len(table):,} rows)")

    sugg = table[table["Unified Model p-value"]
                 < cfg.param("suggestive_threshold")].copy()
    sugg = sugg.sort_values("Unified Model p-value")
    s1_cols = ["Gene Name", "Group", "G4E", "# of Variants",
               "Constraint p-value", "GERP p-value", "Pathogenicity p-value",
               "pLoF p-value", "Missense p-value", "Epi25 DM p-value",
               "Epi25 PTV p-value", "Unified Model p-value"]
    s1 = sugg[[c for c in s1_cols if c in sugg.columns]]
    out = cfg.result(f"Table_S1_suggestive_genes{suffix}.csv")
    s1.to_csv(out, index=False)
    log(f"      wrote {out.name} ({len(s1):,} rows; submitted Table S1 had 146)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-suffix", default="_FIXED")
    ap.add_argument("--label", default="FIXED")
    args = ap.parse_args()

    cfg = get_config()
    log(cfg.describe())
    log("")

    path = cfg.result(
        f"constraint_gerp_am_epi25_variants{args.input_suffix}.tsv.gz")
    unified = fit(prepare(path, cfg))
    unified.to_csv(
        cfg.result(f"epilepsy_unified_model_pvalues{args.input_suffix}.tsv"),
        index=False, sep="\t")

    merged = join_reference_sets(unified, cfg)
    merged.to_csv(
        cfg.result(f"unified_epi25_g4e_merged_n{cfg.param('min_variants')}"
                   f"{args.input_suffix}.tsv.gz"),
        sep="\t", index=False, compression="gzip")

    genes = report(merged, cfg, args.label)
    genes.to_csv(
        cfg.result(f"gene_level_significance{args.input_suffix}.tsv"),
        sep="\t", index=False)

    export_tables(merged, cfg, args.input_suffix)
    log("done.")


if __name__ == "__main__":
    main()
