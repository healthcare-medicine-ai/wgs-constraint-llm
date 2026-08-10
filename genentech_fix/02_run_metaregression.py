#!/usr/bin/env python
"""
Step 2: run the unified meta-regression and build the gene-level comparison
table.

Faithful port of `Epilepsy Analysis.ipynb` cells 18-24. No methodological
changes: same clipping, same mean imputation, same log transforms, same WLS
specification, same min_variants and significance thresholds. The only reason
numbers move is the two upstream corrections in 01_build_regression_input.py.

Reports the gene-level significant set so the corrected run can be compared
against the submitted one directly.
"""

import argparse
import glob

import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant
from tqdm import tqdm

DATA = "/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/"
RESULTS = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
           "osthoag/wgs-constraint-llm/results/")

GENES4EPILEPSY = DATA + "EpilepsyGenes_v2025-09.tsv"

MODERATORS = ["log_constraint", "GERP_RS", "log_pathogenicity",
              "pLoF_ind", "missense_ind"]
MIN_VARIANTS = 25
CAND_THR = 3.4e-7
EPSILON = 1e-2


def log(msg):
    print(msg, flush=True)


def prepare_input(path):
    log(f"[1/4] loading {path} ...")
    df = pd.read_csv(path, sep="\t")

    df["prob_0"] = np.clip(df["prob_0"], EPSILON, 1 - EPSILON)
    df["am_pathogenicity"] = np.clip(df["am_pathogenicity"], EPSILON, 1 - EPSILON)

    cols = ["prob_0", "GERP_RS", "am_pathogenicity"]
    df[cols] = df[cols].fillna(df[cols].mean())
    df["pLoF_ind"] = df["pLoF_ind"].fillna(0)
    df["missense_ind"] = df["missense_ind"].fillna(0)

    df[["log_constraint", "log_pathogenicity"]] = -np.log1p(
        -(df[["prob_0", "am_pathogenicity"]]))

    df = df.dropna(subset=["effect_size", "var_effect_size"])
    df = df[df["var_effect_size"] != 0]

    log(f"      regression rows: {len(df):,}")
    return df


def run_regression(input_df):
    log("[2/4] fitting per-gene WLS models ...")
    grouped = input_df.groupby(["gene_id", "gene_name", "group"])
    results = []

    for gene_key, gene_data in tqdm(grouped, desc="genes", unit="gene"):
        gene_id, gene_name, group = gene_key

        check = MODERATORS + ["effect_size", "var_effect_size"]
        if gene_data[check].isnull().any().any():
            continue

        X = add_constant(gene_data[MODERATORS])
        y = gene_data["effect_size"]
        weights = 1 / gene_data["var_effect_size"]

        try:
            model = WLS(y, X, weights=weights, missing="drop").fit()
            results.append({
                "gene_id": gene_id,
                "gene_name": gene_name,
                "group": group,
                "n_variants": len(gene_data),
                "p_constraint": model.pvalues["log_constraint"],
                "p_gerp": model.pvalues["GERP_RS"],
                "p_pathogenicity": model.pvalues["log_pathogenicity"],
                "p_pLoF": model.pvalues["pLoF_ind"],
                "p_missense": model.pvalues["missense_ind"],
                "p_unified": model.f_pvalue,
            })
        except Exception:
            pass

    out = pd.DataFrame(results)
    log(f"      models fitted: {len(out):,}")
    return out


def merge_epi25_and_g4e(unified_model_df):
    log("[3/4] merging Genes4Epilepsy and Epi25 per-group results ...")

    g4e = pd.read_csv(GENES4EPILEPSY, sep="\t", dtype=str)
    candidate = [c for c in g4e.columns
                 if any(k in c.lower() for k in ["gene", "symbol", "hgnc"])]
    preferred = [c for c in candidate if c.lower() in
                 ["gene", "gene_name", "gene symbol", "genesymbol",
                  "hgnc", "hgnc_symbol", "symbol"]]
    gene_col = preferred[0] if preferred else candidate[0]
    g4e_set = set(g4e[gene_col].dropna().astype(str).str.strip().str.upper())
    log(f"      G4E symbols: {len(g4e_set):,} (column '{gene_col}')")

    df = unified_model_df[unified_model_df["n_variants"] >= MIN_VARIANTS].copy()
    log(f"      gene-groups passing n_variants >= {MIN_VARIANTS}: {len(df):,}")

    df["p_unified"] = pd.to_numeric(df["p_unified"], errors="coerce")
    df["gene_name_upper"] = df["gene_name"].astype(str).str.strip().str.upper()
    df["G4E"] = df["gene_name_upper"].isin(g4e_set)

    merged_list = []
    for group in df["group"].dropna().unique():
        files = glob.glob(DATA + f"{group}_results_2023_12_30_14_36_*.csv")
        if not files:
            log(f"      WARNING: no Epi25 file for group {group}")
            continue
        epi25_df = pd.read_csv(files[0], sep=",")
        for c in ["Damaging Missense p‑val", "PTV p‑val"]:
            if c in epi25_df.columns:
                epi25_df[c] = pd.to_numeric(epi25_df[c], errors="coerce")

        merged = pd.merge(df[df["group"] == group], epi25_df,
                          left_on="gene_id", right_on="Gene", how="inner").copy()
        dm = next((c for c in merged.columns
                   if "Damaging Missense" in c and "p" in c), None)
        ptv = next((c for c in merged.columns
                    if "PTV" in c and "p" in c), None)
        merged["epi25_min_p"] = merged[[dm, ptv]].min(axis=1, skipna=True)
        merged_list.append(merged)

    out = pd.concat(merged_list, ignore_index=True)
    out.replace([np.inf, -np.inf], np.nan, inplace=True)
    log(f"      merged rows: {len(out):,}")
    return out


def summarise(merged_df, label):
    log("[4/4] gene-level significance summary")

    gene_df = (merged_df
               .assign(_g4e=merged_df["G4E"].astype(bool))
               .groupby("gene_id", as_index=False)
               .agg({"gene_name": "first", "_g4e": "max",
                     "p_unified": "min", "epi25_min_p": "min"}))

    gene_df["meta_sig"] = gene_df["p_unified"] < CAND_THR
    gene_df["epi25_sig"] = gene_df["epi25_min_p"] < CAND_THR

    n_meta = int(gene_df["meta_sig"].sum())
    n_epi = int(gene_df["epi25_sig"].sum())
    n_both = int((gene_df["meta_sig"] & gene_df["epi25_sig"]).sum())

    print()
    print("=" * 70)
    print(f"SIGNIFICANCE SUMMARY [{label}]  threshold p < {CAND_THR:g}")
    print("=" * 70)
    print(f"genes tested                     : {len(gene_df):,}")
    print(f"significant, unified model       : {n_meta}")
    print(f"significant, Epi25 (DM or PTV)   : {n_epi}")
    print(f"significant in both              : {n_both}")
    print(f"  ...of the unified hits, in G4E : "
          f"{int((gene_df['meta_sig'] & gene_df['_g4e']).sum())}")

    sig = gene_df[gene_df["meta_sig"]].sort_values("p_unified")
    print("\nUnified-model significant genes:")
    for _, r in sig.iterrows():
        print(f"  {r['gene_name']:<15} {r['gene_id']:<20} "
              f"p={r['p_unified']:.3e}  G4E={bool(r['_g4e'])}")
    print("=" * 70)

    return gene_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-suffix", default="_FIXED")
    ap.add_argument("--label", default="FIXED")
    args = ap.parse_args()

    in_path = (RESULTS + "constraint_gerp_am_epi25_variants"
               + args.input_suffix + ".tsv.gz")
    tag = args.input_suffix

    input_df = prepare_input(in_path)
    unified = run_regression(input_df)
    unified.to_csv(RESULTS + f"epilepsy_unified_model_pvalues{tag}.tsv",
                   index=False, sep="\t")

    merged = merge_epi25_and_g4e(unified)
    merged.to_csv(
        RESULTS + f"unified_epi25_g4e_merged_n{MIN_VARIANTS}{tag}.tsv.gz",
        sep="\t", index=False, compression="gzip")

    gene_df = summarise(merged, args.label)
    gene_df.to_csv(RESULTS + f"gene_level_significance{tag}.tsv",
                   sep="\t", index=False)
    log("\ndone.")


if __name__ == "__main__":
    main()
