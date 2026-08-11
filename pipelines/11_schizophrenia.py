#!/usr/bin/env python
"""Stage 11 -- the schizophrenia (SCHEMA) unified meta-regression.

Extracted from `Schizophrenia Analysis.ipynb` cells 15, 17-21. Cells 17, 18 and
19 were character-for-character identical to the epilepsy notebook's, so the
model itself comes from `wgs_constraint.metareg`; only what is genuinely
cohort-specific lives here.

What differs from the epilepsy analysis:

  * SCHEMA supplies a VEP-style consequence vocabulary, so pLoF and missense
    indicators are list membership rather than a direct comparison, and a set
    of consequences is excluded outright.
  * The reference set is SCHEMA's own gene results (`P meta`) rather than the
    Epi25 per-group files.
  * There is no Genes4Epilepsy equivalent.

Everything else -- the AC <= 5 cut, Haldane-Anscombe effect sizes, the five
moderators, the inclusion floor, the whole-model F-test -- is shared.

  python pipelines/11_schizophrenia.py                # corrected
  python pipelines/11_schizophrenia.py --no-fix       # reproduce as published
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import (  # noqa: E402
    annotate_with_gerp, fit_per_gene, get_config, haldane_effect_sizes,
    load_collapsed, prepare_regression_input,
)
from wgs_constraint.alphamissense import KEY_COLUMNS, load_uncollapsed  # noqa: E402
from wgs_constraint.epi25 import load_gene_annotation  # noqa: E402

# SCHEMA consequence vocabulary.
#
# Cell 15 of the notebook lists ten exclusions, six of them commented out. Taking
# the file at face value -- only the four uncommented entries -- regenerates
# 178,821 rows the September 2025 cached merge does not contain, every one of
# them neither pLoF nor missense.
#
# Those six were recovered empirically: joining the cached merge back to the
# SCHEMA source shows intron, splice_region, 3'/5' UTR and upstream/downstream
# variants absent from it entirely, while every coding consequence is fully
# present (see genentech_fix/diag_scz_consequences.py). So the submitted analysis
# ran with all ten active, and the comment characters were added afterwards. The
# notebook's current state does not describe the published results.
#
# Restoring all ten reproduces the cached merge exactly. This is the analysis
# definition, not one of the two corrections, so it applies on both paths.
#
# Only these four actually occur in SCHEMA: coding_sequence_variant,
# synonymous_variant, intron_variant, splice_region_variant, 3_prime_UTR_variant,
# 5_prime_UTR_variant, upstream_gene_variant, downstream_gene_variant. The rest
# are inert and kept only to mirror the notebook.
EXCLUDED_CONSEQUENCES = [
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "coding_sequence_variant",
    "downstream_gene_variant",
    "intron_variant",
    "mature_miRNA_variant",
    "null",
    "splice_region_variant",
    "synonymous_variant",
    "upstream_gene_variant",
]
PLOF_CONSEQUENCES = [
    "stop_gained", "splice_acceptor_variant", "splice_donor_variant",
    "frameshift_variant",
]
MISSENSE_CONSEQUENCES = [
    "inframe_insertion", "inframe_deletion", "stop_retained_variant",
    "stop_lost", "missense_variant_mpc_<2", "protein_altering_variant",
    "missense_variant_mpc_2-3", "missense_variant_mpc_>=3",
]

VARIANT_COLUMNS = [
    "chr", "pos", "ref", "alt", "gene_id", "gene_name", "group",
    "ac_case", "an_case", "ac_ctrl", "an_ctrl",
    "pLoF_ind", "missense_ind", "effect_size", "var_effect_size",
]

FIGURE_NAME = "Figure A1: p-value comparison for schizophrenia"


def log(m):
    print(m, flush=True)


def load_variants(cfg, genes):
    path = cfg.data_dir / "SCHEMA_variant_results_hg38.tsv.gz"
    df = pd.read_csv(path, sep="\t", compression="gzip")

    # Chained != rather than .isin(), matching cell 15 literally. Kept for
    # fidelity, but be clear about what it did NOT do: switching to this form
    # left effect_size bit-identical to the cached merge on exactly 80.38% of
    # rows, the same figure as with .isin(). So the remaining difference is not
    # this filter, and the memory-layout explanation that worked for epi25.py
    # does not transfer here. Cause still unidentified -- see
    # docs/OVERHAUL_PLAN.md. The consequence filter below stays .isin() because
    # that is what cell 15 uses.
    chrom_mask = pd.Series(True, index=df.index)
    for excluded in cfg.param("excluded_chromosomes"):
        chrom_mask &= (df["chr"] != excluded)
    df = df[chrom_mask]
    df = df[~df["consequence"].isin(EXCLUDED_CONSEQUENCES)]
    df = df[~df["consequence"].isna()]

    df = pd.merge(df, genes[["std_gene_id", "gene_name"]].drop_duplicates(),
                  left_on="gene_id", right_on="std_gene_id",
                  how="left").drop("std_gene_id", axis=1)

    df = haldane_effect_sizes(df, max_allele_count=cfg.param("max_allele_count"))
    df["pLoF_ind"] = df["consequence"].isin(PLOF_CONSEQUENCES).astype("int32")
    df["missense_ind"] = df["consequence"].isin(
        MISSENSE_CONSEQUENCES).astype("int32")
    return df[VARIANT_COLUMNS]


def figure_a1(merged, cfg, suffix):
    """Cell 20: unified vs SCHEMA p-values, one panel per group."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    d = merged.copy()
    # Coerce before taking logs: p_unified arrives as object dtype whenever any
    # model returned a non-float, and numpy cannot vectorise log10 over objects.
    for col in ("p_unified", "P meta"):
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["minus_log_pval"] = -np.log10(d["p_unified"])
    d["minus_log_pub_pval"] = -np.log10(d["P meta"])
    d.replace([np.inf, -np.inf], np.nan, inplace=True)

    written = []
    for group, gdf in d.groupby("group"):
        clean = gdf.dropna(subset=["minus_log_pub_pval", "minus_log_pval"])
        if clean.empty:
            continue
        x0, x1 = clean["minus_log_pub_pval"].min(), clean["minus_log_pub_pval"].max()
        y0, y1 = clean["minus_log_pval"].min(), clean["minus_log_pval"].max()
        nb = 10
        xe = np.linspace(x0, x1 + 1, nb + 1)
        ye = np.linspace(y0, y1 + 1, nb + 1)

        plt.figure(figsize=(8, 6))
        plt.hist2d(clean["minus_log_pub_pval"], clean["minus_log_pval"],
                   bins=[xe, ye], cmap="Blues", norm=LogNorm())
        top = max(x1, y1) + 1
        plt.plot([0, top], [0, top], color="red", linestyle="--")

        annot = clean[
            (np.abs(clean["minus_log_pub_pval"] - clean["minus_log_pval"]) > 2)
            & (clean[["minus_log_pub_pval", "minus_log_pval"]].max(axis=1) > 4)
        ].copy()
        xc = 0.5 * (xe[:-1] + xe[1:])
        yc = 0.5 * (ye[:-1] + ye[1:])
        annot["xi"] = np.digitize(annot["minus_log_pub_pval"], xe) - 1
        annot["yi"] = np.digitize(annot["minus_log_pval"], ye) - 1
        offset = (y1 + 1 - y0) * 0.02
        for (xi, yi), sub in annot.groupby(["xi", "yi"]):
            if not (0 <= xi < len(xc) and 0 <= yi < len(yc)):
                continue
            total = (len(sub) - 1) * offset
            for k, (_, row) in enumerate(
                    sub.sort_values("minus_log_pval").iterrows()):
                plt.text(xc[xi], yc[yi] - total / 2 + k * offset,
                         row["gene_name"], fontsize=8, ha="center")

        plt.xlim(x0, x1 + 1)
        plt.ylim(y0, y1 + 1)
        plt.xlabel(r"SCHEMA $-log_{10}(p)$", fontsize=12)
        plt.ylabel(r"Unified Model $-log_{10}(p)$", fontsize=12)
        plt.title("Unified Model vs SCHEMA p-values for schizophrenia",
                  fontsize=14, fontweight="bold")
        plt.colorbar(label="Log-Scaled Count")
        out = cfg.result(f"{FIGURE_NAME}{suffix}.png")
        plt.savefig(out)
        plt.close()
        written.append(out)
        log(f"      {out.name}")
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fix", action="store_true")
    ap.add_argument("--chroms", default=None)
    args = ap.parse_args()

    fix = not args.no_fix
    suffix = "_FIXED" if fix else "_REPRO"
    cfg = get_config()
    min_variants = cfg.param("min_variants")
    log(cfg.describe())
    log(f"corrections applied: {fix}\n")

    log("[1/6] GENCODE annotation ...")
    genes = load_gene_annotation(cfg.data("gene_annotation"))

    log("[2/6] SCHEMA variants ...")
    variants = load_variants(cfg, genes)
    log(f"      qualifying variant rows: {len(variants):,}")

    log("[3/6] GERP annotation ...")
    pred = pd.read_csv(cfg.derived("hmm_predictions"), sep="\t",
                       dtype={"chr": "string"})
    if "position" in pred.columns and "pos" not in pred.columns:
        pred = pred.rename(columns={"position": "pos"})
    pred["pos"] = pred["pos"].astype("int64")
    pred["chr"] = pred["chr"].astype("string")
    constraint = annotate_with_gerp(
        pred, cfg.data("gerp_bigwig"), position_base=1 if fix else 0,
        only_chromosomes=args.chroms.split(",") if args.chroms else None)
    log(f"      positions with GERP: {len(constraint):,}")

    log("[4/6] AlphaMissense ...")
    wanted = pd.MultiIndex.from_arrays(
        [variants[c] for c in KEY_COLUMNS]).unique()
    if fix:
        am = load_collapsed(cfg.data("alphamissense"), restrict_to=wanted,
                            how=cfg.param("alphamissense_collapse"))
    else:
        am = load_uncollapsed(cfg.data("alphamissense"), restrict_to=wanted)

    merged = pd.merge(
        constraint[["chr", "pos", "prob_0", "GERP_RS"]],
        pd.merge(variants, am, on=KEY_COLUMNS, how="left"),
        on=["chr", "pos"], how="inner")
    out_variants = cfg.result(f"constraint_gerp_am_scz_variants{suffix}.tsv.gz")
    merged.to_csv(out_variants, index=False, compression="gzip", sep="\t")
    log(f"      merged rows: {len(merged):,} -> {out_variants.name}")

    log("[5/6] meta-regression ...")
    unified = fit_per_gene(prepare_regression_input(
        merged, clip_epsilon=cfg.param("clip_epsilon")))
    log(f"      models fitted: {len(unified):,}")
    unified.to_csv(
        cfg.result(f"schizophrenia_variants_unified_gerp_model_pvalues{suffix}.tsv"),
        index=False, sep="\t")

    log("[6/6] SCHEMA comparison, figure and table ...")
    unified = unified[unified["n_variants"] >= min_variants]
    schema = pd.read_csv(cfg.data_dir / "SCHEMA_gene_results.tsv.bgz",
                         sep="\t", compression="gzip")
    merged_pub = pd.merge(unified, schema, on=["gene_id", "group"])
    log(f"      gene-groups after merge: {len(merged_pub):,}")

    figure_a1(merged_pub, cfg, suffix)

    comparison = merged_pub[[
        "gene_id", "gene_name", "group", "n_variants", "p_constraint",
        "p_gerp", "p_pathogenicity", "p_pLoF", "p_missense", "p_unified",
        "P meta"]].rename(columns={
            "gene_id": "Gene ID", "gene_name": "Gene Name", "group": "Group",
            "n_variants": "# of Variants",
            "p_constraint": "Constraint p-value", "p_gerp": "GERP p-value",
            "p_pathogenicity": "Pathogenicity p-value",
            "p_pLoF": "pLoF p-value", "p_missense": "Missense p-value",
            "p_unified": "Unified Model p-value", "P meta": "SCHEMA p-value"})
    table = cfg.result(
        f"table_A2_schizophrenia_full_p_value_comparison_n{min_variants}{suffix}.csv")
    comparison.to_csv(table, index=False)
    log(f"      {table.name} ({len(comparison):,} rows)")

    thr = cfg.param("exome_wide_threshold")
    sig = comparison[comparison["Unified Model p-value"] < thr]
    print()
    print("=" * 64)
    print(f"SCHIZOPHRENIA SUMMARY [{suffix.lstrip('_')}]  p < {thr:g}")
    print("=" * 64)
    print(f"  gene-groups tested         : {len(comparison):,}")
    print(f"  unified-model significant  : {len(sig)}")
    # Column names contain spaces, so itertuples renames them positionally.
    # Select explicitly instead of relying on tuple position.
    for _, r in sig.sort_values("Unified Model p-value")[
            ["Gene Name", "Unified Model p-value", "SCHEMA p-value"]].iterrows():
        print(f"    {str(r['Gene Name']):<14} unified={r['Unified Model p-value']:.3e}"
              f"  SCHEMA={r['SCHEMA p-value']:.3e}")
    print("=" * 64)


if __name__ == "__main__":
    main()
