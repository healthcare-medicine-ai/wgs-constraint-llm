#!/usr/bin/env python
"""Stage 6 -- manuscript Tables 1 and 2, and the Figure 5 group heatmap.

Extracted from `Epilepsy Analysis.ipynb` cells 21 and 23. Those cells rendered
to notebook output rather than to files, which is why Tables 1 and 2 were
transcribed into the manuscript by hand. They are written to disk here.

Definitions follow the manuscript captions verbatim:

  Table 1  gene-group pairs reaching exome-wide significance (p < 3.4e-7) under
           EITHER the Epi25 analysis or the unified model.
  Table 2  suggestive pairs (3.4e-7 <= p < 1e-4) under the unified model that
           Genes4Epilepsy lists as epilepsy-associated.
  Figure 5 per-group 2-D histogram of unified -log10(p) against Epi25
           -log10(p), with outlying genes labelled.

Gate: run with --input-suffix _REPRO and the Table 1 counts must come to 28
pairs / 24 unique genes, with STX1B the only Epi25-only gene, matching the
submitted manuscript.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

DISPLAY = ["Gene Name", "Group", "# of Variants", "Constraint p-value",
           "GERP p-value", "Pathogenicity p-value", "pLoF p-value",
           "Missense p-value", "Epi25 DM p-value", "Epi25 PTV p-value",
           "Unified Model p-value"]


def log(m):
    print(m, flush=True)


def load_comparison(cfg, suffix):
    path = cfg.result(
        f"table_A1_epilepsy_full_p_value_comparison_"
        f"n{cfg.param('min_variants')}{suffix}.csv")
    df = pd.read_csv(path)
    for c in df.columns:
        if "p-value" in c:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    log(f"      loaded {path.name}: {len(df):,} rows")
    return df


def build_tables(df, cfg):
    cand = cfg.param("exome_wide_threshold")
    sugg = cfg.param("suggestive_threshold")

    unified_sig = df["Unified Model p-value"] < cand
    epi25_sig = df["Epi25 min p-value"] < cand

    table1 = (df[unified_sig | epi25_sig]
              .sort_values("Unified Model p-value")[DISPLAY])

    in_g4e = df["G4E"].astype(bool)
    table2 = (df[(df["Unified Model p-value"] >= cand)
                 & (df["Unified Model p-value"] < sugg) & in_g4e]
              .sort_values("Unified Model p-value")[DISPLAY])

    only_unified = df[unified_sig & ~epi25_sig]
    only_epi25 = df[epi25_sig & ~unified_sig]

    stats = {
        "table1_rows": len(table1),
        "table2_rows": len(table2),
        "unified_only_pairs": len(only_unified),
        "unified_only_genes": only_unified["Gene Name"].nunique(),
        "epi25_only": sorted(only_epi25["Gene Name"].dropna().unique()),
    }
    return table1, table2, stats


def figure5(df, cfg, out_path):
    """Cell 21: 2x2 grid of per-group hist2d, unified vs Epi25 -log10(p)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    eps, num_bins = 1e-300, 10
    d = df.dropna(subset=["Unified Model p-value", "Epi25 min p-value",
                          "Group"]).copy()
    d["mlog_unified"] = -np.log10(
        np.clip(d["Unified Model p-value"].to_numpy(float), eps, 1.0))
    d["mlog_epi25"] = -np.log10(
        np.clip(d["Epi25 min p-value"].to_numpy(float), eps, 1.0))
    d.replace([np.inf, -np.inf], np.nan, inplace=True)
    d = d.dropna(subset=["mlog_unified", "mlog_epi25"])

    groups = sorted(d["Group"].unique())
    limit = max(d["mlog_epi25"].max(), d["mlog_unified"].max()) + 1
    edges = np.linspace(0, limit, num_bins + 1)

    fig = plt.figure(figsize=(12, 10))
    L, R, B, T = 0.08, 0.88, 0.08, 0.92
    side = min((R - L) / 2, (T - B) / 2)
    gl = L + ((R - L) - 2 * side) / 2
    gb = B + ((T - B) - 2 * side) / 2
    positions = [[gl, gb + side, side, side], [gl + side, gb + side, side, side],
                 [gl, gb, side, side], [gl + side, gb, side, side]]
    axes = [fig.add_axes(p) for p in positions]
    mappable = None

    for ax, group in zip(axes, groups):
        g = d[d["Group"] == group]
        h = ax.hist2d(g["mlog_epi25"], g["mlog_unified"],
                      bins=[edges, edges], cmap="Blues", norm=LogNorm())
        mappable = h[3]
        ax.plot([0, limit], [0, limit], color="red", linestyle="--", linewidth=1)

        annot = g[(np.abs(g["mlog_epi25"] - g["mlog_unified"]) > 2)
                  & (g[["mlog_epi25", "mlog_unified"]].max(axis=1) > 6)].copy()
        centers = 0.5 * (edges[:-1] + edges[1:])
        annot["xi"] = np.digitize(annot["mlog_epi25"], edges) - 1
        annot["yi"] = np.digitize(annot["mlog_unified"], edges) - 1
        for (xi, yi), sub in annot.groupby(["xi", "yi"]):
            if not (0 <= xi < len(centers) and 0 <= yi < len(centers)):
                continue
            offset = limit * 0.02
            total = (len(sub) - 1) * offset
            for k, (_, row) in enumerate(
                    sub.sort_values("mlog_unified").iterrows()):
                ax.text(centers[xi], centers[yi] - total / 2 + k * offset,
                        str(row["Gene Name"]), fontsize=6, ha="center")

        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.text(0.97 * limit, 0.03 * limit, f"Group: {group}", fontsize=12,
                fontweight="bold", ha="right", va="bottom")

    axes[0].set_xticklabels([])
    axes[1].set_xticklabels([])
    axes[1].set_yticklabels([])
    axes[3].set_yticklabels([])

    fig.suptitle("Unified model vs Epi25 p-values across epilepsy groups",
                 fontsize=16, fontweight="bold", y=0.97)
    fig.supxlabel(r"Epi25 $-log_{10}(p)$: min(Missense, PTV) ", fontsize=14)
    fig.supylabel(r"Unified model $-log_{10}(p)$", fontsize=14)
    cax = fig.add_axes([0.91, gb, 0.02, 2 * side])
    fig.colorbar(mappable, cax=cax).set_label("Log-scaled count")

    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # The figure's provenance is its data, not its pixels: matplotlib versions
    # render differently but the underlying values must be reproducible.
    return d[["Gene Name", "Group", "mlog_unified", "mlog_epi25"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-suffix", default="_FIXED")
    args = ap.parse_args()
    cfg = get_config()
    suffix = args.input_suffix

    log("[1/3] loading comparison table ...")
    df = load_comparison(cfg, suffix)

    log("[2/3] building Tables 1 and 2 ...")
    table1, table2, stats = build_tables(df, cfg)
    t1 = cfg.result(f"Table_1_significant_gene_groups{suffix}.csv")
    t2 = cfg.result(f"Table_2_suggestive_G4E_genes{suffix}.csv")
    table1.to_csv(t1, index=False)
    table2.to_csv(t2, index=False)
    log(f"      Table 1: {stats['table1_rows']} rows -> {t1.name}")
    log(f"      Table 2: {stats['table2_rows']} rows -> {t2.name}")

    log("[3/3] rendering Figure 5 ...")
    fig_path = cfg.result(f"Figure_5_epi25_vs_unified_by_group{suffix}.png")
    fig_data = figure5(df, cfg, fig_path)
    data_path = cfg.result(f"Figure_5_underlying_data{suffix}.tsv.gz")
    fig_data.to_csv(data_path, sep="\t", index=False, compression="gzip")
    log(f"      {fig_path.name}  ({len(fig_data):,} points)")

    print()
    print("=" * 66)
    print(f"MANUSCRIPT CLAIMS [{suffix.lstrip('_')}]")
    print("=" * 66)
    print(f"  Table 1 pairs                        : {stats['table1_rows']}")
    print(f"  significant, Epi25 does not identify : "
          f"{stats['unified_only_pairs']}   (submitted: 28)")
    print(f"  ...unique genes                      : "
          f"{stats['unified_only_genes']}   (submitted: 24)")
    print(f"  Epi25-only genes                     : "
          f"{stats['epi25_only']}   (submitted: ['STX1B'])")
    print(f"  Table 2 rows (suggestive and in G4E) : {stats['table2_rows']}")
    print("=" * 66)


if __name__ == "__main__":
    main()
