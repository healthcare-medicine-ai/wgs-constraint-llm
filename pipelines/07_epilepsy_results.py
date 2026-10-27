#!/usr/bin/env python
"""Epilepsy results: Tables 1, 2 and S1, and Figure 4.

Input (results dir)
  epilepsy_gene_results.tsv  stage 03; one row per gene and epilepsy group with at
                             least 25 variants. Columns used: gene_name, group,
                             n_variants, p_constraint, p_gerp, p_pathogenicity,
                             p_unified, p_epi25_dm, p_epi25_ptv (blank where Epi25
                             has no result) and genes4epilepsy.
Outputs (results dir)
  table_1.csv   gene-groups with p < 3.4e-7 under the unified model or Epi25
  table_2.csv   gene-groups with 3.4e-7 <= p < 1e-4 under the unified model that are
                in Genes4Epilepsy
  table_s1.csv  all gene-groups with p < 1e-4 under the unified model
  figure_4.png  unified model against Epi25 -log10(p), one panel per epilepsy group

The Epi25 p-value of a gene-group is the smaller of its damaging-missense and PTV
p-values. Tables are ordered by unified model p-value.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

COLUMNS = {
    "gene_name": "Gene Name", "group": "Group", "genes4epilepsy": "G4E",
    "n_variants": "# of Variants", "p_constraint": "Constraint p-value",
    "p_gerp": "GERP p-value", "p_pathogenicity": "Pathogenicity p-value",
    "p_epi25_dm": "Epi25 DM p-value", "p_epi25_ptv": "Epi25 PTV p-value",
    "p_unified": "Unified Model p-value",
}


def table(rows, with_g4e=False):
    cols = [c for c in COLUMNS if with_g4e or c != "genes4epilepsy"]
    return (rows.sort_values("p_unified", kind="stable")[cols]
            .rename(columns=COLUMNS))


def figure_4(df, out_path):
    """2x2 grid of per-group 2-D histograms; outlying genes labelled by bin."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    eps, num_bins = 1e-300, 10
    d = df.dropna(subset=["p_unified", "p_epi25", "group"]).copy()
    d["mlog_unified"] = -np.log10(np.clip(d["p_unified"].to_numpy(float), eps, 1.0))
    d["mlog_epi25"] = -np.log10(np.clip(d["p_epi25"].to_numpy(float), eps, 1.0))

    groups = sorted(d["group"].unique())
    limit = max(d["mlog_epi25"].max(), d["mlog_unified"].max()) + 1
    edges = np.linspace(0, limit, num_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

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
        g = d[d["group"] == group]
        h = ax.hist2d(g["mlog_epi25"], g["mlog_unified"],
                      bins=[edges, edges], cmap="Blues", norm=LogNorm())
        mappable = h[3]
        ax.plot([0, limit], [0, limit], color="red", linestyle="--", linewidth=1)

        # Label genes whose two -log10(p) differ by more than 2 and where either
        # exceeds 6, at the centre of their bin, stacked within a bin.
        annot = g[(np.abs(g["mlog_epi25"] - g["mlog_unified"]) > 2)
                  & (g[["mlog_epi25", "mlog_unified"]].max(axis=1) > 6)].copy()
        annot["xi"] = np.digitize(annot["mlog_epi25"], edges) - 1
        annot["yi"] = np.digitize(annot["mlog_unified"], edges) - 1
        for (xi, yi), sub in annot.groupby(["xi", "yi"]):
            if not (0 <= xi < num_bins and 0 <= yi < num_bins):
                continue
            offset = limit * 0.02
            total = (len(sub) - 1) * offset
            for k, (_, row) in enumerate(sub.sort_values("mlog_unified").iterrows()):
                ax.text(centers[xi], centers[yi] - total / 2 + k * offset,
                        str(row["gene_name"]), fontsize=6, fontstyle="italic",
                        ha="center")

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

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    cfg = get_config()
    significant = cfg.param("exome_wide_threshold")
    suggestive = cfg.param("suggestive_threshold")

    df = pd.read_csv(cfg.result("epilepsy_gene_results.tsv"), sep="\t",
                     float_precision="round_trip")
    df["p_epi25"] = df[["p_epi25_dm", "p_epi25_ptv"]].min(axis=1)
    print(f"{len(df):,} gene-groups")

    sig = (df["p_unified"] < significant) | (df["p_epi25"] < significant)
    sugg = (df["p_unified"] >= significant) & (df["p_unified"] < suggestive)
    outputs = {
        "table_1.csv": table(df[sig]),
        "table_2.csv": table(df[sugg & df["genes4epilepsy"]]),
        "table_s1.csv": table(df[df["p_unified"] < suggestive], with_g4e=True),
    }
    for name, rows in outputs.items():
        rows.to_csv(cfg.result(name), index=False)
        print(f"{name}: {len(rows)} rows")

    figure_4(df, cfg.result("figure_4.png"))
    print("figure_4.png")


if __name__ == "__main__":
    main()
