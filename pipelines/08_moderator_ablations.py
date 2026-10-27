#!/usr/bin/env python
"""Moderator ablations: the unified model refitted with one moderator (or the
pLoF/missense pair) removed, and Figures S1-S4.

Input (results dir)
  epilepsy_regression_input.tsv.gz   stage 02
Outputs (results dir)
  ablation_minus_constraint.tsv      without HMM constraint               (S1)
  ablation_minus_gerp.tsv            without GERP RS                      (S2)
  ablation_minus_pathogenicity.tsv   without AlphaMissense pathogenicity  (S3)
  ablation_minus_plof_missense.tsv   without pLoF and missense            (S4)
  figure_s1.png .. figure_s4.png     full model against each ablation, -log10(p),
                                     one panel per epilepsy group

Figures use gene-groups with at least 25 variants in both models. A gene is
labelled when either model gives p < 1e-6 and the two -log10(p) differ by at
least 1; labels sit at the centre of the point's bin, stacked within a bin.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import (  # noqa: E402
    MODERATORS, fit_per_gene, get_config, prepare_regression_input,
)

KEY_COLS = ["gene_id", "gene_name", "group"]
NUM_BINS = 10
P_THRESHOLD = 1e-6
DELTA_THRESHOLD = 1
EPS = 1e-300

# (output stem, moderators dropped, title, x-axis label)
ABLATIONS = [
    ("minus_constraint", ["log_constraint"],
     "Unified Model p-values: with vs without HMM Constraint",
     r"$-log_{10}(p)$ without HMM constraint"),
    ("minus_gerp", ["GERP_RS"],
     "Unified Model p-values: with vs without GERP RS score",
     r"$-log_{10}(p)$ without GERP RS"),
    ("minus_pathogenicity", ["log_pathogenicity"],
     "Unified Model p-values: with vs without AM Pathogenicity",
     r"$-log_{10}(p)$ without AM pathogenicity"),
    ("minus_plof_missense", ["pLoF_ind", "missense_ind"],
     "Unified Model p-values: with vs without pLoF/Missense Annotations",
     r"$-log_{10}(p)$ without pLoF/missense"),
]


def compare(full, other, min_variants):
    """Join full and ablated p-values, keep n >= min_variants in both, take -log10."""
    cols = KEY_COLS + ["n_variants", "p_unified"]
    cmp = (full[cols].rename(columns={"n_variants": "n_with", "p_unified": "p_with"})
           .merge(other[cols].rename(columns={"n_variants": "n_without",
                                              "p_unified": "p_without"}),
                  on=KEY_COLS, how="inner"))
    # p_unified can arrive as an object column of 0-d arrays; make it float.
    for c in ("p_with", "p_without"):
        cmp[c] = cmp[c].to_numpy(dtype=float)
    cmp = cmp[(cmp["n_with"] >= min_variants) & (cmp["n_without"] >= min_variants)]
    cmp = cmp.dropna(subset=["p_with", "p_without"]).reset_index(drop=True)
    cmp["log10_with"] = -np.log10(np.clip(cmp["p_with"], EPS, 1.0))
    cmp["log10_without"] = -np.log10(np.clip(cmp["p_without"], EPS, 1.0))
    cmp["delta"] = cmp["log10_with"] - cmp["log10_without"]
    return cmp


def stack_labels(fig, labels, edges):
    """Move each panel's labels to their bin centre, stacked one label height apart."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    centres = 0.5 * (edges[:-1] + edges[1:])
    limit = edges[-1]
    for ax, texts in labels:
        if not texts:
            continue
        tallest = max(t.get_window_extent(renderer).height for t in texts)
        step = 1.10 * tallest * limit / ax.bbox.height
        cells = {}
        for t in texts:
            x, y = t.get_position()
            xi = int(np.clip(np.digitize(x, edges) - 1, 0, NUM_BINS - 1))
            yi = int(np.clip(np.digitize(y, edges) - 1, 0, NUM_BINS - 1))
            cells.setdefault((xi, yi), []).append(t)
        for (xi, yi), members in cells.items():
            members.sort(key=lambda t: (t.get_position()[1], t.get_text()))
            total = (len(members) - 1) * step
            for k, t in enumerate(members):
                t.set_position((centres[xi], centres[yi] - total / 2 + k * step))


def render(cmp, title, x_label, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    groups = sorted(cmp["group"].dropna().unique())
    limit = max(cmp["log10_with"].max(), cmp["log10_without"].max()) + 1
    edges = np.linspace(0, limit, NUM_BINS + 1)

    fig = plt.figure(figsize=(12, 10))
    left, right, bottom, top = 0.08, 0.88, 0.08, 0.92
    side = min((right - left) / 2, (top - bottom) / 2)
    grid_l = left + ((right - left) - 2 * side) / 2
    grid_b = bottom + ((top - bottom) - 2 * side) / 2
    positions = [[grid_l, grid_b + side, side, side],
                 [grid_l + side, grid_b + side, side, side],
                 [grid_l, grid_b, side, side],
                 [grid_l + side, grid_b, side, side]]
    axes = [fig.add_axes(p) for p in positions]
    mappable = None
    labels = []

    for ax, group in zip(axes, groups):
        df = cmp[cmp["group"] == group]
        hist = ax.hist2d(df["log10_without"], df["log10_with"],
                         bins=[edges, edges], cmap="Blues", norm=LogNorm())
        mappable = hist[3]
        ax.plot([0, limit], [0, limit], linestyle="--", color="red",
                linewidth=1.0, alpha=0.6)

        annot = df[((df["p_with"] < P_THRESHOLD) | (df["p_without"] < P_THRESHOLD))
                   & (np.abs(df["delta"]) >= DELTA_THRESHOLD)]
        labels.append((ax, [
            ax.text(row["log10_without"], row["log10_with"], row["gene_name"],
                    fontsize=6, ha="center", va="center")
            for _, row in annot.iterrows()]))

        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.text(0.97 * limit, 0.03 * limit, f"Group: {group}", fontsize=12,
                fontweight="bold", ha="right", va="bottom")

    axes[0].set_xticklabels([])
    axes[1].set_xticklabels([])
    axes[1].set_yticklabels([])
    axes[3].set_yticklabels([])

    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.97)
    fig.supxlabel(x_label, fontsize=14)
    fig.supylabel(r"$-log_{10}(p)$ full model", fontsize=14)
    cbar_ax = fig.add_axes([0.91, grid_b, 0.02, 2 * side])
    fig.colorbar(mappable, cax=cbar_ax).set_label("Log-scaled count")

    stack_labels(fig, labels, edges)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    cfg = get_config()
    min_variants = cfg.param("min_variants")

    prepared = prepare_regression_input(
        pd.read_csv(cfg.result("epilepsy_regression_input.tsv.gz"), sep="\t"),
        clip_epsilon=cfg.param("clip_epsilon"))
    print(f"{len(prepared):,} variants")

    full = fit_per_gene(prepared)
    for k, (stem, dropped, title, x_label) in enumerate(ABLATIONS, start=1):
        fit = fit_per_gene(prepared, moderators=[m for m in MODERATORS
                                                 if m not in dropped])
        fit.to_csv(cfg.result(f"ablation_{stem}.tsv"), sep="\t", index=False)
        render(compare(full, fit, min_variants), title, x_label,
               cfg.result(f"figure_s{k}.png"))
        print(f"ablation_{stem}.tsv, figure_s{k}.png")


if __name__ == "__main__":
    main()
