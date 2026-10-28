#!/usr/bin/env python
"""Figure 3: HMM constraint against published constraint metrics.

  A  gnomAD v4 missense z-score per transcript
  B  missense tolerance ratio (MTR) per transcript
     x for A and B: fraction of coding bases with P(0) > 0.5, computed per CDS
     interval and averaged (unweighted) over each transcript's intervals
  C  percent difference, observed vs expected joint distribution of HMM P(0)
     (bins of 0.1) and GERP RS (bins of 1), expectation under independence
  D  chi-square contribution of each bin of the same table

Inputs
  data/    gencode.v44.basic.annotation.gtf.gz, gnomad.v4.0.constraint_metrics.tsv
  results/ HMM_rgc_ALL_RS_merged_predictions.tsv.gz (constraint + GERP, stage 02)
Output (results/)
  figure_3.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import mpl_scatter_density  # noqa: E402,F401  (registers the projection)
import seaborn as sns  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402
from wgs_constraint.epi25 import load_gene_annotation  # noqa: E402
from wgs_constraint.gene_constraint import (  # noqa: E402
    add_constraint_proportions, load_gene_constraint, ols_r_squared,
)

OUT = "figure_3.png"
COL = "proportion_over_50"
PROB_BINS = np.arange(0, 1.1, 0.1)
GERP_BINS = np.arange(-12, 7, 1)
LABEL_FS, TITLE_FS, TICK_FS, PANEL_FS = 13, 13, 10, 20
LABEL_BOX = dict(facecolor="white", edgecolor="none", pad=2)
WHITE_VIRIDIS = LinearSegmentedColormap.from_list("white_viridis", [
    (0, "#ffffff"), (1e-20, "#440053"), (0.2, "#404388"),
    (0.4, "#2a788e"), (0.6, "#21a784"), (0.8, "#78d151"), (1, "#fde624"),
], N=256)


def log(msg=""):
    print(msg, flush=True)


def per_transcript(intervals):
    """Unweighted mean of the interval fractions; z-score and MTR are
    transcript-level, repeated on every interval."""
    return (intervals.groupby("transcript", as_index=False)
            .agg(**{COL: (COL, "mean"), "mis.z_score": ("mis.z_score", "first"),
                    "MTR": ("MTR", "first")}))


def joint_distribution(table):
    """Percent difference and chi-square of the binned P(0) x GERP table
    against independence; empty cells are 0."""
    df = table.dropna(subset=["prob_0", "GERP_RS"])
    prob = pd.cut(df["prob_0"], bins=PROB_BINS, precision=2)
    gerp = pd.cut(df["GERP_RS"], bins=GERP_BINS, precision=2)
    marginal_prob = prob.value_counts(normalize=True).sort_index()
    marginal_gerp = gerp.value_counts(normalize=True).sort_index()
    expected = np.outer(marginal_prob, marginal_gerp) * len(df)

    counts = pd.DataFrame({"prob": prob, "gerp": gerp}).value_counts().sort_index()
    actual = np.array([[counts.get((p, g), 0) for g in marginal_gerp.index]
                       for p in marginal_prob.index])
    with np.errstate(divide="ignore", invalid="ignore"):
        percent_diff = np.nan_to_num((actual - expected) / expected)
        chi_sqr = np.nan_to_num((expected - actual) ** 2 / expected)
    return percent_diff, chi_sqr, list(marginal_prob.index), list(marginal_gerp.index)


def scatter_panel(fig, spec, x, y, *, ylim, title, ylabel, r2, r2_at):
    ax = fig.add_subplot(spec, projection="scatter_density")
    ax.scatter_density(x.to_numpy(float), y.to_numpy(float), cmap=WHITE_VIRIDIS)
    ax.set(xlim=(0, 0.7), ylim=ylim)
    ax.set_title(title, fontsize=TITLE_FS)
    ax.set_xlabel(r"Fraction of coding bases with $\mathbb{P}(0) > 0.5$", fontsize=LABEL_FS)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS)
    ax.text(*r2_at, rf"$R^2$ = {r2:.3f}", transform=ax.transAxes, fontsize=12,
            zorder=5, bbox=LABEL_BOX)
    return ax


def figure_3(tx, percent_diff, chi_sqr, prob_labels, gerp_labels, out_path):
    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.24)

    dA = tx.dropna(subset=[COL, "mis.z_score"])
    axA = scatter_panel(
        fig, gs[0, 0], dA[COL], dA["mis.z_score"], ylim=(-10, 10),
        title="HMM Constraint vs gnomAD v4\nMissense Z-score per transcript",
        ylabel="Missense Z-score", r2=ols_r_squared(tx, COL, "mis.z_score"),
        r2_at=(0.66, 0.06))
    dB = tx.dropna(subset=[COL, "MTR"])
    axB = scatter_panel(
        fig, gs[0, 1], dB[COL], dB["MTR"], ylim=(0.70, 1.10),
        title="HMM Constraint vs\nMissense Tolerance Ratio per transcript",
        ylabel="Missense Tolerance Ratio", r2=ols_r_squared(tx, COL, "MTR"),
        r2_at=(0.66, 0.90))

    axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    for ax, values, kw, title in [
        (axC, percent_diff, {"cmap": "coolwarm", "center": 0},
         "Percent Difference,\nObserved vs Expected Joint Distribution"),
        (axD, chi_sqr, {"cmap": "OrRd"},
         "$\\chi^2$ Statistic,\nObserved vs Expected Joint Distribution"),
    ]:
        sns.heatmap(values, xticklabels=gerp_labels, yticklabels=prob_labels,
                    cbar=True, ax=ax, cbar_kws={"pad": 0.02}, **kw)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=TITLE_FS)
        ax.set_xlabel("GERP RS Score", fontsize=LABEL_FS)
        ax.set_ylabel("HMM Constraint Probability", fontsize=LABEL_FS)
        plt.setp(ax.get_xticklabels(), rotation=90)
        plt.setp(ax.get_yticklabels(), rotation=0)

    for ax, letter in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")]:
        ax.tick_params(labelsize=TICK_FS)
        box = ax.get_position()
        fig.text(box.x0 - 0.045, box.y1 + 0.015, letter, fontsize=PANEL_FS,
                 fontweight="bold", va="bottom", ha="left")

    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    cfg = get_config()
    log(cfg.describe())

    log("loading constraint + GERP table ...")
    table = pd.read_csv(cfg.derived("constraint_gerp"), sep="\t",
                        usecols=["chr", "pos", "prob_0", "GERP_RS"])
    log(f"  {len(table):,} positions")

    log("constrained fraction per CDS interval and transcript ...")
    cds = load_gene_annotation(cfg.data("gene_annotation"))
    intervals = load_gene_constraint(
        cds[["chr", "start", "end", "gene_name", "transcript"]],
        cfg.data("gnomad_constraint"))
    intervals = add_constraint_proportions(intervals, table, thresholds=(0.5,))
    tx = per_transcript(intervals)
    log(f"  {len(intervals):,} CDS intervals, {len(tx):,} transcripts")

    log("HMM x GERP joint distribution ...")
    percent_diff, chi_sqr, prob_labels, gerp_labels = joint_distribution(table)

    out = cfg.result(OUT)
    figure_3(tx, percent_diff, chi_sqr, prob_labels, gerp_labels, out)
    log(f"R2 per transcript: z-score {ols_r_squared(tx, COL, 'mis.z_score'):.4f}, "
        f"MTR {ols_r_squared(tx, COL, 'MTR'):.4f}")
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
