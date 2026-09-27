#!/usr/bin/env python
"""Stage 15 -- merged four-panel Figure 3 for Revision 3.

The editor asked for at most eight display items. This merges the old Figure 3
(HMM vs missense z-score, HMM vs MTR) and the old Figure 4 (HMM vs GERP percent
difference and chi-squared) into one figure, and makes panel B what the
Revision 2 response letter said it was:

  A  fraction of bases with P(0) > 0.5 vs gnomAD v4 missense z, per CDS interval
     (as submitted; --panel-a-per-transcript aggregates it like panel B)
  B  fraction of bases with P(0) > 0.5 vs MTR, per transcript: unweighted mean
     of the per-interval fractions, which reproduces the letter's R2 of 0.300
     (0.3058 measured; see stage 09's docstring)
  C  percent difference, observed vs expected HMM x GERP joint distribution
  D  chi-squared contribution for the same table

C and D come from stage 08's corrected output. The per-interval fractions come
from stage 09. They depend only on prob_0, which the GERP/AlphaMissense
corrections do not touch, so the ASPUBLISHED table is the right input for both
the published and the corrected figure.

Ported from `Submissions/4 - Revision 3 (2026-08)/Reviewer 1 analysis/
R3.1 - Merged Figure 3 (four panels).py`.

  python pipelines/15_figure3_merged.py --panel-a-per-transcript   # Revision 3 figure

The Revision 3 figure uses --panel-a-per-transcript, so both scatter panels are
per transcript (R2 = 0.514 and 0.306).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402
from wgs_constraint.gene_constraint import ols_r_squared  # noqa: E402

OUT = "Figure 3 (merged four panel)_FIXED.png"
COL = "proportion_over_50"
LABEL_FS, TITLE_FS, TICK_FS, PANEL_FS = 13, 13, 10, 20
# White box behind the R2 labels so no density point can sit on them.
LABEL_BOX = dict(facecolor="white", edgecolor="none", pad=2)


def log(m=""):
    print(m, flush=True)


def white_viridis():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("white_viridis", [
        (0, "#ffffff"), (1e-20, "#440053"), (0.2, "#404388"),
        (0.4, "#2a788e"), (0.6, "#21a784"), (0.8, "#78d151"), (1, "#fde624"),
    ], N=256)


def per_transcript(intervals):
    """One row per transcript: unweighted mean of the interval fractions.

    mis.z_score and MTR are transcript-level values repeated on every interval,
    so taking the first is exact.
    """
    return (intervals.groupby("transcript", as_index=False)
            .agg(gene_name=("gene_name", "first"), n_intervals=(COL, "size"),
                 **{COL: (COL, "mean")},
                 **{"mis.z_score": ("mis.z_score", "first"), "MTR": ("MTR", "first")}))


def bin_labels(bins):
    """The tick labels stage 08 used: pd.cut categories at precision 2."""
    return list(pd.cut(pd.Series([], dtype=float), bins=bins, precision=2).cat.categories)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-a-per-transcript", action="store_true",
                    help="aggregate panel A to transcripts as well as panel B")
    args = ap.parse_args()
    cfg = get_config()

    log("[1/3] loading inputs ...")
    iv = pd.read_csv(cfg.result("gene_constraint_proportions_ASPUBLISHED.tsv.gz"), sep="\t")
    tx = per_transcript(iv)
    # A transcript's MTR and z-score must be one value; check rather than assume.
    spread = iv.groupby("transcript")[["mis.z_score", "MTR"]].nunique().max()
    log(f"      {len(iv):,} CDS intervals, {len(tx):,} transcripts; "
        f"max distinct values per transcript: {spread.to_dict()}")
    jd = np.load(cfg.result("hmm_gerp_joint_distribution_FIXED.npz"))
    prob_lab, gerp_lab = bin_labels(jd["prob_bins"]), bin_labels(jd["gerp_bins"])

    log("[2/3] R2 at each granularity (P(0) > 0.5):")
    r2 = {(g, y): ols_r_squared(d, COL, y)
          for g, d in (("interval", iv), ("transcript", tx))
          for y in ("mis.z_score", "MTR")}
    for (g, y), v in r2.items():
        log(f"      {y:12s} per {g:10s} {v:.4f}")

    a_gran = "transcript" if args.panel_a_per_transcript else "interval"
    dA = tx if args.panel_a_per_transcript else iv
    dA = dA.dropna(subset=[COL, "mis.z_score"])
    dB = tx.dropna(subset=[COL, "MTR"])

    log("[3/3] rendering ...")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    import mpl_scatter_density  # noqa: F401  (registers the projection)
    import seaborn as sns

    fig = plt.figure(figsize=(14, 12))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.24)
    xlab = r"Fraction of coding bases with $\mathbb{P}(0) > 0.5$"
    unit = {"interval": "coding interval", "transcript": "transcript"}

    axA = fig.add_subplot(gs[0, 0], projection="scatter_density")
    axA.scatter_density(dA[COL].to_numpy(float), dA["mis.z_score"].to_numpy(float),
                        cmap=white_viridis())
    axA.set(xlim=(0, 0.7), ylim=(-10, 10))
    axA.set_title(f"HMM Constraint vs gnomAD v4\nMissense Z-score per {unit[a_gran]}",
                  fontsize=TITLE_FS)
    axA.set_xlabel(xlab, fontsize=LABEL_FS)
    axA.set_ylabel("Missense Z-score", fontsize=LABEL_FS)
    axA.text(0.66, 0.06, rf"$R^2$ = {r2[(a_gran, 'mis.z_score')]:.3f}",
             transform=axA.transAxes, fontsize=12, zorder=5, bbox=LABEL_BOX)

    axB = fig.add_subplot(gs[0, 1], projection="scatter_density")
    axB.scatter_density(dB[COL].to_numpy(float), dB["MTR"].to_numpy(float),
                        cmap=white_viridis())
    axB.set(xlim=(0, 0.7), ylim=(0.70, 1.10))
    axB.set_title("HMM Constraint vs\nMissense Tolerance Ratio per transcript",
                  fontsize=TITLE_FS)
    axB.set_xlabel(xlab, fontsize=LABEL_FS)
    axB.set_ylabel("Missense Tolerance Ratio", fontsize=LABEL_FS)
    axB.text(0.66, 0.90, rf"$R^2$ = {r2[('transcript', 'MTR')]:.3f}",
             transform=axB.transAxes, fontsize=12, zorder=5, bbox=LABEL_BOX)

    axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    for ax, arr, kw, title in [
        (axC, jd["percent_diff"], {"cmap": "coolwarm", "center": 0},
         "Percent Difference,\nObserved vs Expected Joint Distribution"),
        (axD, jd["chi_sqr"], {"cmap": "OrRd"},
         "$\\chi^2$ Statistic,\nObserved vs Expected Joint Distribution"),
    ]:
        sns.heatmap(arr, xticklabels=gerp_lab, yticklabels=prob_lab, cbar=True,
                    ax=ax, cbar_kws={"pad": 0.02}, **kw)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=TITLE_FS)
        ax.set_xlabel("GERP RS Score", fontsize=LABEL_FS)
        ax.set_ylabel("HMM Constraint Probability", fontsize=LABEL_FS)
        plt.setp(ax.get_xticklabels(), rotation=90)
        plt.setp(ax.get_yticklabels(), rotation=0)

    for ax, letter in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")]:
        ax.tick_params(labelsize=TICK_FS)
        b = ax.get_position()
        fig.text(b.x0 - 0.045, b.y1 + 0.015, letter, fontsize=PANEL_FS,
                 fontweight="bold", va="bottom", ha="left")

    out = cfg.result(OUT)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"      wrote {out.name}")


if __name__ == "__main__":
    main()
