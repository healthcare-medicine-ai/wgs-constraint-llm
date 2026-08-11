#!/usr/bin/env python
"""Stage 10 -- Figures 2a and 2b, RGC against AoU constraint predictions.

Extracted from `RGC + AoU Predictions.ipynb` cells 8 and 9.

Provenance note. `Constraint Measures Comparison.ipynb` cell 18 contains a
byte-identical copy of the plotting code and writes the same two filenames, but
in that notebook the variables it plots (`actual_joint`, `chi_sqr`) are left
over from cell 16, which builds the HMM-versus-GERP distribution on a 10x18
grid. Running that notebook top to bottom therefore overwrites Figures 2a and 2b
with mislabelled GERP data on axes captioned "RGC" and "AoU". The figures in the
manuscript come from this notebook; the copy in the other one is a latent
hazard and should be deleted.

Chromosome 2 is excluded because the HMM was trained on it.

  python pipelines/10_rgc_aou_joint.py            # gate against published
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

JOINT = "HMM_rgc_0.9_over20_chr2_joint_predictions_aou_rgc_wes.tsv.gz"
FIG_2A = "Figure 2a: observed joint distribution of RGC vs AoU WES constraint predictions"
FIG_2B = "Figure 2b: chi-square statistics of RGC vs AoU WES constraint predictions"
BINS = np.arange(0, 1.1, 0.1)
TRAINING_CHROMOSOME = "chr2"


def log(m):
    print(m, flush=True)


def joint_distribution(merged):
    """Cell 8. Bins with no observations become NaN rather than zero, which is
    what makes the published panels show gaps."""
    df = merged[merged["chr"] != TRAINING_CHROMOSOME].copy()
    df["prob_0_aou_bin"] = pd.cut(df["prob_0_aou"], bins=BINS, precision=2)
    df["prob_0_rgc_bin"] = pd.cut(df["prob_0_rgc"], bins=BINS, precision=2)

    marginal_aou = df["prob_0_aou_bin"].value_counts(normalize=True).sort_index()
    marginal_rgc = df["prob_0_rgc_bin"].value_counts(normalize=True).sort_index()

    expected = np.outer(marginal_aou, marginal_rgc) * len(df)

    counts = df[["prob_0_aou_bin", "prob_0_rgc_bin"]].value_counts().sort_index()
    actual = np.array([[counts.get((a, r), 0) for r in marginal_rgc.index]
                       for a in marginal_aou.index])

    no_data = actual == 0
    with np.errstate(divide="ignore", invalid="ignore"):
        percent_diff = (actual - expected) / expected
        chi_sqr = (expected - actual) ** 2 / expected
    percent_diff[no_data] = np.nan
    chi_sqr[no_data] = np.nan

    return {"marginal_aou": marginal_aou, "marginal_rgc": marginal_rgc,
            "expected": expected, "actual": actual,
            "percent_diff": percent_diff, "chi_sqr": chi_sqr,
            "n": len(df), "n_empty_bins": int(no_data.sum())}


def render(d, cfg, suffix):
    """Cell 9, preserved including its 11-tick labelling of a 10-bin axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LogNorm

    ticklabels = [f"{i:.1f}" for i in np.arange(0, 1.1, 0.1)]

    fig, ax1 = plt.subplots(1, figsize=(7, 6))
    sns.heatmap(d["actual"], cmap="Blues", annot=False, fmt=".3f", cbar=True,
                norm=LogNorm(), ax=ax1)
    ax1.set_title("Observed Joint Distribution", fontsize=14)
    ax1.set_xlabel("RGC Constraint Probability", fontsize=12)
    ax1.set_ylabel("AoU Constraint Probability", fontsize=12)
    ax1.invert_yaxis()
    ax1.set_xticks(range(11)); ax1.set_xticklabels(ticklabels)
    ax1.set_yticks(range(11)); ax1.set_yticklabels(ticklabels)
    plt.tight_layout()
    a = cfg.result(f"{FIG_2A}{suffix}.png")
    plt.savefig(a, dpi=300); plt.close(fig)
    log(f"      {a.name}")

    fig, ax2 = plt.subplots(1, figsize=(7, 6))
    sns.heatmap(d["chi_sqr"], cmap="Blues", annot=False, fmt=".3f", cbar=True,
                ax=ax2)
    ax2.set_title(r"$\chi^2$ for Observed vs Expected Joint Distribution",
                  y=1, x=0.55, fontsize=14)
    ax2.set_xlabel("RGC Constraint Probability", fontsize=12)
    ax2.set_ylabel("AoU Constraint Probability", fontsize=12)
    ax2.invert_yaxis()
    ax2.set_xticks(range(11)); ax2.set_xticklabels(ticklabels)
    ax2.set_yticks(range(11)); ax2.set_yticklabels(ticklabels)
    plt.tight_layout()
    b = cfg.result(f"{FIG_2B}{suffix}.png")
    plt.savefig(b, dpi=300); plt.close(fig)
    log(f"      {b.name}")
    return a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="_REGEN")
    args = ap.parse_args()
    cfg = get_config()
    log(cfg.describe() + "\n")

    source = cfg.result(JOINT)
    log(f"[1/3] loading {source.name} ...")
    merged = pd.read_csv(source, sep="\t",
                         usecols=["chr", "pos", "prob_0_aou", "prob_0_rgc"])
    log(f"      {len(merged):,} positions")

    log(f"[2/3] joint distribution ({TRAINING_CHROMOSOME} excluded) ...")
    d = joint_distribution(merged)
    log(f"      {d['n']:,} positions after excluding {TRAINING_CHROMOSOME}")
    log(f"      grid {d['actual'].shape}, {d['n_empty_bins']} empty bins")

    npz = cfg.result(f"rgc_aou_joint_distribution{args.suffix}.npz")
    np.savez_compressed(npz, expected=d["expected"], actual=d["actual"],
                        percent_diff=d["percent_diff"], chi_sqr=d["chi_sqr"],
                        bins=BINS)
    log(f"      wrote {npz.name}")

    log("[3/3] rendering ...")
    render(d, cfg, args.suffix)

    print()
    print("=" * 64)
    print("RGC vs AoU JOINT DISTRIBUTION")
    print("=" * 64)
    print(f"  positions used     : {d['n']:,}")
    print(f"  grid               : {d['actual'].shape}")
    print(f"  empty bins (NaN)   : {d['n_empty_bins']}")
    print(f"  total chi-squared  : {np.nansum(d['chi_sqr']):.6e}")
    print("=" * 64)


if __name__ == "__main__":
    main()
