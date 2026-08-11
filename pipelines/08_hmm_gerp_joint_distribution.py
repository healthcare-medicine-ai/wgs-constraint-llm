#!/usr/bin/env python
"""Stage 8 -- HMM constraint vs GERP RS joint distribution (Figures 4a and 4b).

Extracted from `Constraint Measures Comparison.ipynb` cells 16 and 17. These
become Figure 3C and 3D in Revision 3, and they are the figures affected by the
GERP coordinate defect: they are computed directly from the GERP-annotated
constraint table.

The notebook rendered these to display only, so no numeric reference was ever
saved and the published PNGs are the sole record. This stage writes the
underlying arrays alongside the figures so future comparisons are numeric
rather than visual.

  python pipelines/08_hmm_gerp_joint_distribution.py                 # corrected
  python pipelines/08_hmm_gerp_joint_distribution.py --as-published  # gate
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

PROB_BINS = np.arange(0, 1.1, 0.1)      # HMM constraint, width 0.1
GERP_BINS = np.arange(-12, 7, 1)        # GERP RS, width 1

# Published filenames, preserved so a corrected run is directly comparable.
FIG_A = "Figure 4a: percent difference distribution of GERP RS vs HMM constraint.png"
FIG_B = "Figure 4b: chi-square statistics of GERP RS vs HMM constraint.png"


def log(m):
    print(m, flush=True)


def joint_distribution(predictions):
    """Cell 16: binned joint distribution, expectation under independence,
    percent difference and chi-squared contribution."""
    df = predictions.dropna(subset=["prob_0", "GERP_RS"]).copy()
    df["prob_0_bin"] = pd.cut(df["prob_0"], bins=PROB_BINS, precision=2)
    df["gerp_bin"] = pd.cut(df["GERP_RS"], bins=GERP_BINS, precision=2)

    marginal_prob = df["prob_0_bin"].value_counts(normalize=True).sort_index()
    marginal_gerp = df["gerp_bin"].value_counts(normalize=True).sort_index()

    expected = np.outer(marginal_prob, marginal_gerp) * len(df)

    counts = df[["prob_0_bin", "gerp_bin"]].value_counts().sort_index()
    actual = np.array([[counts.get((p, g), 0) for g in marginal_gerp.index]
                       for p in marginal_prob.index])

    with np.errstate(divide="ignore", invalid="ignore"):
        percent_diff = np.nan_to_num((actual - expected) / expected)
        chi_sqr = np.nan_to_num((expected - actual) ** 2 / expected)

    return {
        "marginal_prob": marginal_prob, "marginal_gerp": marginal_gerp,
        "expected": expected, "actual": actual,
        "percent_diff": percent_diff, "chi_sqr": chi_sqr,
        "n": len(df),
    }


def render(d, cfg, suffix):
    """Cell 17: the two published heatmaps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    xt = d["marginal_gerp"].index
    yt = d["marginal_prob"].index
    out = []

    for array, cmap, center, title, fname in [
        (d["percent_diff"], "coolwarm", 0,
         "Percent Difference for Observed vs Expected Joint Distribution", FIG_A),
        (d["chi_sqr"], "OrRd", None,
         "Chi-Squared Statistic for Observed vs Expected Joint Distribution", FIG_B),
    ]:
        plt.figure(figsize=(8, 6))
        kw = {"center": center} if center is not None else {}
        sns.heatmap(array, xticklabels=xt, yticklabels=yt, cmap=cmap,
                    cbar=True, **kw)
        plt.xlabel("GERP RS Score", fontsize=12)
        plt.ylabel("HMM Constraint Probability", fontsize=12)
        plt.title(title, fontsize=14, fontweight="bold")
        plt.gca().invert_yaxis()
        plt.tight_layout()
        stem, ext = fname.rsplit(".", 1)
        path = cfg.result(f"{stem}{suffix}.{ext}")
        # bbox_inches="tight" crops to content. The published PNGs are
        # 2266x1773 and 2297x1773 -- different widths for the same figsize,
        # which only happens with this crop. Cell 17 ended in plt.show(), so
        # the saved files have no recorded provenance; this reproduces them.
        plt.savefig(path, dpi=300, bbox_inches="tight")
        plt.close()
        out.append(path)
        log(f"      {path.name}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--as-published", action="store_true",
                    help="use the pre-correction constraint+GERP table")
    args = ap.parse_args()
    cfg = get_config()

    if args.as_published:
        source = cfg.result("HMM_rgc_ALL_RS_merged_predictions_ASPUBLISHED.tsv.gz")
        suffix = "_ASPUBLISHED"
    else:
        source = cfg.result("HMM_rgc_ALL_RS_merged_predictions.tsv.gz")
        suffix = "_FIXED"

    log(f"[1/3] loading {source.name} ...")
    predictions = pd.read_csv(source, sep="\t",
                              usecols=["prob_0", "GERP_RS"])
    log(f"      {len(predictions):,} positions")

    log("[2/3] computing joint distribution ...")
    d = joint_distribution(predictions)
    log(f"      {d['n']:,} positions binned into "
        f"{d['actual'].shape[0]}x{d['actual'].shape[1]}")

    npz = cfg.result(f"hmm_gerp_joint_distribution{suffix}.npz")
    np.savez_compressed(
        npz, expected=d["expected"], actual=d["actual"],
        percent_diff=d["percent_diff"], chi_sqr=d["chi_sqr"],
        prob_bins=PROB_BINS, gerp_bins=GERP_BINS)
    log(f"      wrote {npz.name}")

    log("[3/3] rendering figures ...")
    render(d, cfg, suffix)

    print()
    print("=" * 64)
    print(f"JOINT DISTRIBUTION SUMMARY [{suffix.lstrip('_')}]")
    print("=" * 64)
    print(f"  positions                 : {d['n']:,}")
    print(f"  total chi-squared         : {d['chi_sqr'].sum():.6e}")
    print(f"  max |percent difference|  : {np.abs(d['percent_diff']).max():.6f}")
    print(f"  mean |percent difference| : {np.abs(d['percent_diff']).mean():.6f}")
    print("=" * 64)


if __name__ == "__main__":
    main()
