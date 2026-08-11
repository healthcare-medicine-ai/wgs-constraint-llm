#!/usr/bin/env python
"""Stage 9 -- Figures 3a and 3b, HMM constraint against published metrics.

Extracted from `Constraint Measures Comparison.ipynb`. Which cells produced the
published panels was ambiguous -- the notebook writes several competing versions
of "Figure 3a" and "Figure 3b" and the last one executed wins. Resolved by
comparing what the AJHG-R2-submitted tag kept against what it deleted:

    Figure 3a  cell 32  fraction of bases with P(0) > 0.5  vs gnomAD missense z
    Figure 3b  cell 28  proportion with P(0) > 0.6         vs MTR

Two quirks of the originals are preserved so the gate can pass, and flagged so
they can be decided on deliberately:

  * The panels use different thresholds, 0.5 and 0.6. Pass --match-thresholds
    to put both on 0.5, as the Revision 3 plan proposes.
  * Cell 28 annotated panel B with a hardcoded "$R^2=0.152$" rather than the
    computed value. This stage always computes it, and reports the difference.

  python pipelines/09_gene_constraint_figures.py --as-published   # gate
  python pipelines/09_gene_constraint_figures.py                  # corrected
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402
from wgs_constraint.gene_constraint import (  # noqa: E402
    add_constraint_proportions, load_gene_constraint, ols_r_squared,
)

FIG_3A = "Figure 3a: HMM vs gnomAD missense z-score per gene"
FIG_3B = "Figure 3b: HMM vs MTR per gene"

# The literal the notebook printed on panel B.
HARDCODED_R2_PANEL_B = 0.152


def log(m):
    print(m, flush=True)


def white_viridis():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("white_viridis", [
        (0, "#ffffff"), (1e-20, "#440053"), (0.2, "#404388"),
        (0.4, "#2a788e"), (0.6, "#21a784"), (0.8, "#78d151"), (1, "#fde624"),
    ], N=256)


def scatter_panel(df, xcol, ycol, *, xlabel, ylabel, title, xlim, ylim,
                  r2, r2_at, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mpl_scatter_density  # noqa: F401  (registers the projection)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection="scatter_density")
    ax.scatter_density(df[xcol], df[ycol], cmap=white_viridis())
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.text(r2_at[0], r2_at[1], f"$R^2={round(r2, 3)}$", fontsize=10)
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    log(f"      {Path(out_path).name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--as-published", action="store_true",
                    help="use the pre-correction constraint+GERP predictions")
    ap.add_argument("--match-thresholds", action="store_true",
                    help="use P(0) > 0.5 for both panels instead of 0.5 and 0.6")
    args = ap.parse_args()

    cfg = get_config()
    if args.as_published:
        source = cfg.result("HMM_rgc_ALL_RS_merged_predictions_ASPUBLISHED.tsv.gz")
        suffix = "_ASPUBLISHED"
    else:
        source = cfg.result("HMM_rgc_ALL_RS_merged_predictions.tsv.gz")
        suffix = "_FIXED"
    if args.match_thresholds:
        suffix += "_MATCHED"

    log(cfg.describe())
    log(f"source: {source.name}\n")

    log("[1/4] GENCODE CDS intervals and gnomAD constraint ...")
    # The shared loader in wgs_constraint.epi25 does not carry a transcript
    # column, which this join needs, so the GTF is parsed here instead.
    gtf = pd.read_csv(cfg.data("gene_annotation"), sep="\t", comment="#",
                      header=None,
                      names=["chr", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attribute"],
                      dtype={"start": int, "end": int})
    gtf["gene_type"] = gtf["attribute"].str.extract(r'gene_type "(.*?)"')
    gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "(.*?)"')
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "(.*?)"')
    gtf["transcript"] = gtf["transcript_id"].str.split(".").str[0]
    gtf = gtf[(gtf["gene_type"] == "protein_coding") & (gtf["feature"] == "CDS")]
    gtf = gtf.drop(columns=["attribute"])
    log(f"      protein-coding CDS intervals: {len(gtf):,}")

    gene_constraint = load_gene_constraint(gtf, cfg.data_dir /
                                           "gnomad.v4.0.constraint_metrics.tsv")
    log(f"      joined to gnomAD constraint: {len(gene_constraint):,} rows")

    log("[2/4] loading predictions ...")
    predictions = pd.read_csv(source, sep="\t", usecols=["chr", "pos", "prob_0"])
    log(f"      {len(predictions):,} positions")

    log("[3/4] computing constrained fractions per CDS interval ...")
    thresholds = (0.5,) if args.match_thresholds else (0.5, 0.6)
    gene_constraint = add_constraint_proportions(
        gene_constraint, predictions, thresholds=thresholds)
    col_a = "proportion_over_50"
    col_b = "proportion_over_50" if args.match_thresholds else "proportion_over_60"

    r2_a = ols_r_squared(gene_constraint, col_a, "mis.z_score")
    r2_b = ols_r_squared(gene_constraint, col_b, "MTR")
    log(f"      panel A R2 = {r2_a:.4f}")
    log(f"      panel B R2 = {r2_b:.4f}  "
        f"(notebook hardcoded {HARDCODED_R2_PANEL_B})")

    log("[4/4] rendering ...")
    scatter_panel(
        gene_constraint, col_a, "mis.z_score",
        xlabel=r"Fraction of coding bases with $\mathbb{P}(0) > 0.5$",
        ylabel="Missense Z-score",
        title="HMM Constraint Proportion vs\n gnomAD v4 Missense Z-score per gene",
        xlim=(0, 0.7), ylim=(-10, 10), r2=r2_a, r2_at=(0.4, -7),
        out_path=cfg.result(f"{FIG_3A}{suffix}.png"))

    thr_label = "0.5" if args.match_thresholds else "0.6"
    scatter_panel(
        gene_constraint, col_b, "MTR",
        xlabel=rf"Proportion of gene with $\mathbb{{P}}(0) > {thr_label}$",
        ylabel="Missense Tolerance Ratio",
        title="HMM Constraint vs MTR per gene",
        xlim=(0, 0.7), ylim=(0.7, 1.1), r2=r2_b, r2_at=(0.5, 1.05),
        out_path=cfg.result(f"{FIG_3B}{suffix}.png"))

    keep = ["chr", "start", "end", "length", "gene_name", "transcript",
            "mis.z_score", "MTR"] + [c for c in gene_constraint.columns
                                     if c.startswith("proportion_over_")]
    out = cfg.result(f"gene_constraint_proportions{suffix}.tsv.gz")
    gene_constraint[keep].to_csv(out, sep="\t", index=False, compression="gzip")

    print()
    print("=" * 64)
    print(f"FIGURE 3 SUMMARY [{suffix.lstrip('_')}]")
    print("=" * 64)
    print(f"  CDS intervals            : {len(gene_constraint):,}")
    print(f"  panel A threshold        : P(0) > 0.5")
    print(f"  panel B threshold        : P(0) > {thr_label}")
    print(f"  panel A R2 (computed)    : {r2_a:.4f}")
    print(f"  panel B R2 (computed)    : {r2_b:.4f}")
    print(f"  panel B R2 (as printed)  : {HARDCODED_R2_PANEL_B}")
    delta = abs(r2_b - HARDCODED_R2_PANEL_B)
    print(f"  difference               : {delta:.4f}"
          f"{'   <-- published label is wrong' if delta > 0.001 else ''}")
    print("=" * 64)


if __name__ == "__main__":
    main()
