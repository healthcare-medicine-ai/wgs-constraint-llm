#!/usr/bin/env python
"""Stage 12 -- Figure 1, observed versus predicted mutations across SCN1A.

Extracted from `HMM Mutation Predictions.ipynb` cells 21 and 22.

Two things about the original are worth stating rather than burying:

  * The plotted window is not the gene. It is the gene's index range shifted by
    +750 at the start and -4700 at the end -- hand-tuned offsets that crop to a
    visually interesting region. They are preserved here as named constants so
    the figure reproduces, but they are arbitrary and the caption should not
    imply the panel spans SCN1A in full.
  * `plot_subsequence` is defined twice in that notebook, at cell 6 and again at
    cell 21. Cell 22, which saves the published figure, runs after the second
    definition, so that is the one implemented here: three panels, with a
    centred moving average of the inverted observation track in the middle.

  python pipelines/12_figure1_scn1a.py            # regenerate and gate
  python pipelines/12_figure1_scn1a.py --full-gene   # no crop, for inspection
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

PREDICTIONS = "HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz"
FIGURE = "Figure 1: observed vs predicted mutations for SCN1A"

GENE = "SCN1A"
MA_WINDOW = 3
START_OFFSET = 750     # cell 22, hand-tuned
END_OFFSET = -4700     # cell 22, hand-tuned


def log(m):
    print(m, flush=True)


def gene_bounds(gtf_path, gene_name):
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None,
                      names=["chr", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attribute"],
                      dtype={"start": int, "end": int})
    gtf["gene_type"] = gtf["attribute"].str.extract(r'gene_type "(.*?)"')
    gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "(.*?)"')
    gtf = gtf[(gtf["gene_type"] == "protein_coding") & (gtf["feature"] == "CDS")]
    g = gtf[gtf["gene_name"] == gene_name]
    if g.empty:
        raise SystemExit(f"{gene_name} not found among protein-coding CDS")
    return max(g["chr"]), int(g["start"].min()), int(g["end"].max())


def plot_subsequence(observations, predictions, start_idx, end_idx, *,
                     gene_name, ma_window, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 6), sharex=True)
    plt.suptitle(f"Observed Mutations vs MA({ma_window}) vs HMM Predictions for "
                 + gene_name, y=0.94)

    axes[0].bar(range(start_idx, end_idx), observations[start_idx:end_idx],
                width=1.0, color="black")
    axes[0].set_ylabel("Observation")
    axes[0].margins(x=0)

    flipped = pd.Series(1 - np.asarray(observations))
    moving_avg = flipped.rolling(window=ma_window, min_periods=1,
                                 center=True).mean()
    axes[1].bar(range(start_idx, end_idx), moving_avg[start_idx:end_idx],
                width=1.0, color="orange")
    axes[1].set_ylabel(f"MA({ma_window})")
    axes[1].margins(x=0)

    axes[2].bar(range(start_idx, end_idx), predictions[start_idx:end_idx, 0],
                width=1.0, color="#8C1515")
    axes[2].set_xlabel("Position")
    axes[2].set_ylabel("Probability of 0")
    axes[2].margins(x=0)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-gene", action="store_true",
                    help="ignore the hand-tuned crop offsets")
    ap.add_argument("--suffix", default="_REGEN")
    args = ap.parse_args()
    cfg = get_config()
    log(cfg.describe() + "\n")

    log("[1/3] loading HMM predictions ...")
    predictions_df = pd.read_csv(cfg.result(PREDICTIONS), sep="\t")
    log(f"      {len(predictions_df):,} positions")

    log(f"[2/3] locating {GENE} ...")
    chrom, start_pos, end_pos = gene_bounds(cfg.data("gene_annotation"), GENE)
    log(f"      {chrom}:{start_pos:,}-{end_pos:,}")

    on_chrom = predictions_df["chr"] == chrom
    start_idx = predictions_df[on_chrom
                               & (predictions_df["pos"] >= start_pos)]["pos"].idxmin()
    end_idx = predictions_df[on_chrom
                             & (predictions_df["pos"] <= end_pos)]["pos"].idxmax()
    log(f"      index range {start_idx:,}-{end_idx:,} "
        f"({end_idx - start_idx:,} positions)")

    if args.full_gene:
        lo, hi, tag = start_idx, end_idx, args.suffix + "_FULLGENE"
    else:
        lo, hi = start_idx + START_OFFSET, end_idx + END_OFFSET
        tag = args.suffix
        log(f"      cropped to {lo:,}-{hi:,} "
            f"(offsets {START_OFFSET:+}, {END_OFFSET:+})")

    log("[3/3] rendering ...")
    out = cfg.result(f"{FIGURE}{tag}.png")
    plot_subsequence(
        predictions_df["observation"],
        predictions_df[["prob_0", "prob_1"]].to_numpy(),
        lo, hi, gene_name=GENE, ma_window=MA_WINDOW, out_path=out)
    log(f"      {out.name}")

    print()
    print("=" * 64)
    print("FIGURE 1")
    print("=" * 64)
    print(f"  gene            : {GENE} at {chrom}:{start_pos:,}-{end_pos:,}")
    print(f"  full index span : {start_idx:,}-{end_idx:,}")
    print(f"  plotted span    : {lo:,}-{hi:,}"
          f"{'  (uncropped)' if args.full_gene else ''}")
    print(f"  fraction shown  : "
          f"{100*(hi-lo)/max(end_idx-start_idx, 1):.1f}% of the gene's positions")
    print(f"  moving average  : {MA_WINDOW}")
    print("=" * 64)


if __name__ == "__main__":
    main()
