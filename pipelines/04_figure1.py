#!/usr/bin/env python
"""Figure 1: observed variants and HMM constraint probability along SCN1A.

Plots a stretch of SCN1A's coding sequence: on top the RGC-ME observation
(1 if any variant is seen at the position), below the HMM probability of
observing no variant. The window is the gene's index range in the chr2
predictions table, trimmed by START_OFFSET and END_OFFSET. The x axis is the
row index in that table, not a genomic coordinate.

Inputs
  data/    gencode.v44.basic.annotation.gtf.gz
  results/ HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz
Output (results/)
  figure_1.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

GENE = "SCN1A"
START_OFFSET = 750     # rows trimmed from the start of the gene's range
END_OFFSET = -4700     # rows trimmed from the end
OUT = "figure_1.png"


def gene_rows(predictions, gtf_path, gene_name):
    """First and last row of the predictions table inside the gene's CDS span."""
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None,
                      names=["chr", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attribute"],
                      dtype={"start": int, "end": int})
    gtf["gene_type"] = gtf["attribute"].str.extract(r'gene_type "(.*?)"')
    gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "(.*?)"')
    gene = gtf[(gtf["gene_type"] == "protein_coding") & (gtf["feature"] == "CDS")
               & (gtf["gene_name"] == gene_name)]
    chrom, start, end = max(gene["chr"]), gene["start"].min(), gene["end"].max()

    on_chrom = predictions["chr"] == chrom
    first = predictions[on_chrom & (predictions["pos"] >= start)]["pos"].idxmin()
    last = predictions[on_chrom & (predictions["pos"] <= end)]["pos"].idxmax()
    return first, last


def title_with_italic_gene(fig, prefix, gene, y):
    """Suptitle "<prefix> <gene>" with the gene symbol in italics, centred.

    Drawn as two text pieces. The full title is also set, invisible, so that
    tight_layout reserves the same space for it.
    """
    style = dict(fontsize=matplotlib.rcParams["figure.titlesize"],
                 fontweight=matplotlib.rcParams["figure.titleweight"],
                 va="top", ha="left")
    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()

    def width(s, **extra):
        t = fig.text(0, y, s, **{**style, **extra})
        w = to_fig.transform_bbox(t.get_window_extent(renderer)).width
        t.remove()
        return w

    space = width("a a") - width("aa")
    w_prefix, w_gene = width(prefix), width(gene, fontstyle="italic")
    left = 0.5 - (w_prefix + space + w_gene) / 2
    fig.text(left, y, prefix, **style)
    fig.text(left + w_prefix + space, y, gene, fontstyle="italic", **style)
    fig.suptitle(f"{prefix} {gene}", y=y, color="none")


def plot(observations, prob_0, lo, hi, out_path):
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 4), sharex=True)
    title_with_italic_gene(fig, "Observed Mutations vs HMM Predictions for", GENE, y=0.9)
    x = range(lo, hi)
    axes[0].bar(x, observations[lo:hi], width=1.0, color="black")
    axes[0].set_ylabel("Observation")
    axes[0].margins(x=0)
    axes[1].bar(x, prob_0[lo:hi], width=1.0, color="#8C1515")
    axes[1].set_xlabel("Position index in the chromosome 2 sequence")
    axes[1].set_ylabel("Probability of 0")
    axes[1].margins(x=0)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    cfg = get_config()
    predictions = pd.read_csv(cfg.derived("hmm_predictions"), sep="\t",
                              usecols=["chr", "pos", "observation", "prob_0"])
    first, last = gene_rows(predictions, cfg.data("gene_annotation"), GENE)
    lo, hi = first + START_OFFSET, last + END_OFFSET
    print(f"{GENE}: rows {first:,}-{last:,}, plotted {lo:,}-{hi:,}")

    out = cfg.result(OUT)
    plot(predictions["observation"], predictions["prob_0"].to_numpy(), lo, hi, out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
