"""Per-interval HMM constraint compared with published constraint metrics (Figure 3).

For each protein-coding CDS interval: the fraction of its bases the HMM calls
constrained, set against gnomAD's missense z-score and the missense tolerance ratio.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["compute_mtr", "load_gene_constraint", "add_constraint_proportions",
           "fraction_in_intervals", "ols_r_squared"]


def compute_mtr(constraint: pd.DataFrame) -> pd.Series:
    """Missense tolerance ratio.

    Observed missense share of observed coding variation, divided by the same
    quantity expected.
    """
    return ((constraint["mis.obs"] / (constraint["mis.obs"] + constraint["syn.obs"]))
            / (constraint["mis.exp"] / (constraint["mis.exp"] + constraint["syn.exp"])))


def load_gene_constraint(gene_annotation: pd.DataFrame,
                         constraint_path) -> pd.DataFrame:
    """Join GENCODE CDS intervals to gnomAD constraint metrics.

    Joined on (gene_name, transcript) against (gene, transcript), inner, so only
    transcripts present in both survive. ``length`` is the inclusive interval
    width, ``abs(end - start) + 1``.
    """
    constraint = pd.read_csv(constraint_path, sep="\t")
    constraint["MTR"] = compute_mtr(constraint)

    merged = pd.merge(gene_annotation, constraint,
                      left_on=["gene_name", "transcript"],
                      right_on=["gene", "transcript"], how="inner")
    merged["length"] = (merged["end"] - merged["start"]).abs() + 1
    return merged


def fraction_in_intervals(intervals: pd.DataFrame,
                          positions_by_chrom: dict[str, np.ndarray]
                          ) -> np.ndarray:
    """Fraction of each interval's bases that appear in the position set.

    ``intervals`` needs ``chr``, ``start``, ``end`` and ``length``. Positions
    are inclusive of both endpoints. Uses two binary searches per interval
    rather than a scan of every position.
    """
    sorted_positions = {c: np.sort(np.asarray(p))
                        for c, p in positions_by_chrom.items()}

    out = np.zeros(len(intervals), dtype=float)
    chroms = intervals["chr"].to_numpy()
    starts = intervals["start"].to_numpy()
    ends = intervals["end"].to_numpy()
    lengths = intervals["length"].to_numpy()

    for chrom in np.unique(chroms):
        positions = sorted_positions.get(chrom)
        mask = chroms == chrom
        if positions is None or positions.size == 0:
            continue
        lo = np.searchsorted(positions, starts[mask], side="left")
        hi = np.searchsorted(positions, ends[mask], side="right")
        out[mask] = (hi - lo) / lengths[mask]
    return out


def add_constraint_proportions(gene_constraint: pd.DataFrame,
                               predictions: pd.DataFrame,
                               thresholds=(0.5, 0.6, 0.8)) -> pd.DataFrame:
    """Attach ``proportion_over_<t>`` columns for each threshold."""
    out = gene_constraint.copy()
    for t in thresholds:
        subset = predictions.loc[predictions["prob_0"] > t, ["chr", "pos"]]
        by_chrom = {c: g["pos"].to_numpy()
                    for c, g in subset.groupby("chr", observed=True)}
        label = f"proportion_over_{int(round(t * 100))}"
        out[label] = fraction_in_intervals(out, by_chrom)
    return out


def ols_r_squared(frame: pd.DataFrame, x: str, y: str) -> float:
    """R-squared of an ordinary least squares fit of y on x, NaNs dropped."""
    import statsmodels.api as sm
    d = frame.dropna(subset=[x, y])
    model = sm.OLS(d[y], sm.add_constant(d[x])).fit()
    return float(model.rsquared)
