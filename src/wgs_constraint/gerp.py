"""GERP RS annotation of HMM constraint predictions.

This module exists because the same twenty lines were copy-pasted into
`Epilepsy Analysis.ipynb`, `Constraint Measures Comparison.ipynb` and
`Schizophrenia Analysis.ipynb`. A single coordinate error therefore appeared in
three places and was fixed in none of them for a year. There is now one
implementation; the notebooks call it.

Coordinate conventions -- the whole point of this module
--------------------------------------------------------
`pyBigWig.values(chrom, 0, chr_len)` returns a vector whose element ``i`` holds
the score at **0-based** coordinate ``i``, i.e. 1-based position ``i + 1``.

HMM prediction positions are **1-based**. This is not an assumption: in
`HMM Mutation Predictions.ipynb`, ``get_sequence()`` builds its masks by
indexing arrays directly with 1-based coordinates taken from a GENCODE GTF
(1-based inclusive) and a VCF (1-based)::

    protein_coding_mask[start:end+1] = 1
    coverage_mask[chr_coverage_df['pos'].to_numpy()] = 1
    sequence[chr_variants_df['pos'].to_numpy()] = 1
    positions = np.where(protein_coding_mask & coverage_mask)[0]

so array index ``i`` corresponds to 1-based position ``i``, and ``positions``
carries that convention into the ``pos`` column.

The correct lookup for 1-based position ``p`` is therefore ``vector[p - 1]``.

The original code used ``vector[p]``, shifting every score one base downstream.
Because GERP within coding sequence varies sharply between adjacent bases (the
wobble position of a codon is far less constrained than the first two), this is
closer to randomising the track than to blurring it: old and corrected values
correlate at r = 0.173.

The bounds guard has to move with the index. With a ``-1`` index, the original
``pos >= 0`` would let position 0 wrap to the chromosome's final base via
negative indexing, and ``pos < chr_len`` silently discarded the last base.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["annotate_with_gerp", "chromosome_sort_key", "GERP_COLUMN"]

GERP_COLUMN = "GERP_RS"


def chromosome_sort_key(name) -> tuple[int, int]:
    """Sort chromosomes 1..22, then X, Y, M/MT, then anything else."""
    stripped = str(name).replace("chr", "")
    if stripped.isdigit():
        return (0, int(stripped))
    return (1, {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(stripped, 1000))


def annotate_with_gerp(
    predictions: pd.DataFrame,
    bigwig_path,
    *,
    position_base: int = 1,
    drop_missing: bool = True,
    only_chromosomes=None,
    progress: bool = True,
) -> pd.DataFrame:
    """Attach a GERP RS score to each row of ``predictions``.

    Parameters
    ----------
    predictions
        Must contain ``chr`` and ``pos``. ``pos`` is interpreted according to
        ``position_base``.
    bigwig_path
        A GERP RS bigWig on the same assembly as ``predictions``.
    position_base
        ``1`` (default) if ``pos`` is 1-based, which is the convention
        throughout this project. ``0`` is accepted for callers working with
        genuinely 0-based coordinates. Passing the wrong value is the original
        defect, so it is explicit and validated rather than assumed.
    drop_missing
        Drop rows with no GERP value, matching the published behaviour.
    only_chromosomes
        Restrict to these chromosomes. Intended for tests and smoke runs.

    Returns
    -------
    A copy of ``predictions`` with a ``GERP_RS`` column added.
    """
    import pyBigWig

    if position_base not in (0, 1):
        raise ValueError(f"position_base must be 0 or 1, got {position_base!r}")
    for column in ("chr", "pos"):
        if column not in predictions.columns:
            raise KeyError(f"predictions is missing required column {column!r}")

    frame = predictions
    if only_chromosomes is not None:
        frame = frame[frame["chr"].isin(set(only_chromosomes))]

    iterator = sorted(frame["chr"].unique().tolist(), key=chromosome_sort_key)
    if progress:
        try:
            from tqdm import tqdm
            iterator = tqdm(iterator, desc="GERP merge")
        except ImportError:
            pass

    bw = pyBigWig.open(str(bigwig_path))
    try:
        available = bw.chroms()
        chunks = []
        n_out_of_bounds = 0

        for chrom in iterator:
            if chrom not in available:
                continue
            chr_len = int(available[chrom])

            sub = frame.loc[frame["chr"] == chrom].copy()
            if sub.empty:
                continue

            # Valid 1-based positions are [1, chr_len]; valid 0-based are
            # [0, chr_len - 1]. Both map to vector indices [0, chr_len - 1].
            if position_base == 1:
                in_range = (sub["pos"] >= 1) & (sub["pos"] <= chr_len)
            else:
                in_range = (sub["pos"] >= 0) & (sub["pos"] < chr_len)

            if not in_range.all():
                n_out_of_bounds += int((~in_range).sum())
                sub = sub.loc[in_range]
            if sub.empty:
                continue

            # pyBigWig hands back float64; float32 halves peak memory and is
            # well beyond GERP's precision.
            vector = np.array(
                bw.values(chrom, 0, chr_len, numpy=True), dtype=np.float32)

            index = sub["pos"].to_numpy(dtype=np.int64) - position_base
            sub[GERP_COLUMN] = vector[index]

            if drop_missing:
                sub = sub.dropna(subset=[GERP_COLUMN])

            chunks.append(sub)
            del vector
    finally:
        bw.close()

    if not chunks:
        raise ValueError("no chromosomes matched between predictions and bigWig")

    merged = pd.concat(chunks, axis=0, ignore_index=True)
    merged.attrs["gerp_out_of_bounds_dropped"] = n_out_of_bounds
    merged.attrs["gerp_position_base"] = position_base
    return merged
