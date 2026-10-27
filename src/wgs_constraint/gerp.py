"""GERP RS annotation of HMM constraint predictions.

Coordinate convention
---------------------
HMM prediction positions (``pos``) are 1-based. ``get_sequence()`` in
`HMM Mutation Predictions.ipynb` builds its masks by indexing arrays directly
with 1-based coordinates from a GENCODE GTF and a VCF, so array index ``i`` is
1-based position ``i`` and that convention carries into ``pos``.

``pyBigWig.values(chrom, 0, chr_len)`` returns a vector whose element ``i`` is
the score at 0-based coordinate ``i``. The score for 1-based position ``p`` is
therefore ``vector[p - 1]``, and valid positions are ``1 <= p <= chr_len``.
Positions outside that range are dropped rather than indexed, so ``p == 0``
cannot wrap round to the last base of the chromosome.
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


def annotate_with_gerp(predictions: pd.DataFrame, bigwig_path, *,
                       only_chromosomes=None, progress: bool = True
                       ) -> pd.DataFrame:
    """Attach a GERP RS score to each row of ``predictions``.

    ``predictions`` must contain ``chr`` and a 1-based ``pos``. Rows with no
    GERP value, on a chromosome absent from the bigWig, or outside the
    chromosome are dropped; the number dropped as out of range is recorded in
    ``attrs["gerp_out_of_bounds_dropped"]``. ``only_chromosomes`` restricts the
    work to a subset, for tests.

    Returns a copy of ``predictions`` with a ``GERP_RS`` column added.
    """
    import pyBigWig

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
            in_range = (sub["pos"] >= 1) & (sub["pos"] <= chr_len)
            if not in_range.all():
                n_out_of_bounds += int((~in_range).sum())
                sub = sub.loc[in_range]
            if sub.empty:
                continue

            # float32 halves peak memory and is beyond GERP's precision.
            vector = np.array(
                bw.values(chrom, 0, chr_len, numpy=True), dtype=np.float32)
            sub[GERP_COLUMN] = vector[sub["pos"].to_numpy(dtype=np.int64) - 1]
            chunks.append(sub.dropna(subset=[GERP_COLUMN]))
            del vector
    finally:
        bw.close()

    if not chunks:
        raise ValueError("no chromosomes matched between predictions and bigWig")

    merged = pd.concat(chunks, axis=0, ignore_index=True)
    merged.attrs["gerp_out_of_bounds_dropped"] = n_out_of_bounds
    return merged
