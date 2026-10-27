"""Two-state HMM of nucleotide-level constraint from variant observations.

For one chromosome, the observation sequence is variant presence (1/0) at every
position that lies in a protein-coding CDS and has gnomAD v4.0 exome coverage
``over_20 > 0.9``, in genomic order. Positions are 1-based: array index ``i`` is
1-based position ``i`` (GTF and VCF coordinates index the arrays directly).

Each position is encoded with the next one as a count of zeros and ones over a
window of ``ORDER = 2`` (so every row sums to 2), and a ``MultinomialHMM`` with
two hidden states is fitted to those counts. The last ``ORDER`` positions of a
chromosome have no full window and are dropped from the output.

Inputs: GENCODE GTF, gnomAD exome coverage summary, a variant table with
``chr`` and ``pos`` columns. Output: per-position ``prob_0``/``prob_1`` (state
posteriors) and the ``observation`` at that position.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "ORDER", "COVERAGE_COLUMN", "COVERAGE_THRESHOLD", "RANDOM_STATE",
    "normalize_chrom", "load_cds_intervals", "load_covered_positions",
    "load_variant_positions", "build_sequence", "encode_counts",
    "fit_hmm", "predict_chromosome",
]

ORDER = 2
COVERAGE_COLUMN = "over_20"
COVERAGE_THRESHOLD = 0.9
RANDOM_STATE = 10

_CHUNK = 20_000_000


def normalize_chrom(name) -> str:
    """``1``, ``"1"`` and ``"chr1"`` all become ``"chr1"``."""
    name = str(name)
    return name if name.startswith("chr") else f"chr{name}"


def load_cds_intervals(gtf_path) -> pd.DataFrame:
    """Protein-coding CDS intervals (``chr``, ``start``, ``end``; 1-based, inclusive)."""
    from .epi25 import load_gene_annotation
    return load_gene_annotation(gtf_path)[["chr", "start", "end"]].reset_index(drop=True)


def load_covered_positions(coverage_path, chroms=None, *,
                           column: str = COVERAGE_COLUMN,
                           threshold: float = COVERAGE_THRESHOLD
                           ) -> dict[str, np.ndarray]:
    """Positions with ``column > threshold`` in a gnomAD coverage summary, by chromosome."""
    keep = None if chroms is None else {normalize_chrom(c) for c in chroms}
    parts: dict[str, list] = {}
    reader = pd.read_csv(coverage_path, sep="\t", compression="gzip",
                         usecols=["locus", column], chunksize=_CHUNK)
    for chunk in reader:
        chunk = chunk.loc[chunk[column] > threshold, "locus"]
        locus = chunk.str.split(":", expand=True)
        if locus.empty:
            continue
        for chrom, pos in locus.groupby(0)[1]:
            if keep is None or chrom in keep:
                parts.setdefault(chrom, []).append(pos.astype(np.int64).to_numpy())
    return {c: np.concatenate(p) for c, p in parts.items()}


def load_variant_positions(variants_path, chroms=None) -> dict[str, np.ndarray]:
    """Variant positions by chromosome from a table with ``chr`` and ``pos`` columns.

    Chromosome names are normalised, so RGC-ME (``1``) and All of Us
    (``chr1``) are matched the same way against the annotation.
    """
    keep = None if chroms is None else {normalize_chrom(c) for c in chroms}
    parts: dict[str, list] = {}
    reader = pd.read_csv(variants_path, sep="\t", usecols=["chr", "pos"],
                         dtype={"chr": str, "pos": np.int64}, chunksize=_CHUNK)
    for chunk in reader:
        for chrom, pos in chunk.groupby("chr")["pos"]:
            chrom = normalize_chrom(chrom)
            if keep is None or chrom in keep:
                parts.setdefault(chrom, []).append(pos.to_numpy())
    return {c: np.concatenate(p) for c, p in parts.items()}


def build_sequence(cds: pd.DataFrame, covered: np.ndarray,
                   variants: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Observation sequence for one chromosome.

    ``cds`` holds that chromosome's ``start``/``end`` (1-based, inclusive);
    ``covered`` and ``variants`` are 1-based positions. Returns
    ``(observations, positions)``: float 0/1 variant presence at, and the
    1-based coordinates of, every CDS position that is also covered.
    """
    starts = cds["start"].to_numpy(dtype=np.int64)
    ends = cds["end"].to_numpy(dtype=np.int64)
    length = 1 + max([int(a.max()) for a in (ends, covered, variants) if len(a)],
                     default=0)

    coding = np.zeros(length, dtype=bool)
    for start, end in zip(starts, ends):
        coding[start:end + 1] = True

    is_covered = np.zeros(length, dtype=bool)
    is_covered[covered] = True

    has_variant = np.zeros(length, dtype=bool)
    has_variant[variants] = True

    mask = coding & is_covered
    return has_variant[mask].astype(np.float64), np.flatnonzero(mask)


def encode_counts(observations: np.ndarray, order: int = ORDER) -> np.ndarray:
    """Count zeros and ones over each window of ``order`` consecutive observations.

    Row ``j`` covers observations ``j .. j + order - 1``; there are
    ``len(observations) - order`` rows, so the final window is not used.
    """
    windows = np.stack([observations[i:i - order] for i in range(order)], axis=1)
    return np.column_stack([(windows == value).sum(axis=1) for value in range(2)])


def fit_hmm(observations: np.ndarray, order: int = ORDER):
    """Fit the two-state multinomial HMM to one chromosome's observations."""
    from hmmlearn import hmm
    model = hmm.MultinomialHMM(n_components=2, random_state=RANDOM_STATE)
    model.fit(encode_counts(observations, order))
    return model


def predict_chromosome(model, chrom: str, observations: np.ndarray,
                       positions: np.ndarray, order: int = ORDER) -> pd.DataFrame:
    """Per-position state posteriors for one chromosome.

    Columns ``chr, pos, prob_0, prob_1, observation``; the last ``order``
    positions are dropped (see :func:`encode_counts`).
    """
    probabilities = model.predict_proba(encode_counts(observations, order))
    return pd.DataFrame({
        "chr": normalize_chrom(chrom),
        "pos": positions[:-order],
        "prob_0": probabilities[:, 0],
        "prob_1": probabilities[:, 1],
        "observation": observations[:-order],
    })
