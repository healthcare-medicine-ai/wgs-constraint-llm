"""AlphaMissense loading, collapsed to one score per variant.

`AlphaMissense_hg38.tsv.gz` is the isoform-level release: 71,697,556 rows
covering 71,236,463 distinct variants, so a variant annotated on N transcripts
appears N times. Joined un-collapsed onto a variant table, each such variant
would become N rows and carry N times its weight in the weighted least-squares
fit. Scores are therefore collapsed to one row per (chr, pos, ref, alt).

The collapse takes the maximum: "could this substitution be damaging in any
annotated isoform". Unlike the mean, it does not depend on how many isoforms
happen to be annotated at a locus. Only 1,854 of 454,089 duplicated variants
span more than one UniProt id, so duplicates are overwhelmingly isoforms of a
single protein.

The file has three comment lines before a header line that itself begins with
``#``, hence ``header=3``.
"""

from __future__ import annotations

import pandas as pd

__all__ = ["KEY_COLUMNS", "SCORE_COLUMN", "load_collapsed", "collapse"]

KEY_COLUMNS = ["chr", "pos", "ref", "alt"]
SCORE_COLUMN = "am_pathogenicity"

_HEADER_ROW = 3
_RENAME = {"#CHROM": "chr", "POS": "pos", "REF": "ref", "ALT": "alt"}
_DEFAULT_CHUNK = 5_000_000


def collapse(frame: pd.DataFrame, how: str = "max") -> pd.DataFrame:
    """Reduce to one row per variant, by ``max`` (default) or ``mean``."""
    if how not in ("max", "mean"):
        raise ValueError(f"unsupported collapse {how!r}")
    grouped = frame.groupby(KEY_COLUMNS, sort=False, as_index=False)[SCORE_COLUMN]
    return getattr(grouped, how)()


def load_collapsed(path, *, restrict_to=None, how: str = "max",
                   chunksize: int = _DEFAULT_CHUNK, verbose: bool = True
                   ) -> pd.DataFrame:
    """Load AlphaMissense as one row per variant.

    ``restrict_to`` is an optional ``pandas.MultiIndex`` of (chr, pos, ref,
    alt) to keep. A caller that will left-join onto a known variant set should
    pass it: rows that cannot match never contribute, so filtering first gives
    the same result while keeping the groupby small.
    """
    kept = []
    n_scanned = 0
    reader = pd.read_csv(
        path, sep="\t", header=_HEADER_ROW,
        usecols=list(_RENAME) + [SCORE_COLUMN], chunksize=chunksize)
    for chunk in reader:
        n_scanned += len(chunk)
        chunk = chunk.rename(columns=_RENAME)
        if restrict_to is not None:
            index = pd.MultiIndex.from_arrays([chunk[c] for c in KEY_COLUMNS])
            chunk = chunk[index.isin(restrict_to)]
        if len(chunk):
            kept.append(chunk)
    if not kept:
        raise ValueError("no AlphaMissense rows matched the requested variants")

    frame = pd.concat(kept, ignore_index=True)
    n_before = len(frame)
    frame = collapse(frame, how=how)
    if verbose:
        print(f"      AlphaMissense: scanned {n_scanned:,} rows; "
              f"{n_before:,} matched; collapsed to {len(frame):,}", flush=True)
    return frame
