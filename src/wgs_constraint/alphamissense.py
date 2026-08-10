"""AlphaMissense loading, collapsed to one score per variant.

`AlphaMissense_hg38.tsv.gz` is the **isoform-level** release: 71,697,556 rows
covering 71,236,463 distinct variants. A variant annotated on N transcripts
appears N times.

Merging it un-collapsed into a variant table multiplies those rows. That is not
merely a cosmetic count problem -- in a weighted least squares fit the duplicated
rows also carry duplicated **weight**, so a locus with many annotated transcripts
accumulates apparent evidence it has not earned. Measured on the Epi25 analysis:
0.77% of rows duplicated, but eight genes changed exome-wide significance,
including all five protocadherins.

Collapsing by maximum is the project default. Duplicates are overwhelmingly
alternative isoforms of a single protein -- only 1,854 of 454,089 duplicated
variants span more than one UniProt id -- so a per-variant summary is well
defined. The maximum answers "could this substitution be damaging in any
expressed isoform" and, unlike the mean, does not depend on how many isoforms
happen to be annotated at that locus, which is an artefact of annotation depth.

The file ships three comment lines before a header line that itself begins with
``#``, hence ``skiprows``/``header=3``.
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
    """Reduce to one row per variant.

    ``how`` is one of ``max`` (default), ``mean`` or ``first``. ``first`` exists
    only to reproduce arbitrary-pick behaviour for comparison; it is not a
    defensible analysis choice because row order is an artefact of the file.
    """
    if how not in ("max", "mean", "first"):
        raise ValueError(f"unsupported collapse {how!r}")
    grouped = frame.groupby(KEY_COLUMNS, sort=False, as_index=False)[SCORE_COLUMN]
    return getattr(grouped, how)()


def load_collapsed(
    path,
    *,
    restrict_to=None,
    how: str = "max",
    chunksize: int = _DEFAULT_CHUNK,
    verbose: bool = True,
) -> pd.DataFrame:
    """Load AlphaMissense as one row per variant.

    Parameters
    ----------
    restrict_to
        Optional ``pandas.MultiIndex`` of (chr, pos, ref, alt) to keep. Callers
        that will left-join onto a known variant set should pass it: rows that
        cannot match never contribute, so filtering first is exactly equivalent
        and turns a 71M-row groupby into one sized to the caller's data.
    how
        Collapse strategy, see :func:`collapse`.
    """
    kept = []
    n_scanned = 0

    reader = pd.read_csv(
        path, sep="\t", header=_HEADER_ROW,
        usecols=list(_RENAME) + [SCORE_COLUMN],
        chunksize=chunksize,
    )
    for chunk in reader:
        n_scanned += len(chunk)
        chunk = chunk.rename(columns=_RENAME)
        if restrict_to is not None:
            index = pd.MultiIndex.from_arrays(
                [chunk[c] for c in KEY_COLUMNS])
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
              f"{n_before:,} matched; collapsed to {len(frame):,} "
              f"({n_before - len(frame):,} transcript duplicates removed)",
              flush=True)

    frame.attrs["alphamissense_rows_scanned"] = n_scanned
    frame.attrs["alphamissense_duplicates_removed"] = n_before - len(frame)
    return frame
