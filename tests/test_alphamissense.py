"""Tests for the AlphaMissense per-variant collapse.

The defect these guard against is subtle: the row count inflation was under 1%,
but because duplicated rows also duplicate their weight in a WLS fit, eight
genes changed exome-wide significance. The last test states that mechanism
directly.
"""

import numpy as np
import pandas as pd
import pytest

from wgs_constraint.alphamissense import KEY_COLUMNS, SCORE_COLUMN, collapse


def frame(rows):
    return pd.DataFrame(rows, columns=KEY_COLUMNS + [SCORE_COLUMN])


def test_isoforms_of_one_variant_collapse_to_one_row():
    df = frame([
        ("chr1", 100, "A", "T", 0.10),
        ("chr1", 100, "A", "T", 0.90),   # same variant, another transcript
        ("chr1", 100, "A", "T", 0.50),
        ("chr1", 200, "G", "C", 0.30),
    ])

    out = collapse(df, how="max")

    assert len(out) == 2
    got = out.set_index(KEY_COLUMNS)[SCORE_COLUMN]
    assert got[("chr1", 100, "A", "T")] == 0.90
    assert got[("chr1", 200, "G", "C")] == 0.30


def test_multiallelic_sites_are_kept_separate():
    """Different ALT alleles at one position are different variants."""
    df = frame([
        ("chr1", 100, "A", "T", 0.10),
        ("chr1", 100, "A", "G", 0.80),
    ])

    assert len(collapse(df)) == 2


def test_mean_and_max_differ_as_expected():
    df = frame([
        ("chr1", 100, "A", "T", 0.20),
        ("chr1", 100, "A", "T", 0.80),
    ])

    assert collapse(df, how="max")[SCORE_COLUMN].iloc[0] == pytest.approx(0.80)
    assert collapse(df, how="mean")[SCORE_COLUMN].iloc[0] == pytest.approx(0.50)


def test_already_unique_input_is_unchanged():
    df = frame([
        ("chr1", 100, "A", "T", 0.10),
        ("chr1", 200, "G", "C", 0.30),
    ])

    out = collapse(df).sort_values(KEY_COLUMNS).reset_index(drop=True)
    expected = df.sort_values(KEY_COLUMNS).reset_index(drop=True)
    pd.testing.assert_frame_equal(out, expected)


def test_unsupported_strategy_is_rejected():
    with pytest.raises(ValueError, match="unsupported collapse"):
        collapse(frame([("chr1", 1, "A", "T", 0.1)]), how="median")


def test_uncollapsed_merge_inflates_both_count_and_weight():
    """Why this matters: the WLS weight is what actually moves p-values.

    One variant on three transcripts contributes three identical rows, so its
    inverse-variance weight is counted three times.
    """
    variants = pd.DataFrame({
        "chr": ["chr1", "chr1"],
        "pos": [100, 200],
        "ref": ["A", "G"],
        "alt": ["T", "C"],
        "var_effect_size": [0.25, 0.25],
    })
    am = frame([
        ("chr1", 100, "A", "T", 0.10),
        ("chr1", 100, "A", "T", 0.50),
        ("chr1", 100, "A", "T", 0.90),   # three transcripts
        ("chr1", 200, "G", "C", 0.30),
    ])

    naive = variants.merge(am, on=KEY_COLUMNS, how="left")
    fixed = variants.merge(collapse(am), on=KEY_COLUMNS, how="left")

    assert len(naive) == 4 and len(fixed) == 2

    weight = lambda d: (1 / d["var_effect_size"]).sum()
    # The first variant's weight is tripled relative to the second.
    assert weight(naive) == pytest.approx(2 * weight(fixed))
