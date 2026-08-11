"""Tests for gene-level constraint aggregation.

The interval-overlap function replaces a linear scan with binary search. That is
the only place in the overhaul where an algorithm changed rather than moving,
so it is checked against a literal transcription of the original.
"""

import numpy as np
import pandas as pd
import pytest

from wgs_constraint.gene_constraint import (
    compute_mtr, fraction_in_intervals, ols_r_squared,
)


def notebook_overlap(intervals, positions_by_chrom):
    """Verbatim transcription of the notebook's implementation."""
    def calculate_overlap(row, chr_pos_dict):
        positions = chr_pos_dict[row["chr"]]
        return np.sum((row["start"] <= positions)
                      & (positions <= row["end"])) / row["length"]
    return intervals.apply(
        lambda r: calculate_overlap(r, positions_by_chrom), axis=1).to_numpy()


def make_intervals(rows):
    df = pd.DataFrame(rows, columns=["chr", "start", "end"])
    df["length"] = (df["end"] - df["start"]).abs() + 1
    return df


def test_matches_the_notebook_on_a_worked_example():
    intervals = make_intervals([("chr1", 10, 19), ("chr1", 100, 109)])
    positions = {"chr1": np.array([9, 10, 15, 19, 20, 105])}

    fast = fraction_in_intervals(intervals, positions)
    slow = notebook_overlap(intervals, positions)

    assert fast == pytest.approx(slow)
    # 10, 15 and 19 fall inside [10, 19]; length 10
    assert fast[0] == pytest.approx(0.3)
    assert fast[1] == pytest.approx(0.1)


def test_endpoints_are_inclusive_at_both_ends():
    intervals = make_intervals([("chr1", 5, 7)])
    positions = {"chr1": np.array([5, 7])}
    assert fraction_in_intervals(intervals, positions)[0] == pytest.approx(2 / 3)


def test_unsorted_positions_are_handled():
    """The notebook scanned linearly so order never mattered; binary search
    requires sorting, which the function must do itself."""
    intervals = make_intervals([("chr1", 10, 20)])
    shuffled = {"chr1": np.array([20, 11, 5, 19, 100, 10])}
    ordered = {"chr1": np.array([5, 10, 11, 19, 20, 100])}
    assert (fraction_in_intervals(intervals, shuffled)
            == pytest.approx(fraction_in_intervals(intervals, ordered)))


def test_agrees_with_the_notebook_on_random_data():
    rng = np.random.default_rng(20260811)
    for _ in range(25):
        positions = {
            "chr1": rng.integers(0, 5000, size=rng.integers(1, 400)),
            "chr2": rng.integers(0, 5000, size=rng.integers(1, 400)),
        }
        rows = []
        for _ in range(40):
            c = rng.choice(["chr1", "chr2"])
            s = int(rng.integers(0, 4900))
            rows.append((c, s, s + int(rng.integers(0, 100))))
        intervals = make_intervals(rows)

        assert fraction_in_intervals(intervals, positions) == pytest.approx(
            notebook_overlap(intervals, positions))


def test_chromosome_absent_from_predictions_gives_zero():
    intervals = make_intervals([("chrZ", 1, 10)])
    assert fraction_in_intervals(intervals, {"chr1": np.array([1])})[0] == 0.0


def test_empty_position_array_gives_zero():
    intervals = make_intervals([("chr1", 1, 10)])
    assert fraction_in_intervals(intervals, {"chr1": np.array([])})[0] == 0.0


def test_mtr_definition():
    df = pd.DataFrame({"mis.obs": [30.0], "syn.obs": [70.0],
                       "mis.exp": [25.0], "syn.exp": [75.0]})
    # (30/100) / (25/100) = 1.2
    assert compute_mtr(df).iloc[0] == pytest.approx(1.2)


def test_r_squared_is_computed_not_assumed():
    rng = np.random.default_rng(1)
    x = rng.normal(size=500)
    df = pd.DataFrame({"x": x, "y": 3 * x + rng.normal(scale=0.1, size=500)})
    assert ols_r_squared(df, "x", "y") > 0.99
