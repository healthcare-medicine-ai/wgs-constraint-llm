"""Coordinate-convention tests for the GERP annotation.

These are the tests whose absence let a one-base shift survive from August 2025
into a submitted manuscript. They build a small bigWig with a known, distinct
value at every base, so any off-by-one is unambiguous rather than plausible.

Run:  pytest tests/ -v
"""

import numpy as np
import pandas as pd
import pyBigWig
import pytest

from wgs_constraint.gerp import GERP_COLUMN, annotate_with_gerp, chromosome_sort_key

CHROM = "chr1"
LENGTH = 50


@pytest.fixture
def bigwig(tmp_path):
    """A bigWig where 0-based coordinate i holds value i * 10.

    So 1-based position p holds (p - 1) * 10, which makes the expected value
    for any lookup trivially checkable by eye.
    """
    path = tmp_path / "test.bw"
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(CHROM, LENGTH)])
    bw.addEntries(
        [CHROM] * LENGTH,
        list(range(LENGTH)),
        ends=list(range(1, LENGTH + 1)),
        values=[float(i * 10) for i in range(LENGTH)],
    )
    bw.close()
    return path


def test_one_based_positions_read_the_correct_base(bigwig):
    """1-based position p must return the value at 0-based p-1."""
    positions = [1, 2, 10, 25, 50]
    pred = pd.DataFrame({"chr": [CHROM] * len(positions), "pos": positions})

    out = annotate_with_gerp(pred, bigwig, position_base=1, progress=False)

    expected = [(p - 1) * 10.0 for p in positions]
    assert out[GERP_COLUMN].tolist() == expected


def test_zero_based_positions_read_the_correct_base(bigwig):
    positions = [0, 1, 9, 24, 49]
    pred = pd.DataFrame({"chr": [CHROM] * len(positions), "pos": positions})

    out = annotate_with_gerp(pred, bigwig, position_base=0, progress=False)

    assert out[GERP_COLUMN].tolist() == [float(p * 10) for p in positions]


def test_the_two_conventions_differ_by_exactly_one_base(bigwig):
    """The original defect, stated as a test.

    Reading 1-based data as though it were 0-based returns the next base along.
    """
    pred = pd.DataFrame({"chr": [CHROM] * 5, "pos": [10, 11, 12, 13, 14]})

    correct = annotate_with_gerp(pred, bigwig, position_base=1, progress=False)
    shifted = annotate_with_gerp(pred, bigwig, position_base=0, progress=False)

    assert (shifted[GERP_COLUMN] - correct[GERP_COLUMN] == 10.0).all()


def test_final_base_of_chromosome_is_retained(bigwig):
    """The old guard (pos < chr_len) silently discarded the last base."""
    pred = pd.DataFrame({"chr": [CHROM], "pos": [LENGTH]})

    out = annotate_with_gerp(pred, bigwig, position_base=1, progress=False)

    assert len(out) == 1
    assert out[GERP_COLUMN].iloc[0] == (LENGTH - 1) * 10.0


def test_position_zero_is_rejected_not_wrapped(bigwig):
    """With a -1 index, pos==0 would wrap to the chromosome's last base.

    It must be dropped as out of range instead.
    """
    pred = pd.DataFrame({"chr": [CHROM, CHROM], "pos": [0, 5]})

    out = annotate_with_gerp(pred, bigwig, position_base=1, progress=False)

    assert out["pos"].tolist() == [5]
    assert out.attrs["gerp_out_of_bounds_dropped"] == 1
    assert (LENGTH - 1) * 10.0 not in out[GERP_COLUMN].tolist()


def test_position_beyond_chromosome_is_dropped(bigwig):
    pred = pd.DataFrame({"chr": [CHROM, CHROM], "pos": [LENGTH + 1, 5]})

    out = annotate_with_gerp(pred, bigwig, position_base=1, progress=False)

    assert out["pos"].tolist() == [5]


def test_invalid_position_base_is_rejected(bigwig):
    pred = pd.DataFrame({"chr": [CHROM], "pos": [1]})
    with pytest.raises(ValueError, match="position_base"):
        annotate_with_gerp(pred, bigwig, position_base=2, progress=False)


def test_missing_required_column_is_rejected(bigwig):
    with pytest.raises(KeyError, match="pos"):
        annotate_with_gerp(pd.DataFrame({"chr": [CHROM]}), bigwig,
                           progress=False)


def test_chromosomes_absent_from_the_bigwig_are_skipped(bigwig):
    pred = pd.DataFrame({"chr": [CHROM, "chrZZ"], "pos": [5, 5]})

    out = annotate_with_gerp(pred, bigwig, progress=False)

    assert out["chr"].tolist() == [CHROM]


def test_chromosome_sort_order():
    names = ["chr10", "chrM", "chr2", "chrX", "chr1", "chrY"]
    assert sorted(names, key=chromosome_sort_key) == [
        "chr1", "chr2", "chr10", "chrX", "chrY", "chrM"]
