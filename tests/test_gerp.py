"""Coordinate-convention tests for the GERP annotation.

They build a small bigWig with a known, distinct value at every base, so any
off-by-one is unambiguous.

Run:  pytest tests/ -v
"""

import pandas as pd
import pyBigWig
import pytest

from wgs_constraint.gerp import GERP_COLUMN, annotate_with_gerp, chromosome_sort_key

CHROM = "chr1"
LENGTH = 50


@pytest.fixture
def bigwig(tmp_path):
    """A bigWig where 0-based coordinate i holds value i * 10.

    So 1-based position p holds (p - 1) * 10.
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
    """1-based position p must return the value at 0-based p-1, not p."""
    positions = [1, 2, 10, 25, 50]
    pred = pd.DataFrame({"chr": [CHROM] * len(positions), "pos": positions})

    out = annotate_with_gerp(pred, bigwig, progress=False)

    assert out[GERP_COLUMN].tolist() == [(p - 1) * 10.0 for p in positions]


def test_final_base_of_chromosome_is_retained(bigwig):
    pred = pd.DataFrame({"chr": [CHROM], "pos": [LENGTH]})

    out = annotate_with_gerp(pred, bigwig, progress=False)

    assert len(out) == 1
    assert out[GERP_COLUMN].iloc[0] == (LENGTH - 1) * 10.0


def test_position_zero_is_rejected_not_wrapped(bigwig):
    """pos == 0 would index -1, the chromosome's last base; it must be dropped."""
    pred = pd.DataFrame({"chr": [CHROM, CHROM], "pos": [0, 5]})

    out = annotate_with_gerp(pred, bigwig, progress=False)

    assert out["pos"].tolist() == [5]
    assert out.attrs["gerp_out_of_bounds_dropped"] == 1
    assert (LENGTH - 1) * 10.0 not in out[GERP_COLUMN].tolist()


def test_position_beyond_chromosome_is_dropped(bigwig):
    pred = pd.DataFrame({"chr": [CHROM, CHROM], "pos": [LENGTH + 1, 5]})

    out = annotate_with_gerp(pred, bigwig, progress=False)

    assert out["pos"].tolist() == [5]


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
