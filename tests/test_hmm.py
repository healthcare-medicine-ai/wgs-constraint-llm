"""Observation-sequence construction and order-2 encoding for the HMM, on a toy chromosome."""

import gzip

import numpy as np
import pandas as pd
import pytest

from wgs_constraint import hmm

# CDS 3-7 and 10-12 (1-based inclusive, overlapping a second record at 6-8);
# coverage everywhere except 5 and 11; variants at 4, 8, 10 and 15 (15 is
# outside any CDS, so it must not appear).
CDS = pd.DataFrame({"chr": ["chr1"] * 3, "start": [3, 6, 10], "end": [7, 8, 12]})
COVERED = np.array([p for p in range(1, 16) if p not in (5, 11)])
VARIANTS = np.array([4, 8, 10, 15])

POSITIONS = [3, 4, 6, 7, 8, 10, 12]
OBSERVATIONS = [0, 1, 0, 0, 1, 1, 0]


def original_get_sequence(chr_gene_df, chr_coverage_pos, chr_variant_pos):
    """The original get_sequence(), unchanged apart from its inputs."""
    sequence_length = max(chr_gene_df["end"].max(), chr_coverage_pos.max(), chr_variant_pos.max())
    protein_coding_mask = np.zeros(sequence_length + 1, dtype=bool)
    for start, end in zip(chr_gene_df["start"].to_numpy(), chr_gene_df["end"].to_numpy()):
        protein_coding_mask[start:end + 1] = 1
    coverage_mask = np.zeros(sequence_length + 1, dtype=bool)
    coverage_mask[chr_coverage_pos] = 1
    sequence = np.zeros(sequence_length + 1)
    sequence[chr_variant_pos] = 1
    observations = np.array(sequence[protein_coding_mask & coverage_mask])
    positions = np.where(protein_coding_mask & coverage_mask)[0]
    return observations, positions


def test_sequence_on_toy_chromosome():
    observations, positions = hmm.build_sequence(CDS, COVERED, VARIANTS)
    assert positions.tolist() == POSITIONS
    assert observations.tolist() == OBSERVATIONS
    assert observations.dtype == np.float64


def test_sequence_matches_original_on_random_chromosome():
    rng = np.random.default_rng(0)
    starts = rng.integers(1, 5000, 200)
    cds = pd.DataFrame({"start": starts, "end": starts + rng.integers(0, 300, 200)})
    covered = np.sort(rng.choice(np.arange(1, 6000), 4000, replace=False))
    variants = rng.choice(np.arange(1, 6500), 1500)

    ours = hmm.build_sequence(cds, covered, variants)
    theirs = original_get_sequence(cds, covered, variants)
    np.testing.assert_array_equal(ours[0], theirs[0])
    np.testing.assert_array_equal(ours[1], theirs[1])


def test_counts_use_each_position_and_the_next_and_drop_the_last_window():
    counts = hmm.encode_counts(np.array(OBSERVATIONS, dtype=float))
    # windows (0,1) (1,0) (0,0) (0,1) (1,1); the final (1,0) is not used
    assert counts.tolist() == [[1, 1], [1, 1], [2, 0], [1, 1], [0, 2]]
    assert (counts.sum(axis=1) == hmm.ORDER).all()


def test_predictions_are_trimmed_to_the_encoded_rows():
    pytest.importorskip("hmmlearn")
    rng = np.random.default_rng(1)
    observations = (rng.random(400) < 0.3).astype(float)
    positions = np.arange(100, 500)
    model = hmm.fit_hmm(observations)

    out = hmm.predict_chromosome(model, "1", observations, positions)
    assert list(out.columns) == ["chr", "pos", "prob_0", "prob_1", "observation"]
    assert len(out) == len(observations) - hmm.ORDER
    assert out["pos"].tolist() == positions[:-hmm.ORDER].tolist()
    assert out["observation"].tolist() == observations[:-hmm.ORDER].tolist()
    assert (out["chr"] == "chr1").all()
    np.testing.assert_allclose(out["prob_0"] + out["prob_1"], 1.0)


def test_loaders_threshold_coverage_and_normalise_chromosomes(tmp_path):
    coverage = tmp_path / "coverage.tsv.bgz"
    with gzip.open(coverage, "wt") as f:
        f.write("locus\tover_10\tover_20\n")
        f.write("chr1:3\t1.0\t0.95\nchr1:4\t1.0\t0.9\nchr1:5\t1.0\t0.91\nchr2:7\t1.0\t0.99\n")
    covered = hmm.load_covered_positions(coverage)
    assert {c: v.tolist() for c, v in covered.items()} == {"chr1": [3, 5], "chr2": [7]}
    assert list(hmm.load_covered_positions(coverage, ["2"])) == ["chr2"]

    rgc = tmp_path / "rgc.tsv"
    rgc.write_text("chr\tpos\tvariant\n1\t4\ta\n1\t4\tb\n2\t9\tc\n")
    aou = tmp_path / "aou.tsv"
    aou.write_text("chr\tpos\nchr1\t4\nchr2\t9\n")
    for path in (rgc, aou):
        variants = hmm.load_variant_positions(path, ["chr1"])
        assert list(variants) == ["chr1"]
        assert set(variants["chr1"].tolist()) == {4}
