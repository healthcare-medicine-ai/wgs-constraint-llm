"""Epi25 variant loading, filtering and effect-size construction.

Kept in one place so the notebooks and the pipeline scripts cannot drift apart
on the filters, which are the part of this analysis that most affects the
result and the part external groups most often get wrong.

The rare-variant restriction is an **allele-count** cut, ``ac_case + ac_ctrl <=
5``, not a MAF cut. Across Epi25's ~108,800 alleles that is a MAF of roughly
5e-5, some three orders of magnitude stricter than the ``MAF <= 0.05`` used by
the separate phenome-wide project. Applying the latter admits far more common
variants, which carry smaller variances and therefore much larger weights in
the WLS, and inflates per-gene counts so more genes clear the inclusion floor.

Effect sizes are computed directly from allele counts as the Haldane--Anscombe
corrected log odds ratio with its Woolf variance. No Firth BETA/SE are used or
required.
"""

from __future__ import annotations

import pandas as pd

from .metareg import haldane_effect_sizes

__all__ = ["load_gene_annotation", "load_epi25_variants", "add_effect_sizes"]

VARIANT_KEY = ["chr", "pos", "ref", "alt"]


def load_gene_annotation(gtf_path) -> pd.DataFrame:
    """Protein-coding CDS records from a GENCODE GTF, with a versionless id."""
    genes = pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None,
        names=["chr", "source", "feature", "start", "end", "score",
               "strand", "frame", "attribute"],
        dtype={"start": int, "end": int},
    )
    for field in ("gene_id", "gene_type", "gene_name"):
        genes[field] = genes["attribute"].str.extract(rf'{field} "(.*?)"')
    genes = genes.drop("attribute", axis=1)
    genes = genes[(genes["gene_type"] == "protein_coding")
                  & (genes["feature"] == "CDS")]
    genes["std_gene_id"] = genes["gene_id"].str.split(".").str[0]
    return genes


def load_epi25_variants(
    variants_path,
    gene_annotation: pd.DataFrame,
    *,
    excluded_chromosomes=("chrX", "chrY", "chrMT"),
    excluded_consequences=("synonymous", "non_coding", "NA"),
) -> pd.DataFrame:
    """Load Epi25 results and apply the published inclusion filters."""
    df = pd.read_csv(variants_path, sep="\t")
    df[VARIANT_KEY] = df["variant_id"].str.split(":", expand=True)
    df["pos"] = df["pos"].astype(int)

    df = df[~df["chr"].isin(set(excluded_chromosomes))]
    df = df[~df["consequence"].isin(set(excluded_consequences))]

    return pd.merge(
        df,
        gene_annotation[["std_gene_id", "gene_name"]].drop_duplicates(),
        left_on="gene_id", right_on="std_gene_id", how="left",
    ).drop("std_gene_id", axis=1)


def add_effect_sizes(df: pd.DataFrame, *, max_allele_count: int = 5) -> pd.DataFrame:
    """Apply the rare-variant cut, attach effect sizes, and set indicators.

    The effect-size computation is shared with the schizophrenia analysis and
    lives in :mod:`wgs_constraint.metareg`; only the consequence vocabulary is
    cohort-specific. Epi25 uses a small controlled set, so the mapping is a
    direct comparison here rather than the list membership SCHEMA needs.
    """
    out = haldane_effect_sizes(df, max_allele_count=max_allele_count)
    out["pLoF_ind"] = (out["consequence"] == "pLoF").astype("int32")
    out["missense_ind"] = out["consequence"].isin(
        ["damaging_missense", "other_missense"]).astype("int32")
    return out
