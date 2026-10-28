"""Epi25 variant loading, filtering and effect-size construction.

The rare-variant restriction is an allele-count cut, ``ac_case + ac_ctrl <=
5`` (a MAF of roughly 5e-5 across Epi25's ~108,800 alleles), not a MAF cut.
Effect sizes are the Haldane-Anscombe corrected log odds ratio from allele
counts with its Woolf variance; no Firth BETA/SE are used.
"""

from __future__ import annotations

import pandas as pd

from .metareg import haldane_effect_sizes

__all__ = ["load_gene_annotation", "load_epi25_variants", "add_effect_sizes"]

VARIANT_KEY = ["chr", "pos", "ref", "alt"]


def load_gene_annotation(gtf_path) -> pd.DataFrame:
    """Protein-coding CDS records from a GENCODE GTF, with versionless gene
    (``std_gene_id``) and transcript (``transcript``) ids."""
    genes = pd.read_csv(
        gtf_path, sep="\t", comment="#", header=None,
        names=["chr", "source", "feature", "start", "end", "score",
               "strand", "frame", "attribute"],
        dtype={"start": int, "end": int},
    )
    for field in ("gene_id", "gene_type", "gene_name", "transcript_id"):
        genes[field] = genes["attribute"].str.extract(rf'{field} "(.*?)"')
    genes = genes.drop("attribute", axis=1)
    genes = genes[(genes["gene_type"] == "protein_coding")
                  & (genes["feature"] == "CDS")]
    genes["std_gene_id"] = genes["gene_id"].str.split(".").str[0]
    genes["transcript"] = genes["transcript_id"].str.split(".").str[0]
    return genes


def load_epi25_variants(
    variants_path,
    gene_annotation: pd.DataFrame,
    *,
    excluded_chromosomes=("chrX", "chrY", "chrMT"),
    excluded_consequences=("synonymous", "non_coding", "NA"),
) -> pd.DataFrame:
    """Load Epi25 variant results, drop excluded chromosomes and consequences,
    and attach gene names from the annotation."""
    df = pd.read_csv(variants_path, sep="\t")
    df[VARIANT_KEY] = df["variant_id"].str.split(":", expand=True)
    df["pos"] = df["pos"].astype(int)

    # Chained != rather than .isin(). Both select the same rows, but the
    # resulting frame's memory layout differs, which changes the last bit of
    # effect_size downstream. This form reproduces the reference results
    # bit for bit.
    chrom_mask = pd.Series(True, index=df.index)
    for excluded in excluded_chromosomes:
        chrom_mask &= (df["chr"] != excluded)
    df = df[chrom_mask]

    consequence_mask = pd.Series(True, index=df.index)
    for excluded in excluded_consequences:
        consequence_mask &= (df["consequence"] != excluded)
    df = df[consequence_mask]

    return pd.merge(
        df,
        gene_annotation[["std_gene_id", "gene_name"]].drop_duplicates(),
        left_on="gene_id", right_on="std_gene_id", how="left",
    ).drop("std_gene_id", axis=1)


def add_effect_sizes(df: pd.DataFrame, *, max_allele_count: int = 5) -> pd.DataFrame:
    """Apply the rare-variant cut, attach effect sizes, and set the pLoF and
    missense indicators from Epi25's consequence vocabulary."""
    out = haldane_effect_sizes(df, max_allele_count=max_allele_count)
    out["pLoF_ind"] = (out["consequence"] == "pLoF").astype("int32")
    out["missense_ind"] = out["consequence"].isin(
        ["damaging_missense", "other_missense"]).astype("int32")
    return out
