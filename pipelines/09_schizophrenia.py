#!/usr/bin/env python
"""Schizophrenia (SCHEMA): unified meta-regression, Table 3 and Figure 5.

Runs the epilepsy model (wgs_constraint.metareg) on SCHEMA's variant results,
after lifting them from hg19 to hg38. What is specific to SCHEMA is the
consequence vocabulary: the non-coding and
synonymous classes below are excluded, and pLoF and missense are set by list
membership. SCHEMA's own gene-level meta-analysis p-value is the comparator.

Inputs
  data/    gencode.v44.basic.annotation.gtf.gz, AlphaMissense_hg38.tsv.gz,
           SCHEMA_variant_results.tsv.bgz, SCHEMA_gene_results.tsv.bgz,
           hg19ToHg38.over.chain.gz
  results/ HMM_rgc_ALL_RS_merged_predictions.tsv.gz (constraint + GERP, stage 02)
Outputs (results/)
  schizophrenia_regression_input.tsv.gz  one row per (variant, gene, group)
  schizophrenia_gene_results.tsv         unified model, every gene-group
  schizophrenia_schema_comparison.csv    gene-groups with >= 25 variants, with
                                         the SCHEMA p-value
  table_3.csv                            those with unified p <= 1e-4, by p
  figure_5.png                           unified vs SCHEMA -log10(p)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from pyliftover import LiftOver  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import (  # noqa: E402
    fit_per_gene, get_config, haldane_effect_sizes,
    load_collapsed, prepare_regression_input,
)
from wgs_constraint.alphamissense import KEY_COLUMNS  # noqa: E402
from wgs_constraint.epi25 import load_gene_annotation  # noqa: E402


EXCLUDED_CONSEQUENCES = [
    "3_prime_UTR_variant", "5_prime_UTR_variant", "coding_sequence_variant",
    "downstream_gene_variant", "intron_variant", "mature_miRNA_variant",
    "null", "splice_region_variant", "synonymous_variant",
    "upstream_gene_variant",
]
PLOF_CONSEQUENCES = [
    "stop_gained", "splice_acceptor_variant", "splice_donor_variant",
    "frameshift_variant",
]
MISSENSE_CONSEQUENCES = [
    "inframe_insertion", "inframe_deletion", "stop_retained_variant",
    "stop_lost", "missense_variant_mpc_<2", "protein_altering_variant",
    "missense_variant_mpc_2-3", "missense_variant_mpc_>=3",
]
VARIANT_COLUMNS = [
    "chr", "pos", "ref", "alt", "gene_id", "gene_name", "group",
    "ac_case", "an_case", "ac_ctrl", "an_ctrl",
    "pLoF_ind", "missense_ind", "effect_size", "var_effect_size",
]
COMPARISON_COLUMNS = {
    "gene_id": "Gene ID", "gene_name": "Gene Name", "group": "Group",
    "n_variants": "# of Variants", "p_constraint": "Constraint p-value",
    "p_gerp": "GERP p-value", "p_pathogenicity": "Pathogenicity p-value",
    "p_pLoF": "pLoF p-value", "p_missense": "Missense p-value",
    "p_unified": "Unified Model p-value", "P meta": "SCHEMA p-value",
}
TABLE_3_COLUMNS = [
    "Gene Name", "Group", "# of Variants", "Constraint p-value", "GERP p-value",
    "Pathogenicity p-value", "SCHEMA p-value", "Unified Model p-value",
]


def log(msg):
    print(msg, flush=True)


def lift_to_hg38(df, chain):
    """Replace SCHEMA's hg19 ``locus`` with hg38 ``chr`` and ``pos``, and split
    ``alleles`` into ``ref`` and ``alt``. Loci that do not map are dropped."""
    locus = df["locus"].str.split(":", n=1, expand=True)
    chroms = ("chr" + locus[0]).to_numpy()
    positions = locus[1].astype("int64").to_numpy()
    lift = LiftOver(str(chain))
    mapped = {}
    for chrom, pos in set(zip(chroms, positions)):
        hits = lift.convert_coordinate(chrom, int(pos) - 1)  # pyliftover is 0-based
        mapped[chrom, pos] = (hits[0][0], hits[0][1] + 1) if hits else (None, -1)
    lifted = [mapped[key] for key in zip(chroms, positions)]
    alleles = df["alleles"].map(json.loads)
    df = df.assign(chr=[c for c, _ in lifted], pos=[p for _, p in lifted],
                   ref=alleles.str[0], alt=alleles.str[1])
    return df[df["pos"] > 0]


def load_variants(cfg, genes):
    """SCHEMA variants on hg38 with effect sizes and pLoF/missense indicators."""
    df = pd.read_csv(cfg.data("schema_variants"), sep="\t", compression="gzip")
    df = lift_to_hg38(df, cfg.data("liftover_chain"))
    # Chained != (not .isin) keeps effect_size bit-identical; see metareg.
    chrom_mask = pd.Series(True, index=df.index)
    for excluded in cfg.param("excluded_chromosomes"):
        chrom_mask &= (df["chr"] != excluded)
    df = df[chrom_mask]
    df = df[~df["consequence"].isin(EXCLUDED_CONSEQUENCES)]
    df = df[~df["consequence"].isna()]

    df = pd.merge(df, genes[["std_gene_id", "gene_name"]].drop_duplicates(),
                  left_on="gene_id", right_on="std_gene_id",
                  how="left").drop("std_gene_id", axis=1)

    df = haldane_effect_sizes(df, max_allele_count=cfg.param("max_allele_count"))
    df["pLoF_ind"] = df["consequence"].isin(PLOF_CONSEQUENCES).astype("int32")
    df["missense_ind"] = df["consequence"].isin(MISSENSE_CONSEQUENCES).astype("int32")
    return df[VARIANT_COLUMNS]


def figure_5(comparison, out_path):
    """2D histogram of unified vs SCHEMA -log10(p). Genes are labelled where the
    two differ by more than 2 and either exceeds 4, stacked within each bin."""
    d = pd.DataFrame({
        "gene_name": comparison["Gene Name"],
        "x": -np.log10(comparison["SCHEMA p-value"]),
        "y": -np.log10(comparison["Unified Model p-value"]),
    }).replace([np.inf, -np.inf], np.nan).dropna()

    x0, x1 = d["x"].min(), d["x"].max()
    y0, y1 = d["y"].min(), d["y"].max()
    x_edges = np.linspace(x0, x1 + 1, 11)
    y_edges = np.linspace(y0, y1 + 1, 11)

    plt.figure(figsize=(8, 6))
    plt.hist2d(d["x"], d["y"], bins=[x_edges, y_edges], cmap="Blues", norm=LogNorm())
    top = max(x1, y1) + 1
    plt.plot([0, top], [0, top], color="red", linestyle="--")

    labelled = d[(abs(d["x"] - d["y"]) > 2) & (d[["x", "y"]].max(axis=1) > 4)].copy()
    labelled["xi"] = np.digitize(labelled["x"], bins=x_edges) - 1
    labelled["yi"] = np.digitize(labelled["y"], bins=y_edges) - 1
    x_centres = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centres = 0.5 * (y_edges[:-1] + y_edges[1:])
    step = (y1 + 1 - y0) * 0.02
    for (xi, yi), group in labelled.groupby(["xi", "yi"]):
        if 0 <= xi < len(x_centres) and 0 <= yi < len(y_centres):
            height = (len(group) - 1) * step
            for k, (_, row) in enumerate(group.sort_values(by="y").iterrows()):
                plt.text(x_centres[xi], y_centres[yi] - height / 2 + k * step,
                         row["gene_name"], fontsize=8, ha="center", fontstyle="italic")

    plt.xlim(x0, x1 + 1)
    plt.ylim(y0, y1 + 1)
    plt.xlabel(r"SCHEMA $-log_{10}(p)$", fontsize=12)
    plt.ylabel(r"Unified Model $-log_{10}(p)$", fontsize=12)
    plt.title("Unified Model vs SCHEMA p-values for schizophrenia",
              fontsize=14, fontweight="bold")
    plt.colorbar(label="Log-Scaled Count")
    plt.savefig(out_path, dpi=300)
    plt.close()


def main():
    cfg = get_config()
    log(cfg.describe())

    log("[1/5] SCHEMA variants ...")
    genes = load_gene_annotation(cfg.data("gene_annotation"))
    variants = load_variants(cfg, genes)
    log(f"      qualifying variant rows: {len(variants):,}")

    log("[2/5] constraint + GERP ...")
    # GERP was annotated as float32 (see gerp.py); reading it back as float32,
    # and prob_0 with round_trip, recovers exactly the values stage 02 wrote.
    constraint = pd.read_csv(cfg.derived("constraint_gerp"), sep="\t",
                             usecols=["chr", "pos", "prob_0", "GERP_RS"],
                             dtype={"chr": "string", "pos": "int64", "GERP_RS": "float32"},
                             float_precision="round_trip")
    log(f"      positions with GERP: {len(constraint):,}")

    log("[3/5] AlphaMissense ...")
    wanted = pd.MultiIndex.from_arrays([variants[c] for c in KEY_COLUMNS]).unique()
    alphamissense = load_collapsed(cfg.data("alphamissense"), restrict_to=wanted,
                                   how=cfg.param("alphamissense_collapse"))
    merged = pd.merge(
        constraint[["chr", "pos", "prob_0", "GERP_RS"]],
        pd.merge(variants, alphamissense, on=KEY_COLUMNS, how="left"),
        on=["chr", "pos"], how="inner")
    merged.to_csv(cfg.result("schizophrenia_regression_input.tsv.gz"),
                  index=False, compression="gzip", sep="\t")
    log(f"      regression input rows: {len(merged):,}")

    log("[4/5] meta-regression ...")
    results = fit_per_gene(prepare_regression_input(
        merged, clip_epsilon=cfg.param("clip_epsilon")))
    results.to_csv(cfg.result("schizophrenia_gene_results.tsv"), index=False, sep="\t")
    results["p_unified"] = results["p_unified"].astype(float)
    log(f"      gene-groups fitted: {len(results):,}")

    log("[5/5] SCHEMA comparison, Table 3, Figure 5 ...")
    results = results[results["n_variants"] >= cfg.param("min_variants")]
    schema = pd.read_csv(cfg.data("schema_genes"), sep="\t", compression="gzip")
    comparison = (pd.merge(results, schema, on=["gene_id", "group"])
                  [list(COMPARISON_COLUMNS)].rename(columns=COMPARISON_COLUMNS))
    comparison.to_csv(cfg.result("schizophrenia_schema_comparison.csv"), index=False)

    table_3 = (comparison[comparison["Unified Model p-value"]
                          <= cfg.param("suggestive_threshold")]
               .sort_values("Unified Model p-value")[TABLE_3_COLUMNS])
    table_3.to_csv(cfg.result("table_3.csv"), index=False)
    figure_5(comparison, cfg.result("figure_5.png"))
    log(f"      compared: {len(comparison):,}; Table 3 rows: {len(table_3)}")


if __name__ == "__main__":
    main()
