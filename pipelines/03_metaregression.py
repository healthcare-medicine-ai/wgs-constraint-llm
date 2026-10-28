#!/usr/bin/env python
"""Per gene-group meta-regression on the Epi25 variants.

For each (gene, group) the variant effect sizes are regressed on the five
moderators (log constraint, GERP RS, log AlphaMissense pathogenicity, pLoF
and missense indicators) by inverse-variance weighted least squares; the
unified p-value is the whole-model F-test (src/wgs_constraint/metareg.py).
Gene-groups with fewer than min_variants variants (config.yaml) are dropped,
and each remaining one is annotated with the Epi25 gene-level p-values for the
same group and with Genes4Epilepsy membership.

Inputs
  results/epilepsy_regression_input.tsv.gz    from stage 02
  data/<group>_results_2023_12_30_14_36_*.csv  Epi25 gene results, per group (config
                                               keys epi25_genes_<group>)
  data/EpilepsyGenes_v2025-09.tsv              Genes4Epilepsy

Output
  results/epilepsy_gene_results.tsv, one row per gene-group:
    gene_id, gene_name, group, n_variants,
    p_constraint, p_gerp, p_pathogenicity, p_pLoF, p_missense  coefficient p
    p_unified                                                    F-test p
    in_epi25        the gene has a row in that group's Epi25 results
    p_epi25_dm, p_epi25_ptv
                    Epi25 damaging-missense and PTV p-values (blank where
                    Epi25 has no result)
    genes4epilepsy  listed in Genes4Epilepsy

The manuscript's comparison with Epi25 (Tables 1, 2 and S1) uses the
rows with in_epi25 true.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import (  # noqa: E402
    fit_per_gene, get_config, prepare_regression_input,
)

REGRESSION_INPUT = "epilepsy_regression_input.tsv.gz"
GENE_RESULTS = "epilepsy_gene_results.tsv"


def log(msg):
    print(msg, flush=True)


def epi25_gene_results(cfg, groups):
    """Epi25 DM and PTV p-values for every gene in each group, stacked."""
    parts = []
    for group in groups:
        epi25 = pd.read_csv(cfg.data(f"epi25_genes_{group}"))
        dm = next(c for c in epi25.columns if c.startswith("Damaging Missense p"))
        ptv = next(c for c in epi25.columns if c.startswith("PTV p"))
        parts.append(pd.DataFrame({
            "gene_id": epi25["Gene"],
            "group": group,
            "in_epi25": True,
            "p_epi25_dm": pd.to_numeric(epi25[dm], errors="coerce"),
            "p_epi25_ptv": pd.to_numeric(epi25[ptv], errors="coerce"),
        }))
    return pd.concat(parts, ignore_index=True)


def genes4epilepsy(cfg):
    g4e = pd.read_csv(cfg.data("genes4epilepsy"), sep="\t", dtype=str)
    return set(g4e["Gene"].dropna().str.strip().str.upper())


def main():
    cfg = get_config()
    log(cfg.describe())
    min_variants = cfg.param("min_variants")
    threshold = cfg.param("exome_wide_threshold")

    log("[1/3] loading regression input ...")
    variants = prepare_regression_input(
        pd.read_csv(cfg.result(REGRESSION_INPUT), sep="\t"),
        clip_epsilon=cfg.param("clip_epsilon"))
    log(f"      variant rows: {len(variants):,}")

    log("[2/3] fitting per gene-group WLS ...")
    fits = fit_per_gene(variants)
    log(f"      models fitted: {len(fits):,}")

    log("[3/3] annotating with Epi25 and Genes4Epilepsy ...")
    genes = fits[fits["n_variants"] >= min_variants].copy()
    genes["p_unified"] = pd.to_numeric(genes["p_unified"], errors="coerce")
    epi25 = epi25_gene_results(cfg, sorted(genes["group"].unique()))
    genes = genes.merge(epi25, on=["gene_id", "group"], how="left")
    genes["in_epi25"] = genes["in_epi25"].fillna(False).astype(bool)
    names = genes["gene_name"].astype(str).str.strip().str.upper()
    genes["genes4epilepsy"] = names.isin(genes4epilepsy(cfg))
    genes = genes.replace([np.inf, -np.inf], np.nan)

    out = cfg.result(GENE_RESULTS)
    genes.to_csv(out, sep="\t", index=False)
    log(f"      wrote {out} ({len(genes):,} gene-groups with "
        f">= {min_variants} variants; {int(genes['in_epi25'].sum()):,} "
        f"in Epi25)")

    significant = genes[genes["p_unified"] < threshold]
    log(f"\nexome-wide significant (p < {threshold:g}): {len(significant)} "
        f"gene-groups, {significant['gene_name'].nunique()} genes")
    for r in significant.sort_values("p_unified").itertuples():
        log(f"  {r.gene_name:<12}{r.group:<6}{r.p_unified:.3e}  "
            f"n={r.n_variants}  in_epi25={r.in_epi25}  g4e={r.genes4epilepsy}")


if __name__ == "__main__":
    main()
