"""The unified meta-regression, shared by the epilepsy and schizophrenia analyses.

Both cohorts run the identical model: Haldane-Anscombe effect sizes from allele
counts, inverse-variance weighted least squares per gene-group, five moderators,
whole-model F-test as the reported p-value. In the notebooks this was two
verbatim copies that differed only in file paths -- the same duplication that
turned one GERP defect into three.

What is genuinely cohort-specific is only the mapping from consequence
annotation to the pLoF and missense indicators: Epi25 uses a small controlled
vocabulary, SCHEMA a larger VEP-style one. That mapping stays with the caller.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant

__all__ = ["MODERATORS", "GROUP_KEY", "haldane_effect_sizes",
           "prepare_regression_input", "fit_per_gene"]

MODERATORS = ["log_constraint", "GERP_RS", "log_pathogenicity",
              "pLoF_ind", "missense_ind"]
GROUP_KEY = ["gene_id", "gene_name", "group"]


def haldane_effect_sizes(df: pd.DataFrame, *, max_allele_count: int = 5
                         ) -> pd.DataFrame:
    """Apply the rare-variant cut and attach effect sizes.

    ``effect_size`` is the Haldane--Anscombe corrected log odds ratio and
    ``var_effect_size`` its Woolf variance, both computed directly from allele
    counts with a 0.5 continuity correction. No Firth statistics are used.

    The cut is on **allele count**, not frequency: across Epi25's ~108,800
    alleles ``ac_case + ac_ctrl <= 5`` is a MAF of roughly 5e-5.
    """
    # No .copy() here, deliberately. Boolean masking already returns a new
    # frame; adding a second copy changes the array's memory layout, which
    # changes which vectorised path numpy selects for np.log, which perturbs
    # effect_size by one unit in the last place. That is numerically
    # irrelevant but it breaks bit-for-bit reproduction of the submitted
    # results, and an exact gate is worth more than a tidier line.
    out = df[(df["ac_ctrl"] + df["ac_case"] <= max_allele_count)
             & (df["an_case"] > 0) & (df["an_ctrl"] > 0)]

    alt_case, alt_ctrl = out["ac_case"], out["ac_ctrl"]
    ref_case = out["an_case"] - out["ac_case"]
    ref_ctrl = out["an_ctrl"] - out["ac_ctrl"]

    out["effect_size"] = np.log(
        ((0.5 + alt_case) * (0.5 + ref_ctrl))
        / ((0.5 + ref_case) * (0.5 + alt_ctrl)))
    out["var_effect_size"] = (1 / (0.5 + ref_case) + 1 / (0.5 + ref_ctrl)
                              + 1 / (0.5 + alt_case) + 1 / (0.5 + alt_ctrl))
    return out


def prepare_regression_input(df: pd.DataFrame, *, clip_epsilon: float = 0.01
                             ) -> pd.DataFrame:
    """Clip, impute, transform and drop, exactly as the notebooks did.

    Note that GERP enters the model untransformed; only the two probabilities
    are put through ``-log1p(-x)``. Missing moderator values are filled with
    the column mean over the whole table, not within gene.
    """
    out = df.copy()
    eps = clip_epsilon
    out["prob_0"] = np.clip(out["prob_0"], eps, 1 - eps)
    out["am_pathogenicity"] = np.clip(out["am_pathogenicity"], eps, 1 - eps)

    cols = ["prob_0", "GERP_RS", "am_pathogenicity"]
    out[cols] = out[cols].fillna(out[cols].mean())
    out["pLoF_ind"] = out["pLoF_ind"].fillna(0)
    out["missense_ind"] = out["missense_ind"].fillna(0)

    out[["log_constraint", "log_pathogenicity"]] = -np.log1p(
        -(out[["prob_0", "am_pathogenicity"]]))

    out = out.dropna(subset=["effect_size", "var_effect_size"])
    return out[out["var_effect_size"] != 0]


def fit_per_gene(input_df: pd.DataFrame, *, progress: bool = True
                 ) -> pd.DataFrame:
    """Weighted least squares per (gene_id, gene_name, group).

    Genes with any missing moderator are skipped, matching the notebooks. The
    reported ``p_unified`` is ``model.f_pvalue`` -- the F-test for the whole
    regression, not a coefficient p-value. External groups have repeatedly
    misread this as a coefficient test.
    """
    grouped = input_df.groupby(GROUP_KEY)
    iterator = grouped
    if progress:
        try:
            from tqdm import tqdm
            iterator = tqdm(grouped, desc="genes", unit="gene")
        except ImportError:
            pass

    required = MODERATORS + ["effect_size", "var_effect_size"]
    results = []

    for (gene_id, gene_name, group), data in iterator:
        if data[required].isnull().any().any():
            continue
        try:
            model = WLS(data["effect_size"], add_constant(data[MODERATORS]),
                        weights=1 / data["var_effect_size"],
                        missing="drop").fit()
        except Exception:
            continue
        results.append({
            "gene_id": gene_id, "gene_name": gene_name, "group": group,
            "n_variants": len(data),
            "p_constraint": model.pvalues["log_constraint"],
            "p_gerp": model.pvalues["GERP_RS"],
            "p_pathogenicity": model.pvalues["log_pathogenicity"],
            "p_pLoF": model.pvalues["pLoF_ind"],
            "p_missense": model.pvalues["missense_ind"],
            "p_unified": model.f_pvalue,
        })

    return pd.DataFrame(results)
