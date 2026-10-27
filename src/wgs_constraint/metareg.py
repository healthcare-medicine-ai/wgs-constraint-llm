"""The unified meta-regression, shared by the epilepsy and schizophrenia analyses.

Both cohorts run the same model: Haldane-Anscombe effect sizes from allele
counts, inverse-variance weighted least squares per gene-group on five
moderators, and the whole-model F-test as the reported p-value. Only the
mapping from consequence annotation to the pLoF and missense indicators is
cohort-specific, and that stays with the caller.
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

    ``effect_size`` is the Haldane-Anscombe corrected log odds ratio and
    ``var_effect_size`` its Woolf variance, both from allele counts with a 0.5
    continuity correction. The cut is on allele count, ``ac_case + ac_ctrl <=
    max_allele_count``, not on frequency.

    How ``df`` was filtered before this call can change the last bit of
    ``effect_size``: ``np.log`` picks its vectorised loop from the array's
    memory layout, and boolean masking, ``.isin()`` and ``.copy()`` each build
    that array differently. Values agree to about fifteen figures either way,
    but to stay bit-identical with the reference results, filter as
    ``epi25.load_epi25_variants`` does (chained ``!=``), and do not add a
    ``.copy()`` here.
    """
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
    """Clip, impute and transform the moderators.

    The two probabilities (``prob_0``, ``am_pathogenicity``) are clipped to
    ``[eps, 1 - eps]`` and transformed as ``-log1p(-x)``; GERP enters the model
    untransformed. Missing moderator values are filled with the column mean
    over the whole table, not within gene. Rows with a missing or zero
    variance are dropped.
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


def fit_per_gene(input_df: pd.DataFrame, *, moderators: list[str] | None = None,
                 progress: bool = True) -> pd.DataFrame:
    """Weighted least squares per (gene_id, gene_name, group).

    Weights are ``1 / var_effect_size``. ``p_unified`` is ``model.f_pvalue``,
    the F-test for the whole regression, not a coefficient p-value. The full
    model also reports each moderator's coefficient p-value.

    Gene-groups with any missing moderator, or whose fit fails, are skipped.

    ``moderators`` fits a subset of :data:`MODERATORS` (the leave-one-out
    ablations) and then reports only ``n_variants`` and ``p_unified``. The skip
    rule still tests all five moderators, so an ablation covers exactly the
    gene-groups of the full model.

    ``p_unified`` can come back as an object column, because
    ``scipy.stats.f.sf`` returns a 0-d array for some models; use
    ``.to_numpy(dtype=float)`` before taking logs.
    """
    subset = moderators is not None
    mods = list(MODERATORS) if moderators is None else list(moderators)
    unknown = set(mods) - set(MODERATORS)
    if unknown:
        raise ValueError(f"unknown moderator(s): {sorted(unknown)}")
    if not mods:
        raise ValueError("at least one moderator is required")

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
            model = WLS(data["effect_size"], add_constant(data[mods]),
                        weights=1 / data["var_effect_size"],
                        missing="drop").fit()
        except Exception:
            continue
        row = {"gene_id": gene_id, "gene_name": gene_name, "group": group,
               "n_variants": len(data)}
        if not subset:
            row.update({
                "p_constraint": model.pvalues["log_constraint"],
                "p_gerp": model.pvalues["GERP_RS"],
                "p_pathogenicity": model.pvalues["log_pathogenicity"],
                "p_pLoF": model.pvalues["pLoF_ind"],
                "p_missense": model.pvalues["missense_ind"],
            })
        row["p_unified"] = model.f_pvalue
        results.append(row)

    return pd.DataFrame(results)
