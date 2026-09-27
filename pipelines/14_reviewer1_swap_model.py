#!/usr/bin/env python
"""Stage 14 -- Reviewer 1, point 1b: RGC variant status in place of HMM constraint.

Revision 3. Reviewer 1: "A more comprehensive approach could be to fit two
models, one with RGC binary status in place of HMM constraint, and compare
their performance to quantify the relative contribution of each."

Fits the SWAP model -- the unified model with binary RGC variant status
replacing log_constraint, all other moderators unchanged -- on the corrected
regression input, and compares it with the full model and with the model that
simply drops HMM constraint.

Read the within-gene variance diagnostic first. The meta-regression is fit
within each gene, so a moderator carries information only if it varies across
that gene's variants. Every Epi25 site is a site where a variant was seen in a
human cohort, so RGC variant status may be nearly constant within most genes.
If it is, the swap model is not merely weaker but unidentifiable for those
genes, and that is the substantive answer: binary variant status cannot act as
a within-gene moderator, which is why a continuous constraint estimate is
needed.

Ported from `Submissions/4 - Revision 3 (2026-08)/Reviewer 1 analysis/
R1.2 - Competing meta-regression models.py`. Its own WLS loop is kept local
rather than extending `fit_per_gene`, because `rgc_binary` is not one of the
five moderators and the shared path is under a bit-identity gate.

  python pipelines/14_reviewer1_swap_model.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config, prepare_regression_input  # noqa: E402

KEY = ["gene_id", "gene_name", "group"]
SWAP = ["rgc_binary", "GERP_RS", "log_pathogenicity", "pLoF_ind", "missense_ind"]
# The genes named in the Revision 3 abstract (Manny Rivas, 21 Sept 2026).
REPORTED = ["KCNQ2", "STXBP1", "CACNA1A", "SLC6A1", "DYRK1A", "KCNB1",
            "NPRL3", "KANSL1", "ANKRD11", "TBL1XR1", "SP4"]


def log(message=""):
    print(message, flush=True)


def fit(df, predictors):
    """WLS per gene-group, whole-model F-test, as in metareg.fit_per_gene."""
    rows = []
    need = predictors + ["effect_size", "var_effect_size"]
    for (gid, gname, grp), data in df.groupby(KEY):
        if data[need].isnull().any().any():
            continue
        try:
            model = WLS(data["effect_size"], add_constant(data[predictors]),
                        weights=1 / data["var_effect_size"], missing="drop").fit()
        except Exception:
            continue
        rows.append(dict(gene_id=gid, gene_name=gname, group=grp,
                         n_variants=len(data), p_unified=float(model.f_pvalue)))
    return pd.DataFrame(rows)


def main():
    cfg = get_config()
    floor, thr = cfg.param("min_variants"), cfg.param("exome_wide_threshold")
    log(cfg.describe())

    src = cfg.result("constraint_gerp_am_epi25_variants_FIXED.tsv.gz")
    log(f"\n[1/4] loading corrected input {src.name} ...")
    df = prepare_regression_input(pd.read_csv(src, sep="\t"),
                                  clip_epsilon=cfg.param("clip_epsilon"))
    df["chr"] = df["chr"].astype(str)
    df["pos"] = df["pos"].astype(np.int64)

    pred_path = cfg.derived("hmm_predictions")
    head = pd.read_csv(pred_path, sep="\t", nrows=2)
    pos_col = "pos" if "pos" in head.columns else "position"
    if "observation" not in head.columns:
        sys.exit(f"no 'observation' column in {pred_path.name}: {list(head.columns)}")
    rgc = pd.read_csv(pred_path, sep="\t", usecols=["chr", pos_col, "observation"],
                      dtype={"chr": str}).rename(columns={pos_col: "pos",
                                                          "observation": "rgc_binary"})
    rgc["pos"] = rgc["pos"].astype(np.int64)
    before = len(df)
    df = df.merge(rgc, on=["chr", "pos"], how="left")
    log(f"      rows {before:,}; matched to an RGC site {int(df.rgc_binary.notna().sum()):,}")
    df = df.dropna(subset=["rgc_binary"])
    df["rgc_binary"] = df["rgc_binary"].astype(int)
    log(f"      fraction of variant sites that are variant in RGC: {df.rgc_binary.mean():.4f}")

    log("\n[2/4] DIAGNOSTIC -- within-gene variance of rgc_binary")
    gv = (df.groupby(KEY).agg(n_variants=("rgc_binary", "size"),
                              n_unique=("rgc_binary", "nunique"),
                              mean_rgc=("rgc_binary", "mean")).reset_index())
    el = gv[gv.n_variants >= floor]
    const = int((el.n_unique < 2).sum())
    log(f"      gene-groups with >= {floor} variants: {len(el):,}")
    log(f"      ... rgc_binary CONSTANT within gene-group: {const:,} "
        f"({const / max(len(el), 1):.1%})")
    log("      distribution of within-gene mean(rgc_binary):")
    log(el["mean_rgc"].describe().to_string())
    gv.to_csv(cfg.result("R1_rgc_binary_within_gene_variance.csv"), index=False)

    log("\n[3/4] fitting SWAP model (rgc_binary replaces log_constraint) ...")
    swap = fit(df, SWAP)
    swap.to_csv(cfg.result("epilepsy_unified_model_pvalues_SWAP_rgc_binary_FIXED.tsv"),
                sep="\t", index=False)
    log(f"      {len(swap):,} gene-groups fitted")

    log("\n[4/4] comparison: FULL vs SWAP vs MINUS-constraint")

    def load(path, label):
        d = pd.read_csv(path, sep="\t")
        d["p_unified"] = pd.to_numeric(d["p_unified"], errors="coerce")
        d = d[d.n_variants >= floor].dropna(subset=["p_unified"])
        return d[KEY + ["p_unified"]].rename(columns={"p_unified": f"p_{label}"})
    full = load(cfg.result("epilepsy_unified_model_pvalues_FIXED.tsv"), "full")
    minus = load(cfg.result("epilepsy_unified_model_pvalues_MINUS_log_constraint_FIXED.tsv"),
                 "minus")
    sw = swap[swap.n_variants >= floor][KEY + ["p_unified"]].rename(columns={"p_unified": "p_swap"})
    cmp = full.merge(sw, on=KEY).merge(minus, on=KEY)
    log(f"      gene-groups compared on identical footing: {len(cmp):,}")
    rows = []
    for m in ("full", "swap", "minus"):
        sig = cmp[f"p_{m}"] < thr
        found = set(cmp.loc[sig, "gene_name"].str.upper())
        rows.append(dict(model=m, n_exome_wide=int(sig.sum()),
                         n_of_11_reported=sum(g in found for g in REPORTED),
                         median_mlog10p=float(
                             -np.log10(np.clip(cmp[f"p_{m}"], 1e-300, 1)).median())))
    summary = pd.DataFrame(rows)
    log(summary.to_string(index=False))
    detail = cmp[cmp.gene_name.str.upper().isin(REPORTED)].sort_values("p_full")
    log("\n      reported genes under each model:")
    log(detail.to_string(index=False, float_format=lambda v: f"{v:.2e}"))
    only = cmp[(cmp.p_full < thr) != (cmp.p_swap < thr)]
    log(f"\n      gene-groups whose significance differs between FULL and SWAP: {len(only)}")
    if len(only):
        log(only.sort_values("p_full").to_string(index=False, float_format=lambda v: f"{v:.2e}"))
    cmp.to_csv(cfg.result("R1_model_comparison_full_vs_swap.csv"), index=False)
    summary.to_csv(cfg.result("R1_model_comparison_summary.csv"), index=False)
    log("\nwrote R1_model_comparison_full_vs_swap.csv, R1_model_comparison_summary.csv")


if __name__ == "__main__":
    main()
