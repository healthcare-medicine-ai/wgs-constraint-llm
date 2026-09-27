#!/usr/bin/env python
"""Stage 13 -- Reviewer 1, point 1a: RGC variant status vs HMM constraint.

Revision 3. Reviewer 1 objected that the invariant-sites chi-squared test
conditions on observation_rgc == 0, which pins the raw training signal to a
constant and so cannot compare it with the HMM. This compares both predictors
of AoU variant status across ALL sites, with effect sizes and confidence
intervals, and decomposes their contributions with nested logistic models.

Ported from `Submissions/4 - Revision 3 (2026-08)/Reviewer 1 analysis/
R1.1 - RGC vs HMM predictive comparison.py`, which was written to paste into
`RGC + AoU Predictions.ipynb`. The analysis is unchanged; only file access is.
It does not touch GERP or AlphaMissense, so the August 2026 corrections do not
affect it.

Direction convention: prob_0 is the probability a site is CONSTRAINED, so the
raw signal is coded as "invariant in RGC". Both odds ratios are expected < 1.

  python pipelines/13_reviewer1_rgc_vs_hmm.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2 as chi2_dist
from scipy.stats import chi2_contingency

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

JOINT = "HMM_rgc_0.9_over20_chr2_joint_predictions_aou_rgc_wes.tsv.gz"
HMM_CUT = 0.5          # matches the > 0.5 rule used for gene-level aggregation
EXCLUDED_CHR = "chr2"  # the HMM was trained on chromosome 2
CI_Z = 1.959963984540054


def log(message=""):
    print(message, flush=True)


def two_by_two(pred, outcome, label):
    """pred, outcome: 0/1 Series; pred == 1 means 'constrained'."""
    a = int(((pred == 1) & (outcome == 1)).sum())
    b = int(((pred == 1) & (outcome == 0)).sum())
    c = int(((pred == 0) & (outcome == 1)).sum())
    d = int(((pred == 0) & (outcome == 0)).sum())
    if min(a, b, c, d) == 0:
        log(f"  [{label}] zero cell (a={a}, b={b}, c={c}, d={d}); "
            "Haldane 0.5 correction applied to the CI only")
        a_, b_, c_, d_ = a + .5, b + .5, c + .5, d + .5
    else:
        a_, b_, c_, d_ = a, b, c, d
    chi2, p, _, _ = chi2_contingency(np.array([[a, b], [c, d]]), correction=False)
    odds = (a_ * d_) / (b_ * c_)
    se_log = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    lo, hi = np.exp(np.log(odds) - CI_Z * se_log), np.exp(np.log(odds) + CI_Z * se_log)
    r1_, r0_ = a_ / (a_ + b_), c_ / (c_ + d_)
    rd = r1_ - r0_
    se_rd = np.sqrt(r1_ * (1 - r1_) / (a_ + b_) + r0_ * (1 - r0_) / (c_ + d_))
    n = a + b + c + d
    phi = np.sqrt(chi2 / n)
    log(f"  [{label}]")
    log(f"    AoU variant | constrained   : {a:>12,} / {a + b:>12,}  = {r1_:.4%}")
    log(f"    AoU variant | unconstrained : {c:>12,} / {c + d:>12,}  = {r0_:.4%}")
    log(f"    chi2 = {chi2:,.1f}   (p = {p:.3g})")
    log(f"    OR   = {odds:.4f}  95% CI [{lo:.4f}, {hi:.4f}]")
    log(f"    RD   = {rd:+.4%}   95% CI [{rd - CI_Z * se_rd:+.4%}, {rd + CI_Z * se_rd:+.4%}]")
    log(f"    phi  = {phi:.5f}\n")
    return dict(label=label, a=a, b=b, c=c, d=d, n=n, chi2=chi2, p=p,
                risk_constrained=r1_, risk_unconstrained=r0_,
                odds_ratio=odds, or_lo=lo, or_hi=hi, risk_diff=rd, phi=phi)


def bernoulli_ll(k, n, p):
    p = np.clip(np.asarray(p, dtype=float), 1e-15, 1 - 1e-15)
    return float(np.sum(k * np.log(p) + (n - k) * np.log(1 - p)))


def fit_grouped(cells, predictors, ll_null):
    endog = np.column_stack([cells["k"], cells["n"] - cells["k"]])
    if predictors:
        X = sm.add_constant(cells[predictors].astype(float), has_constant="add")
    else:
        X = pd.DataFrame({"const": np.ones(len(cells))})
    res = sm.GLM(endog, X, family=sm.families.Binomial()).fit()
    ll = bernoulli_ll(cells["k"], cells["n"], res.fittedvalues)
    return dict(res=res, ll=ll, mcfadden=1 - ll / ll_null)


def main():
    cfg = get_config()
    log(cfg.describe())
    path = cfg.result(JOINT)
    log(f"\nloading {path.name} ...")
    df = pd.read_csv(path, sep="\t",
                     usecols=["chr", "prob_0_rgc", "observation_rgc", "observation_aou"])
    df = df[df["chr"] != EXCLUDED_CHR]
    df = df.dropna(subset=["prob_0_rgc", "observation_rgc", "observation_aou"])
    df["rgc_constrained"] = (df["observation_rgc"].astype(int) == 0).astype(int)
    df["hmm_constrained"] = (df["prob_0_rgc"] > HMM_CUT).astype(int)
    df["aou_var"] = df["observation_aou"].astype(int)

    log(f"Sites analysed: {len(df):,}  ({EXCLUDED_CHR} excluded)")
    log(f"AoU variant rate overall: {df['aou_var'].mean():.4%}")
    log(f"RGC-invariant fraction:   {df['rgc_constrained'].mean():.4%}")
    log(f"HMM-constrained fraction (p>{HMM_CUT}): {df['hmm_constrained'].mean():.4%}\n")

    out = {}
    log("=" * 78)
    log("R1.1a  Marginal association with AoU variant status, ALL sites")
    log("=" * 78)
    out["rgc"] = two_by_two(df["rgc_constrained"], df["aou_var"],
                            "binary RGC variant status (invariant = constrained)")
    out["hmm"] = two_by_two(df["hmm_constrained"], df["aou_var"],
                            f"binarized HMM constraint (prob_0 > {HMM_CUT})")

    log("=" * 78)
    log("R1.1b  HMM effect WITHIN each stratum of RGC variant status")
    log("=" * 78)
    for val, name in [(0, "RGC-VARIANT sites"), (1, "RGC-INVARIANT sites")]:
        sub = df[df["rgc_constrained"] == val]
        log(f"  stratum: {name}   n = {len(sub):,}")
        if sub["hmm_constrained"].nunique() < 2 or sub["aou_var"].nunique() < 2:
            log("    degenerate stratum; skipped\n")
            continue
        out[f"hmm_within_{val}"] = two_by_two(sub["hmm_constrained"], sub["aou_var"],
                                              f"HMM | {name}")

    log("=" * 78)
    log("R1.1c  Nested logistic models -- decomposition of predictive power")
    log("=" * 78)
    cells = (df.groupby(["rgc_constrained", "hmm_constrained"])["aou_var"]
               .agg(k="sum", n="size").reset_index().astype(float))
    log(cells.to_string(index=False) + "\n")
    ll_null = bernoulli_ll(cells["k"], cells["n"], cells["k"].sum() / cells["n"].sum())
    m_rgc = fit_grouped(cells, ["rgc_constrained"], ll_null)
    m_hmm = fit_grouped(cells, ["hmm_constrained"], ll_null)
    m_both = fit_grouped(cells, ["rgc_constrained", "hmm_constrained"], ll_null)
    for name, m in (("RGC binary only", m_rgc), ("HMM only", m_hmm), ("both", m_both)):
        log(f"  {name:<18} logLik = {m['ll']:>18,.1f}   McFadden R2 = {m['mcfadden']:.6f}")
    inc_hmm = m_both["mcfadden"] - m_rgc["mcfadden"]
    inc_rgc = m_both["mcfadden"] - m_hmm["mcfadden"]
    joint = m_both["mcfadden"]
    log(f"\n  Increment of HMM over RGC alone : {inc_hmm:.6f}   "
        f"({inc_hmm / joint:.1%} of joint R2)")
    log(f"  Increment of RGC over HMM alone : {inc_rgc:.6f}   ({inc_rgc / joint:.1%} of joint R2)")
    log(f"  Shared / non-separable          : {(joint - inc_hmm - inc_rgc) / joint:.1%}")
    res = m_both["res"]
    conf = res.conf_int()
    log("\n  Adjusted odds ratios (joint model):")
    for term in ["rgc_constrained", "hmm_constrained"]:
        lo, hi = conf.loc[term]
        log(f"    {term:<18} OR = {np.exp(res.params[term]):.4f}  "
            f"95% CI [{np.exp(lo):.4f}, {np.exp(hi):.4f}]")
    lr = 2 * (m_both["ll"] - m_rgc["ll"])
    log(f"\n  LRT, adding HMM to RGC-only model: chi2(1) = {lr:,.1f}, "
        f"p = {chi2_dist.sf(lr, 1):.3g}")

    log("\n" + "=" * 78)
    log("R1.1d  Sensitivity: binarization threshold, and continuous HMM")
    log("=" * 78)
    rows = []
    for cut in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        hb = (df["prob_0_rgc"] > cut).astype(int)
        a = int(((hb == 1) & (df["aou_var"] == 1)).sum())
        b = int(((hb == 1) & (df["aou_var"] == 0)).sum())
        c = int(((hb == 0) & (df["aou_var"] == 1)).sum())
        d = int(((hb == 0) & (df["aou_var"] == 0)).sum())
        if min(a, b, c, d) == 0:
            continue
        odds = (a * d) / (b * c)
        se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        rows.append(dict(threshold=cut, n_constrained=a + b, odds_ratio=odds,
                         ci_lo=np.exp(np.log(odds) - CI_Z * se),
                         ci_hi=np.exp(np.log(odds) + CI_Z * se)))
    sens = pd.DataFrame(rows)
    log(sens.to_string(index=False, float_format=lambda v: f"{v:.4f}") + "\n")

    cont = (df.assign(p_round=df["prob_0_rgc"].round(3))
              .groupby(["rgc_constrained", "p_round"])["aou_var"]
              .agg(k="sum", n="size").reset_index().astype(float))
    log(f"  continuous-HMM aggregation: {len(cont):,} cells")
    ll_null_c = bernoulli_ll(cont["k"], cont["n"], cont["k"].sum() / cont["n"].sum())
    c_rgc = fit_grouped(cont, ["rgc_constrained"], ll_null_c)
    c_hmm = fit_grouped(cont, ["p_round"], ll_null_c)
    c_both = fit_grouped(cont, ["rgc_constrained", "p_round"], ll_null_c)
    for name, m in (("RGC binary only", c_rgc), ("HMM continuous only", c_hmm), ("both", c_both)):
        log(f"  {name:<22} McFadden R2 = {m['mcfadden']:.6f}")
    inc_c = c_both["mcfadden"] - c_rgc["mcfadden"]
    log(f"  Increment of continuous HMM over RGC alone: {inc_c:.6f}   "
        f"({inc_c / c_both['mcfadden']:.1%} of joint R2)")

    summary = pd.DataFrame([
        dict(analysis="RGC binary, marginal", odds_ratio=out["rgc"]["odds_ratio"],
             ci_lo=out["rgc"]["or_lo"], ci_hi=out["rgc"]["or_hi"], phi=out["rgc"]["phi"]),
        dict(analysis="HMM binarized, marginal", odds_ratio=out["hmm"]["odds_ratio"],
             ci_lo=out["hmm"]["or_lo"], ci_hi=out["hmm"]["or_hi"], phi=out["hmm"]["phi"]),
        dict(analysis="nested: McFadden RGC / HMM / both",
             odds_ratio=m_rgc["mcfadden"], ci_lo=m_hmm["mcfadden"], ci_hi=m_both["mcfadden"],
             phi=inc_hmm / joint),
    ])
    summary.to_csv(cfg.result("R1_rgc_vs_hmm_effect_sizes.csv"), index=False)
    sens.to_csv(cfg.result("R1_hmm_threshold_sensitivity.csv"), index=False)
    log("\nwrote R1_rgc_vs_hmm_effect_sizes.csv, R1_hmm_threshold_sensitivity.csv")


if __name__ == "__main__":
    main()
