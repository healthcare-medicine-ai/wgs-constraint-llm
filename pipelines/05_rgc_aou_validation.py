#!/usr/bin/env python
"""Validation of the HMM in All of Us: Figure 2 and the RGC vs HMM comparison.

Uses every site outside chromosome 2 (the HMM's training chromosome) with
predictions in both RGC-ME and AoU.

Figure 2. (A) Observed joint distribution of the RGC and AoU constraint
probabilities, 10 x 10 bins, log colour scale. (B) Chi-square contribution of
each bin against the expectation under independence. Empty bins are left blank.

RGC vs HMM. How well binary RGC variant status (invariant = constrained) and
the binarized HMM (P(0) > 0.5) each predict variant presence in AoU:
  * 2x2 odds ratios with 95% CIs, risk differences, phi, for each predictor;
  * the HMM's odds ratio within RGC-invariant sites (at RGC-variant sites
    the HMM call does not vary);
  * nested logistic models (RGC, HMM, both): McFadden R2 and the likelihood
    ratio test for adding the HMM;
  * sensitivity to the binarization threshold, and the continuous HMM.

Input (results/)
  HMM_rgc_0.9_over20_chr2_joint_predictions_aou_rgc_wes.tsv.gz
Outputs (results/)
  figure_2.png
  rgc_vs_hmm_effect_sizes.csv           marginal odds ratios; McFadden R2 of
                                        the RGC, HMM and joint models
  rgc_vs_hmm_threshold_sensitivity.csv  HMM odds ratio at thresholds 0.2-0.9
  rgc_vs_hmm_by_rgc_status.csv          HMM 2x2 within RGC-invariant sites
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2 as chi2_dist
from scipy.stats import chi2_contingency

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

JOINT = "HMM_rgc_0.9_over20_chr2_joint_predictions_aou_rgc_wes.tsv.gz"
TRAINING_CHROMOSOME = "chr2"
BINS = np.arange(0, 1.1, 0.1)
HMM_CUT = 0.5
CI_Z = 1.959963984540054


def log(msg=""):
    print(msg, flush=True)


# -- Figure 2 ---------------------------------------------------------------

def joint_distribution(df):
    """Binned joint counts and per-bin chi-square; empty bins become NaN."""
    aou = pd.cut(df["prob_0_aou"], bins=BINS, precision=2)
    rgc = pd.cut(df["prob_0_rgc"], bins=BINS, precision=2)
    marginal_aou = aou.value_counts(normalize=True).sort_index()
    marginal_rgc = rgc.value_counts(normalize=True).sort_index()
    expected = np.outer(marginal_aou, marginal_rgc) * len(df)

    counts = pd.DataFrame({"aou": aou, "rgc": rgc}).value_counts().sort_index()
    actual = np.array([[counts.get((a, r), 0) for r in marginal_rgc.index]
                       for a in marginal_aou.index])
    with np.errstate(divide="ignore", invalid="ignore"):
        chi_sqr = (expected - actual) ** 2 / expected
    chi_sqr[actual == 0] = np.nan
    return actual, chi_sqr


def heatmap_panel(values, title, title_kw, heatmap_kw):
    """One 7 x 6 inch panel at 300 dpi, returned as an RGB image."""
    ticklabels = [f"{i:.1f}" for i in BINS]
    fig, ax = plt.subplots(1, figsize=(7, 6))
    sns.heatmap(values, cmap="Blues", annot=False, cbar=True, ax=ax, **heatmap_kw)
    ax.set_title(title, fontsize=14, **title_kw)
    ax.set_xlabel("RGC Constraint Probability", fontsize=12)
    ax.set_ylabel("AoU Constraint Probability", fontsize=12)
    ax.invert_yaxis()
    ax.set_xticks(range(11))
    ax.set_xticklabels(ticklabels)
    ax.set_yticks(range(11))
    ax.set_yticklabels(ticklabels)
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, dpi=300)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def figure_2(actual, chi_sqr, out_path):
    """Panels side by side with a 40 px gap, labelled A and B in Arial Bold."""
    a = heatmap_panel(actual, "Observed Joint Distribution", {}, {"norm": LogNorm()})
    b = heatmap_panel(chi_sqr, r"$\chi^2$ for Observed vs Expected Joint Distribution",
                      {"y": 1, "x": 0.55}, {})
    gap = 40
    fig = Image.new("RGB", (a.width + gap + b.width, a.height), "white")
    fig.paste(a, (0, 0))
    fig.paste(b, (a.width + gap, 0))

    arial_bold = font_manager.findfont(
        font_manager.FontProperties(family="Arial", weight="bold"),
        fallback_to_default=False)
    font = ImageFont.truetype(arial_bold, 99)
    draw = ImageDraw.Draw(fig)
    draw.text((12, 8), "A", fill="black", font=font)
    draw.text((a.width + gap + 12, 8), "B", fill="black", font=font)
    fig.save(out_path, dpi=(300, 300))


# -- RGC vs HMM -------------------------------------------------------------

def two_by_two(pred, outcome, label):
    """Association of a 0/1 predictor (1 = constrained) with AoU variant status."""
    a = int(((pred == 1) & (outcome == 1)).sum())
    b = int(((pred == 1) & (outcome == 0)).sum())
    c = int(((pred == 0) & (outcome == 1)).sum())
    d = int(((pred == 0) & (outcome == 0)).sum())
    # Haldane 0.5 correction for the estimates only if a cell is empty.
    a_, b_, c_, d_ = (a, b, c, d) if min(a, b, c, d) else (a + .5, b + .5, c + .5, d + .5)
    chi2, p, _, _ = chi2_contingency(np.array([[a, b], [c, d]]), correction=False)
    odds = (a_ * d_) / (b_ * c_)
    se_log = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    lo, hi = np.exp(np.log(odds) - CI_Z * se_log), np.exp(np.log(odds) + CI_Z * se_log)
    r1, r0 = a_ / (a_ + b_), c_ / (c_ + d_)
    rd = r1 - r0
    se_rd = np.sqrt(r1 * (1 - r1) / (a_ + b_) + r0 * (1 - r0) / (c_ + d_))
    n = a + b + c + d
    phi = np.sqrt(chi2 / n)
    log(f"  {label}")
    log(f"    AoU variant | constrained   : {a:>12,} / {a + b:>12,} = {r1:.4%}")
    log(f"    AoU variant | unconstrained : {c:>12,} / {c + d:>12,} = {r0:.4%}")
    log(f"    chi2 = {chi2:,.1f} (p = {p:.3g})   phi = {phi:.5f}")
    log(f"    OR = {odds:.4f} [{lo:.4f}, {hi:.4f}]   "
        f"RD = {rd:+.4%} [{rd - CI_Z * se_rd:+.4%}, {rd + CI_Z * se_rd:+.4%}]")
    return dict(label=label, a=a, b=b, c=c, d=d, n=n, chi2=chi2, p=p,
                risk_constrained=r1, risk_unconstrained=r0,
                odds_ratio=odds, or_lo=lo, or_hi=hi, risk_diff=rd, phi=phi)


def bernoulli_ll(k, n, p):
    p = np.clip(np.asarray(p, dtype=float), 1e-15, 1 - 1e-15)
    return float(np.sum(k * np.log(p) + (n - k) * np.log(1 - p)))


def fit_grouped(cells, predictors, ll_null):
    """Binomial GLM on aggregated cells (k variants out of n sites)."""
    endog = np.column_stack([cells["k"], cells["n"] - cells["k"]])
    X = sm.add_constant(cells[predictors].astype(float), has_constant="add")
    res = sm.GLM(endog, X, family=sm.families.Binomial()).fit()
    ll = bernoulli_ll(cells["k"], cells["n"], res.fittedvalues)
    return dict(res=res, ll=ll, mcfadden=1 - ll / ll_null)


def nested_models(df, x_cols):
    """Fit RGC-only, HMM-only and joint models on cells grouped by x_cols."""
    cells = (df.groupby(x_cols)["aou_var"].agg(k="sum", n="size")
               .reset_index().astype(float))
    ll_null = bernoulli_ll(cells["k"], cells["n"], cells["k"].sum() / cells["n"].sum())
    rgc, hmm = x_cols
    return (fit_grouped(cells, [rgc], ll_null), fit_grouped(cells, [hmm], ll_null),
            fit_grouped(cells, [rgc, hmm], ll_null))


def threshold_sensitivity(df):
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
    return pd.DataFrame(rows)


def rgc_vs_hmm(df, cfg):
    df = df.dropna(subset=["prob_0_rgc", "observation_rgc", "observation_aou"])
    df = df.assign(rgc_constrained=(df["observation_rgc"].astype(int) == 0).astype(int),
                   hmm_constrained=(df["prob_0_rgc"] > HMM_CUT).astype(int),
                   aou_var=df["observation_aou"].astype(int))
    log(f"sites: {len(df):,}   AoU variant rate {df['aou_var'].mean():.4%}   "
        f"RGC-invariant {df['rgc_constrained'].mean():.4%}   "
        f"HMM-constrained {df['hmm_constrained'].mean():.4%}")

    log("\nMarginal association with AoU variant status")
    rgc = two_by_two(df["rgc_constrained"], df["aou_var"], "RGC variant status (invariant = constrained)")
    hmm = two_by_two(df["hmm_constrained"], df["aou_var"], f"HMM, P(0) > {HMM_CUT}")

    log("\nHMM within each RGC stratum")
    strata = []
    for value, name in [(0, "RGC-variant sites"), (1, "RGC-invariant sites")]:
        sub = df[df["rgc_constrained"] == value]
        if sub["hmm_constrained"].nunique() < 2:
            log(f"  HMM | {name}: n = {len(sub):,}, HMM call constant; skipped")
            continue
        strata.append(two_by_two(sub["hmm_constrained"], sub["aou_var"], f"HMM | {name}"))

    log("\nNested logistic models")
    m_rgc, m_hmm, m_both = nested_models(df, ["rgc_constrained", "hmm_constrained"])
    joint = m_both["mcfadden"]
    inc_hmm = joint - m_rgc["mcfadden"]
    inc_rgc = joint - m_hmm["mcfadden"]
    for name, m in (("RGC only", m_rgc), ("HMM only", m_hmm), ("both", m_both)):
        log(f"  {name:<9} logLik = {m['ll']:>18,.1f}   McFadden R2 = {m['mcfadden']:.6f}")
    log(f"  HMM over RGC: +{inc_hmm:.6f} ({inc_hmm / joint:.1%} of joint R2)")
    log(f"  RGC over HMM: +{inc_rgc:.6f} ({inc_rgc / joint:.1%} of joint R2)")
    res, conf = m_both["res"], m_both["res"].conf_int()
    for term in ["rgc_constrained", "hmm_constrained"]:
        lo, hi = conf.loc[term]
        log(f"  adjusted OR {term:<16} {np.exp(res.params[term]):.4f} "
            f"[{np.exp(lo):.4f}, {np.exp(hi):.4f}]")
    lr = 2 * (m_both["ll"] - m_rgc["ll"])
    log(f"  LRT adding HMM to RGC: chi2(1) = {lr:,.1f}, p = {chi2_dist.sf(lr, 1):.3g}")

    log("\nThreshold sensitivity")
    sens = threshold_sensitivity(df)
    log(sens.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    log("\nContinuous HMM (P(0) rounded to 3 decimals)")
    c_rgc, c_hmm, c_both = nested_models(
        df.assign(p_round=df["prob_0_rgc"].round(3)), ["rgc_constrained", "p_round"])
    for name, m in (("RGC only", c_rgc), ("HMM only", c_hmm), ("both", c_both)):
        log(f"  {name:<9} McFadden R2 = {m['mcfadden']:.6f}")
    inc_c = c_both["mcfadden"] - c_rgc["mcfadden"]
    log(f"  HMM over RGC: +{inc_c:.6f} ({inc_c / c_both['mcfadden']:.1%} of joint R2)")

    summary = pd.DataFrame([
        dict(analysis="RGC binary, marginal", odds_ratio=rgc["odds_ratio"],
             ci_lo=rgc["or_lo"], ci_hi=rgc["or_hi"], phi=rgc["phi"]),
        dict(analysis="HMM binarized, marginal", odds_ratio=hmm["odds_ratio"],
             ci_lo=hmm["or_lo"], ci_hi=hmm["or_hi"], phi=hmm["phi"]),
        dict(analysis="nested: McFadden RGC / HMM / both",
             odds_ratio=m_rgc["mcfadden"], ci_lo=m_hmm["mcfadden"], ci_hi=joint,
             phi=inc_hmm / joint),
    ])
    summary.to_csv(cfg.result("rgc_vs_hmm_effect_sizes.csv"), index=False)
    sens.to_csv(cfg.result("rgc_vs_hmm_threshold_sensitivity.csv"), index=False)
    pd.DataFrame(strata).to_csv(cfg.result("rgc_vs_hmm_by_rgc_status.csv"), index=False)


def main():
    cfg = get_config()
    log(cfg.describe() + "\n")
    df = pd.read_csv(cfg.result(JOINT), sep="\t",
                     usecols=["chr", "prob_0_aou", "prob_0_rgc",
                              "observation_rgc", "observation_aou"])
    df = df[df["chr"] != TRAINING_CHROMOSOME]
    log(f"{len(df):,} sites outside {TRAINING_CHROMOSOME}")

    actual, chi_sqr = joint_distribution(df)
    figure_2(actual, chi_sqr, cfg.result("figure_2.png"))
    log(f"Figure 2: {int((actual == 0).sum())} empty bins, "
        f"total chi-square {np.nansum(chi_sqr):.6e}\n")

    rgc_vs_hmm(df, cfg)
    log("\nwrote figure_2.png, rgc_vs_hmm_effect_sizes.csv, "
        "rgc_vs_hmm_threshold_sensitivity.csv, rgc_vs_hmm_by_rgc_status.csv")


if __name__ == "__main__":
    main()
