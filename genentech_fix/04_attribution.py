#!/usr/bin/env python
"""
Step 4: attribute each gene's movement to a specific defect.

Compares four configurations of the same pipeline:

    submitted   both defects present   (= the R2 baseline)
    gerponly    GERP off-by-one fixed, AlphaMissense fan-out still present
    amonly      AlphaMissense collapsed, GERP off-by-one still present
    fixed       both corrected

For every gene that is significant in either the submitted or the fully
corrected run, reports p under all four and assigns a cause. Written so the
response letter can state which correction moved which gene as a measured
fact rather than an inference from n_variants.

Gene-level output only.
"""

import os

import pandas as pd

RESULTS = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
           "osthoag/wgs-constraint-llm/results/")
BASELINE = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
            "osthoag/wgs-constraint-llm/genentech_fix/baseline_R2/results/")

RUNS = {
    "submitted": BASELINE + "epilepsy_unified_model_pvalues.tsv",
    "gerponly": RESULTS + "epilepsy_unified_model_pvalues_GERPONLY.tsv",
    "amonly": RESULTS + "epilepsy_unified_model_pvalues_AMONLY.tsv",
    "fixed": RESULTS + "epilepsy_unified_model_pvalues_FIXED.tsv",
}

MIN_VARIANTS = 25
CAND_THR = 3.4e-7

# The manuscript's named result (Results section).
ELEVEN = ["KCNQ2", "SCN2A", "STXBP1", "CACNA1A", "SLC6A1", "DYRK1A",
          "KCNB1", "SATB1", "PCDHAC2", "SP4", "RYR2"]


def gene_level(path):
    df = pd.read_csv(path, sep="\t")
    df = df[df["n_variants"] >= MIN_VARIANTS]
    return df.groupby("gene_id", as_index=False).agg(
        {"gene_name": "first", "p_unified": "min", "n_variants": "max"})


def classify(p_sub, p_gerp, p_am, p_fix):
    """Which correction is responsible for the change in significance?"""
    sig = lambda p: pd.notna(p) and p < CAND_THR
    if sig(p_sub) == sig(p_fix):
        return "unchanged"
    # Significance flipped. Which single fix reproduces the final state?
    gerp_flips = sig(p_gerp) == sig(p_fix)
    am_flips = sig(p_am) == sig(p_fix)
    if gerp_flips and not am_flips:
        return "GERP"
    if am_flips and not gerp_flips:
        return "AlphaMissense"
    if gerp_flips and am_flips:
        return "either"
    return "interaction"


def main():
    missing = [k for k, v in RUNS.items() if not os.path.exists(v)]
    if missing:
        print(f"missing run outputs: {missing}")
        print("(attribution jobs may still be running)")
        return

    frames = {}
    for name, path in RUNS.items():
        frames[name] = gene_level(path)
        print(f"{name:<10} genes tested: {len(frames[name]):,}")

    df = frames["submitted"][["gene_id", "gene_name", "p_unified", "n_variants"]]
    df = df.rename(columns={"p_unified": "p_submitted",
                            "n_variants": "n_submitted"})
    for name in ["gerponly", "amonly", "fixed"]:
        col = frames[name][["gene_id", "p_unified", "n_variants"]].rename(
            columns={"p_unified": f"p_{name}", "n_variants": f"n_{name}"})
        df = df.merge(col, on="gene_id", how="outer")

    df["cause"] = [
        classify(r.p_submitted, r.p_gerponly, r.p_amonly, r.p_fixed)
        for r in df.itertuples()
    ]

    changed = df[df["cause"] != "unchanged"].copy()
    changed["sort_p"] = changed[["p_submitted", "p_fixed"]].min(axis=1)
    changed = changed.sort_values("sort_p")

    def show(sub, title):
        print("\n" + "=" * 96)
        print(title)
        print("=" * 96)
        print(f"{'GENE':<16}{'p_submitted':>12}{'p_gerponly':>12}"
              f"{'p_amonly':>12}{'p_fixed':>12}   {'n_sub':>6}{'n_fix':>7}  CAUSE")
        for r in sub.itertuples():
            f = lambda v: f"{v:.2e}" if pd.notna(v) else "-"
            ns = int(r.n_submitted) if pd.notna(r.n_submitted) else 0
            nf = int(r.n_fixed) if pd.notna(r.n_fixed) else 0
            print(f"{str(r.gene_name):<16}{f(r.p_submitted):>12}"
                  f"{f(r.p_gerponly):>12}{f(r.p_amonly):>12}{f(r.p_fixed):>12}"
                  f"   {ns:>6}{nf:>7}  {r.cause}")

    eleven = df[df["gene_name"].isin(ELEVEN)].copy()
    eleven["order"] = eleven["gene_name"].map({g: i for i, g in enumerate(ELEVEN)})
    show(eleven.sort_values("order"),
         "THE 11 NAMED GENES -- p under each configuration")

    show(changed, "ALL GENES WHOSE SIGNIFICANCE CHANGED")

    print("\n" + "=" * 96)
    print("ATTRIBUTION SUMMARY")
    print("=" * 96)
    for cause, n in changed["cause"].value_counts().items():
        print(f"  {cause:<16} {n}")
    print("\n  GERP          = corrected coordinate alone reproduces the outcome")
    print("  AlphaMissense = transcript collapse alone reproduces the outcome")
    print("  either        = both single fixes independently reach the outcome")
    print("  interaction   = neither alone; only the combination changes it")

    out = RESULTS + "attribution_by_gene.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nwritten: {out}")


if __name__ == "__main__":
    main()
