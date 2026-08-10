#!/usr/bin/env python
"""
Step 3: compare runs against the snapshotted R2 submission.

Two questions, in order:

  A. FIDELITY -- does the control run (--no-fix) reproduce the submitted
     numbers? If not, the port is not faithful and the fixed run means
     nothing. This gate must pass first.

  B. IMPACT -- how far does the corrected run move from the submitted one,
     and specifically: does the set of significant genes change?

Gene-level output only. No variant-level data is printed.
"""

import numpy as np
import pandas as pd

RESULTS = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
           "osthoag/wgs-constraint-llm/results/")
BASELINE = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
            "osthoag/wgs-constraint-llm/genentech_fix/baseline_R2/results/")

SUBMITTED = BASELINE + "epilepsy_unified_model_pvalues.tsv"
REPRO = RESULTS + "epilepsy_unified_model_pvalues_REPRO.tsv"
FIXED = RESULTS + "epilepsy_unified_model_pvalues_FIXED.tsv"

KEY = ["gene_id", "group"]
PCOLS = ["p_constraint", "p_gerp", "p_pathogenicity",
         "p_pLoF", "p_missense", "p_unified"]
MIN_VARIANTS = 25
CAND_THR = 3.4e-7


def rule(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def load(path, tag):
    df = pd.read_csv(path, sep="\t")
    print(f"{tag:<12} rows={len(df):,}  file={path.split('/')[-1]}")
    return df


def compare(a, b, label_a, label_b):
    m = a.merge(b, on=KEY, suffixes=("_a", "_b"), how="outer", indicator=True)

    only_a = int((m["_merge"] == "left_only").sum())
    only_b = int((m["_merge"] == "right_only").sum())
    both = int((m["_merge"] == "both").sum())

    print(f"\ngene-groups in {label_a} only : {only_a:,}")
    print(f"gene-groups in {label_b} only : {only_b:,}")
    print(f"gene-groups in both           : {both:,}")

    m = m[m["_merge"] == "both"].copy()

    nv_same = int((m["n_variants_a"] == m["n_variants_b"]).sum())
    print(f"\nn_variants identical          : {nv_same:,} / {len(m):,} "
          f"({100 * nv_same / max(len(m), 1):.2f}%)")
    if nv_same < len(m):
        d = (m["n_variants_b"] - m["n_variants_a"])
        d = d[d != 0]
        print(f"  n_variants delta: min={d.min():+d} max={d.max():+d} "
              f"mean={d.mean():+.3f}")

    print("\np-value agreement:")
    for c in PCOLS:
        va, vb = m[c + "_a"].to_numpy(float), m[c + "_b"].to_numpy(float)
        ok = np.isfinite(va) & np.isfinite(vb)
        if ok.sum() == 0:
            print(f"  {c:<20} no comparable values")
            continue
        exact = np.isclose(va[ok], vb[ok], rtol=1e-9, atol=0)
        maxrel = np.max(np.abs(va[ok] - vb[ok])
                        / np.maximum(np.abs(va[ok]), 1e-300))
        print(f"  {c:<20} identical={100 * exact.mean():6.2f}%   "
              f"max rel diff={maxrel:.3e}")
    return m


def significant_set(df):
    d = df[df["n_variants"] >= MIN_VARIANTS]
    g = d.groupby("gene_id", as_index=False).agg(
        {"gene_name": "first", "p_unified": "min"})
    return g[g["p_unified"] < CAND_THR].sort_values("p_unified")


def main():
    rule("INPUTS")
    submitted = load(SUBMITTED, "submitted")
    try:
        repro = load(REPRO, "control")
    except FileNotFoundError:
        repro = None
        print("control      NOT FOUND (job may still be running)")
    try:
        fixed = load(FIXED, "fixed")
    except FileNotFoundError:
        fixed = None
        print("fixed        NOT FOUND (job may still be running)")

    if repro is not None:
        rule("A. FIDELITY GATE -- control run vs submitted")
        m = compare(submitted, repro, "submitted", "control")
        pu = m[["p_unified_a", "p_unified_b"]].to_numpy(float)
        ok = np.isfinite(pu).all(axis=1)
        identical = np.isclose(pu[ok, 0], pu[ok, 1], rtol=1e-9, atol=0).all()
        print("\n" + ("PASS: control reproduces the submitted numbers."
                      if identical else
                      "FAIL: control does NOT reproduce the submitted numbers."
                      "\n      Do not interpret the fixed run until this is"
                      " resolved."))

    if fixed is not None:
        rule("B. IMPACT -- fixed run vs submitted")
        compare(submitted, fixed, "submitted", "fixed")

        sub_sig = significant_set(submitted)
        fix_sig = significant_set(fixed)

        rule(f"SIGNIFICANT GENES (p_unified < {CAND_THR:g}, "
             f"n_variants >= {MIN_VARIANTS})")
        print(f"submitted : {len(sub_sig)} genes")
        print(f"fixed     : {len(fix_sig)} genes")

        s_ids, f_ids = set(sub_sig["gene_id"]), set(fix_sig["gene_id"])
        retained = s_ids & f_ids
        lost = s_ids - f_ids
        gained = f_ids - s_ids

        print(f"\nretained  : {len(retained)}")
        print(f"lost      : {len(lost)}")
        print(f"gained    : {len(gained)}")

        name = dict(zip(sub_sig["gene_id"], sub_sig["gene_name"]))
        name.update(dict(zip(fix_sig["gene_id"], fix_sig["gene_name"])))
        sp = dict(zip(sub_sig["gene_id"], sub_sig["p_unified"]))
        fp = dict(zip(fix_sig["gene_id"], fix_sig["p_unified"]))

        print("\n%-14s %-18s %-12s %-12s %s"
              % ("GENE", "GENE_ID", "p_submitted", "p_fixed", "STATUS"))
        for gid in sorted(s_ids | f_ids,
                          key=lambda g: fp.get(g, sp.get(g, 1))):
            status = ("retained" if gid in retained
                      else "LOST" if gid in lost else "GAINED")
            a = f"{sp[gid]:.3e}" if gid in sp else "-"
            b = f"{fp[gid]:.3e}" if gid in fp else "-"
            print(f"{name.get(gid, '?'):<14} {gid:<18} {a:<12} {b:<12} {status}")

    print()


if __name__ == "__main__":
    main()
