"""Are the associations the HMM adds enriched for known epilepsy genes?

Revision 3, Reviewer 1. Removing HMM constraint from the unified model drops the
number of suggestive gene-group pairs (p < 1e-4) from 136 to 101 at unchanged
calibration. This asks whether the pairs the HMM gains are more often
Genes4Epilepsy (G4E) genes than expected.

Plan, fixed before looking:
  primary   G4E rate among pairs suggestive with the HMM but not without it
            ("gained"), against the G4E rate among all gene-group pairs tested
            (Fisher exact, one-sided for enrichment).
  context   the same for pairs suggestive without the HMM but not with it
            ("lost"), for the full model's and the no-HMM model's suggestive
            sets, and for the RGC-binary swap model; the same at the
            exome-wide threshold.

Inputs are stage 14's comparison table (the three models on identical
gene-groups) and the G4E release named in config.yaml, matched on Ensembl ID.

  python genentech_fix/diag_hmm_gained_g4e.py
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

SUGGESTIVE = 1e-4


def main():
    cfg = get_config()
    thr = cfg.param("exome_wide_threshold")
    cmp = pd.read_csv(cfg.result("R1_model_comparison_full_vs_swap.csv"))
    g4e = pd.read_csv(cfg.data_dir / "EpilepsyGenes_v2025-09.tsv", sep="\t")
    ids = set(g4e["Ensemble_ID"].dropna().astype(str))
    cmp["G4E"] = cmp["gene_id"].astype(str).isin(ids)
    n, k = len(cmp), int(cmp["G4E"].sum())
    base = k / n
    print(f"gene-group pairs tested: {n:,}; in G4E: {k:,} ({base:.2%})  [G4E genes: {len(ids):,}]")

    def report(label, mask):
        m = int(mask.sum())
        g = int(cmp.loc[mask, "G4E"].sum())
        # 2x2: in set vs not, G4E vs not
        table = [[g, m - g], [k - g, (n - k) - (m - g)]]
        odds, p = fisher_exact(table, alternative="greater") if m else (float("nan"), float("nan"))
        genes = sorted(set(cmp.loc[mask & cmp["G4E"], "gene_name"]))
        print(f"  {label:48s} n={m:4d}  G4E={g:3d} ({(g / m if m else 0):6.1%})  "
              f"enrichment OR={odds:5.2f}  one-sided p={p:.2g}")
        return genes

    for name, t in (("suggestive (p < 1e-4)", SUGGESTIVE), (f"exome-wide (p < {thr:g})", thr)):
        full, minus, swap = (cmp[f"p_{m}"] < t for m in ("full", "minus", "swap"))
        print(f"\n== {name}")
        report("full model (with HMM)", full)
        report("without HMM", minus)
        report("swap (RGC binary in place of HMM)", swap)
        gained = report("GAINED with HMM (full, not without-HMM)", full & ~minus)
        lost = report("LOST with HMM (without-HMM, not full)", minus & ~full)
        report("swap only (swap, not full)", swap & ~full)
        print(f"     G4E genes gained: {gained}")
        print(f"     G4E genes lost:   {lost}")


if __name__ == "__main__":
    main()
