"""Which aggregation produces the response letter's R2 = 0.300?

Established already: per CDS interval at P(0) > 0.5 gives 0.1458, matching the
letter's "was 0.146" exactly. So the "before" is pinned. This searches the
plausible "after" variants for 0.300.
"""
import sys
sys.path.insert(0, "src")
import pandas as pd
from wgs_constraint import get_config
from wgs_constraint.gene_constraint import (
    add_constraint_proportions, load_gene_constraint, ols_r_squared)

cfg = get_config()
gtf = pd.read_csv(cfg.data("gene_annotation"), sep="\t", comment="#", header=None,
                  names=["chr", "source", "feature", "start", "end", "score",
                         "strand", "frame", "attribute"], dtype={"start": int, "end": int})
gtf["gene_type"] = gtf["attribute"].str.extract(r'gene_type "(.*?)"')
gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "(.*?)"')
gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "(.*?)"')
gtf["transcript"] = gtf["transcript_id"].str.split(".").str[0]
gtf = gtf[(gtf["gene_type"] == "protein_coding") & (gtf["feature"] == "CDS")]
gtf = gtf.drop(columns=["attribute"])

gc = load_gene_constraint(gtf, cfg.data_dir / "gnomad.v4.0.constraint_metrics.tsv")
pred = pd.read_csv(cfg.result("HMM_rgc_ALL_RS_merged_predictions_ASPUBLISHED.tsv.gz"),
                   sep="\t", usecols=["chr", "pos", "prob_0"])
gc = add_constraint_proportions(gc, pred, thresholds=(0.5, 0.6))
print(f"CDS intervals: {len(gc):,}\n", flush=True)

for t in (50, 60):
    col = f"proportion_over_{t}"
    thr = f"0.{t//10}"
    gc["_n"] = gc[col] * gc["length"]

    lw = gc.groupby("transcript", as_index=False).agg(
        _n=("_n", "sum"), length=("length", "sum"), MTR=("MTR", "first"))
    lw[col] = lw["_n"] / lw["length"]

    um = gc.groupby("transcript", as_index=False).agg(
        **{col: (col, "mean"), "MTR": ("MTR", "first")})

    gene_lw = gc.groupby("gene_name", as_index=False).agg(
        _n=("_n", "sum"), length=("length", "sum"), MTR=("MTR", "first"))
    gene_lw[col] = gene_lw["_n"] / gene_lw["length"]

    gene_um = gc.groupby("gene_name", as_index=False).agg(
        **{col: (col, "mean"), "MTR": ("MTR", "first")})

    print(f"P(0) > {thr}")
    print(f"  per CDS interval                    R2 = {ols_r_squared(gc, col, 'MTR'):.4f}  n={len(gc):,}")
    print(f"  per transcript, length-weighted     R2 = {ols_r_squared(lw, col, 'MTR'):.4f}  n={len(lw):,}")
    print(f"  per transcript, unweighted mean     R2 = {ols_r_squared(um, col, 'MTR'):.4f}  n={len(um):,}")
    print(f"  per gene, length-weighted           R2 = {ols_r_squared(gene_lw, col, 'MTR'):.4f}  n={len(gene_lw):,}")
    print(f"  per gene, unweighted mean           R2 = {ols_r_squared(gene_um, col, 'MTR'):.4f}  n={len(gene_um):,}")
    print()
print("target: letter says 0.300 (from 0.146); figure label says 0.152")
