"""Is the R2 response letter's "R2 improved from 0.146 to 0.300" reproducible?

The letter tells the reviewer Figure 3b was changed from per-gene to
per-transcript. The submitted figure is still titled "per gene" at P(0) > 0.6
with a hardcoded R2 = 0.152. This computes every granularity so the discrepancy
can be settled with a number rather than an argument.
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
                         "strand", "frame", "attribute"],
                  dtype={"start": int, "end": int})
gtf["gene_type"] = gtf["attribute"].str.extract(r'gene_type "(.*?)"')
gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "(.*?)"')
gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "(.*?)"')
gtf["transcript"] = gtf["transcript_id"].str.split(".").str[0]
gtf = gtf[(gtf["gene_type"] == "protein_coding") & (gtf["feature"] == "CDS")]
gtf = gtf.drop(columns=["attribute"])

gc = load_gene_constraint(gtf, cfg.data_dir / "gnomad.v4.0.constraint_metrics.tsv")
print(f"CDS intervals joined to gnomAD: {len(gc):,}", flush=True)

pred = pd.read_csv(cfg.result("HMM_rgc_ALL_RS_merged_predictions_ASPUBLISHED.tsv.gz"),
                   sep="\t", usecols=["chr", "pos", "prob_0"])
print(f"positions: {len(pred):,}", flush=True)
gc = add_constraint_proportions(gc, pred, thresholds=(0.5, 0.6))

print()
print("=" * 62)
print("PER CDS INTERVAL  (what the submitted figure actually plots)")
print("=" * 62)
for t in (50, 60):
    print(f"  P(0) > 0.{t//10}   R2 vs MTR = {ols_r_squared(gc, f'proportion_over_{t}', 'MTR'):.4f}"
          f"   (n={gc[['proportion_over_'+str(t),'MTR']].dropna().shape[0]:,})")
print(f"  P(0) > 0.5   R2 vs mis.z_score = {ols_r_squared(gc, 'proportion_over_50', 'mis.z_score'):.4f}")

# Aggregate every CDS interval of a transcript into one row: constrained bases
# over total coding bases.
print()
print("=" * 62)
print("PER TRANSCRIPT  (what the caption and response letter describe)")
print("=" * 62)
gc = gc.copy()
for t in (50, 60):
    gc[f"nconstr_{t}"] = gc[f"proportion_over_{t}"] * gc["length"]
agg = gc.groupby("transcript", as_index=False).agg(
    nconstr_50=("nconstr_50", "sum"), nconstr_60=("nconstr_60", "sum"),
    length=("length", "sum"), MTR=("MTR", "first"),
    mis_z=("mis.z_score", "first"))
for t in (50, 60):
    agg[f"proportion_over_{t}"] = agg[f"nconstr_{t}"] / agg["length"]
print(f"  transcripts: {len(agg):,}")
for t in (50, 60):
    print(f"  P(0) > 0.{t//10}   R2 vs MTR = {ols_r_squared(agg, f'proportion_over_{t}', 'MTR'):.4f}"
          f"   (n={agg[['proportion_over_'+str(t),'MTR']].dropna().shape[0]:,})")
print(f"  P(0) > 0.5   R2 vs mis.z_score = {ols_r_squared(agg, 'proportion_over_50', 'mis_z'):.4f}")
print()
print("  submitted figure says: per gene, P(0) > 0.6, R2 = 0.152 (hardcoded)")
print("  response letter says : per transcript, R2 = 0.300 (was 0.146)")
