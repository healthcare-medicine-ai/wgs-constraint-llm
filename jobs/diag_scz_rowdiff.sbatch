#!/bin/bash
#SBATCH --job-name=sczdiff
#SBATCH --partition=mrivas
#SBATCH --time=02:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
export PYTHONUNBUFFERED=1
/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python - <<'PY'
import pandas as pd
R = "results/"
key = ["chr","pos","ref","alt","gene_id","group"]
old = pd.read_csv(R+"constraint_gerp_am_scz_variants.tsv.gz", sep="\t", usecols=key+["prob_0"])
new = pd.read_csv(R+"constraint_gerp_am_scz_variants_REPRO.tsv.gz", sep="\t", usecols=key+["prob_0"])
print(f"cached (2025-09-15) : {len(old):,} rows")
print(f"regenerated         : {len(new):,} rows   ({len(new)-len(old):+,})")

m = old[key].merge(new[key], on=key, how="outer", indicator=True)
vc = m["_merge"].value_counts()
print(f"\n  in both        : {vc.get('both',0):,}")
print(f"  cached only    : {vc.get('left_only',0):,}")
print(f"  new only       : {vc.get('right_only',0):,}")

extra = m[m["_merge"]=="right_only"]
if len(extra):
    print("\n  what are the extra rows?")
    print("   by chromosome (top 6):")
    print(extra["chr"].value_counts().head(6).to_string())
    print("   by group:")
    print(extra["group"].value_counts().to_string())
    print(f"   distinct genes affected: {extra['gene_id'].nunique():,}")
missing = m[m["_merge"]=="left_only"]
if len(missing):
    print(f"\n  rows present in the cache but not regenerated: {len(missing):,}")
    print(missing["chr"].value_counts().head(5).to_string())
PY
