#!/bin/bash
#SBATCH --job-name=viafig
#SBATCH --partition=mrivas
#SBATCH --time=08:00:00
#SBATCH --mem=192G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Proves the claim made to Genentech: that the published constraint+GERP table
# plus the public inputs rebuilds the prepared variant table exactly, so they
# never need the controlled-access file or the bigWig.
set -uo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/01_build_regression_input.py \
    --constraint-gerp results/HMM_rgc_ALL_RS_merged_predictions.tsv.gz \
    --out-suffix _VIAFIG || exit 1
echo
echo "=== decompressed md5: via published table vs via bigWig ==="
A=$(zcat results/constraint_gerp_am_epi25_variants_FIXED.tsv.gz | md5sum | cut -d' ' -f1)
B=$(zcat results/constraint_gerp_am_epi25_variants_VIAFIG.tsv.gz | md5sum | cut -d' ' -f1)
echo "  via bigWig        : $A"
echo "  via published tbl : $B"
if [ "$A" = "$B" ]; then echo; echo "IDENTICAL -- the published route reproduces the requested file"; exit 0; fi
echo; echo "DIFFER -- investigating"
$PY - <<'PY'
import pandas as pd
a = pd.read_csv("results/constraint_gerp_am_epi25_variants_FIXED.tsv.gz", sep="\t")
b = pd.read_csv("results/constraint_gerp_am_epi25_variants_VIAFIG.tsv.gz", sep="\t")
print(f"  rows {len(a):,} vs {len(b):,}   cols equal: {list(a.columns)==list(b.columns)}")
PY
exit 1
