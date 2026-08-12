#!/bin/bash
#SBATCH --job-name=gate01
#SBATCH --partition=mrivas
#SBATCH --time=04:00:00
#SBATCH --mem=192G
# Pin the CPU vendor. numpy dispatches a different SIMD kernel for np.log on
# AVX-512 (Intel Skylake-X) than on AVX2 (AMD Rome), and the two differ in the
# last bit. effect_size then differs on ~18% of rows at up to 3.3e-15, so a
# bit-identity gate passes or fails depending purely on which node Slurm picked.
# The reference artifacts were built on Intel.
#SBATCH --constraint=CPU_MNF:INTEL
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Regenerates the epilepsy regression input through the current code in
# reproduce mode and compares DECOMPRESSED content against the run that was
# already verified bit-identical to the Revision 2 submission.
#
# Compare decompressed, never the .gz: gzip stores a timestamp in its header,
# so identical data yields different compressed bytes.
set -euo pipefail
cd /oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm
PY=/oak/stanford/groups/mrivas/users/osthoag/anaconda3/bin/python
export PYTHONUNBUFFERED=1
$PY pipelines/01_build_regression_input.py --no-fix --out-suffix _REGATE
echo
A=$(zcat results/constraint_gerp_am_epi25_variants_REPRO.tsv.gz  | md5sum | cut -d' ' -f1)
B=$(zcat results/constraint_gerp_am_epi25_variants_REGATE.tsv.gz | md5sum | cut -d' ' -f1)
echo "  reference (_REPRO) : $A"
echo "  current   (_REGATE): $B"
if [ "$A" = "$B" ]; then
  echo "GATE 01: PASS"
else
  echo "GATE 01: FAIL"
  echo "  first differing lines:"
  diff <(zcat results/constraint_gerp_am_epi25_variants_REPRO.tsv.gz) \
       <(zcat results/constraint_gerp_am_epi25_variants_REGATE.tsv.gz) | head -4
  exit 1
fi
