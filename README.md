# wgs-constraint-llm

Code for **"A unified meta-regression model identifies genes associated with
epilepsy"** (Aguilar, Rivas & Rivas), submitted to the *American Journal of
Human Genetics* as AJHG-D-24-00613.

An HMM is applied to whole-exome sequencing data (RGC-ME, All of Us) to predict,
for each position in the genome, the probability of observing no variant across
the population. Those constraint predictions are then combined with GERP RS,
AlphaMissense pathogenicity, and pLoF/missense annotations in a weighted
least-squares meta-regression that tests rare-variant burden per gene.

---

## Layout

```
config.yaml            all filesystem paths and analysis parameters
requirements.txt       exact pinned versions
src/wgs_constraint/    importable logic shared by notebooks and pipelines
  config.py              path resolution, env-var overrides
  gerp.py                GERP RS annotation  <- coordinate conventions live here
  alphamissense.py       AlphaMissense loading and per-variant collapse
  epi25.py               Epi25 loading, filters, effect sizes
  metareg.py             shared WLS meta-regression
  gene_constraint.py     per-gene constrained fractions
pipelines/             end-to-end stages, runnable non-interactively
jobs/                  Slurm wrappers, one per stage or gate
tests/                 pytest suite; no cluster data required
notebooks/             ARCHIVED implementations -- read, do not run
docs/                  overhaul plan, input manifest with checksums
genentech_fix/         August 2026 correction workspace and diagnostics
results/               gene-level outputs and figures
pyproject.toml         packaging, pytest and ruff configuration
```

Anything that produces a number in the manuscript runs from `pipelines/`.
Notebooks import from `src/` rather than carrying their own copies — the GERP
annotation was previously pasted into three of them, which is how one
coordinate error became three.

## Reproducing the epilepsy analysis

```bash
pip install -r requirements.txt

# point at your own data if you are not in the Rivas Lab OAK space
export WGS_DATA_DIR=/path/to/data
export WGS_RESULTS_DIR=/path/to/results

python pipelines/01_build_regression_input.py     # variant-level input
python pipelines/02_run_metaregression.py         # per-gene WLS + tables
```

On Slurm:

```bash
sbatch --job-name=fixed jobs/run_pipeline.sbatch fixed
```

`genentech_fix/03_compare_runs.py` diffs any run against another;
`genentech_fix/04_attribution.py` isolates the effect of each correction. Both
need the Revision 2 baseline snapshot, which stays on the cluster.

### Every stage

`pipelines/` is numbered by dependency order, not contiguously -- 03, 04 and 07
were folded into other stages or into `genentech_fix/` and the numbers were left
alone so existing references keep resolving.

| Stage | Produces |
|---|---|
| `00_input_manifest.py` | `docs/INPUT_MANIFEST.md`; `--verify` rechecks every input |
| `01_build_regression_input.py` | Epi25 variant-level regression input |
| `02_run_metaregression.py` | per-gene WLS p-values |
| `05_build_constraint_gerp_predictions.py` | constraint + GERP prediction table |
| `06_epilepsy_tables_and_figures.py` | Tables 1 and 2, Figure 5 |
| `08_hmm_gerp_joint_distribution.py` | Figures 4a, 4b |
| `09_gene_constraint_figures.py` | Figures 3a, 3b |
| `10_rgc_aou_joint.py` | Figures 2a, 2b |
| `11_schizophrenia.py` | SCHEMA replication, Table A2, Figure A1 |
| `12_figure1_scn1a.py` | Figure 1 |

Every stage takes `--no-fix` or an equivalent published-fidelity flag, so the
Revision 2 output can be regenerated for comparison. `jobs/gate_*.sbatch` runs
each against its published artifact; all seven published figures reproduce
pixel-identically.

### Reproducing the *submitted* Revision 2 numbers

```bash
python pipelines/01_build_regression_input.py --no-fix --out-suffix _REPRO
python pipelines/02_run_metaregression.py --input-suffix _REPRO
```

This reproduces the Revision 2 output bit-for-bit across all 72,386 gene×group
models. It exists so that any difference in a corrected run is attributable to
the corrections rather than to refactoring. Tag `AJHG-R2-submitted` marks the
tree as submitted.

## Configuration

No path is hardcoded. Resolution order:

1. `WGS_DATA_DIR` / `WGS_RESULTS_DIR`
2. the file named by `WGS_CONFIG`
3. `config.yaml` at the repository root

`config.yaml` also holds the analysis parameters — the inclusion floor, the
significance thresholds, the allele-count cut — so they are version-controlled
rather than retyped per notebook cell.

## Analysis choices worth knowing

These are the points external groups have most often got wrong.

| | |
|---|---|
| Rare-variant cut | **allele count**, `ac_case + ac_ctrl <= 5` — *not* a MAF filter |
| Effect size | Haldane–Anscombe corrected log odds ratio from allele counts |
| Variance | Woolf variance; weights are `1 / var_effect_size` |
| Firth BETA/SE | **not used and not required** |
| Moderators | five: HMM constraint, GERP RS, AlphaMissense, pLoF, missense |
| Transforms | `-log1p(-x)` on constraint and pathogenicity; **GERP enters untransformed** |
| Gene inclusion | `n_variants >= 25` |
| Reported p-value | `model.f_pvalue`, the **whole-model F-test** — not a coefficient p-value |
| Exome-wide significance | `p < 3.4e-7` |

**This repository is not `rivas-lab/phenome-wide-unified-model`.** That project
applies a related four-moderator model (no GERP) with `MAF <= 0.05` and
pre-computed `BETA`/`SE` columns to phenome-wide biobank data. Its
`unified_reg_MAF.05.py` will not reproduce the epilepsy results, and applying it
to Epi25 yields substantially more significant genes because `MAF <= 0.05`
admits variants roughly three orders of magnitude more common than `AC <= 5`.

## Corrections, August 2026

Two defects were found in the Revision 2 pipeline during an external
reproduction attempt and corrected. Both are documented in
`genentech_fix/README.md`, with per-gene attribution.

1. **GERP coordinate off-by-one.** `pyBigWig.values()` returns a 0-based vector;
   HMM positions are 1-based. Indexing one with the other shifted every GERP
   score one base downstream. Old and corrected values correlate at r = 0.173 —
   within coding sequence, adjacent bases differ sharply in constraint, so the
   shift randomises rather than blurs. Fixed in `src/wgs_constraint/gerp.py`,
   where the coordinate convention is now an explicit, validated argument.

2. **AlphaMissense transcript fan-out.** `AlphaMissense_hg38.tsv.gz` is the
   isoform-level release. Merged without de-duplication, a variant annotated on
   N transcripts became N rows — inflating `n_variants` *and* its weight in the
   WLS. Fixed in `src/wgs_constraint/alphamissense.py` by collapsing to one row
   per variant by maximum pathogenicity.

## Tests

```bash
pytest -q
```

Runs in seconds against synthetic fixtures; no cluster data needed --
`pyproject.toml` puts `src/` on the path. GitHub Actions runs the suite and
`ruff check` on every push and pull request (`.github/workflows/ci.yml`). CI
cannot verify the analysis itself: every real input is either controlled-access
or hundreds of gigabytes. Reproduction gates run on Sherlock. The GERP
tests build a small bigWig whose value at each base is known, so a one-base
shift is arithmetic rather than judgement — `test_the_two_conventions_differ_by_
exactly_one_base` fails against the pre-correction code.

## Data

Inputs are not redistributed here.

| Input | Source |
|---|---|
| Epi25 variant results | https://epi25.broadinstitute.org/results |
| GENCODE v44 basic annotation | GENCODE release 44 |
| AlphaMissense (hg38, isoform-level) | Zenodo record 10813168 |
| Genes4Epilepsy | https://github.com/bahlolab/Genes4Epilepsy |
| WES constraint predictions (AoU, RGC-ME) | https://doi.org/10.6084/m9.figshare.27184245.v1 |
| GERP RS, hg38 | lifted from the hg19 track with UCSC `hg19ToHg38.over.chain.gz`; see `genentech_fix/README.md` |

Variant-level intermediates derived from Epi25 are controlled-access and are
excluded by `.gitignore`.
