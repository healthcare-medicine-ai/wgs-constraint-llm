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

`pipelines/` is numbered by dependency order, not contiguously -- 03 and 04 live
in `genentech_fix/` as correction-specific diagnostics, and the numbers were left
alone so existing references keep resolving.

| Stage | Produces |
|---|---|
| `00_input_manifest.py` | `docs/INPUT_MANIFEST.md`; `--verify` rechecks every input |
| `01_build_regression_input.py` | Epi25 variant-level regression input |
| `02_run_metaregression.py` | per-gene WLS p-values |
| `05_build_constraint_gerp_predictions.py` | constraint + GERP prediction table |
| `06_epilepsy_tables_and_figures.py` | Tables 1 and 2, Figure 5 |
| `07_moderator_ablations.py` | leave-one-out ablations, Supplementary Figures S1-S4 |
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

Inputs are not redistributed here. Every one is checksummed in
[`docs/INPUT_MANIFEST.md`](docs/INPUT_MANIFEST.md); verify a local copy with
`python pipelines/00_input_manifest.py --verify`.

| Input | File | Source |
|---|---|---|
| Epi25 variant results | — | https://epi25.broadinstitute.org/results |
| GENCODE v44, basic annotation | `gencode.v44.basic.annotation.gtf.gz` | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz |
| AlphaMissense, hg38, isoform-level | `AlphaMissense_hg38.tsv.gz` | https://console.cloud.google.com/storage/browser/dm_alphamissense |
| Genes4Epilepsy | — | https://github.com/bahlolab/Genes4Epilepsy |
| gnomAD v4.0 constraint metrics | `gnomad.v4.0.constraint_metrics.tsv` | https://gnomad.broadinstitute.org/downloads |
| SCHEMA variant and gene results | `SCHEMA_variant_results_hg38.tsv.gz`, `SCHEMA_gene_results.tsv.bgz` | https://schema.broadinstitute.org/ |
| GERP RS, hg19 (liftover source) | `All_hg19_RS.bw` | https://genome-asia.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=compGeno&hgta_track=allHg19RS_BW&hgta_table=allHg19RS_BW |
| UCSC chain, hg19 to hg38 | `hg19ToHg38.over.chain.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz |
| **GERP RS, hg38 (lifted here)** | `All_hg38_RS.bw` | https://doi.org/10.6084/m9.figshare.33201549.v1 |
| WES constraint predictions, AoU and RGC-ME | — | https://doi.org/10.6084/m9.figshare.27184245.v1 |

**Which GENCODE file matters.** Release 44 ships several GTFs and they give
different CDS interval sets. This analysis uses the *basic* annotation,
29,570,410 bytes, md5 `7450ef42cf9cb3d29625320b22d4bb45` — confirmed by byte
match against the URL above. The `primary_assembly` file at the same release is
49,730,393 bytes and is **not** what was used, and some notebook comments cite a
`chr_patch_hapl_scaff` URL, which is a third, different file.

**The hg38 GERP track was produced here, not downloaded.** UCSC publishes the
hg19 track only. `All_hg38_RS.bw` (16.4 GB, md5
`10a03ee7969d3000ffd2e6e6f84f2453`) was lifted from `All_hg19_RS.bw` with the
chain above, using `bigWigLiftOver` with no options. Its positional correctness
was verified by sampling: of 1,500 informative positions, 1,498 carry the hg19
value at the chain-mapped hg38 coordinate at offset 0, and none at plus or minus
one or two bases (`genentech_fix/verify_liftover.py`).

Variant-level intermediates derived from Epi25 are controlled access and are
excluded by `.gitignore`. They cannot be redistributed; `pipelines/01` rebuilds
them from the sources above.
