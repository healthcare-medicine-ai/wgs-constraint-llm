# wgs-constraint-llm

Code for "A unified meta-regression model identifies genes associated with epilepsy"
(Aguilar, Rivas and Rivas; *American Journal of Human Genetics*, manuscript AJHG-D-24-00613).

A two-state hidden Markov model (HMM) is fit to rare-variant presence in the RGC Million Exome
data to estimate, for each coding position, the probability that it is constrained. These
predictions are combined with GERP RS scores, AlphaMissense pathogenicity and pLoF/missense
annotations in a weighted least-squares meta-regression of Epi25 variant effect sizes, giving
one test per gene and epilepsy group.

## Layout

```
config.yaml          paths and analysis parameters
DATA.md              where each input comes from, with checksums
src/wgs_constraint/  shared code used by the pipeline stages
pipelines/           one script per step of the analysis, run in order
tests/               unit tests on synthetic data
results/             gene-level results and the figures and tables in the paper
dashboard/           Streamlit browser for the gene-level results
```

## Setup

```bash
pip install -r requirements.txt
```

Inputs are not redistributed here. Download them as described in `DATA.md` into `data/`, or
point `WGS_DATA_DIR` at a directory that has them. Outputs go to `results/` (override with
`WGS_RESULTS_DIR`).

## Running

```bash
bash run_all.sh
```

or one stage at a time:

| Stage | Produces |
|---|---|
| `01_hmm_constraint.py` | HMM trained on RGC-ME chromosome 2; constraint predictions for RGC-ME and All of Us |
| `02_regression_input.py` | constraint and GERP for every coding position; Epi25 variants annotated with every moderator |
| `03_metaregression.py` | p-values for every gene and epilepsy group |
| `04_figure1.py` | Figure 1 |
| `05_rgc_aou_validation.py` | Figure 2; RGC variant status vs HMM constraint as predictors of All of Us variants |
| `06_constraint_comparison.py` | Figure 3 |
| `07_epilepsy_results.py` | Tables 1, 2 and S1; Figure 4 |
| `08_moderator_ablations.py` | Figures S1–S4 |
| `09_schizophrenia.py` | Table 3 and Figure 5 |

Stages 01, 02 + 03, 08 and 09 each take 40–45 minutes on one node; the others take a few
minutes. Stage 02 needs about 128 GB of memory, the rest 32–64 GB. Stage 01 loads the trained
model from `results/` if it is there. Training gives bit-identical parameters only on the same
CPU generation, so use the saved model to reproduce the published predictions exactly.

## Tests

```bash
pytest -q
```

## Earlier versions

Tag `ajhg-r4-submitted` is the code as linked from the submitted manuscript, including the
original notebooks and the scripts used to check the August 2026 corrections. Tag
`AJHG-R2-submitted` is the second revision.

## Licence

MIT. See `LICENSE`.
