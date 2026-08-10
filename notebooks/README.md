# Archived notebooks

These are the original implementations, kept for the record. **They are not the
source of truth.** Anything that produces a number or figure in the manuscript
is being moved to `pipelines/`, verified against the artifact it replaces. See
`docs/OVERHAUL_PLAN.md` for the inventory and progress.

Read them to understand how something was originally done. Do not run them to
produce a result that will be published.

## Why they are not authoritative

Two defects found in August 2026 illustrate the problem these notebooks create:

- The **GERP coordinate off-by-one** existed because the same bigWig lookup was
  pasted into three notebooks. Fixing one would not have fixed the others.
  All three now call `wgs_constraint.annotate_with_gerp` instead.
- The **AlphaMissense transcript fan-out** is fixed in
  `wgs_constraint.alphamissense`, which `pipelines/01` uses. The merge inside
  `Epilepsy Analysis.ipynb` cell 16 is **not** fixed.

Cell outputs have been cleared. They were rendered by the defective code, and
leaving them would present superseded numbers as authoritative.

## Per-notebook status

| Notebook | Superseded by | Still authoritative for |
|---|---|---|
| `Epilepsy Analysis.ipynb` | `pipelines/01`, `02`, `06` — fully superseded | nothing |
| `Constraint Measures Comparison.ipynb` | `pipelines/05` for the constraint+GERP table | Figures 3a/3b, 4a/4b |
| `Schizophrenia Analysis.ipynb` | — | everything; extraction pending |
| `HMM Mutation Predictions.ipynb` | — | Figure 1, and the HMM predictions themselves |
| `RGC + AoU Predictions.ipynb` | — | Figures 2a/2b |
| `AoU Mutation Predictions.ipynb` | — | AoU predictions |
| `Constraint + AM Analysis.ipynb` | — | exploratory |
| `SCZ Liftover.ipynb` | — | SCZ hg19→hg38 liftover |

⚠️ **`Epilepsy Analysis.ipynb` cells 15–22 duplicate `pipelines/01` and `02`.**
The notebook path carries the uncorrected AlphaMissense merge, so it and the
scripts can disagree. Use the scripts.

## Running one anyway

They expect to be run from this directory, and insert `../src` on the path so
they can import the shared modules.

```bash
cd notebooks
jupyter lab      # or via https://ondemand.sherlock.stanford.edu
```

The state exactly as submitted for Revision 2 — before any correction, with
original outputs intact — is tagged `AJHG-R2-submitted`:

```bash
git show AJHG-R2-submitted:"Epilepsy Analysis.ipynb" > /tmp/as_submitted.ipynb
```
