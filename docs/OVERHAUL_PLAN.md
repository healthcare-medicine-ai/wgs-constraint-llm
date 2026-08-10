# Codebase overhaul — plan and progress

**Goal.** Every artifact the manuscript depends on is produced by a documented,
runnable, version-controlled pipeline stage, and every stage is verified to
reproduce the artifact it replaces. Notebooks are retained as archived
implementations, not as the source of truth.

Started 2026-08-09, during the Revision 3 correction work.

---

## The validation rule

No stage is "done" until it has passed a **fidelity gate**: the extracted code
reproduces the existing artifact exactly, before any correction is applied.
Corrections are then applied on top, so any change in output is attributable to
the correction rather than to the extraction.

Two forms of gate, depending on cost:

**Full gate** — re-run the stage end to end and compare every value. Used where
the stage is cheap enough (the epilepsy meta-regression: ~30 min).

**Slice gate** — re-run the stage over a deterministic slice (a chromosome, a
fixed row range) and compare against the *same slice* of the archived output.
Used where full regeneration is expensive. A slice gate is a real proof of
equivalence for deterministic, row-independent transformations; it is *not*
sufficient for stages with global state (a fit over all rows, a normalisation,
a sort-dependent operation). Record which kind was used.

**Pinned artifact** — for stages that are stochastic or take hours (HMM
training), the output is treated as a fixed input: checksummed, its provenance
documented, and never re-derived. Only the deterministic transformation
downstream of it is gated. Do not regenerate the figshare predictions.

Every gate is a test under `tests/` or a script under `pipelines/`, not a
one-off command in someone's history.

---

## Inventory

Manuscript artifacts and the code that produces them. Status as of the last
edit to this file.

| Artifact | Source notebook | Target stage | Gate | Status |
|---|---|---|---|---|
| Gene-level p-values | Epilepsy Analysis | `pipelines/01`+`02` | full | **done** — bit-identical, 72,386 models |
| Table S1 | Epilepsy Analysis | `pipelines/02` | full | **done** — 146/146 rows |
| Table A1 comparison | Epilepsy Analysis | `pipelines/02` | full | **done** |
| Tables 1 and 2 | Epilepsy Analysis | `pipelines/06` | full | **done** — gate: 28 pairs / 24 genes / STX1B |
| Figure 5 (→ Fig 4) | Epilepsy Analysis | `pipelines/06` | data | **done** — underlying data exported alongside |
| Constraint + GERP table | Constraint Measures Comparison | `pipelines/05` | full | **done** — md5 of decompressed output identical to published |
| Figure 3a/3b | Constraint Measures Comparison | `pipelines/08` | slice | todo |
| Figure 4a/4b (→ Fig 3C/3D) | Constraint Measures Comparison | `pipelines/08` | slice | todo |
| Figure 2a/2b | RGC + AoU Predictions | `pipelines/09` | slice | todo |
| Figure 1 | HMM Mutation Predictions | `pipelines/10` | slice | todo |
| HMM constraint predictions | HMM Mutation Predictions | pinned artifact | checksum | todo |
| Schizophrenia results + Figure 6 | Schizophrenia Analysis | `pipelines/11` | full | todo |
| SCZ hg38 liftover | SCZ Liftover | `pipelines/12` | slice | todo |
| AoU predictions | AoU Mutation Predictions | pinned artifact | checksum | todo |
| Constraint + AlphaMissense | Constraint + AM Analysis | fold into `src/` | slice | todo |

### Known blockers and defects

**`pipelines/05` is gated.** The `--no-fix` rebuild reproduces the published
`HMM_rgc_ALL_RS_merged_predictions.tsv.gz` exactly: identical md5 of the
decompressed content (c56b4f5d376cfc5d29975a9828db6848), 28,922,258 rows, identical
header. An earlier note in this file claimed a 538 KB discrepancy and that the
published file had been overwritten; both were wrong. The published file lives at the
repository root (the notebook wrote it there via a relative path) while the pipeline
writes to `results/`, so nothing was overwritten, and the 538 KB figure came from
comparing against an unrelated file.

**`Epilepsy Analysis.ipynb` still reimplements the pipeline** in cells 15, 16,
19, 20 and 22, and its cell 16 still carries the *uncorrected* AlphaMissense
merge. Until those cells are removed or superseded, the repository can produce
two different answers depending on whether the notebook or the script is run.
Highest-priority item.

**Correction status.** The GERP off-by-one is fixed in `src/wgs_constraint/
gerp.py` and all three notebooks call it. The AlphaMissense fan-out is fixed in
`src/wgs_constraint/alphamissense.py` and used by `pipelines/01`, but the
notebook merge paths are unfixed.

---

## Conventions

- Analysis parameters live in `config.yaml`, never inline.
- No filesystem path is hardcoded; resolve through `wgs_constraint.get_config()`.
- Logic used by more than one stage goes in `src/wgs_constraint/`.
- Coordinate conventions are explicit arguments, never assumptions in comments.
  The GERP defect existed because a comment asserted "both datasets are 0-based"
  and nothing checked it.
- Anything expensive runs via `jobs/*.sbatch`, not an interactive session.
- Outputs are suffixed by run configuration; nothing overwrites a published
  artifact in place. Checksum before overwriting anything that cannot be
  regenerated cheaply.

## Notebooks

Retained under `notebooks/` as archived implementations. They document how the
analysis was originally done and are not the source of truth. Once a stage is
extracted and gated, the corresponding notebook cells are superseded — see
`notebooks/README.md` for per-notebook status.

The exact state behind the Revision 2 submission is tagged `AJHG-R2-submitted`.

## Pinned artifact checksums

Decompressed-content md5, so gzip metadata differences do not confuse comparisons.

| Artifact | md5 (decompressed) |
|---|---|
| `HMM_rgc_ALL_RS_merged_predictions.tsv.gz` (as published, repo root) | `c56b4f5d376cfc5d29975a9828db6848` |
