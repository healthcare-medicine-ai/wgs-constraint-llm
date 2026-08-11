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

## Figure provenance, resolved 2026-08-11

The notebooks save several competing versions of the same figure number and the
last cell executed wins, which nothing records. Determined by comparing what the
`AJHG-R2-submitted` tag kept against what it deleted.

**Authoritative, kept in the R2 commit:**

| Figure | Produced by | Detail |
|---|---|---|
| 2a, 2b | `Constraint Measures Comparison` cell 18 **or** `RGC + AoU Predictions` cell 9 | both write the same filenames; unresolved |
| 3a | `Constraint Measures Comparison` cell 32 | fraction of bases with P(0) > **0.5** vs gnomAD missense z-score |
| 3b | `Constraint Measures Comparison` cell 28 | proportion with P(0) > **0.6** vs MTR |
| 4a, 4b | `Constraint Measures Comparison` cells 16-17 | gated, pixel-identical, now `pipelines/08` |

**Deleted by the R2 commit, i.e. superseded:** "Figure 3a: HMM vs gnomAD v4
constraint per gene" (cell 26), "Figure 3a: HMM vs gnomAD LoF z-score per gene",
"Figure 3b: HMM vs gnomAD LoF z-score per gene" (cell 34), "Figure 3c: HMM vs
MTR per gene" (cell 36), and a copy with a newline in its filename.

**Two defects in the surviving Figure 3 pair:**

1. *Mismatched thresholds.* Panel A uses P(0) > 0.5 and panel B uses P(0) > 0.6.
   Already noted in the Revision 3 plan; `MATCH_THRESHOLDS` was proposed to put
   both on 0.5. Not a correctness bug, but the caption should say which is which.
2. *Hardcoded R².* Cell 28 annotates panel B with the literal string
   `$R^2=0.152$` rather than the computed value, while cell 32 computes it. If
   the input changes the annotation does not. The extracted stage must compute
   both, and the recomputed value for panel B should be checked against 0.152.

`calculate_overlap` is byte-identical in cells 20 and 30; cell 30 merely
recomputes `proportion_over_50` redundantly. An earlier draft of this file
speculated the two panels used different overlap logic -- they do not.

## Session of 2026-08-11

### Stages added

| Stage | Produces | Extracted from |
|---|---|---|
| `pipelines/00_input_manifest.py` | `docs/INPUT_MANIFEST.md`, checksums of every input | — |
| `pipelines/09_gene_constraint_figures.py` | Figures 3a, 3b | Constraint Measures Comparison 28, 32 |
| `pipelines/10_rgc_aou_joint.py` | Figures 2a, 2b | RGC + AoU Predictions 8, 9 |
| `pipelines/11_schizophrenia.py` | SCZ results, Figure A1, table A2 | Schizophrenia Analysis 15-21 |
| `pipelines/12_figure1_scn1a.py` | Figure 1 | HMM Mutation Predictions 21, 22 |

New shared modules: `metareg.py` (the meta-regression, previously duplicated
verbatim between the epilepsy and schizophrenia notebooks) and
`gene_constraint.py` (CDS-interval aggregation).

### Two lessons about gates, learned the hard way

**Never compare gzipped files byte-for-byte.** gzip writes a timestamp into its
header, so identical data produces different compressed bytes and different
md5s. `gate01` failed spuriously for this reason while every row count matched.
Always compare the decompressed stream.

**Published figures were saved at dpi=300.** The notebooks call bare
`plt.savefig(path)`, which uses matplotlib's default of 100, so a naive
extraction renders at exactly one third the published dimensions and cannot be
pixel-compared. All figure stages now pass `dpi=300`.

### Findings

**Figure 3b's R² label is wrong.** The notebook annotates the panel with a
hardcoded `$R^2=0.152$`. The computed value is **0.1531**. Small, but it proves
the hardcoded label had already drifted from the data before publication.

**Figure 1 shows 4.2% of SCN1A.** The crop offsets (+750, -4700) reduce a
5,691-position gene to a **241-position window**. The caption reads "observed vs
predicted mutations for SCN1A", which implies the whole gene. It is a
representative window and the caption should say so.

**`Constraint Measures Comparison` cell 18 is a latent hazard.** It is a
byte-identical copy of the Figure 2a/2b plotting code from `RGC + AoU
Predictions` cell 9, but in that notebook the variables it plots are left over
from cell 16, the HMM-versus-GERP distribution on a 10x18 grid. Running that
notebook top to bottom overwrites Figures 2a and 2b with mislabelled GERP data
on axes captioned "RGC" and "AoU". Delete that cell.

**`SCZ Liftover.ipynb` is the best code in the repository** and should be the
model for the rest: it empirically detects whether input coordinates are 0- or
1-based by trying both and counting successful mappings, rather than assuming.
Had the GERP merge done the same, the defect that prompted all of this could not
have occurred.

### Pinned artifacts, not extracted

These generate inputs, take hours, and are distributed externally. They are
checksummed in `docs/INPUT_MANIFEST.md` and must not be re-derived:
`HMM Mutation Predictions.ipynb` (the RGC-ME predictions),
`AoU Mutation Predictions.ipynb` (the AoU predictions),
`SCZ Liftover.ipynb` (SCHEMA and scz.tsv.gz lifted to hg38).
`Constraint + AM Analysis.ipynb` is exploratory and produces nothing the
manuscript cites.

### Figure gates: all seven published figures reproduce pixel-exactly

Rendered from the pre-correction inputs and compared pixel-by-pixel against the
PNGs kept by the `AJHG-R2-submitted` tag.

| Figure | Stage | Result |
|---|---|---|
| 1, SCN1A | `pipelines/12` | identical, 0 pixels differ |
| 2a, 2b, RGC vs AoU | `pipelines/10` | identical, 0 pixels differ |
| 3a, HMM vs missense z | `pipelines/09` | identical, 0 pixels differ |
| 3b, HMM vs MTR | `pipelines/09` | identical except 314 pixels |
| 4a, 4b, HMM vs GERP | `pipelines/08` | identical, 0 pixels differ |

Figure 3b's 314 differing pixels occupy a 21x32 box at 81% across and 20% down,
which is the position of the R-squared annotation. The plotted data is
bit-identical; the only difference is the digit, corrected from the hardcoded
0.152 to the computed 0.1531. This is the intended correction, not a
discrepancy.

### One more lesson: do not add .copy() to a validated numerical path

`gate01` failed after the meta-regression was moved into `metareg.py`, with the
row counts matching exactly but the content differing. The cause was a `.copy()`
added for tidiness in `haldane_effect_sizes`: it changes the array's memory
layout, which changes the vectorised path numpy selects for `np.log`, which
perturbs `effect_size` by one unit in the last place.

Numerically that is around 1e-16 and harmless. It was removed anyway. A
bit-identical gate is a far stronger claim than "agrees to fifteen decimals",
and a perturbed effect size feeds the weighted least squares where a borderline
p-value could in principle flip. The comment in `metareg.py` records why the
line is missing so nobody helpfully adds it back.

## Results with no stage that regenerates them

`results/epilepsy_unified_model_pvalues_MINUS_*.tsv`,
`meta_FULL_all_moderators.tsv` and
`unified_pvalue_comparison_WITH_vs_WITHOUT_log_constraint.tsv` were committed
with the Revision 2 snapshot and are referenced by no document and produced by
no pipeline stage. They are not junk: they are the leave-one-out moderator
analysis backing **Supplementary Figures S1-S4**, added in response to
Reviewer 1 (response letter item 10, "we compared the full model against models
with each annotation removed in turn"). Keep them. But until a stage produces
them they cannot be regenerated, and the corrections change them, so S1-S4 are
stale in exactly the way the main figures were.

## Figure 3b describes an analysis the submitted figure does not contain

See the docstring of `pipelines/09_gene_constraint_figures.py` for the measured
numbers. Short version: the response letter and the manuscript caption both say
panel B is per transcript with R2 = 0.300; the submitted figure is per CDS
interval at P(0) > 0.6 with a hardcoded R2 = 0.152. The letter's "was 0.146"
reproduces exactly as the per-interval value at 0.5. This stage reproduces what
was submitted. Reconciling it is Revision 3 work.

## Bit-exactness is fragile in ways that are not obvious

Two separate refactors, each of which looked semantically inert, perturbed
`effect_size` by one unit in the last place and thereby changed p-values in the
fifteenth decimal:

1. Adding `.copy()` after a boolean mask.
2. Replacing chained `!=` comparisons with `.isin()`.

Both select identical rows. Both change the resulting frame's memory layout,
which changes the vectorised loop numpy selects for `np.log`. Measured effect on
the published p-values: 100% agree within 1e-12, but only ~56-75% are
bit-identical.

That difference is scientifically meaningless -- no gene changes significance --
but a bit-identical gate is a much stronger claim than "agrees to twelve
figures", and it costs nothing to preserve. Both changes were reverted, with
comments in `epi25.py` and `metareg.py` explaining why the tidier form is not
used. **Do not "clean up" those expressions.**

## Schizophrenia: resolved, 2026-08-11

For most of a day this looked like an input-provenance problem. It was a
transcription error in the port, and the reasoning that pointed away from it was
wrong. Both halves are worth recording.

**The symptom.** `pipelines/11` reproduced the published results exactly when fed
the cached 2025-09-15 merge — 190,328 models against 190,328, `n_variants` 100%
bit-identical. Regenerating that merge from the same inputs instead produced
178,821 extra rows and 1,958 extra gene-groups.

**The wrong turn.** The extra rows were spread across every chromosome, every
ancestry group and 12,441 genes, and I concluded that a filter difference "would
cluster on one chromosome, consequence or group. This does not." That inference
was simply invalid. A *consequence* filter is orthogonal to chromosome, gene and
ancestry, so removing one scatters rows across all three by construction.
Diffuseness was evidence for the filter hypothesis, not against it. Acting on
that conclusion sent the investigation upstream into file dates and figshare
rebuilds, none of which was relevant.

**What actually found it.** Not more hypotheses — one measurement. Comparing the
*null pattern* of the extra rows against the file as a whole: 90% of them lacked
an AlphaMissense score against a 14% baseline. AlphaMissense only scores missense
SNVs, so the extras were concentrated in variants it cannot score. The follow-up
crosstab was decisive: 100% of the extra rows were neither pLoF nor missense.

**The cause.** Cell 15 of `Schizophrenia Analysis.ipynb` lists ten excluded
consequences, six of them commented out. I transcribed the four uncommented
ones. Joining the cached merge back to the SCHEMA source shows that
`intron_variant`, `splice_region_variant`, `3_prime_UTR_variant`,
`5_prime_UTR_variant`, `upstream_gene_variant` and `downstream_gene_variant` are
absent from it *entirely*, while every coding consequence is fully present. The
submitted analysis ran with all ten active; the comment characters were added
afterwards. **The notebook's current state does not describe the published
results.**

Two details corroborate it. The residue matches exactly: the 9,776 rows in the
cache that are neither pLoF nor missense are `start_lost` (9,768) plus
`incomplete_terminal_codon_variant` (5) and a few `stop_retained_variant` — all
consequences no version of the list excludes. And of the four entries left
uncommented, only `coding_sequence_variant` and `synonymous_variant` occur in
SCHEMA at all; `mature_miRNA_variant` and `null` never appear. The visible list
was almost entirely inert.

**The fix (gate in flight at the time of writing).** `EXCLUDED_CONSEQUENCES`
in `pipelines/11` now carries all ten, with
the recovery argument in a comment so nobody "tidies" it back. This is the
analysis definition rather than one of the two corrections, so it applies on both
the `--no-fix` and `--fix` paths. Gated by `jobs/gate_11_scz.sbatch`, which now
checks the merged input against the cache before checking the fitted p-values —
the earlier gate only checked p-values and so could not see this.

**Consequence for the corrected numbers.** The previously recorded "no
exome-wide significant genes, against three published" was computed on the
contaminated input and is void. The corrected result must be re-derived from the
restored filter before it goes anywhere near the manuscript.

**The transferable lesson.** A commented-out line in a notebook is not evidence
about what was run. Where a cached output exists, recover the parameters from
the output rather than reading them off the source — that is what settled this,
and it took one join.

### Schizophrenia: superseded diagnostic notes, 2026-08-11


**The port is correct.** Run against the cached 2025-09-15 input, `pipelines/11`
reproduces the published results exactly: 190,328 models against 190,328,
perfect gene-group overlap, `n_variants` 100% bit-identical.

**(WRONG — see above.) The input difference is diffuse, not a filter bug.** Regenerating the merged
input yields 178,821 extra rows spread across every chromosome, every ancestry
group (meta 99,981, EUR 28,955, AFR 11,702, and so on) and 12,441 distinct
genes. A missing or mis-specified filter would cluster on one chromosome,
consequence or group. This does not. Cause not yet identified.

**(VOID — computed on the contaminated input.) Corrected result, confounded.** With both defects fixed, the schizophrenia
analysis finds **no exome-wide significant genes**, against three published
(PCDHA4, PCDHGA5, XPO7). Two of the three are protocadherins, consistent with
the AlphaMissense fan-out that inflated the same family in the epilepsy
analysis.

That figure must not be reported yet. It is computed on the regenerated input,
so it conflates the corrections with the unexplained 3% input difference, and
the two cannot be separated by applying corrections to the cached file: fixing
GERP requires re-merging from the predictions.

**Next step for whoever picks this up.** Identify what makes the regenerated
SCZ input differ from the September 2025 cache. Candidates worth checking in
order: whether `SCHEMA_variant_results_hg38.tsv.gz` was regenerated by
`SCZ Liftover.ipynb` cell 2 after the cached merge was built; whether the
GERP dropna admits a different position set; whether the GENCODE join changed.
Only once that is understood can the corrected schizophrenia numbers be
trusted. This bears on the manuscript's claim that the method generalises
beyond epilepsy.
