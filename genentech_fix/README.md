# Epi25 unified meta-regression — pipeline, provenance, and corrections

Companion material for **AJHG-D-24-00613**, *A unified meta-regression model
identifies genes associated with epilepsy* (Aguilar, Rivas & Rivas).

Prepared in response to the reproduction attempt by Jihoon Choi (Roche) and
Tushar Bhangale (Genentech), which surfaced two defects in the published
pipeline. Both are described in [Corrections](#corrections) below.

---

## 1. What is here

| File | Purpose |
|---|---|
| `01_build_regression_input.py` | GENCODE + Epi25 + HMM constraint + GERP + AlphaMissense → variant-level regression input |
| `02_run_metaregression.py` | Per-gene WLS meta-regression → gene-level p-values and comparison table |
| `03_compare_runs.py` | Compares any run against the submitted R2 baseline |
| `04_attribution.py` | Attributes each gene's change to a specific defect |
| `diag_duplication.py` | Quantifies row duplication in the inputs |
| `job_full.sbatch`, `job_attrib.sbatch`, `job_diag.sbatch` | Slurm wrappers |
| `baseline_R2/` | Snapshot of the exact notebooks and outputs behind the R2 submission, with MD5s |
| `environment-versions.txt` | Pinned versions of everything that touched the results |

These scripts are a direct port of `Epilepsy Analysis.ipynb` (cells 4–24). The
port was validated by running it with both corrections **disabled** and
confirming it reproduces the submitted output bit-for-bit across all 72,386
gene×group models (max relative difference `0.000e+00` on every p-value
column). That control is the reason any difference in the corrected run can be
attributed to the corrections rather than to the rewrite.

---

## 2. Inputs and their provenance

| Input | Source |
|---|---|
| `gencode.v44.basic.annotation.gtf.gz` | GENCODE release 44, `Gencode_human/release_44/gencode.v44.chr_patch_hapl_scaff.basic.annotation.gtf.gz` |
| `epi25_variant_results.tsv.gz` | Epi25 Collaborative, https://epi25.broadinstitute.org/results |
| `AlphaMissense_hg38.tsv.gz` | Zenodo record 10813168 (**isoform-level release** — see Correction 2) |
| `EpilepsyGenes_v2025-09.tsv` | Genes4Epilepsy, bahlolab/Genes4Epilepsy @ 3a8b2eb |
| HMM constraint predictions | https://doi.org/10.6084/m9.figshare.27184245.v1 |
| `All_hg38_RS.bw` | Derived — see below |

### The GERP hg38 track

GERP RS is published on hg19. The hg38 track used here was produced locally by
lifting it over. The exact commands, recovered from shell history:

```bash
cd $OAK/projects/wgs-constraint-llm/data
wget -O hg19ToHg38.over.chain.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget -O hg38.chrom.sizes \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
bigWigLiftOver All_hg19_RS.bw hg19ToHg38.over.chain.gz All_hg38_RS.bw unmapped.hg19.bed
```

**Caveat, stated plainly:** `bigWigLiftOver` is no longer present on our
systems and we cannot determine its version or origin, and `unmapped.hg19.bed`
was not retained. We are therefore shipping `All_hg38_RS.bw` itself rather than
a recipe, so that reproduction does not depend on recovering that tool. A
re-derived track from a different tool would not be the file that produced the
published results, and we would rather you have the actual input.

Header comparison of source and derived tracks:

| | `All_hg19_RS.bw` | `All_hg38_RS.bw` |
|---|---|---|
| bases covered | 3,095,677,412 | 2,864,325,273 |
| min / max | −12 / 6 | −12 / 6 |
| sum of squares | 11,135,093,840 | 11,127,259,590 |
| chromosome naming | UCSC (`chr1`…), with hg38 alt contigs in the derived file |

Coverage falls 7.5% across the liftover while the sum of squares is preserved
to within 0.07%, i.e. the dropped bases carry GERP ≈ 0 — unplaced and
repetitive sequence, not conserved coding sequence. Within the HMM position set
the loss was 0.0392%.

---

## 3. Pipeline

```
GENCODE v44 ──┐
              ├─→ Epi25 variants, filtered ──┐
Epi25 ────────┘                              │
                                             ├─→ regression input ─→ WLS ─→ gene p-values
HMM predictions ─→ + GERP (bigWig) ──────────┤
                                             │
AlphaMissense ─→ collapsed per variant ──────┘
```

### Filters and choices, in the order applied

| Step | Rule |
|---|---|
| Chromosomes | autosomes only; `chrX`, `chrY`, `chrMT` excluded |
| Consequence | `synonymous`, `non_coding`, `NA` excluded |
| Allele count | `ac_case + ac_ctrl <= 5` |
| Coverage | `an_case > 0` and `an_ctrl > 0` |
| Effect size | Haldane–Anscombe corrected log odds ratio, 0.5 continuity correction |
| Variance | Woolf variance of the same |
| Probability clipping | `prob_0`, `am_pathogenicity` clipped to `[0.01, 0.99]` |
| Missing values | `prob_0`, `GERP_RS`, `am_pathogenicity` imputed with the column mean; `pLoF_ind`, `missense_ind` filled with 0 |
| Transform | `-log1p(-x)` applied to `prob_0` and `am_pathogenicity`; **GERP enters the model untransformed** |
| Dropped rows | missing or zero `var_effect_size` |
| Model | WLS, weights `1 / var_effect_size`, grouped by (`gene_id`, `gene_name`, `group`) |
| Gene inclusion | `n_variants >= 25` |
| Reported p-value | `model.f_pvalue`, the **F-test for the whole regression** — not a coefficient p-value |
| Exome-wide significance | `p < 3.4e-7` |
| Suggestive | `p < 1.0e-4` |

### Two points that have caused confusion

**There is no MAF filter.** The public script is named
`unified_reg_MAF.05.py`, which implies MAF ≤ 0.05, but in the analysis notebook
that line is commented out:

```python
# frequency threshold
# input_df = input_df[input_df['AF'] <= 0.05]
```

The rare-variant restriction is the `ac_case + ac_ctrl <= 5` cut applied
upstream. Across Epi25's ~108,800 alleles that is a MAF of roughly 5e-5 —
about three orders of magnitude stricter than MAF ≤ 0.05. Applying a MAF filter
instead admits far more common variants, which carry smaller variances, hence
much larger WLS weights, and also inflates per-gene variant counts so more
genes clear the 25-variant floor.

**Firth statistics are not used.** The published pipeline computes effect sizes
directly from allele counts using the Haldane–Anscombe corrected log odds ratio
and its Woolf variance. Reconstructing effect sizes from allele counts with a
0.5 continuity correction reproduces our estimator exactly; no Firth BETA/SE
are required.

---

## 4. Corrections

Both defects were present in the submitted Revision 2 and were found through
the Genentech/Roche reproduction.

### Correction 1 — GERP coordinate off-by-one

`pyBigWig`'s `bw.values(chrom, 0, chr_len)` returns a vector whose element *i*
is the score at **0-based** coordinate *i*. The HMM position column is
**1-based** — `get_sequence()` in `HMM Mutation Predictions.ipynb` indexes its
masks directly with 1-based GTF and VCF coordinates. Indexing the 0-based
vector with a 1-based position returns the score for the next base along.

```diff
- mask = (sub["pos"] >= 0) & (sub["pos"] < chr_len)
+ mask = (sub["pos"] >= 1) & (sub["pos"] <= chr_len)

- gerp_vals = gerp_vec[pos_idx]
+ gerp_vals = gerp_vec[pos_idx - 1]
```

The bounds guard changes with it: with a `-1` index the old `pos >= 0` would
wrap position 0 to the last base of the chromosome through negative indexing,
and `pos < chr_len` silently discarded the final base.

The assumption is recorded in the notebook, cell 11:
`# (Assuming both datasets are 0-based already per your note)`.

**Magnitude.** 81.68% of GERP values change; old and corrected values correlate
at **r = 0.173**. A one-base shift does not merely blur this track: adjacent
bases within a codon differ sharply in constraint, so the shift is closer to
randomising GERP than to attenuating it. The 18% unchanged sit inside
constant-value runs in the bigWig.

The same defect exists in `Constraint Measures Comparison.ipynb` (which
produces the HMM-vs-GERP figure panels) and in `Schizophrenia Analysis.ipynb`.

### Correction 2 — AlphaMissense transcript fan-out

`AlphaMissense_hg38.tsv.gz` is the isoform-level release: 71,697,556 rows for
71,236,463 distinct variants. Merged without de-duplication, a variant
annotated on *N* transcripts became *N* rows — counted *N* times in
`n_variants` **and** weighted *N*× in the WLS fit.

```diff
+ am = am.groupby(["chr","pos","ref","alt"], sort=False,
+                 as_index=False)["am_pathogenicity"].max()
```

Collapsed by **maximum** across transcripts. Duplicates are almost entirely
alternative isoforms of one protein — only 1,854 of 454,089 duplicated variants
span more than one UniProt id — so a per-variant collapse is well defined, and
the maximum is insensitive to how many isoforms happen to be annotated at a
locus.

**Magnitude.** 74,378 duplicate rows, 0.77% of the regression input. 738
gene×group models carried an inflated `n_variants`; 48 gene-groups passed the
`n_variants >= 25` filter only on duplicate rows.

### Effect on results

The exome-wide significant set moves from 29 genes to 18: 16 retained, 13 lost,
2 gained. Per-gene attribution is in `attribution_by_gene.tsv`.

---

## 5. Environment

Versions that produced both the submitted and the corrected results are in
`environment-versions.txt`:

```
python 3.11.4
numpy==1.24.3      pandas==1.5.3       scipy==1.11.4
statsmodels==0.14.2 pyBigWig==0.3.24   pyliftover==0.4
tqdm==4.65.0       matplotlib==3.7.1
```

## 6. Running it

```bash
# reproduce the submitted numbers (fidelity control)
sbatch --job-name=repro job_full.sbatch control

# corrected
sbatch --job-name=fixed job_full.sbatch fixed

# isolate one defect at a time
sbatch --job-name=gerponly job_attrib.sbatch gerponly
sbatch --job-name=amonly   job_attrib.sbatch amonly

# compare and attribute
python 03_compare_runs.py
python 04_attribution.py
```

Outputs are suffixed `_REPRO`, `_FIXED`, `_GERPONLY`, `_AMONLY`. Nothing
overwrites the submitted files.

## 7. Known gaps

- The Figure 5 heatmap (notebook cell 21) and the Table S1 export (cell 22)
  are not yet ported to scripts.
- `bigWigLiftOver` and `unmapped.hg19.bed` are not recoverable, as above.
- The analysis behind Revision 2 was never committed to version control; the
  `baseline_R2/` snapshot exists to make that state permanently recoverable.
