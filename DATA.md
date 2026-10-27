# Data

None of the inputs are redistributed here. Put them in `data/` (or set `WGS_DATA_DIR`).
File names must match `config.yaml`. md5 is of the file as downloaded.

| File | Source | Bytes | md5 |
|---|---|---|---|
| `rgc_me_variant_frequencies_20231004.vcf.gz` | RGC Million Exome Variant Browser, https://rgc-research.regeneron.com/me/home (download after accepting the terms of use) | 267,585,341 | `c8278678167697ffdeb0d1fc58ff8142` |
| `aou_variants.tsv.gz` | All of Us Researcher Workbench, controlled tier (see note) | 2,358,692,639 | `41a7b4b0995b10e8dbe19c58edd495d4` |
| `gnomad.exomes.v4.0.coverage.summary.tsv.bgz` | gnomAD v4.0 exome coverage summary, https://gnomad.broadinstitute.org/downloads | 4,053,048,598 | `a6955332c9cccae7efb9c95581282a73` |
| `gencode.v44.basic.annotation.gtf.gz` | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz | 29,570,410 | `7450ef42cf9cb3d29625320b22d4bb45` |
| `epi25_variant_results.tsv.gz` | Epi25 Collaborative, https://epi25.broadinstitute.org/results | 252,492,252 | `96eb632897360f4e851e576c9a4d1381` |
| `AlphaMissense_hg38.tsv.gz` | https://doi.org/10.5281/zenodo.8208688 | 642,961,469 | `9fd167735f16a1b87da6eb3e4c25fcb5` |
| `All_hg38_RS.bw` | GERP RS lifted to hg38 (see note), https://doi.org/10.6084/m9.figshare.33201549.v1 | 16,358,613,092 | `10a03ee7969d3000ffd2e6e6f84f2453` |
| `gnomad.v4.0.constraint_metrics.tsv` | https://gnomad.broadinstitute.org/downloads#v4-constraint | 85,888,818 | `612f3899a85a1b7c6513a9dd65b00485` |
| `EpilepsyGenes_v2025-09.tsv` | Genes4Epilepsy, https://github.com/bahlolab/Genes4Epilepsy | 56,126 | `499af8f86f603813538fac3dc9ec03f6` |
| `SCHEMA_variant_results_hg38.tsv.gz` | SCHEMA variant results (hg19), https://schema.broadinstitute.org/, lifted to hg38 (see note) | 386,775,391 | `9d60c244f05f9d41825dc7b38fccbb9f` |
| `SCHEMA_gene_results.tsv.bgz` | SCHEMA gene results, https://schema.broadinstitute.org/ | 516,894 | `97877e41d294ce835bb6c69f346a8f76` |

The RGC-ME constraint predictions written by stage 01 (`HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz`,
538,175,139 bytes, md5 `d0b013f019b01f50f1f522b760da8e2c`) are also published at
https://doi.org/10.6084/m9.figshare.27184245.v1. With that file in `results/`, stages 02–04 and 06–09
run without the RGC-ME, All of Us and gnomAD coverage inputs. *[Oscar: the All of Us file in that
figshare record has no observed variants; see the handover notes.]*

## Notes

**GENCODE.** Use the *basic* annotation of release 44. Other GTFs from the same release give
different coding intervals.

**GERP on hg38.** UCSC publishes GERP RS for hg19 only
(https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=allHg19RS_BW). We lifted `All_hg19_RS.bw`
(7,580,335,760 bytes, md5 `f73544d231178bfb75ae30a431b008b0`) to hg38 with `bigWigLiftOver` and
the UCSC chain `hg19ToHg38.over.chain.gz`, with no options. The lifted file is on figshare.

**All of Us.** Variant presence across the All of Us exome data can only be produced inside the
Researcher Workbench and cannot be redistributed. *[Oscar: describe the export — dataset
version, region and the columns kept.]*

**SCHEMA on hg38.** SCHEMA publishes variant results on hg19 (`SCHEMA_variant_results.tsv.bgz`).
We lifted them to hg38 with pyliftover (UCSC hg19-to-hg38 chain); the code is in
`notebooks/SCZ Liftover.ipynb` at tag `ajhg-r4-submitted`. SCHEMA loci are 1-based and pyliftover
works in 0-based coordinates; the notebook passes the locus unchanged and keeps the 0-based result,
so the two offsets cancel. Checked against the hg38 reference: 29 of 30 sampled SNVs carry the
SCHEMA reference allele at the lifted position; the exception lies in the 1q21 segmental
duplication, where the hg38 reference base is the SCHEMA alternate allele.
*[To decide: move this step into the pipeline.]*

**Epi25 and SCHEMA.** Only the published summary results are used. Variant-level tables derived
from them are not committed.
