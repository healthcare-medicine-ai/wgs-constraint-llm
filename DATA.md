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
| `DEE_results_2023_12_30_14_36_36.csv` | Epi25 gene results for DEE, exported from https://epi25.broadinstitute.org/results | 1,534,357 | `bc5021e614379248afc2a357a3290a7c` |
| `EPI_results_2023_12_30_14_36_32.csv` | Epi25 gene results for EPI, as above | 1,572,565 | `f6c4acbba5868ff8c2dd9a90e711eb71` |
| `GGE_results_2023_12_30_14_36_38.csv` | Epi25 gene results for GGE, as above | 1,539,599 | `329709d5db3416d4591364bed5403d9f` |
| `NAFE_results_2023_12_30_14_36_41.csv` | Epi25 gene results for NAFE, as above (one value edited, see note) | 1,544,077 | `e6b7303350e854da2ad73917113a069a` |
| `AlphaMissense_hg38.tsv.gz` | https://doi.org/10.5281/zenodo.8208688 | 642,961,469 | `9fd167735f16a1b87da6eb3e4c25fcb5` |
| `All_hg38_RS.bw` | GERP RS lifted to hg38 (see note), https://doi.org/10.6084/m9.figshare.33201549.v1 | 16,358,613,092 | `10a03ee7969d3000ffd2e6e6f84f2453` |
| `gnomad.v4.0.constraint_metrics.tsv` | https://gnomad.broadinstitute.org/downloads#v4-constraint | 85,888,818 | `612f3899a85a1b7c6513a9dd65b00485` |
| `EpilepsyGenes_v2025-09.tsv` | Genes4Epilepsy, https://github.com/bahlolab/Genes4Epilepsy | 56,126 | `499af8f86f603813538fac3dc9ec03f6` |
| `SCHEMA_variant_results.tsv.bgz` | SCHEMA variant results (hg19), https://schema.broadinstitute.org/ | 414,642,856 | `e72dada66d9b70457b7dc94b4e16c7e4` |
| `SCHEMA_gene_results.tsv.bgz` | SCHEMA gene results, https://schema.broadinstitute.org/ | 516,894 | `97877e41d294ce835bb6c69f346a8f76` |
| `hg19ToHg38.over.chain.gz` | UCSC, https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | 227,698 | `35887f73fe5e2231656504d1f6430900` |

The RGC-ME constraint predictions written by stage 01 (`HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz`,
538,175,139 bytes, md5 `d0b013f019b01f50f1f522b760da8e2c`) are also published at
https://doi.org/10.6084/m9.figshare.27184245.v1. With that file in `results/`, stages 02–04 and 06–09
run without the RGC-ME, All of Us and gnomAD coverage inputs.

## Notes

**GENCODE.** Use the *basic* annotation of release 44. Other GTFs from the same release give
different coding intervals.

**GERP on hg38.** UCSC publishes GERP RS for hg19 only
(https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=allHg19RS_BW). We lifted `All_hg19_RS.bw`
(7,580,335,760 bytes, md5 `f73544d231178bfb75ae30a431b008b0`) to hg38 with `bigWigLiftOver` and
the UCSC chain above, with no options. The lifted file is on figshare.

**All of Us.** `aou_variants.tsv.gz` has one row per variant site (`chr`, `pos`) in the All of Us
short-read genome callset *[confirm: CDR v7, 245,388 participants; all sites or only PASS sites]*.
It was exported from the controlled-tier Researcher Workbench as eight files of variant IDs
(`chr-pos-ref-alt`), which were concatenated and reduced to unique positions. It cannot be
redistributed. Only stage 01 reads it, to make the All of Us predictions that stage 05 uses.

**SCHEMA on hg38.** SCHEMA publishes variant results on hg19. Stage 09 lifts each locus to hg38
with pyliftover and the UCSC chain above, and drops the 0.01% of loci that do not map.

**Epi25 gene results.** Stage 03 reads one CSV per epilepsy group for the Epi25 damaging-missense
and PTV p-values. As exported, the NAFE file gives DEPDC5's PTV p-value as `0`; we replaced it
with `2.2e-16` (md5 of the unedited file: `15eb6952a7251e4bd631f3b2d247f467`). That is the value
in Table 1. No other value was changed.

**Epi25 and SCHEMA.** Only the published summary results are used. Variant-level tables derived
from them are not committed.
