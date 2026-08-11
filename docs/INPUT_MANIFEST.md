# Input manifest

Generated 2026-08-11 03:57 UTC by `pipelines/00_input_manifest.py`.

Verify with `python pipelines/00_input_manifest.py --verify`.

`kind` is **raw** for a downloaded input, **pinned** for something that cannot practically be regenerated and must be treated as fixed, and **derived** for something this repository produces.

For gzipped files the decompressed md5 is the meaningful one; gzip headers differ between writes of identical content.

| File | Kind | Bytes | md5 (file) | md5 (decompressed) |
|---|---|---|---|---|
| `gencode.v44.basic.annotation.gtf.gz` | raw | 29,570,410 | `7450ef42cf9cb3d29625320b22d4bb45` | `a724f79aa2e8d586ff9f84bef4628b0c` |
| `epi25_variant_results.tsv.gz` | raw | 252,492,252 | `96eb632897360f4e851e576c9a4d1381` | `a40a3fcd4ffb49dc35d2e84a13c403d0` |
| `AlphaMissense_hg38.tsv.gz` | raw | 642,961,469 | `9fd167735f16a1b87da6eb3e4c25fcb5` | `178b56ce32b5529f0d28a46c03301b20` |
| `EpilepsyGenes_v2025-09.tsv` | raw | 56,126 | `499af8f86f603813538fac3dc9ec03f6` | `—` |
| `All_hg38_RS.bw` | pinned | 16,358,613,092 | `10a03ee7969d3000ffd2e6e6f84f2453` | `—` |
| `All_hg19_RS.bw` | raw | 7,580,335,760 | `f73544d231178bfb75ae30a431b008b0` | `—` |
| `hg19ToHg38.over.chain.gz` | raw | 227,698 | `35887f73fe5e2231656504d1f6430900` | `9267c9fef79b54962da8efadd0ddf6b6` |
| `gnomad.v4.0.constraint_metrics.tsv` | raw | 85,888,818 | `612f3899a85a1b7c6513a9dd65b00485` | `—` |
| `SCHEMA_variant_results_hg38.tsv.gz` | raw | 386,775,391 | `9d60c244f05f9d41825dc7b38fccbb9f` | `ac14d58300c88f787a680c6de08a9a02` |
| `SCHEMA_gene_results.tsv.bgz` | raw | 516,894 | `97877e41d294ce835bb6c69f346a8f76` | `705a457328f035d8b20d06d5d1e56e79` |
| `HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz` | pinned | 538,175,139 | `d0b013f019b01f50f1f522b760da8e2c` | `c89590067da0ffcb56e30cfa7f9624bb` |
| `HMM_rgc_ALL_RS_merged_predictions.tsv.gz` | derived | 577,699,523 | `b6cc9daefe3344d2db84350da2683114` | `9adf4de7f1850dbd7d5462ddf5793da3` |

## Notes

- **`gencode.v44.basic.annotation.gtf.gz`** — GENCODE v44 primary assembly. NOTE: notebook comments cite the chr_patch_hapl_scaff URL, which is a different, larger file.
- **`epi25_variant_results.tsv.gz`** — Epi25 summary statistics (controlled access)
- **`AlphaMissense_hg38.tsv.gz`** — AlphaMissense hg38, ISOFORM-level release -- must be collapsed per variant
- **`EpilepsyGenes_v2025-09.tsv`** — Genes4Epilepsy reference gene set
- **`All_hg38_RS.bw`** — GERP RS hg38, lifted from hg19 with bigWigLiftOver (tool no longer available). Positional correctness verified 2026-08-09.
- **`All_hg19_RS.bw`** — GERP RS hg19, source of the liftover
- **`hg19ToHg38.over.chain.gz`** — UCSC chain used for the liftover
- **`gnomad.v4.0.constraint_metrics.tsv`** — gnomAD v4 constraint, source of missense z-score and MTR
- **`SCHEMA_variant_results_hg38.tsv.gz`** — SCHEMA variants, lifted to hg38
- **`SCHEMA_gene_results.tsv.bgz`** — SCHEMA gene-level results
- **`HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz`** — HMM constraint predictions, RGC-ME. Hours to regenerate; distributed at doi:10.6084/m9.figshare.27184245. Do not re-derive.
- **`HMM_rgc_ALL_RS_merged_predictions.tsv.gz`** — Predictions annotated with corrected GERP (pipelines/05)

## Locations

- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/gencode.v44.basic.annotation.gtf.gz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/epi25_variant_results.tsv.gz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/AlphaMissense_hg38.tsv.gz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/EpilepsyGenes_v2025-09.tsv`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/All_hg38_RS.bw`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/All_hg19_RS.bw`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/hg19ToHg38.over.chain.gz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/gnomad.v4.0.constraint_metrics.tsv`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/SCHEMA_variant_results_hg38.tsv.gz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/SCHEMA_gene_results.tsv.bgz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/results/HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz`
- `/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/osthoag/wgs-constraint-llm/results/HMM_rgc_ALL_RS_merged_predictions.tsv.gz`
