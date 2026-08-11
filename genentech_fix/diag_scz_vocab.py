import sys; sys.path.insert(0, "src")
import pandas as pd
from wgs_constraint import get_config
cfg = get_config()
c = pd.read_csv(cfg.data_dir / "SCHEMA_variant_results_hg38.tsv.gz", sep="\t",
                compression="gzip", usecols=["consequence"])
vc = c["consequence"].value_counts(dropna=False)
print("full SCHEMA consequence vocabulary:")
print(vc.to_string())
print()
print("the three still-commented entries:")
for c3 in ("intergenic_variant", "non_coding_transcript_exon_variant",
           "non_coding_transcript_variant"):
    print(f"  {c3:<38} {int(vc.get(c3, 0)):,}")
