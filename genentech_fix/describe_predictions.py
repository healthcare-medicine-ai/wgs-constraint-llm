import pandas as pd
F = "results/HMM_rgc_ALL_RS_merged_predictions.tsv.gz"
df = pd.read_csv(F, sep="\t")
print(f"rows        : {len(df):,}")
print(f"columns     : {list(df.columns)}")
print(f"chromosomes : {df['chr'].nunique()}")
print(f"GERP present: {df['GERP_RS'].notna().sum():,} "
      f"({100*df['GERP_RS'].notna().mean():.4f}%)")
print(f"GERP range  : {df['GERP_RS'].min():.2f} to {df['GERP_RS'].max():.2f}")
print(f"prob_0 range: {df['prob_0'].min():.3g} to {df['prob_0'].max():.3g}")
print(f"observation : {sorted(df['observation'].dropna().unique())[:5]}")
