"""What does p_unified test? For named genes, compare the unified F-test (do variant
effect sizes track the moderators?) with a test of the weighted mean effect size
(are variants, on average, enriched in cases? -- a burden-like question)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from statsmodels.regression.linear_model import WLS
from statsmodels.tools.tools import add_constant
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config, prepare_regression_input, MODERATORS
cfg = get_config()
want = [("STX1B","DEE"),("DEPDC5","EPI"),("GABRG2","EPI"),("DEPDC5","NAFE"),
        ("SP4","GGE"),("KCNQ2","DEE"),("STXBP1","DEE"),("SCN1A","DEE"),("SCN2A","DEE")]
df = pd.read_csv(cfg.result("constraint_gerp_am_epi25_variants_FIXED.tsv.gz"), sep="\t")
df = df[df.set_index(["gene_name","group"]).index.isin(want)]
df = prepare_regression_input(df, clip_epsilon=cfg.param("clip_epsilon"))
print(f"{'gene':7s} {'group':5s} {'n':>4s}  {'mean effect':>11s} {'p(mean)':>9s}  {'p_unified (F)':>13s}")
for (g, grp), d in df.groupby(["gene_name","group"]):
    w = 1 / d["var_effect_size"]
    m0 = WLS(d["effect_size"], np.ones(len(d)), weights=w).fit()          # intercept only
    m1 = WLS(d["effect_size"], add_constant(d[MODERATORS]), weights=w).fit()
    print(f"{g:7s} {grp:5s} {len(d):4d}  {m0.params.iloc[0]:11.3f} {m0.pvalues.iloc[0]:9.1e}  {m1.f_pvalue:13.1e}")
