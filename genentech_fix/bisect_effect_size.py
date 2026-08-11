"""Fast bisect for the SCZ effect_size last-bit difference.

Loading variants takes a few minutes; the GERP and AlphaMissense stages take
most of an hour. This does only the variant load, so construction variants can
be compared against the cached merge quickly.

Each variant builds the frame a different way, then reports how many rows'
effect_size are bit-identical to the September 2025 cache.
"""
import sys
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from wgs_constraint import get_config
from wgs_constraint.epi25 import load_gene_annotation

sys.path.insert(0, "pipelines")
import importlib.util
spec = importlib.util.spec_from_file_location("s11", "pipelines/11_schizophrenia.py")
s11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s11)

cfg = get_config()
KEY = ["chr", "pos", "ref", "alt", "gene_id", "group"]
MAXAC = cfg.param("max_allele_count")
EXCL_CHR = cfg.param("excluded_chromosomes")

print("loading the cached merge for reference ...", flush=True)
cached = pd.read_csv(cfg.result("constraint_gerp_am_scz_variants.tsv.gz"), sep="\t",
                     usecols=KEY + ["effect_size"]).drop_duplicates(subset=KEY)
print(f"  {len(cached):,} distinct keys", flush=True)

print("loading GENCODE ...", flush=True)
genes = load_gene_annotation(cfg.data("gene_annotation"))
print("loading SCHEMA ...", flush=True)
raw = pd.read_csv(cfg.data_dir / "SCHEMA_variant_results_hg38.tsv.gz", sep="\t",
                  compression="gzip")
print(f"  {len(raw):,} rows\n", flush=True)


def filtered(df):
    """Everything up to but excluding the allele-count cut and effect sizes."""
    df = df[~df["consequence"].isin(s11.EXCLUDED_CONSEQUENCES)]
    df = df[~df["consequence"].isna()]
    return pd.merge(df, genes[["std_gene_id", "gene_name"]].drop_duplicates(),
                    left_on="gene_id", right_on="std_gene_id",
                    how="left").drop("std_gene_id", axis=1)


def chr_chained(df):
    m = pd.Series(True, index=df.index)
    for e in EXCL_CHR:
        m &= (df["chr"] != e)
    return df[m]


def chr_compound(df):
    return df[(df["chr"] != "chrX") & (df["chr"] != "chrY") & (df["chr"] != "chrMT")]


def chr_isin(df):
    return df[~df["chr"].isin(EXCL_CHR)]


def effect(out, *, copy=False, numpy_path=False):
    out = out[(out["ac_ctrl"] + out["ac_case"] <= MAXAC)
              & (out["an_case"] > 0) & (out["an_ctrl"] > 0)]
    if copy:
        out = out.copy()
    a_case, a_ctrl = out["ac_case"], out["ac_ctrl"]
    r_case = out["an_case"] - out["ac_case"]
    r_ctrl = out["an_ctrl"] - out["ac_ctrl"]
    num = (0.5 + a_case) * (0.5 + r_ctrl)
    den = (0.5 + r_case) * (0.5 + a_ctrl)
    if numpy_path:
        out["effect_size"] = np.log(np.asarray(num / den, dtype=np.float64))
    else:
        out["effect_size"] = np.log(num / den)
    return out


VARIANTS = {
    "chained !=  (current)":        lambda: effect(filtered(chr_chained(raw))),
    "compound mask (cell 15 form)": lambda: effect(filtered(chr_compound(raw))),
    ".isin()":                      lambda: effect(filtered(chr_isin(raw))),
    "chained != + .copy()":         lambda: effect(filtered(chr_chained(raw)), copy=True),
    "chained != + np.asarray":      lambda: effect(filtered(chr_chained(raw)), numpy_path=True),
    "chr filter LAST":              lambda: effect(chr_chained(filtered(raw))),
}

print("=" * 72)
for name, build in VARIANTS.items():
    try:
        out = build()
        j = out[KEY + ["effect_size"]].merge(cached, on=KEY, suffixes=("_n", "_c"))
        same = (j["effect_size_n"] == j["effect_size_c"]).mean()
        x = j["effect_size_c"].to_numpy(float)
        y = j["effect_size_n"].to_numpy(float)
        f = np.isfinite(x) & np.isfinite(y) & (x != 0)
        rel = np.abs(x[f] - y[f]) / np.abs(x[f])
        print(f"  {name:<30} rows {len(out):>9,}  matched {len(j):>9,}  "
              f"bit-identical {100*same:6.2f}%  maxrel {rel.max():.2e}")
    except Exception as exc:
        print(f"  {name:<30} ERROR {type(exc).__name__}: {exc}")
print("=" * 72)
