"""Epilepsy unified meta-regression results: gene-level browser for the AJHG manuscript.

Data: epilepsy_unified_model_results.tsv.gz, one row per gene-group with at least 25 variants
(66,748 rows), unified model and moderator p-values from the corrected pipeline, and Epi25
damaging-missense and PTV gene-level p-values where Epi25 reports them.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

DATA = Path(__file__).parent / "epilepsy_unified_model_results.tsv.gz"
EXOME_WIDE = 3.4e-7
SUGGESTIVE = 1e-4
P_COLS = ["p_unified", "p_constraint", "p_gerp", "p_pathogenicity", "p_epi25_dm", "p_epi25_ptv"]
LABELS = {
    "gene_name": "Gene", "gene_id": "Ensembl ID", "group": "Group", "n_variants": "Variants",
    "p_unified": "Unified model p", "p_constraint": "HMM constraint p", "p_gerp": "GERP p",
    "p_pathogenicity": "Pathogenicity p", "p_epi25_dm": "Epi25 DM p", "p_epi25_ptv": "Epi25 PTV p",
    "genes4epilepsy": "Genes4Epilepsy",
}
GROUPS = {"DEE": "Developmental and epileptic encephalopathy", "GGE": "Genetic generalized epilepsy",
          "NAFE": "Non-acquired focal epilepsy", "EPI": "All epilepsy"}

st.set_page_config(page_title="Epilepsy unified meta-regression results", layout="wide")


@st.cache_data
def load():
    df = pd.read_csv(DATA, sep="\t")
    df["g4e"] = df["genes4epilepsy"].map({True: "yes", False: "no"}).fillna("not assessed")
    df["mlog10_unified"] = -np.log10(df["p_unified"])
    df["mlog10_epi25"] = -np.log10(df[["p_epi25_dm", "p_epi25_ptv"]].min(axis=1))
    return df


df = load()

st.title("Epilepsy unified meta-regression results")
st.markdown(
    "Gene-level results from *A unified meta-regression model identifies genes associated with "
    "epilepsy* (Aguilar, Rivas, and Rivas). The unified model tests whether rare-variant effect "
    "sizes in Epi25 track HMM constraint, GERP conservation, AlphaMissense pathogenicity, and "
    "pLoF/missense annotations. Code and data: "
    "[github.com/healthcare-medicine-ai/wgs-constraint-llm]"
    "(https://github.com/healthcare-medicine-ai/wgs-constraint-llm)."
)

c1, c2, c3 = st.columns([1, 1, 2])
group = c1.selectbox("Epilepsy group", list(GROUPS), format_func=lambda g: f"{g}: {GROUPS[g]}")
cutoff = c2.selectbox("Unified model p below", [EXOME_WIDE, SUGGESTIVE, 1e-3, 1.0],
                      index=1, format_func=lambda p: "all genes" if p == 1.0 else f"{p:.1e}")
query = c3.text_input("Find genes (comma-separated symbols, searches all groups)", "")

if query.strip():
    genes = [g.strip().upper() for g in query.split(",") if g.strip()]
    view = df[df["gene_name"].str.upper().isin(genes)]
    st.subheader(f"{len(view)} gene–group results for {', '.join(genes)}")
else:
    view = df[(df["group"] == group) & (df["p_unified"] < cutoff)]
    st.subheader(f"{GROUPS[group]} ({group}): {len(view):,} genes")

view = view.sort_values("p_unified")
show = view[list(LABELS)[:-1] + ["g4e"]].rename(columns={**LABELS, "g4e": "Genes4Epilepsy"})
fmt = [LABELS[c] for c in P_COLS]
for col in fmt:  # rows are already sorted by unified p, so text formatting loses nothing
    show[col] = show[col].map(lambda x: "–" if pd.isna(x) else f"{x:.2e}")
st.dataframe(show, width="stretch", hide_index=True)
st.download_button("Download this selection (TSV)", view.drop(columns=["g4e", "mlog10_unified", "mlog10_epi25"])
                   .to_csv(sep="\t", index=False), f"unified_model_{group}.tsv", "text/tab-separated-values")

if not query.strip():
    grp = df[df["group"] == group].dropna(subset=["mlog10_epi25"])
    fig = px.scatter(grp, x="mlog10_epi25", y="mlog10_unified", hover_name="gene_name",
                     hover_data={"n_variants": True, "mlog10_epi25": ":.2f", "mlog10_unified": ":.2f"},
                     opacity=0.5, labels={"mlog10_epi25": "Epi25 −log10(p), min(DM, PTV)",
                                          "mlog10_unified": "Unified model −log10(p)"},
                     title=f"Unified model vs Epi25, {group}")
    top = max(grp["mlog10_epi25"].max(), grp["mlog10_unified"].max()) + 1
    fig.add_shape(type="line", x0=0, y0=0, x1=top, y1=top, line=dict(color="red", dash="dash"))
    fig.add_hline(y=-np.log10(EXOME_WIDE), line=dict(color="grey", dash="dot"),
                  annotation_text="exome-wide (3.4e-7)")
    st.plotly_chart(fig, width="stretch")

st.caption(
    "Exome-wide significance p < 3.4 × 10⁻⁷ (the Epi25 threshold); suggestive p < 1 × 10⁻⁴. "
    "Only gene–groups with at least 25 rare variants are modeled. Epi25 p-values are gene-level "
    "damaging-missense (DM) and protein-truncating (PTV) burden results from the "
    "[Epi25 Collaborative](https://epi25.broadinstitute.org); – means Epi25 reports no result. "
    "Genes4Epilepsy is marked “not assessed” for genes without an Epi25 result."
)
