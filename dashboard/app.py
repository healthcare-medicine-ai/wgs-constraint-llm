"""Epilepsy exome sequencing dashboard: gene-level results of the unified meta-regression model.

Layout follows the Rivas Lab dashboard (github.com/marivascruz/epilepsy-dashboard); data are the corrected
results reported in the AJHG manuscript: one row per gene-group with at least 25 rare variants (66,748 rows),
with Epi25 damaging-missense and PTV gene-level p-values where Epi25 reports them.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

DATA = Path(__file__).parent / "epilepsy_unified_model_results.tsv.gz"
EXOME_WIDE = 3.4e-7
GROUPS = {"DEE": "Developmental and epileptic encephalopathy", "GGE": "Genetic generalized epilepsy",
          "NAFE": "Non-acquired focal epilepsy", "EPI": "All epilepsy"}
COLUMNS = {"gene_name": "Gene", "gene_id": "Ensembl ID", "group": "Group", "n_variants": "Variants",
           "p_unified": "Unified model p", "p_constraint": "HMM constraint p", "p_gerp": "GERP p",
           "p_pathogenicity": "Pathogenicity p", "p_epi25_dm": "Epi25 DM p", "p_epi25_ptv": "Epi25 PTV p",
           "g4e": "Genes4Epilepsy"}
P_COLUMNS = [c for c in COLUMNS if c.startswith("p_")]

st.set_page_config(page_title="Epilepsy exome sequencing dashboard", layout="wide")


@st.cache_data
def load():
    df = pd.read_csv(DATA, sep="\t")
    df["g4e"] = df["genes4epilepsy"].map({True: "yes", False: "no"}).fillna("not assessed")
    df["-log10 p"] = -np.log10(df["p_unified"])
    df["-log10 p Epi25"] = -np.log10(df[["p_epi25_dm", "p_epi25_ptv"]].min(axis=1))
    return df.sort_values("p_unified")


def table(rows):
    out = rows[list(COLUMNS)].rename(columns=COLUMNS)
    for c in P_COLUMNS:
        out[COLUMNS[c]] = out[COLUMNS[c]].map(lambda x: "–" if pd.isna(x) else f"{x:.2e}")
    return out


st.markdown(
    """
    <div style="display:flex;align-items:center;gap:16px;">
      <img src="https://upload.wikimedia.org/wikipedia/commons/2/26/Spike-waves.png" style="height:120px;">
      <span style="font-size:3rem;font-weight:700;line-height:1.1;">Epilepsy exome sequencing dashboard</span>
    </div>
    <a href="https://rivaslab.stanford.edu" target="_blank"
       style="display:flex;align-items:center;gap:8px;font-size:0.9em;margin-top:8px;">
      <img src="https://mrivas.su.domains/gbe/wp-content/uploads/2025/01/gbe.png" style="height:36px;"/>
      Built by Rivas Lab: Oscar Aguilar, Mijail Rivas, Manuel A. Rivas</a>
    """,
    unsafe_allow_html=True,
)

df = load()
left, mid, right = st.columns([1, 1, 2])
group = left.selectbox("Select group to visualize", list(GROUPS), format_func=lambda g: f"{g}: {GROUPS[g]}")
cutoff = mid.selectbox("Show genes with unified model p below", [1e-3, 1e-4, EXOME_WIDE, 1.0],
                       format_func=lambda p: "all genes" if p == 1.0 else f"{p:.1e}")
query = right.text_input("Find genes (comma-separated, searches all groups)")

genes = [g.strip().upper() for g in query.split(",") if g.strip()]
if genes:
    rows = df[df["gene_name"].str.upper().isin(genes)]
    st.subheader(f"Results for {', '.join(genes)}")
else:
    rows = df[(df["group"] == group) & (df["p_unified"] < cutoff)]
    st.subheader(f"Epilepsy results for analysis group: {group} ({len(rows):,} genes)")
st.dataframe(table(rows), hide_index=True, width="stretch")
download = rows[list(COLUMNS)[:-1] + ["genes4epilepsy"]].to_csv(sep="\t", index=False)
st.download_button("Download these rows (TSV)", download,
                   f"unified_model_{'search' if genes else group}.tsv", "text/tab-separated-values")

grp = df[df["group"] == group]
top = grp[grp["p_unified"] <= 1e-5]
fig = px.scatter(top, x="gene_name", y="-log10 p", color="g4e", hover_data=["gene_id", "n_variants"],
                 color_discrete_map={"yes": "#00c4b4", "no": "#f5a623", "not assessed": "#9e9e9e"},
                 labels={"gene_name": "Gene", "-log10 p": "-log10(p unified)", "g4e": "Genes4Epilepsy"},
                 title=f"Top genes for group {group} (unified model p ≤ 1e-5)", template="plotly_dark")
fig.add_hline(y=-np.log10(EXOME_WIDE), line=dict(color="red", dash="dash"), annotation_text="exome-wide 3.4e-7")
fig.update_traces(marker=dict(size=11))
order = top["gene_name"].tolist()  # already sorted by p
fig.update_xaxes(categoryorder="array", categoryarray=order, tickmode="array", tickvals=order, tickangle=-60)
st.plotly_chart(fig, width="stretch")

cmp = grp.dropna(subset=["-log10 p Epi25"])
fig2 = px.scatter(cmp, x="-log10 p Epi25", y="-log10 p", hover_name="gene_name", opacity=0.6,
                  labels={"-log10 p Epi25": "Epi25 -log10(p), min(DM, PTV)", "-log10 p": "Unified model -log10(p)"},
                  title=f"Unified model vs Epi25 for group {group}", template="plotly_dark")
edge = max(cmp["-log10 p Epi25"].max(), cmp["-log10 p"].max()) + 1
fig2.add_shape(type="line", x0=0, y0=0, x1=edge, y1=edge, line=dict(color="red", dash="dash"))
st.plotly_chart(fig2, width="stretch")

st.markdown("---")
st.markdown(
    "🧠 **Note:** only gene–groups with at least 25 rare variants were included in the unified meta-regression "
    "model. Exome-wide significance is p < 3.4 × 10⁻⁷, the Epi25 threshold. "
    "Code: [github.com/healthcare-medicine-ai/wgs-constraint-llm](https://github.com/healthcare-medicine-ai/wgs-constraint-llm)."
)
st.caption("Source of single-variant summary statistics and gene-level burden p-values: "
           "[Epi25](https://epi25.broadinstitute.org). – means Epi25 reports no result for that class.")
