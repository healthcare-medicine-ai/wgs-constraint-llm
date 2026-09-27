"""P29: OLS R2 of GERP_RS on prob_0 (notebook cell 13), as published vs corrected;
plus where the corrected percent-difference table is most over/under-represented."""
import numpy as np, pandas as pd
R = "results/"
for tag, f in [("ASPUBLISHED", "HMM_rgc_ALL_RS_merged_predictions_ASPUBLISHED.tsv.gz"),
               ("FIXED", "HMM_rgc_ALL_RS_merged_predictions.tsv.gz")]:
    d = pd.read_csv(R + f, sep="\t", usecols=["prob_0", "GERP_RS"]).dropna()
    r = np.corrcoef(d.prob_0, d.GERP_RS)[0, 1]
    print(f"{tag}: n={len(d):,}  r={r:.4f}  R2={r*r:.6f}", flush=True)
for tag in ("ASPUBLISHED", "FIXED"):
    z = np.load(R + f"hmm_gerp_joint_distribution_{tag}.npz")
    pdf = pd.DataFrame(z["percent_diff"], index=z["prob_bins"][:-1], columns=z["gerp_bins"][:-1])
    print(f"\n{tag} percent_diff (rows prob_0 bin left edge, cols GERP bin left edge)")
    print(pdf.round(2).to_string())
