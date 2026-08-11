"""Quantify how the GERP correction moved the HMM-vs-GERP joint distribution."""
import numpy as np
from pathlib import Path
R = Path("results")
a = np.load(R / "hmm_gerp_joint_distribution_ASPUBLISHED.npz")
b = np.load(R / "hmm_gerp_joint_distribution_FIXED.npz")
for k in ("percent_diff", "chi_sqr"):
    x, y = a[k], b[k]
    print(f"  {k}")
    print(f"    as published : max|{np.abs(x).max():.4f}|   sum {x.sum():.4e}")
    print(f"    corrected    : max|{np.abs(y).max():.4f}|   sum {y.sum():.4e}")
    print(f"    correlation  : {np.corrcoef(x.ravel(), y.ravel())[0, 1]:.6f}")
print(f"  total chi-squared: {a['chi_sqr'].sum():.4e} -> {b['chi_sqr'].sum():.4e} "
      f"({100*(b['chi_sqr'].sum()/a['chi_sqr'].sum()-1):+.2f}%)")
