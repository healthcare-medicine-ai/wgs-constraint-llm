"""Is pandas' numexpr path the source of the run-to-run difference?

pandas hands elementwise ops on large Series to numexpr, which is multithreaded
and may use Intel VML. Thread count is detected at import and varies with what
the node is doing, so the same expression can land on a different code path
between runs.
"""
import os
import sys

import numpy as np
import pandas as pd

print("pandas", pd.__version__, "| numpy", np.__version__)
try:
    import numexpr
    print("numexpr", numexpr.__version__,
          "| threads", numexpr.nthreads,
          "| VML:", getattr(numexpr, "use_vml", False))
except ImportError:
    print("numexpr NOT installed")
    numexpr = None
print("pandas uses numexpr:", pd.get_option("compute.use_numexpr"))
print("NUMEXPR_MAX_THREADS =", os.environ.get("NUMEXPR_MAX_THREADS"))
print()

rng = np.random.default_rng(0)
n = 9_600_000
ac_case = pd.Series(rng.integers(0, 6, n))
ac_ctrl = pd.Series(rng.integers(0, 6, n))
an_case = pd.Series(rng.integers(20000, 40000, n))
an_ctrl = pd.Series(rng.integers(90000, 130000, n))


def effect(use_numexpr):
    pd.set_option("compute.use_numexpr", use_numexpr)
    r_case = an_case - ac_case
    r_ctrl = an_ctrl - ac_ctrl
    ratio = ((0.5 + ac_case) * (0.5 + r_ctrl)) / ((0.5 + r_case) * (0.5 + ac_ctrl))
    return np.log(ratio).to_numpy(dtype=float), ratio.to_numpy(dtype=float)


log_ne, ratio_ne = effect(True)
log_np, ratio_np = effect(False)

same_ratio = (ratio_ne == ratio_np).mean()
same_log = (log_ne == log_np).mean()
rel = np.abs(log_ne - log_np) / np.abs(np.where(log_np == 0, 1, log_np))
print(f"numexpr vs plain numpy, same process:")
print(f"  ratio bit-identical : {100*same_ratio:8.4f}%")
print(f"  log   bit-identical : {100*same_log:8.4f}%   max rel {rel.max():.3e}")
print()

if numexpr is not None:
    for nt in (1, 2, 4, 8):
        if nt > numexpr.MAX_THREADS:
            continue
        numexpr.set_num_threads(nt)
        pd.set_option("compute.use_numexpr", True)
        lg, _ = effect(True)
        s = (lg == log_np).mean()
        r = np.abs(lg - log_np) / np.abs(np.where(log_np == 0, 1, log_np))
        print(f"  numexpr {nt} thread(s) vs numpy: bit-identical {100*s:8.4f}%"
              f"   max rel {r.max():.3e}")
