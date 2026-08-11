"""Is the residual effect_size difference an environment artifact?

Every frame-construction variant gave the same 72.79%, so it is not pandas.
The remaining candidate is np.log itself: numpy dispatches a SIMD kernel chosen
by the CPU's vector extensions, and Sherlock's nodes are heterogeneous. If the
cached file was written on a node with different extensions, no code change can
reproduce it here.

Test: recompute effect_size for the differing rows three ways -- numpy
vectorised, numpy on a length-1 array (scalar path), and libm via math.log --
and see which, if any, matches the cache.
"""
import math
import sys
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
from wgs_constraint import get_config

cfg = get_config()
KEY = ["chr", "pos", "ref", "alt", "gene_id", "group"]
cached = pd.read_csv(cfg.result("constraint_gerp_am_scz_variants.tsv.gz"), sep="\t",
                     usecols=KEY + ["ac_case", "an_case", "ac_ctrl", "an_ctrl",
                                    "effect_size"]).drop_duplicates(subset=KEY)
print(f"cached rows: {len(cached):,}", flush=True)

a_case = cached["ac_case"].to_numpy(np.int64)
a_ctrl = cached["ac_ctrl"].to_numpy(np.int64)
r_case = cached["an_case"].to_numpy(np.int64) - a_case
r_ctrl = cached["an_ctrl"].to_numpy(np.int64) - a_ctrl
ratio = ((0.5 + a_case) * (0.5 + r_ctrl)) / ((0.5 + r_case) * (0.5 + a_ctrl))
ref = cached["effect_size"].to_numpy(float)

vec = np.log(ratio)
scalar = np.array([np.log(np.array([r]))[0] for r in ratio[:200000]])
libm = np.array([math.log(r) for r in ratio[:200000]])

print()
print(f"  numpy vectorised, full array : {100*(vec == ref).mean():6.2f}% match cache")
print(f"  numpy one-element-at-a-time  : {100*(scalar == ref[:200000]).mean():6.2f}% (first 200k)")
print(f"  math.log (libm scalar)       : {100*(libm == ref[:200000]).mean():6.2f}% (first 200k)")
print(f"  vectorised, same 200k slice  : {100*(vec[:200000] == ref[:200000]).mean():6.2f}%")
print()
print("  numpy", np.__version__, "| pandas", pd.__version__)
import subprocess
flags = subprocess.run(["bash", "-c", "grep -m1 ^flags /proc/cpuinfo"],
                       capture_output=True, text=True).stdout
have = [f for f in ("avx2", "avx512f", "avx512dq", "fma") if f" {f} " in flags]
print("  cpu vector extensions:", have)
print("  node:", subprocess.run(["hostname"], capture_output=True, text=True).stdout.strip())
