"""Does np.log depend on array alignment? That would explain everything.

numpy dispatches a SIMD kernel for transcendentals over aligned blocks and
handles the unaligned head and tail with scalar libm. The two can differ by an
ulp. Where a given element falls depends on the array's base address, which the
allocator chooses afresh each run -- so the SAME computation would differ on a
different, essentially random, subset of elements each time.
"""
import numpy as np

rng = np.random.default_rng(0)
n = 5_000_000
base = rng.uniform(0.05, 20.0, size=n + 8)

ref = np.log(base[:n].copy())

print("offset  bit-identical to offset 0   max rel diff")
for off in range(8):
    shifted = base[off:off + n].copy()          # same values, new allocation
    shifted[:] = base[:n]                        # force identical contents
    got = np.log(shifted)
    same = (got == ref).mean()
    rel = np.abs(got - ref) / np.abs(np.where(ref == 0, 1, ref))
    print(f"  {off}      {100*same:8.4f}%                 {rel.max():.3e}")

print()
# A view at a non-zero offset is genuinely unaligned rather than reallocated.
print("unaligned VIEWS of one buffer (no copy):")
buf = np.empty(n + 8, dtype=np.float64)
buf[:] = 0.0
for off in range(4):
    view = buf[off:off + n]
    view[:] = base[:n]
    got = np.log(view)
    same = (got == ref).mean()
    rel = np.abs(got - ref) / np.abs(np.where(ref == 0, 1, ref))
    print(f"  view offset {off}: bit-identical {100*same:8.4f}%  "
          f"max rel {rel.max():.3e}  aligned={view.ctypes.data % 64 == 0}")
