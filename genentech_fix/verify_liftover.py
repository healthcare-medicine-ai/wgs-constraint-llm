#!/usr/bin/env python
"""Verify that All_hg38_RS.bw is positionally correct relative to All_hg19_RS.bw.

We have already found one coordinate-convention defect in this project. If the
liftover that produced the hg38 track introduced a shift of its own, correcting
the lookup index would be fixing one error on top of another and the GERP values
would still be wrong.

Method
------
Sample positions on hg19 where GERP varies sharply between neighbouring bases --
flat regions cannot distinguish a shift. Map each through the UCSC chain with
pyliftover, then compare the hg19 value at the source position against the hg38
value at the mapped position.

Both bigWigs are read with the interval API, `bw.values(chrom, p-1, p)[0]`,
which is 0-based half-open and therefore returns the score at 1-based position
p by definition. That sidesteps the vector-indexing convention entirely, so this
check is independent of the bug under investigation.

To detect a residual shift, each hg38 position is also probed at -2..+2 offsets
and we report which offset matches best. A correct liftover peaks at 0.
"""

import argparse
import random
from collections import Counter

import numpy as np
import pyBigWig
from pyliftover import LiftOver

DATA = "/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/"
HG19 = DATA + "All_hg19_RS.bw"
HG38 = DATA + "All_hg38_RS.bw"
CHAIN = DATA + "hg19ToHg38.over.chain.gz"

# Checked in this order, breaking on the first match. Offset 0 must come first:
# GERP values recur, so a neighbouring offset can match coincidentally, and
# testing it before 0 would attribute a correct position to a shift.
OFFSETS = [0, -1, 1, -2, 2]
TOL = 1e-4


def value_at(bw, chrom, pos_1based):
    """GERP at a 1-based position, via the 0-based half-open interval API."""
    if pos_1based < 1:
        return None
    try:
        vals = bw.values(chrom, pos_1based - 1, pos_1based)
    except (RuntimeError, OverflowError):
        return None
    if not vals:
        return None
    v = vals[0]
    return None if v is None or np.isnan(v) else float(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="chr1,chr7,chr12,chr20",
                    help="comma-separated chromosomes to sample")
    ap.add_argument("--per-chrom", type=int, default=150,
                    help="informative positions to test per chromosome")
    ap.add_argument("--seed", type=int, default=20260809)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    print("opening inputs ...", flush=True)
    lo = LiftOver(CHAIN)
    bw19 = pyBigWig.open(HG19)
    bw38 = pyBigWig.open(HG38)
    chroms19 = bw19.chroms()
    chroms38 = bw38.chroms()

    offset_hits = Counter()
    n_tested = n_unmapped = n_no_hg38 = n_multi = 0
    mismatches = []

    for chrom in args.chroms.split(","):
        if chrom not in chroms19:
            print(f"  {chrom}: absent from hg19 track, skipping")
            continue
        length = int(chroms19[chrom])
        found = 0
        attempts = 0

        # Sample until we have enough positions whose neighbourhood is not flat.
        while found < args.per_chrom and attempts < args.per_chrom * 400:
            attempts += 1
            p = rng.randint(2_000_000, length - 2_000_000)

            here = value_at(bw19, chrom, p)
            nxt = value_at(bw19, chrom, p + 1)
            prv = value_at(bw19, chrom, p - 1)
            if here is None or nxt is None or prv is None:
                continue
            # Require a locally distinctive value, otherwise a one-base shift
            # is invisible and the test proves nothing.
            if abs(here - nxt) < 0.5 or abs(here - prv) < 0.5:
                continue

            mapped = lo.convert_coordinate(chrom, p - 1)  # pyliftover is 0-based
            if not mapped:
                n_unmapped += 1
                continue
            if len(mapped) > 1:
                n_multi += 1
            m_chrom, m_pos0 = mapped[0][0], mapped[0][1]
            if m_chrom not in chroms38:
                n_no_hg38 += 1
                continue
            m_pos1 = m_pos0 + 1  # back to 1-based

            best = None
            for off in OFFSETS:
                v38 = value_at(bw38, m_chrom, m_pos1 + off)
                if v38 is not None and abs(v38 - here) <= TOL:
                    best = off
                    break

            found += 1
            n_tested += 1
            if best is None:
                offset_hits["no match"] += 1
                if len(mismatches) < 8:
                    v0 = value_at(bw38, m_chrom, m_pos1)
                    mismatches.append(
                        (chrom, p, m_chrom, m_pos1, here, v0))
            else:
                offset_hits[best] += 1

        print(f"  {chrom}: {found} informative positions tested "
              f"({attempts} sampled)", flush=True)

    bw19.close()
    bw38.close()

    print("\n" + "=" * 66)
    print("LIFTOVER POSITIONAL VERIFICATION")
    print("=" * 66)
    print(f"positions tested        : {n_tested}")
    print(f"unmapped by chain       : {n_unmapped}")
    print(f"mapped to absent contig : {n_no_hg38}")
    print(f"multi-mapping positions : {n_multi}")
    print("\nhg38 offset at which the hg19 value was found"
          " (offset 0 tested first):")
    for off in sorted(OFFSETS):
        n = offset_hits.get(off, 0)
        marker = "   <-- expected" if off == 0 else ""
        pct = 100 * n / n_tested if n_tested else 0
        print(f"  {off:+d} : {n:5d}  ({pct:5.1f}%){marker}")
    nm = offset_hits.get("no match", 0)
    print(f"  no match within +/-2 : {nm:5d}  "
          f"({100 * nm / n_tested if n_tested else 0:5.1f}%)")

    if mismatches:
        print("\nexamples with no match (hg19 pos -> hg38 pos, hg19 val, hg38 val):")
        for c19, p19, c38, p38, v19, v38 in mismatches:
            shown = f"{v38:.4f}" if v38 is not None else "None"
            print(f"  {c19}:{p19} -> {c38}:{p38}   {v19:.4f} vs {shown}")

    exact = offset_hits.get(0, 0)
    print()
    if n_tested == 0:
        print("VERDICT: no positions tested; cannot conclude.")
    elif exact / n_tested > 0.95:
        print("VERDICT: liftover is positionally correct. The hg38 track carries "
              "the hg19\n         value at the chain-mapped coordinate, offset 0.")
    elif max(offset_hits.get(o, 0) for o in OFFSETS if o != 0) > exact:
        worst = max((o for o in OFFSETS if o != 0),
                    key=lambda o: offset_hits.get(o, 0))
        print(f"VERDICT: SYSTEMATIC SHIFT OF {worst:+d} BASES in the hg38 track.")
        print("         The GERP values are misaligned independently of the "
              "lookup bug.")
    else:
        print("VERDICT: inconclusive -- neither clean agreement nor a consistent "
              "shift.\n         Inspect the examples above.")


if __name__ == "__main__":
    main()
