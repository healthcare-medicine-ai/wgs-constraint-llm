#!/usr/bin/env python
"""
Diagnostic: locate the source of row duplication in the epilepsy meta-regression
input, and quantify how much it inflates n_variants.

Emits AGGREGATE COUNTS ONLY. No variant-level rows are printed.

AlphaMissense layout (confirmed from the file on disk):
    #CHROM POS REF ALT genome uniprot_id transcript_id protein_variant
    am_pathogenicity am_class
The header line begins with '#', so it is consumed as a comment; columns are
addressed positionally and validated against the first data row.

Each section is independent -- a failure in one is reported and the rest
still run.
"""

import gzip
import traceback
from itertools import chain

import pandas as pd

DATA = "/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/data/"
RESULTS = ("/oak/stanford/groups/mrivas/projects/wgs-constraint-llm/"
           "osthoag/wgs-constraint-llm/results/")

MERGED = RESULTS + "constraint_gerp_am_epi25_variants.tsv.gz"
ALPHAMISSENSE = DATA + "AlphaMissense_hg38.tsv.gz"

MIN_VARIANTS = 25

I_CHR, I_POS, I_REF, I_ALT = 0, 1, 2, 3
I_UNIPROT, I_TX = 5, 6
N_EXPECTED_FIELDS = 10


def rule(title):
    print("\n" + "=" * 70, flush=True)
    print(title, flush=True)
    print("=" * 70, flush=True)


# ---------------------------------------------------------------------------
# 1. AlphaMissense: how many rows per variant, and why?
# ---------------------------------------------------------------------------
def section_alphamissense():
    rule("1. AlphaMissense duplication per (chr, pos, ref, alt)")

    with gzip.open(ALPHAMISSENSE, "rt") as fh:
        last_comment = None
        first_data = None
        for line in fh:
            if line.startswith("#"):
                last_comment = line.rstrip("\n")
                continue
            first_data = line.rstrip("\n")
            break

        if first_data is None:
            print("ERROR: no data rows found.")
            return

        print(f"header (last comment line): {last_comment}")
        fields = first_data.split("\t")
        print(f"first data row has {len(fields)} fields")

        if len(fields) != N_EXPECTED_FIELDS or fields[4] != "hg38":
            print(f"WARNING: unexpected layout -> {fields}")
            print("Positional column assumptions may be wrong. Continuing anyway.")

        n_rows = 0
        n_unique = 0
        n_dup_rows = 0
        n_keys_with_dups = 0
        max_run = 1
        out_of_order = 0

        # Of the variants that appear more than once, how many span more than
        # one uniprot_id? That separates "isoforms of one protein" from
        # "overlapping genes" -- they call for different collapse rules.
        n_multi_uniprot = 0

        prev_key = None
        prev_sortable = None
        run = 0
        run_uniprots = set()

        def close_run():
            nonlocal n_keys_with_dups, n_multi_uniprot, max_run
            if run > 1:
                n_keys_with_dups += 1
                if len(run_uniprots) > 1:
                    n_multi_uniprot += 1
                if run > max_run:
                    max_run = run

        # Re-include the first data row, then stream the rest lazily.
        for line in chain([first_data + "\n"], fh):
            f = line.rstrip("\n").split("\t")
            if len(f) < N_EXPECTED_FIELDS:
                continue
            key = (f[I_CHR], f[I_POS], f[I_REF], f[I_ALT])
            n_rows += 1

            if key == prev_key:
                run += 1
                n_dup_rows += 1
                run_uniprots.add(f[I_UNIPROT])
            else:
                close_run()
                sortable = (f[I_CHR], int(f[I_POS]))
                if prev_sortable is not None and sortable < prev_sortable:
                    out_of_order += 1
                prev_sortable = sortable
                prev_key = key
                run = 1
                run_uniprots = {f[I_UNIPROT]}
                n_unique += 1

        close_run()

    print(f"\ntotal rows                       : {n_rows:,}")
    print(f"distinct (chr,pos,ref,alt)       : {n_unique:,}")
    print(f"duplicate rows                   : {n_dup_rows:,}")
    print(f"variants appearing >1 time       : {n_keys_with_dups:,}")
    print(f"  ...of those, spanning >1 uniprot: {n_multi_uniprot:,}")
    print(f"max rows for a single variant    : {max_run}")
    print(f"rows per variant (mean)          : {n_rows / max(n_unique, 1):.4f}")
    print(f"sort-order violations            : {out_of_order:,}"
          f"{'  <-- SCAN UNRELIABLE' if out_of_order else '  (file is sorted)'}")

    if n_dup_rows == 0:
        print("\nVERDICT: one row per variant. AlphaMissense is NOT a fan-out")
        print("         source -- the duplication comes from somewhere else.")
    else:
        pct = 100 * n_keys_with_dups / max(n_unique, 1)
        print(f"\nVERDICT: AlphaMissense IS a fan-out source "
              f"({pct:.2f}% of variants duplicated).")
        if n_multi_uniprot > n_keys_with_dups * 0.5:
            print("         Mostly OVERLAPPING GENES (multiple uniprot ids), not")
            print("         isoforms -- collapsing needs a gene-aware rule.")
        else:
            print("         Mostly ISOFORMS of a single protein -- a simple")
            print("         per-variant collapse (max/mean) is defensible.")


# ---------------------------------------------------------------------------
# 2. The merged table actually consumed by the regression
# ---------------------------------------------------------------------------
def section_merged():
    rule("2. Merged regression input: where duplication actually shows up")

    merged = pd.read_csv(
        MERGED, sep="\t",
        usecols=["chr", "pos", "ref", "alt", "gene_id", "group"],
        low_memory=False,
    )

    variant_key = ["chr", "pos", "ref", "alt"]
    group_key = variant_key + ["gene_id", "group"]

    n_total = len(merged)
    n_variants_distinct = merged[variant_key].drop_duplicates().shape[0]
    n_group_distinct = merged[group_key].drop_duplicates().shape[0]

    print(f"total rows                              : {n_total:,}")
    print(f"distinct (chr,pos,ref,alt)              : {n_variants_distinct:,}")
    print(f"distinct (variant, gene_id, group)      : {n_group_distinct:,}")
    print(f"inflation factor                        : "
          f"{n_total / max(n_group_distinct, 1):.4f}")

    if n_total == n_group_distinct:
        print("\nVERDICT: no duplication -- n_variants is already a true count.")
    else:
        extra = n_total - n_group_distinct
        print(f"\nVERDICT: {extra:,} duplicate rows "
              f"({100 * extra / n_total:.2f}% of the table).")
        print("         Double-counted in n_variants AND double-weighted in WLS.")

    return merged, variant_key, group_key


# ---------------------------------------------------------------------------
# 3. Impact on the min_variants inclusion filter
# ---------------------------------------------------------------------------
def section_filter(merged, group_key):
    rule(f"3. Effect on the min_variants >= {MIN_VARIANTS} inclusion filter")

    as_rows = merged.groupby(["gene_id", "group"], sort=False).size()
    as_distinct = (merged.drop_duplicates(group_key)
                         .groupby(["gene_id", "group"], sort=False).size())

    cmp = pd.DataFrame({"rows": as_rows, "distinct": as_distinct}).fillna(0)

    pass_rows = int((cmp["rows"] >= MIN_VARIANTS).sum())
    pass_distinct = int((cmp["distinct"] >= MIN_VARIANTS).sum())
    inflated = int((cmp["rows"] > cmp["distinct"]).sum())

    print(f"gene x group combinations               : {len(cmp):,}")
    print(f"pass filter, counting ROWS (current)    : {pass_rows:,}")
    print(f"pass filter, counting VARIANTS (fixed)  : {pass_distinct:,}")
    print(f"admitted ONLY by duplicate rows         : {pass_rows - pass_distinct:,}")
    print(f"gene-groups with inflated n_variants    : {inflated:,} "
          f"({100 * inflated / max(len(cmp), 1):.2f}%)")


if __name__ == "__main__":
    try:
        section_alphamissense()
    except Exception:
        print("\nSECTION 1 FAILED:")
        traceback.print_exc()

    merged = group_key = None
    try:
        merged, _, group_key = section_merged()
    except Exception:
        print("\nSECTION 2 FAILED:")
        traceback.print_exc()

    if merged is not None:
        try:
            section_filter(merged, group_key)
        except Exception:
            print("\nSECTION 3 FAILED:")
            traceback.print_exc()

    print("\nDone. All output above is aggregate counts only.", flush=True)
