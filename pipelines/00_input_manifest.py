#!/usr/bin/env python
"""Stage 0 -- checksum every input and pinned artifact the analysis depends on.

Some inputs cannot practically be regenerated: the HMM predictions took hours
and are distributed via figshare, the AoU predictions likewise, and the GERP
liftover used a tool no longer present on the system. Those are treated as
pinned artifacts -- never re-derived, only verified.

This writes `docs/INPUT_MANIFEST.md`, which answers "which version of
AlphaMissense did we use" a year from now without archaeology. Re-run with
--verify to check nothing has changed underneath the analysis.

For gzipped files the md5 of the *decompressed* stream is also recorded, since
gzip metadata differs between writes of identical content and would otherwise
produce spurious mismatches.

  python pipelines/00_input_manifest.py            # write the manifest
  python pipelines/00_input_manifest.py --verify   # check against it
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402

CHUNK = 8 * 1024 * 1024

# (config key or literal path, kind, why it matters)
INPUTS = [
    ("data:gene_annotation", "raw",
     "GENCODE v44 BASIC annotation, 29,570,410 bytes -- confirmed by byte "
     "match against the GENCODE release_44 URL. Not the primary_assembly "
     "file (49,730,393 bytes), and not the chr_patch_hapl_scaff file that "
     "some notebook comments cite. The three give different CDS sets."),
    ("data:epi25_variants", "raw", "Epi25 summary statistics (controlled access)"),
    ("data:alphamissense", "raw",
     "AlphaMissense hg38, ISOFORM-level release -- must be collapsed per variant"),
    ("data:genes4epilepsy", "raw", "Genes4Epilepsy reference gene set"),
    ("data:gerp_bigwig", "pinned",
     "GERP RS hg38, lifted from hg19 with bigWigLiftOver (tool no longer "
     "available). Positional correctness verified 2026-08-09."),
    ("data:All_hg19_RS.bw", "raw", "GERP RS hg19, source of the liftover"),
    ("data:hg19ToHg38.over.chain.gz", "raw", "UCSC chain used for the liftover"),
    ("data:gnomad.v4.0.constraint_metrics.tsv", "raw",
     "gnomAD v4 constraint, source of missense z-score and MTR"),
    ("data:SCHEMA_variant_results_hg38.tsv.gz", "raw",
     "SCHEMA variants, lifted to hg38"),
    ("data:SCHEMA_gene_results.tsv.bgz", "raw", "SCHEMA gene-level results"),
    ("derived:hmm_predictions", "pinned",
     "HMM constraint predictions, RGC-ME. Hours to regenerate; distributed at "
     "doi:10.6084/m9.figshare.27184245. Do not re-derive."),
    ("result:HMM_rgc_ALL_RS_merged_predictions.tsv.gz", "derived",
     "Predictions annotated with corrected GERP (pipelines/05)"),
]


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def md5_decompressed(path: Path) -> str | None:
    if path.suffix not in (".gz", ".bgz"):
        return None
    h = hashlib.md5()
    try:
        with gzip.open(path, "rb") as fh:
            while chunk := fh.read(CHUNK):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def resolve(cfg, spec: str) -> Path:
    kind, _, name = spec.partition(":")
    if kind == "data":
        try:
            return cfg.data(name)
        except KeyError:
            return cfg.data_dir / name
    if kind == "derived":
        return cfg.derived(name)
    return cfg.result(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()
    cfg = get_config()
    manifest_path = Path(__file__).resolve().parents[1] / "docs" / "INPUT_MANIFEST.md"

    rows = []
    for spec, kind, note in INPUTS:
        path = resolve(cfg, spec)
        if not path.exists():
            print(f"  MISSING  {path}", flush=True)
            rows.append((path, kind, note, None, None, None))
            continue
        print(f"  hashing  {path.name} ...", flush=True)
        rows.append((path, kind, note, path.stat().st_size,
                     md5_file(path), md5_decompressed(path)))

    if args.verify:
        if not manifest_path.exists():
            sys.exit("no manifest to verify against")
        recorded = manifest_path.read_text()
        bad = [p.name for p, _, _, _, m, _ in rows if m and m not in recorded]
        print()
        if bad:
            print("CHANGED SINCE THE MANIFEST WAS WRITTEN:")
            for n in bad:
                print(f"  {n}")
            sys.exit(1)
        print("all recorded checksums still match")
        return

    lines = [
        "# Input manifest",
        "",
        f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} "
        f"by `pipelines/00_input_manifest.py`.",
        "",
        "Verify with `python pipelines/00_input_manifest.py --verify`.",
        "",
        "`kind` is **raw** for a downloaded input, **pinned** for something that "
        "cannot practically be regenerated and must be treated as fixed, and "
        "**derived** for something this repository produces.",
        "",
        "For gzipped files the decompressed md5 is the meaningful one; gzip "
        "headers differ between writes of identical content.",
        "",
        "| File | Kind | Bytes | md5 (file) | md5 (decompressed) |",
        "|---|---|---|---|---|",
    ]
    for path, kind, note, size, m, dm in rows:
        size_cell = f"{size:,}" if size else "—"
        lines.append(f"| `{path.name}` | {kind} | {size_cell} | "
                     f"`{m or '—'}` | `{dm or '—'}` |")
    lines += ["", "## Notes", ""]
    for path, kind, note, *_ in rows:
        lines.append(f"- **`{path.name}`** — {note}")
    lines += ["", "## Locations", ""]
    for path, *_ in rows:
        lines.append(f"- `{path}`")
    lines.append("")

    manifest_path.parent.mkdir(exist_ok=True)
    manifest_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nwrote {manifest_path}")


if __name__ == "__main__":
    main()
