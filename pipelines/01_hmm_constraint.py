#!/usr/bin/env python
"""HMM constraint predictions for RGC-ME and All of Us.

Fits the two-state HMM on RGC-ME chromosome 2 (or loads the saved model if
present, saving it otherwise) and predicts every autosome for RGC-ME and, with
the same model, for All of Us, then inner-joins the two on (chr, pos). Both
cohorts share one position set, so the join keeps every row. See
src/wgs_constraint/hmm.py for the method.

Inputs (config.yaml): gene_annotation, gnomad_exome_coverage, rgc_variants,
aou_variants. Outputs, under results_dir:

  HMM_rgc_0.9_over20_chr2_model.joblib
  HMM_rgc_0.9_over20_chr2_predictions_rgc_wes.tsv.gz
  HMM_rgc_0.9_over20_chr2_predictions_aou_wes.tsv.gz
  HMM_rgc_0.9_over20_chr2_joint_predictions_aou_rgc_wes.tsv.gz
      (chr, pos, prob_0_aou, observation_aou, prob_0_rgc, observation_rgc)

  python pipelines/01_hmm_constraint.py
  python pipelines/01_hmm_constraint.py --chroms chr21,chr22   # quick check
"""

from __future__ import annotations

import argparse
import gzip
from contextlib import ExitStack
import sys
import time
from pathlib import Path

import joblib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgs_constraint import get_config  # noqa: E402
from wgs_constraint import hmm  # noqa: E402

MODEL_STEM = "HMM_rgc_0.9_over20_chr2"
TRAIN_CHROM = "chr2"
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
COHORTS = {"rgc": "rgc_variants", "aou": "aou_variants"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sequence(chrom, cds, covered, variants):
    empty = np.zeros(0, dtype=np.int64)
    return hmm.build_sequence(cds[cds["chr"] == chrom],
                              covered.get(chrom, empty),
                              variants.get(chrom, empty))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=None,
                    help="comma-separated subset to predict (default chr1-chr22); "
                         "output names get the subset as a suffix")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace existing prediction files")
    args = ap.parse_args()

    chroms = ([hmm.normalize_chrom(c) for c in args.chroms.split(",")]
              if args.chroms else AUTOSOMES)
    suffix = "" if chroms == AUTOSOMES else "_" + "_".join(chroms)

    cfg = get_config()
    print(cfg.describe(), flush=True)

    model_path = cfg.result(f"{MODEL_STEM}_model.joblib")
    outputs = {c: cfg.result(f"{MODEL_STEM}_predictions_{c}_wes{suffix}.tsv.gz")
               for c in COHORTS}
    outputs["joint"] = cfg.result(f"{MODEL_STEM}_joint_predictions_aou_rgc_wes{suffix}.tsv.gz")
    existing = [str(p) for p in outputs.values() if p.exists()]
    if existing and not args.overwrite:
        sys.exit("refusing to overwrite (pass --overwrite):\n  " + "\n  ".join(existing))

    train = not model_path.exists()
    needed = sorted(set(chroms) | ({TRAIN_CHROM} if train else set()))

    log("reading CDS intervals ...")
    cds = hmm.load_cds_intervals(cfg.data("gene_annotation"))
    log(f"reading exome coverage ({hmm.COVERAGE_COLUMN} > {hmm.COVERAGE_THRESHOLD}) ...")
    covered = hmm.load_covered_positions(cfg.data("gnomad_exome_coverage"), needed)
    log(f"  {sum(len(v) for v in covered.values()):,} covered positions")

    rgc_by_chrom = {}  # RGC prob_0/observation per chromosome, for the joint table
    for cohort, key in COHORTS.items():
        log(f"[{cohort}] reading variants ...")
        variants = hmm.load_variant_positions(
            cfg.data(key), needed if cohort == "rgc" else chroms)
        log(f"[{cohort}]   {sum(len(v) for v in variants.values()):,} variant positions")

        if cohort == "rgc" and train:
            log(f"training on RGC-ME {TRAIN_CHROM} ...")
            observations, _ = sequence(TRAIN_CHROM, cds, covered, variants)
            model = hmm.fit_hmm(observations)
            joblib.dump(model, model_path)
            log(f"  saved {model_path} ({model.monitor_.iter} EM iterations)")
        elif cohort == "rgc":
            model = joblib.load(model_path)
            log(f"loaded {model_path}")

        paths = [outputs[cohort]] + ([outputs["joint"]] if cohort == "aou" else [])
        with ExitStack() as stack:
            handles = [stack.enter_context(gzip.open(p, "wt")) for p in paths]
            for i, chrom in enumerate(chroms):
                observations, positions = sequence(chrom, cds, covered, variants)
                frame = hmm.predict_chromosome(model, chrom, observations, positions)
                frame.to_csv(handles[0], sep="\t", index=False, header=i == 0)
                log(f"[{cohort}] {chrom}: {len(frame):,} positions")
                if cohort == "rgc":
                    rgc_by_chrom[chrom] = frame.drop(columns="prob_1")
                else:
                    joint = frame.drop(columns="prob_1").merge(
                        rgc_by_chrom.pop(chrom), on=["chr", "pos"], how="inner",
                        suffixes=("_aou", "_rgc"))
                    joint.to_csv(handles[1], sep="\t", index=False, header=i == 0)
        for path in paths:
            log(f"[{cohort}] wrote {path}")
        del variants


if __name__ == "__main__":
    main()
