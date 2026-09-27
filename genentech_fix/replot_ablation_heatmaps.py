#!/usr/bin/env python
"""Re-render Supplementary Figures S1-S4 with non-overlapping gene labels.

Only the placement of the gene labels changes. Everything else -- data,
n_variants floor, clipping, binning, colours, axes, titles, the labelling
rule -- comes from stage 07 itself: this script imports stage 07 and calls its
own `make_cmp` and `render` unmodified. It intercepts the final `savefig`,
moves the gene labels, and then saves.

Placement follows stage 06 (Figure 5): labels are grouped by the heatmap bin
their point falls in, placed at the bin centre, and stacked vertically in
ascending order of the full-model -log10(p). The one difference from stage 06
is the stacking step: stage 06 uses 2% of the axis limit, which is smaller
than a 6 pt label, so here the step is the measured label height plus 10%.

No model is refitted. The inputs are the p-value tables stage 02 (full model)
and stage 07 (ablations) already wrote to results/.

Checks, all of which must pass or the script exits non-zero:
  * the labels drawn in each panel are exactly the set given by the stage 07
    rule applied to the data;
  * no two gene-label bounding boxes overlap, no gene label overlaps a
    "Group:" label, and every gene label lies inside its own panel;
  * the new PNG has the same pixel dimensions and DPI as the original;
  * the new figure differs from the same figure rendered with the original
    label positions only inside the old and new label boxes (and that
    original-position render is reported against the published _FIXED PNG).

  python genentech_fix/replot_ablation_heatmaps.py
"""

from __future__ import annotations

import importlib.util
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SUFFIX = "_FIXED"
OUT_TAG = "_LABELFIX"

_spec = importlib.util.spec_from_file_location(
    "stage07", ROOT / "pipelines" / "07_moderator_ablations.py")
s07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s07)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.figure  # noqa: E402
from PIL import Image  # noqa: E402

MIN_VARIANTS = 25   # config.yaml min_variants; stage 07 reads the same value
STACK_PAD = 1.10    # stacking step = tallest label height x this
DIFF_PAD_PX = 4     # antialiasing margin around label boxes in the pixel diff


def gene_texts(ax):
    return [t for t in ax.texts if t.get_fontsize() == s07.ANNOT_FONTSIZE]


def group_text(ax):
    return next(t for t in ax.texts if t.get_text().startswith("Group: "))


def overlaps(a, b):
    return (a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1)


def collisions(fig):
    """Every overlapping pair of labels, and every label outside its panel."""
    renderer = fig.canvas.get_renderer()
    boxes, groups = [], []
    problems = []
    for ax in fig.axes[:4]:
        gbox = group_text(ax).get_window_extent(renderer)
        groups.append(gbox)
        for t in gene_texts(ax):
            box = t.get_window_extent(renderer)
            boxes.append((t.get_text(), box))
            if not (ax.bbox.x0 <= box.x0 and box.x1 <= ax.bbox.x1
                    and ax.bbox.y0 <= box.y0 and box.y1 <= ax.bbox.y1):
                problems.append(f"{t.get_text()} outside its panel")
            if overlaps(box, gbox):
                problems.append(f"{t.get_text()} x {group_text(ax).get_text()}")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if overlaps(boxes[i][1], boxes[j][1]):
                problems.append(f"{boxes[i][0]} x {boxes[j][0]}")
    return problems


def label_boxes_px(fig, dpi):
    """Label boxes in the saved (bbox_inches='tight') image's pixel frame."""
    renderer = fig.canvas.get_renderer()
    tight = fig.get_tightbbox(renderer).padded(
        matplotlib.rcParams["savefig.pad_inches"])
    scale = dpi / fig.dpi
    out = []
    for ax in fig.axes[:4]:
        for t in gene_texts(ax):
            b = t.get_window_extent(renderer)
            x0 = (b.x0 - tight.x0 * fig.dpi) * scale
            x1 = (b.x1 - tight.x0 * fig.dpi) * scale
            y0 = (tight.y1 * fig.dpi - b.y1) * scale
            y1 = (tight.y1 * fig.dpi - b.y0) * scale
            out.append((x0, y0, x1, y1))
    return out


def relayout(fig):
    """Stage 06's placement: bin centre, stacked by full-model -log10(p)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes[:4]:
        texts = gene_texts(ax)
        if not texts:
            continue
        limit = ax.get_xlim()[1]
        assert np.isclose(limit, ax.get_ylim()[1])
        edges = np.linspace(0, limit, s07.NUM_BINS + 1)
        centres = 0.5 * (edges[:-1] + edges[1:])
        tallest_px = max(t.get_window_extent(renderer).height for t in texts)
        step = tallest_px * STACK_PAD * limit / ax.bbox.height

        cells = {}
        for t in texts:
            x, y = t.get_position()
            xi = int(np.clip(np.digitize(x, edges) - 1, 0, s07.NUM_BINS - 1))
            yi = int(np.clip(np.digitize(y, edges) - 1, 0, s07.NUM_BINS - 1))
            cells.setdefault((xi, yi), []).append(t)
        for (xi, yi), members in cells.items():
            members.sort(key=lambda t: (t.get_position()[1], t.get_text()))
            total = (len(members) - 1) * step
            for k, t in enumerate(members):
                t.set_position((centres[xi], centres[yi] - total / 2 + k * step))
    fig.canvas.draw()


def expected_labels(cmp):
    """The stage 07 rule, recomputed from the data independently of render."""
    rule = (((cmp["p_with"] < s07.P_THRESHOLD)
             | (cmp["p_without"] < s07.P_THRESHOLD))
            & (np.abs(cmp["delta"]) >= s07.DELTA_THRESHOLD))
    # Captions say "exceeds 1"; stage 07 uses >= 1. Confirm they coincide.
    strict = (((cmp["p_with"] < s07.P_THRESHOLD)
               | (cmp["p_without"] < s07.P_THRESHOLD))
              & (np.abs(cmp["delta"]) > s07.DELTA_THRESHOLD))
    assert rule.equals(strict), "'>= 1' and '> 1' give different label sets"
    sel = cmp[rule]
    return {g: sorted(sel.loc[sel["group"] == g, "gene_name"].astype(str))
            for g in sorted(cmp["group"].unique())}


def main():
    full = pd.read_csv(RESULTS / f"epilepsy_unified_model_pvalues{SUFFIX}.tsv",
                       sep="\t", usecols=s07.KEY_COLS + ["n_variants", "p_unified"])
    full["p_unified"] = pd.to_numeric(full["p_unified"], errors="coerce")

    original_savefig = matplotlib.figure.Figure.savefig
    failed = False

    for stem, _, figure, title, x_label, png in s07.ABLATIONS:
        if figure is None:
            continue
        other = pd.read_csv(
            RESULTS / f"epilepsy_unified_model_pvalues_{stem}{SUFFIX}.tsv", sep="\t")
        cmp = s07.make_cmp(full, other, MIN_VARIANTS)
        expected = expected_labels(cmp)
        orig_png = RESULTS / f"{png}{SUFFIX}.png"
        new_png = RESULTS / f"{png}{SUFFIX}{OUT_TAG}.png"
        report = {}

        def patched_savefig(fig, fname, *args, **kwargs):
            # 1. Original positions: must reproduce the published PNG.
            fig.canvas.draw()
            report["before"] = collisions(fig)
            old_boxes = label_boxes_px(fig, kwargs["dpi"])
            buf = io.BytesIO()
            original_savefig(fig, buf, *args, **kwargs)
            buf.seek(0)
            report["as_published"] = np.asarray(Image.open(buf).convert("RGBA"))
            # 2. Move the labels, check, save.
            relayout(fig)
            report["after"] = collisions(fig)
            report["boxes"] = old_boxes + label_boxes_px(fig, kwargs["dpi"])
            report["drawn"] = {group_text(ax).get_text()[len("Group: "):]:
                               sorted(t.get_text() for t in gene_texts(ax))
                               for ax in fig.axes[:4]}
            return original_savefig(fig, fname, *args, **kwargs)

        matplotlib.figure.Figure.savefig = patched_savefig
        try:
            annotated = s07.render(cmp, title, x_label, new_png)
        finally:
            matplotlib.figure.Figure.savefig = original_savefig

        orig_img = Image.open(orig_png)
        new_img = Image.open(new_png)
        orig_px = np.asarray(orig_img.convert("RGBA"))
        new_px = np.asarray(new_img.convert("RGBA"))

        same_set = (report["drawn"] == expected
                    == {g: sorted(n) for g, n in annotated.items()})
        same_size = orig_img.size == new_img.size
        same_dpi = (np.round(orig_img.info.get("dpi", (0, 0)))
                    == np.round(new_img.info.get("dpi", (0, 0)))).all()
        reproduces = (report["as_published"].shape == orig_px.shape
                      and np.array_equal(report["as_published"], orig_px))
        # Diff against the in-process render at original positions, so the
        # check isolates the label move even if the published PNG were made
        # under a different matplotlib; `reproduces` ties that render to it.
        if report["as_published"].shape == new_px.shape:
            diff = np.any(report["as_published"] != new_px, axis=2)
            mask = np.zeros_like(diff)
            for x0, y0, x1, y1 in report["boxes"]:
                mask[max(0, int(y0) - DIFF_PAD_PX):int(np.ceil(y1)) + DIFF_PAD_PX,
                     max(0, int(x0) - DIFF_PAD_PX):int(np.ceil(x1)) + DIFF_PAD_PX] = True
            stray = int((diff & ~mask).sum())
        else:
            stray = -1

        print(f"\n{figure}  {new_png.name}")
        for g in sorted(expected):
            print(f"  {g:<5} labels: {', '.join(expected[g]) or '-'}")
        print(f"  label set identical to rule (drawn == data == stage 07): {same_set}")
        print(f"  overlaps before: {len(report['before'])}"
              + (f"  [{'; '.join(report['before'])}]" if report["before"] else ""))
        print(f"  overlaps after : {len(report['after'])}"
              + (f"  [{'; '.join(report['after'])}]" if report["after"] else ""))
        print(f"  size {new_img.size} vs original {orig_img.size}; "
              f"dpi {new_img.info.get('dpi')} vs {orig_img.info.get('dpi')}")
        print(f"  original positions reproduce published PNG pixel-for-pixel: {reproduces}")
        print(f"  changed pixels outside label boxes: {stray}")

        if not (same_set and not report["after"] and same_size and same_dpi
                and stray == 0):
            failed = True

    print("\nRESULT:", "FAIL" if failed else "PASS")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
