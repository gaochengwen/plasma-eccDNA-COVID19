#!/usr/bin/env python3
"""Redraw Figure 3d from the v2 recurrent-interval table.

The submitted panel shows the v1 result: 27 COVID-19-specific exact intervals.
Under the methods-consistent call set the same rule gives 16, which is what the
Results text and the Figure 3 legend already say, so the artwork contradicts
both. The panel is a set of concentric rings, one interval each, so the eleven
dropped intervals cannot be erased individually -- the radial layout depends on
how many there are.

Geometry is measured off `Submit/Figure_3.pdf` and reproduced exactly, so the
output can be placed into the .ai file at the artboard origin with no scaling:

    centre            (176.306, 537.307) pt
    ring annulus      r = 60.826 .. 143.214 pt, outermost ring = most recurrent
    sample cells      39 x 6.4231 deg with 0.5 deg gaps, opening 90.5 deg
    first cell        NS120, at -6.67 .. -0.25 deg, running clockwise
    bar chart         x0 = 183.20 pt, 4.549 pt per sample, ticks at 0, 8, 16
    palette           #C82829 detected, #E5E5E5 not detected

The 16 rings are spread over the same annulus the 27 occupied, and the 16 bars
over the same vertical span, so the panel keeps its original footprint and the
surrounding sample labels stay where they are.

Usage
-----
    python3 scripts/local/rebuild_figure_3d_panel_v2.py [--outdir DIR]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
PANEL_DATA = RUN / "figures" / "illustrator_panel_data"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "pdf.compression": 6,
})

PAGE_3 = (518.739990234375, 709.3330078125)

CENTRE = (176.306, 537.307)
R_INNER, R_OUTER = 60.826, 143.214
CELL_DEG, CELL_GAP_DEG = 6.4231, 0.5
FIRST_CELL_END_DEG = -0.25          # trailing edge of the NS120 cell

BAR_X0 = 183.20
BAR_PT_PER_SAMPLE = 4.549
BAR_TOP = 394.36                    # top of the first bar
BAR_SPAN = 81.85                    # vertical extent the 27 v1 bars occupied
BAR_TICKS = (0, 8, 16)
LABEL_PAD = 2.2                     # gap between a bar and its label
LABEL_PT = 5.0                      # interval-label size, matched to the original
AXIS_PT = 7.0                       # tick-label size, matched to the original

DETECTED = "#C82829"
ABSENT = "#E5E5E5"

# Sample order clockwise around the ring, read off the submitted artwork; it is
# the COVID-19 block of the Figure 1e order.
COVID_ORDER = [
    "NS120", "NS158", "NS167", "NS169", "NS186", "NS200", "NS230", "NS24",
    "NS254", "NS257", "NS271", "NS281", "NS289", "NS300", "NS336", "NS337",
    "NS363", "NS371", "NS382", "NS385", "NS44", "NS46", "NS61", "NS62",
    "NS82", "XS108", "XS11", "XS120", "XS21", "XS24", "XS2", "XS3", "XS46",
    "XS49", "XS57", "XS5", "XS75", "XS82", "XS87",
]


def read_panel() -> tuple[list[dict], list[str]]:
    path = PANEL_DATA / "Figure_3d_recurrent_intervals.tsv"
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    samples = [c for c in rows[0] if c not in
               ("interval", "label", "recurrence_of_39_COVID")]
    if sorted(samples) != sorted(COVID_ORDER):
        raise SystemExit("panel data samples do not match the ring order")
    return rows, samples


def new_page():
    w, h = PAGE_3
    fig = plt.figure(figsize=(w / 72.0, h / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_axis_off()
    fig.patch.set_alpha(0.0)
    return fig, ax


def draw(ax, rows: list[dict]) -> None:
    n = len(rows)
    cx, cy = CENTRE
    pitch = (R_OUTER - R_INNER) / n

    # ---- rings: outermost is the most recurrent interval ------------------
    for i, row in enumerate(rows):
        r_out = R_OUTER - i * pitch
        r_in = r_out - pitch
        for j, sample in enumerate(COVID_ORDER):
            end = FIRST_CELL_END_DEG - j * (CELL_DEG + CELL_GAP_DEG)
            start = end - CELL_DEG
            # the page axes run y-downward, so a wedge drawn with the measured
            # (y-up) angles would come out mirrored; negate to compensate
            ax.add_patch(Wedge(
                (cx, cy), r_out, -end, -start, width=pitch,
                facecolor=DETECTED if row[sample] == "1" else ABSENT,
                edgecolor="none", lw=0, clip_on=False, zorder=2))

    # ---- bar chart --------------------------------------------------------
    bar_pitch = BAR_SPAN / n
    bar_h = bar_pitch * 0.72
    for i, row in enumerate(rows):
        y = BAR_TOP + i * bar_pitch
        width = int(row["recurrence_of_39_COVID"]) * BAR_PT_PER_SAMPLE
        ax.add_patch(plt.Rectangle((BAR_X0, y), width, bar_h,
                                   facecolor=DETECTED, edgecolor="none",
                                   clip_on=False, zorder=3))
        # each label sits immediately to the right of its own bar, as in the
        # submitted panel, so the label column steps with the bar lengths
        ax.text(BAR_X0 + width + LABEL_PAD, y + bar_h / 2, row["label"],
                ha="left", va="center", fontsize=LABEL_PT, color="black",
                clip_on=False, zorder=4)

    axis_y = BAR_TOP + BAR_SPAN + 3.0
    ax.plot([BAR_X0, BAR_X0 + BAR_TICKS[-1] * BAR_PT_PER_SAMPLE],
            [axis_y, axis_y], lw=0.7, color="black", solid_capstyle="butt",
            clip_on=False, zorder=3)
    for t in BAR_TICKS:
        x = BAR_X0 + t * BAR_PT_PER_SAMPLE
        ax.plot([x, x], [axis_y, axis_y + 2.2], lw=0.7, color="black",
                solid_capstyle="butt", clip_on=False, zorder=3)
        ax.text(x, axis_y + 7.5, str(t), ha="center", va="center",
                fontsize=AXIS_PT, color="black", clip_on=False, zorder=4)

    # ---- sample labels ----------------------------------------------------
    for j, sample in enumerate(COVID_ORDER):
        end = FIRST_CELL_END_DEG - j * (CELL_DEG + CELL_GAP_DEG)
        mid = end - CELL_DEG / 2
        rad = math.radians(mid)
        r = R_OUTER + 3.0
        x = cx + r * math.cos(rad)
        y = cy - r * math.sin(rad)
        # labels read outward, flipping on the left half so none are upside down
        flip = math.cos(rad) < 0
        ax.text(x, y, sample, ha="right" if flip else "left", va="center",
                rotation=(mid + 180) if flip else mid,
                rotation_mode="anchor", fontsize=5.0, color="black",
                clip_on=False, zorder=4)


def save(fig, outdir: Path, stem: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg"):
        fig.savefig(outdir / f"{stem}.{ext}", transparent=True)
    fig.savefig(outdir / f"{stem}.png", dpi=600, transparent=False,
                facecolor="white")
    plt.close(fig)
    print(f"  wrote {stem}.pdf/.svg/.png")


def preview(outdir: Path) -> None:
    """Lay the rebuilt panel d onto the figure that already carries v2 a and b.

    The base is `preview_Figure_3_with_v2_panels.pdf`, not the submitted PDF, so
    the result is a Figure 3 whose every data panel is on the v2 call set.
    """
    try:
        import fitz
    except ImportError:
        print("  PyMuPDF not available -- skipping preview")
        return
    base = outdir / "preview_Figure_3_with_v2_panels.pdf"
    if not base.is_file():
        print(f"  {base.name} not found -- skipping preview")
        return
    doc = fitz.open(base)
    page = doc[0]
    page.draw_rect(fitz.Rect(0.0, 385.0, 352.0, 709.34), color=None,
                   fill=(1, 1, 1), overlay=True)
    panel = fitz.open(outdir / "Figure_3d_v2.pdf")
    page.show_pdf_page(page.rect, panel, 0, overlay=True)
    out = outdir / "preview_Figure_3_with_v2_panels_and_3d.pdf"
    doc.save(out)
    fitz.open(out)[0].get_pixmap(dpi=200).save(out.with_suffix(".png"))
    print(f"  wrote {out.name} (+ .png)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path,
                    default=RUN / "figures" / "illustrator_panels_v2")
    args = ap.parse_args()

    rows, _ = read_panel()
    print(f"Figure 3d: {len(rows)} intervals, recurrence "
          f"{rows[-1]['recurrence_of_39_COVID']}-{rows[0]['recurrence_of_39_COVID']}")
    total = sum(sum(int(r[s]) for s in COVID_ORDER) for r in rows)
    print(f"  {total} detections across {len(COVID_ORDER)} COVID-19 samples, 0 in HC")

    fig, ax = new_page()
    draw(ax, rows)
    save(fig, args.outdir, "Figure_3d_v2")
    preview(args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
