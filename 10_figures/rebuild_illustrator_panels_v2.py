#!/usr/bin/env python3
"""Redraw the data panels of the three Illustrator main figures from the v2 tables.

Figures 1, 2 and 3 were assembled by hand in Adobe Illustrator and have no
generator in this repository, so the v2 migration could only hand over numbers
(`figures/illustrator_panel_data/`).  This script closes that gap: it redraws
the seven data panels themselves, at the exact geometry they occupy in the
submitted PDFs, so each output can be dropped straight into the .ai file.

  Figure 1b  eccDNA count, COVID-19 vs HC
  Figure 1c  EPM, COVID-19 vs HC
  Figure 1e  per-sample length-class composition
  Figure 2b  per-chromosome eccDNA density
  Figure 2d  gene-element normalised ratio (O/E)
  Figure 3a  eccGene abundance heat map
  Figure 3b  top-30 HC-zero detection oncoprint

Every panel is emitted on a page the size of the *original* figure, with the
panel drawn at its original coordinates and the rest of the page empty.  Placing
such a PDF into the .ai at the artboard origin therefore lands the panel exactly
on top of the one it replaces; no scaling or nudging is needed.

Geometry (axis positions, box widths, marker sizes, colours, font sizes) was
measured off `Submit/Figure_{1,2,3}.pdf` with PyMuPDF and is transcribed below in
PDF points with the origin at the top-left, which is why every axes here uses an
inverted y-limit: measured coordinates can then be used verbatim.

Statistics
----------
Box geometry uses the median-of-halves hinge (lower/upper half excluding the
median) -- the convention already used by `build_illustrator_panel_data.iqr()`
and the one the submitted panels were drawn with (verified: Figure 1c v1 box
edges reproduce `panel_data_summary.json` -> COVID_v1 epm_iqr [1190.1, 5689.8]).

The Figure 2 asterisks are two-sided Wilcoxon rank-sum with tie-corrected normal
approximation and BH-FDR within family, taken *verbatim* from the v1 analysis
code (`analyse/EPM/scripts/run_epm_reanalysis.py`).  That import is gated: the
script first replays the v1 chromosome family and refuses to run unless all 24
p-values and q-values reproduce `analyse/EPM/results/all_group_comparisons.tsv`
exactly.  Only then are the same functions applied to the v2 tables.

Usage
-----
    python3 scripts/local/rebuild_illustrator_panels_v2.py [--outdir DIR]
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import statistics
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Polygon, Rectangle

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
PANEL_DATA = RUN / "figures" / "illustrator_panel_data"
TABLES = RUN / "tables"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",   # keep text editable in Illustrator
    "pdf.fonttype": 42,       # embed TrueType, not Type 3
    "pdf.compression": 6,
})

# Glyph U+2731 (the asterisk the original figures use) is not in Arial; the
# submitted PDFs set it in a symbol font for exactly this reason.
STAR = "✱"
STAR_FONT = "DejaVu Sans"

# ---------------------------------------------------------------- palette ----
RED_POINT = "#D70000"    # Figure 1b/1c COVID-19 markers
TEAL_POINT = "#034E61"   # Figure 1b/1c HC markers
RED_BOX = "#FF8080"      # Figure 1b/1c COVID-19 box fill
BLUE_BOX = "#8BABD3"     # Figure 1b/1c HC box fill
RED_GROUP = "#C84F50"    # Figure 1e/2 COVID-19
TEAL_GROUP = "#36617B"   # Figure 1e/2 HC
BAR_LT1K = "#0070C0"
BAR_1_2K = "#ED7D31"
BAR_GT2K = "#C00000"
GRID_GREY = "#D9D9D9"
ONCO_BLUE = "#3F73BF"
ONCO_GREY = "#D2D2D2"
RIBBON_RED = "#CF3F42"   # Figure 3 group ribbons
RIBBON_TEAL = "#2F627A"

# Diverging map sampled from the colour bar embedded in Submit/Figure_3.pdf.
HEAT_ANCHORS = [
    "#3D70B3", "#4E7DBA", "#608BC0", "#7097C6", "#81A5CD", "#93B2D4",
    "#AEC6DD", "#D3DFEA", "#F6F6F7", "#F5DBD7", "#F3BDB6", "#F0A299",
    "#ED8980", "#EA7067", "#E7574E", "#E43E34", "#E1251B",
]
HEAT_CMAP = LinearSegmentedColormap.from_list("eccdna_abundance", HEAT_ANCHORS, N=512)

# Page sizes of the submitted figures, in points.
PAGE = {
    1: (518.739990234375, 551.468017578125),
    2: (518.739990234375, 465.9119873046875),
    3: (518.739990234375, 709.3330078125),
}

# Sample order printed along Figure 1e, read out of Submit/Figure_1.pdf.
FIG1E_ORDER = [
    "NS120", "NS158", "NS167", "NS169", "NS186", "NS200", "NS230", "NS24",
    "NS254", "NS257", "NS271", "NS281", "NS289", "NS300", "NS336", "NS337",
    "NS363", "NS371", "NS382", "NS385", "NS44", "NS46", "NS61", "NS62",
    "NS82", "XS108", "XS11", "XS120", "XS21", "XS24", "XS2", "XS3", "XS46",
    "XS49", "XS57", "XS5", "XS75", "XS82", "XS87",
    "SP164", "SP175", "SP179", "SP188", "SP218", "SP220", "SP223", "SP225",
    "SP238", "SP244", "SP251", "SP260", "ZXS108", "ZXS11", "ZXS121",
    "ZXS122", "ZXS133", "ZXS134", "ZXS141", "ZXS142", "ZXS144", "ZXS149",
    "ZXS14", "ZXS16", "ZXS178", "ZXS181", "ZXS193", "ZXS1", "ZXS200",
    "ZXS202", "ZXS26", "ZXS38", "ZXS39", "ZXS50", "ZXS55", "ZXS59",
    "ZXS80", "ZXS82", "ZXS85",
]

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ELEMENTS = ["3UTR", "5UTR", "CpG_islands", "Gene2KbD", "Gene2KbU", "exon", "intron"]
ELEMENT_LABEL = {"3UTR": "3UTR", "5UTR": "5UTR", "CpG_islands": "CpG",
                 "Gene2KbD": "Gene2KbD", "Gene2KbU": "Gene2KbU",
                 "exon": "Exon", "intron": "Intron"}
LENGTH_BINS = ["lt200", "200_399", "400_599", "600_999", "1000_1999", "ge2000"]


# ------------------------------------------------------------------ i/o ------
def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_v1_stats_module():
    """Import the v1 statistics helpers and refuse to continue if they drift."""
    path = PROJECT / "analyse" / "EPM" / "scripts" / "run_epm_reanalysis.py"
    spec = importlib.util.spec_from_file_location("epm_v1", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    ref = {}
    for row in read_tsv(PROJECT / "analyse" / "EPM" / "results" / "all_group_comparisons.tsv"):
        ref[(row["fdr_family"], row["feature"])] = (float(row["p_value"]),
                                                    float(row["p_adj_bh"]))

    rows = read_tsv(PROJECT / "analyse" / "EPM" / "input" / "chromosome_distribution.tsv")

    def canon(name: str) -> str:
        return "chr" + name[3:].upper() if name[3:].lower() in ("x", "y") else name

    ps = []
    for chrom in CHROMS:
        sub = [r for r in rows if canon(r["chromosome"]) == chrom]
        x = [float(r["normalized_fraction_per_mb"]) for r in sub if r["group"].startswith("COVID")]
        y = [float(r["normalized_fraction_per_mb"]) for r in sub if r["group"] == "HC"]
        ps.append(module.mann_whitney_two_sided(x, y)[1])
    qs = module.bh_adjust(ps)

    for chrom, p, q in zip(CHROMS, ps, qs):
        rp, rq = ref[("chromosome_distribution", chrom)]
        if not (math.isclose(p, rp, rel_tol=1e-9) and math.isclose(q, rq, rel_tol=1e-9)):
            raise SystemExit(
                f"fidelity gate failed: v1 chromosome family not reproduced at {chrom} "
                f"(p {p!r} vs {rp!r}, q {q!r} vs {rq!r})")
    print("  fidelity gate: v1 chromosome family reproduced exactly (24/24 p and q)")
    return module


def hinges(values) -> tuple[float, float, float]:
    """Lower hinge, median, upper hinge -- halves exclude the median (odd n).

    This is the convention `build_illustrator_panel_data.iqr()` uses and the one
    the submitted boxes were drawn with.
    """
    s = sorted(float(v) for v in values)
    n = len(s)
    return (statistics.median(s[: n // 2]),
            statistics.median(s),
            statistics.median(s[(n + 1) // 2:]))


def group_test(stats, values_by_feature: dict[str, tuple[list[float], list[float]]]):
    """Two-sided Wilcoxon rank-sum + BH within family, using the v1 helpers."""
    features = list(values_by_feature)
    ps, deltas = [], []
    for feature in features:
        x, y = values_by_feature[feature]
        ps.append(stats.mann_whitney_two_sided(x, y)[1])
        deltas.append(stats.cliffs_delta(x, y))
    qs = stats.bh_adjust(ps)
    return {f: {"p": p, "q": q, "cliffs_delta": d}
            for f, p, q, d in zip(features, ps, qs, deltas)}


# ------------------------------------------------------------- canvas --------
def new_page(fignum: int):
    """A figure the size of the original page, in points, y increasing downward."""
    w, h = PAGE[fignum]
    fig = plt.figure(figsize=(w / 72.0, h / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)          # PDF convention: y grows downward
    ax.set_axis_off()
    ax.set_facecolor("none")
    fig.patch.set_alpha(0.0)
    return fig, ax


def line(ax, x0, y0, x1, y1, lw=0.75, color="black", **kw):
    ax.plot([x0, x1], [y0, y1], lw=lw, color=color, solid_capstyle="butt",
            clip_on=False, **kw)


def save(fig, outdir: Path, stem: str):
    for ext in ("pdf", "svg"):
        fig.savefig(outdir / f"{stem}.{ext}", transparent=True)
    fig.savefig(outdir / f"{stem}.png", dpi=600, transparent=False,
                facecolor="white")
    plt.close(fig)
    print(f"  wrote {stem}.pdf/.svg/.png")


# ------------------------------------------------- shared box-plot drawing ---
def draw_box(ax, cx, values, y_of, half_width, cap_half, face, edge_lw,
             edge="black"):
    """Box = hinges, median line, whiskers to the extremes (the original style)."""
    q1, med, q3 = hinges(values)
    lo, hi = min(values), max(values)
    y_q1, y_q3 = y_of(q1), y_of(q3)
    ax.add_patch(Rectangle((cx - half_width, y_q3), 2 * half_width, y_q1 - y_q3,
                           facecolor=face, edgecolor=edge, lw=edge_lw,
                           joinstyle="miter", clip_on=False, zorder=3))
    line(ax, cx - half_width, y_of(med), cx + half_width, y_of(med),
         lw=edge_lw, zorder=4)
    line(ax, cx, y_q3, cx, y_of(hi), lw=edge_lw, zorder=2)
    line(ax, cx - cap_half, y_of(hi), cx + cap_half, y_of(hi), lw=edge_lw, zorder=2)
    line(ax, cx, y_q1, cx, y_of(lo), lw=edge_lw, zorder=2)
    line(ax, cx - cap_half, y_of(lo), cx + cap_half, y_of(lo), lw=edge_lw, zorder=2)


def draw_points(ax, cx, values, y_of, kind, size, jitter, rng, color):
    """Jittered per-sample markers: triangles for COVID-19, discs for HC.

    The submitted markers are filled *and* stroked in the same colour at 0.75 pt,
    which is what gives them their weight; reproduce both.
    """
    offsets = rng.uniform(-jitter, jitter, size=len(values))
    half = size / 2.0
    for value, dx in zip(values, offsets):
        y = y_of(value)
        x = cx + dx
        if kind == "triangle":
            ax.add_patch(Polygon([[x, y - half], [x + half, y + half],
                                  [x - half, y + half]],
                                 closed=True, facecolor=color, edgecolor=color,
                                 lw=0.75, joinstyle="miter", clip_on=False,
                                 zorder=5))
        else:
            ax.add_patch(Circle((x, y), half, facecolor=color, edgecolor=color,
                                lw=0.75, clip_on=False, zorder=5))


def sig_bracket(ax, x0, x1, y_bar, drop, label, star_size, gap=1.3):
    line(ax, x0, y_bar + drop, x0, y_bar, lw=0.75)
    line(ax, x0, y_bar, x1, y_bar, lw=0.75)
    line(ax, x1, y_bar, x1, y_bar + drop, lw=0.75)
    ax.text((x0 + x1) / 2.0, y_bar - gap, label, ha="center", va="bottom",
            fontsize=star_size, family=STAR_FONT, color="black")


# ============================================================== FIGURE 1 =====
def figure1(outdir: Path, report: dict):
    rows = {r["sample_id"]: r for r in read_tsv(PANEL_DATA / "Figure_1bce_per_sample.tsv")}
    covid = [r for r in rows.values() if r["group"] == "COVID-19"]
    hc = [r for r in rows.values() if r["group"] == "HC"]

    # ---- 1b -----------------------------------------------------------------
    counts_c = [float(r["count_v2"]) for r in covid]
    counts_h = [float(r["count_v2"]) for r in hc]
    y_zero, y_top, top_value = 397.31, 272.42, 1_000_000.0

    def y_of(v):
        return y_zero - v * (y_zero - y_top) / top_value

    fig, ax = new_page(1)
    ax.text(5.7, 259.9, "b", ha="left", va="top", fontsize=8, fontweight="bold")
    line(ax, 59.57, 397.94, 59.57, y_top - 0.62, lw=0.75)          # y spine
    for tick in (0, 250_000, 500_000, 750_000, 1_000_000):
        y = y_of(tick)
        line(ax, 59.57, y, 54.29, y, lw=0.75)
        ax.text(53.2, y, f"{tick:d}", ha="right", va="center", fontsize=6)
    line(ax, 58.95, y_zero, 147.29, y_zero, lw=0.75)               # x spine
    for cx in (83.18, 118.71):
        line(ax, cx, y_zero, cx, 402.60, lw=0.75)
    ax.text(80.45, 404.8, "COVID-19", ha="center", va="top", fontsize=6)
    ax.text(117.35, 404.8, "HC", ha="center", va="top", fontsize=6)
    ax.text(16.4, (y_zero + y_top) / 2, "eccDNA count", ha="center", va="center",
            fontsize=6, rotation=90)

    rng = np.random.default_rng(20260809)
    draw_box(ax, 83.18, counts_c, y_of, 11.805, 5.905, RED_BOX, 0.75)
    draw_box(ax, 118.71, counts_h, y_of, 11.805, 5.905, BLUE_BOX, 0.75)
    draw_points(ax, 83.18, counts_c, y_of, "triangle", 2.49, 10.5, rng, RED_POINT)
    draw_points(ax, 118.71, counts_h, y_of, "circle", 2.49, 10.5, rng, TEAL_POINT)
    stats_b = report["figure_1b"]
    sig_bracket(ax, 83.18, 118.71, 273.97, 7.36, stats_b["stars"], 5.29)
    save(fig, outdir, "Figure_1b_v2")

    # ---- 1c -----------------------------------------------------------------
    epm_c = [float(r["epm_v2"]) for r in covid]
    epm_h = [float(r["epm_v2"]) for r in hc]
    y_zero, y_top, top_value = 397.15, 281.25, 10_000.0

    def y_of(v):
        return y_zero - v * (y_zero - y_top) / top_value

    fig, ax = new_page(1)
    ax.text(160.3, 259.9, "c", ha="left", va="top", fontsize=8, fontweight="bold")
    line(ax, 193.04, 397.72, 193.04, y_top - 0.62, lw=0.75)
    for tick in (0, 2500, 5000, 7500, 10000):
        y = y_of(tick)
        line(ax, 193.04, y, 188.14, y, lw=0.75)
        ax.text(186.9, y, f"{tick:d}", ha="right", va="center", fontsize=6)
    line(ax, 192.46, y_zero, 270.50, y_zero, lw=0.75)
    for cx in (214.95, 247.92):
        line(ax, cx, y_zero, cx, 402.05, lw=0.75)
    ax.text(213.4, 403.7, "COVID-19", ha="center", va="top", fontsize=6)
    ax.text(247.0, 403.7, "HC", ha="center", va="top", fontsize=6)
    ax.text(159.6, (y_zero + y_top) / 2, "EPM", ha="center", va="center",
            fontsize=6, rotation=90)

    rng = np.random.default_rng(20260809)
    draw_box(ax, 214.95, epm_c, y_of, 10.95, 5.5, RED_BOX, 0.75)
    draw_box(ax, 247.92, epm_h, y_of, 10.95, 5.5, BLUE_BOX, 0.75)
    draw_points(ax, 214.95, epm_c, y_of, "triangle", 2.30, 9.8, rng, RED_POINT)
    draw_points(ax, 247.92, epm_h, y_of, "circle", 2.30, 9.8, rng, TEAL_POINT)
    sig_bracket(ax, 214.95, 247.92, 268.37, 6.82, report["figure_1c"]["stars"], 5.29)
    save(fig, outdir, "Figure_1c_v2")

    # ---- 1e -----------------------------------------------------------------
    x_left, slot, bar_w = 31.63, 6.1064, 4.898
    y_100, y_0 = 432.70, 503.70

    def y_pct(p):
        return y_0 - p * (y_0 - y_100) / 100.0

    fig, ax = new_page(1)
    ax.text(5.7, 414.5, "e", ha="left", va="top", fontsize=8, fontweight="bold")
    for pct in (0, 20, 40, 60, 80, 100):
        y = y_pct(pct)
        line(ax, x_left - 0.61, y, 508.92, y, lw=0.65, color=GRID_GREY, zorder=0)
        ax.text(27.0, y, f"{pct:d}", ha="right", va="center", fontsize=5.33)
    ax.text(12.1, (y_0 + y_100) / 2, "Percent of eccDNA ( % )", ha="center",
            va="center", fontsize=6, rotation=90)

    for i, sample in enumerate(FIG1E_ORDER):
        row = rows[sample]
        frac = {k: float(row[f"frac_{k}_v2"]) for k in LENGTH_BINS}
        # <1K / 1K-2K / >2K are exact sums of the six tabulated bins.
        parts = [
            (100.0 * (frac["lt200"] + frac["200_399"] + frac["400_599"] + frac["600_999"]), BAR_LT1K),
            (100.0 * frac["1000_1999"], BAR_1_2K),
            (100.0 * frac["ge2000"], BAR_GT2K),
        ]
        # The six tabulated fractions are rounded to 6 dp, so the three plotted
        # classes sum to 100% only up to that rounding.
        assert abs(sum(p for p, _ in parts) - 100.0) < 1e-2, sample
        x0 = x_left + i * slot
        base = 0.0
        for value, colour in parts:
            ax.add_patch(Rectangle((x0, y_pct(base + value)), bar_w,
                                   y_pct(base) - y_pct(base + value),
                                   facecolor=colour, edgecolor="none",
                                   clip_on=False, zorder=2))
            base += value
        ax.text(x0 + bar_w / 2, 506.8, sample, ha="center", va="top",
                fontsize=5.33, rotation=90)

    for x0, x1, colour, label in (
            (31.63, 268.87, RED_GROUP, "COVID-19"),
            (271.05, 508.33, TEAL_GROUP, "HC")):
        line(ax, x0, 530.5, x1, 530.5, lw=3.0, color=colour)
        ax.text((x0 + x1) / 2, 533.4, label, ha="center", va="top", fontsize=6)

    for x_sw, label in ((423.78, "<1K"), (448.71, "1K-2K"), (480.55, ">2K")):
        colour = {"<1K": BAR_LT1K, "1K-2K": BAR_1_2K, ">2K": BAR_GT2K}[label]
        ax.add_patch(Rectangle((x_sw, 423.23), 4.22, 4.22, facecolor=colour,
                               edgecolor="none", clip_on=False))
        ax.text(x_sw + 5.6, 425.34, label, ha="left", va="center", fontsize=6)
    save(fig, outdir, "Figure_1e_v2")


# ============================================================== FIGURE 2 =====
def figure2(outdir: Path, report: dict):
    # ---- 2b -----------------------------------------------------------------
    chrom_rows = read_tsv(PANEL_DATA / "Figure_2b_chromosome_per_sample.tsv")
    by_chrom = {c: {"COVID-19": [], "HC": []} for c in CHROMS}
    for row in chrom_rows:
        for chrom in CHROMS:
            by_chrom[chrom][row["group"]].append(float(row[chrom]) * 100.0)

    y_zero, y_top, top_value = 427.096, 306.94, 6.0
    x0_covid, step, pair_gap = 29.843, 13.6698, 5.144

    def y_of(v):
        return y_zero - v * (y_zero - y_top) / top_value

    fig, ax = new_page(2)
    ax.text(5.7, 301.3, "b", ha="left", va="top", fontsize=8, fontweight="bold")
    line(ax, 26.47, 427.56, 26.47, y_top, lw=1.0)
    for tick in (0, 2, 4, 6):
        y = y_of(tick)
        line(ax, 26.47, y, 22.55, y, lw=1.0)
        ax.text(21.2, y, f"{tick:d}", ha="right", va="center", fontsize=6.82)
    line(ax, 26.004, y_zero, 353.30, y_zero, lw=1.0)
    ax.text(11.0, (y_zero + y_top) / 2, "Percentage eccDNAs per Mb",
            ha="center", va="center", fontsize=6, rotation=90)

    for i, chrom in enumerate(CHROMS):
        cx_c = x0_covid + i * step
        cx_h = cx_c + pair_gap
        draw_box(ax, cx_c, by_chrom[chrom]["COVID-19"], y_of, 1.689, 0.8445,
                 RED_GROUP, 0.487)
        draw_box(ax, cx_h, by_chrom[chrom]["HC"], y_of, 1.689, 0.8445,
                 TEAL_GROUP, 0.487)
        ax.text(cx_c + pair_gap / 2 - 1.0, 432.4, chrom, ha="right", va="top",
                fontsize=6, rotation=45, rotation_mode="anchor")
        if report["figure_2b"][chrom]["q"] < 0.05:
            top = min(y_of(max(by_chrom[chrom]["COVID-19"])),
                      y_of(max(by_chrom[chrom]["HC"])))
            y_bar = top - 5.26
            line(ax, cx_c, y_bar, cx_h, y_bar, lw=0.487)
            ax.text((cx_c + cx_h) / 2, y_bar - 1.1, STAR, ha="center",
                    va="bottom", fontsize=4.48, family=STAR_FONT)

    for y_sw, colour, label in ((313.20, RED_GROUP, "COVID-19"),
                                (325.17, TEAL_GROUP, "HC")):
        ax.add_patch(Rectangle((314.41, y_sw), 9.37, 4.14, facecolor=colour,
                               edgecolor="black", lw=0.487, clip_on=False))
        ax.text(327.0, y_sw + 2.07, label, ha="left", va="center", fontsize=6)
    save(fig, outdir, "Figure_2b_v2")

    # ---- 2d -----------------------------------------------------------------
    elem_rows = read_tsv(PANEL_DATA / "Figure_2d_element_oe_per_sample.tsv")
    by_elem = {e: {"COVID-19": [], "HC": []} for e in ELEMENTS}
    for row in elem_rows:
        by_elem[row["element"]][row["group"]].append(float(row["oe_ratio"]))

    y_zero, y_top, top_value = 426.61, 349.58, 2.0
    x0_hc, step_e, pair_gap_e = 390.97, 16.74, 6.26

    def y_of(v):
        return y_zero - v * (y_zero - y_top) / top_value

    fig, ax = new_page(2)
    ax.text(370.2, 301.4, "d", ha="left", va="top", fontsize=8, fontweight="bold")
    line(ax, 386.79, 427.00, 386.79, y_top, lw=1.0)
    # v2 O/E spans ~0.6-1.6, so the axis runs 0-2 instead of the v1 0-8.  The
    # wider tick labels ("1.5" vs "8") force the two-line axis title 6 pt further
    # left than in the submitted panel; it still clears panel b.
    for tick in (0.0, 0.5, 1.0, 1.5, 2.0):
        y = y_of(tick)
        line(ax, 386.79, y, 383.53, y, lw=1.0)
        label = f"{tick:g}"
        ax.text(382.5, y, label, ha="right", va="center", fontsize=6.82)
    line(ax, 386.41, y_zero, 502.14, y_zero, lw=1.0)
    ax.text(359.3, (y_zero + y_top) / 2, "Normalized ratio of eccDNA",
            ha="center", va="center", fontsize=6, rotation=90)
    ax.text(366.4, (y_zero + y_top) / 2, "in different elements",
            ha="center", va="center", fontsize=6, rotation=90)

    for i, element in enumerate(ELEMENTS):
        cx_h = x0_hc + i * step_e
        cx_c = cx_h + pair_gap_e
        draw_box(ax, cx_h, by_elem[element]["HC"], y_of, 2.076, 1.038,
                 TEAL_GROUP, 0.487)
        draw_box(ax, cx_c, by_elem[element]["COVID-19"], y_of, 2.076, 1.038,
                 RED_GROUP, 0.487)
        ax.text(cx_h + pair_gap_e / 2 - 0.3, 432.2, ELEMENT_LABEL[element],
                ha="right", va="top", fontsize=6, rotation=45,
                rotation_mode="anchor")
        if report["figure_2d"][element]["q"] < 0.05:
            line(ax, cx_h - 0.03, 345.5, cx_c - 0.03, 345.5, lw=0.974)
            ax.text((cx_h + cx_c) / 2, 344.4, STAR, ha="center", va="bottom",
                    fontsize=3.7, family=STAR_FONT)

    for y_sw, colour, label in ((320.862, RED_GROUP, "COVID-19"),
                                (330.206, TEAL_GROUP, "HC")):
        ax.add_patch(Rectangle((474.01, y_sw), 7.79, 3.71, facecolor=colour,
                               edgecolor="black", lw=0.487, clip_on=False))
        ax.text(484.05, y_sw + 1.86, label, ha="left", va="center", fontsize=6)
    save(fig, outdir, "Figure_2d_v2")


# ============================================================== FIGURE 3 =====
def load_matrix(path: Path):
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        genes, values = [], []
        for raw in fh:
            fields = raw.rstrip("\n").split("\t")
            genes.append(fields[0])
            values.append([float(v) for v in fields[1:]])
    return genes, samples, np.asarray(values, dtype=float)


def figure3(outdir: Path, report: dict, matrix_dir: Path):
    groups = {r["sample_id"]: r["group"]
              for r in read_tsv(PANEL_DATA / "Figure_1bce_per_sample.tsv")}

    # ---- 3a -----------------------------------------------------------------
    sig = read_tsv(PANEL_DATA / "Figure_3a_significant_eccgenes.tsv")
    genes_all, samples, ea = load_matrix(matrix_dir / "junction_EA.tsv")
    index = {g: i for i, g in enumerate(genes_all)}

    # Rows: significant eccGenes, ordered so that ubiquitously detected genes sit
    # at the top and COVID-19-restricted genes at the bottom, reproducing the
    # top-to-bottom gradient of the submitted panel.
    ordered = sorted(sig, key=lambda r: (-int(r["hc_detected"]),
                                         -int(r["covid_detected"]), r["Gene"]))
    keep = [r["Gene"] for r in ordered if r["Gene"] in index]
    missing = len(ordered) - len(keep)

    covid_cols = [i for i, s in enumerate(samples) if groups[s] == "COVID-19"]
    hc_cols = [i for i, s in enumerate(samples) if groups[s] == "HC"]
    order_cols = covid_cols + hc_cols
    block = ea[[index[g] for g in keep]][:, order_cols]

    # log1p before the row z-score: EA is 73% zeros and heavily right-skewed, and
    # on the raw scale the row means are dragged so low that most detected COVID-19
    # cells fall below the row mean.  log1p is the transform the eccGene analysis
    # already carries as its abundance sensitivity
    # (tables/eccGene/junction_abundance_log1p_sensitivity.tsv), and it reproduces
    # the continuous red COVID-19 block of the submitted panel.
    logged = np.log1p(block)
    mean = logged.mean(axis=1, keepdims=True)
    sd = logged.std(axis=1, ddof=0, keepdims=True)
    z = np.divide(logged - mean, sd, out=np.zeros_like(logged), where=sd > 0)

    fig, ax = new_page(3)
    ax.text(12.0, 7.5, "a", ha="left", va="top", fontsize=8, fontweight="bold")
    ax.imshow(z, cmap=HEAT_CMAP, vmin=-2.0, vmax=2.0, aspect="auto",
              interpolation="nearest", origin="upper",
              extent=(20.94, 160.19, 186.22, 30.86), zorder=2)
    gradient = np.linspace(0, 1, 512).reshape(1, -1)
    ax.imshow(gradient, cmap=HEAT_CMAP, aspect="auto", interpolation="bilinear",
              extent=(77.02, 140.71, 22.37, 14.38), zorder=2)
    ax.text(36.0, 11.9, "eccDNA", ha="left", va="top", fontsize=6.24)
    ax.text(36.0, 18.7, "abundance", ha="left", va="top", fontsize=6.24)
    ax.text(77.02, 24.7, "Low", ha="center", va="top", fontsize=6.24)
    ax.text(140.71, 24.7, "High", ha="center", va="top", fontsize=6.24)
    for x0, x1, colour, label in ((20.89, 90.53, RIBBON_RED, "COVID-19"),
                                  (90.53, 160.17, RIBBON_TEAL, "HC")):
        ax.add_patch(Rectangle((x0, 189.54), x1 - x0, 2.48, facecolor=colour,
                               edgecolor="none", clip_on=False))
        ax.text((x0 + x1) / 2, 196.9, label, ha="center", va="top", fontsize=7)
    save(fig, outdir, "Figure_3a_v2")
    report["figure_3a"] = {"genes_plotted": len(keep),
                           "genes_not_in_matrix": missing,
                           "row_order": "hc_detected desc, covid_detected desc, symbol",
                           "column_order": "matrix order, COVID-19 block then HC block",
                           "value": "row z-score of log1p(junction eccDNA abundance EA)",
                           "colour_limits": [-2.0, 2.0]}

    # ---- 3b -----------------------------------------------------------------
    top30 = read_tsv(PANEL_DATA / "Figure_3b_top30_detection.tsv")
    genes_det, samples_det, det = load_matrix(matrix_dir / "junction_detected.tsv")
    index_det = {g: i for i, g in enumerate(genes_det)}

    rows_idx = [index_det[r["Gene"]] for r in top30]
    sub = det[rows_idx]
    for row, gene_row in zip(top30, sub):
        n_covid = int(sum(gene_row[i] for i, s in enumerate(samples_det)
                          if groups[s] == "COVID-19"))
        n_hc = int(sum(gene_row[i] for i, s in enumerate(samples_det)
                       if groups[s] == "HC"))
        assert n_covid == int(row["covid_detected_n"]), row["Gene"]
        assert n_hc == int(row["hc_detected_n"]), row["Gene"]

    per_sample = sub.sum(axis=0)
    covid_idx = [i for i, s in enumerate(samples_det) if groups[s] == "COVID-19"]
    hc_idx = [i for i, s in enumerate(samples_det) if groups[s] == "HC"]
    covid_idx.sort(key=lambda i: (-per_sample[i], samples_det[i]))
    hc_idx.sort(key=lambda i: (-per_sample[i], samples_det[i]))
    col_order = covid_idx + hc_idx
    grid = sub[:, col_order]
    counts = per_sample[col_order]

    x_left, x_right = 197.69, 457.06
    col_w = (x_right - x_left) / 78.0
    y_top_grid, y_bot_grid = 51.28, 187.26
    row_h = (y_bot_grid - y_top_grid) / 30.0
    bar_base, bar_top, bar_max = 46.862, 21.435, 30.0

    fig, ax = new_page(3)
    ax.text(173.1, 7.6, "b", ha="left", va="top", fontsize=8, fontweight="bold")
    line(ax, x_left, bar_top, x_left, bar_base, lw=0.5)
    for tick in (0, 10, 20, 30):
        y = bar_base - tick * (bar_base - bar_top) / bar_max
        line(ax, 197.69, y, 195.69, y, lw=0.5)
        ax.text(194.7, y, f"{tick:d}", ha="right", va="center", fontsize=5)
    for i, value in enumerate(counts):
        if value <= 0:
            continue
        height = value * (bar_base - bar_top) / bar_max
        ax.add_patch(Rectangle((x_left + i * col_w + 0.365, bar_base - height),
                               2.594, height, facecolor=ONCO_BLUE,
                               edgecolor="none", clip_on=False, zorder=2))

    ax.add_patch(Rectangle((x_left, y_top_grid), x_right - x_left,
                           y_bot_grid - y_top_grid, facecolor=ONCO_GREY,
                           edgecolor="none", clip_on=False, zorder=2))
    for j in range(30):
        for i in range(78):
            if grid[j, i]:
                ax.add_patch(Rectangle((x_left + i * col_w, y_top_grid + j * row_h),
                                       col_w, row_h, facecolor=ONCO_BLUE,
                                       edgecolor="none", clip_on=False, zorder=3))
    for i in range(79):
        x = x_left + i * col_w
        line(ax, x, y_top_grid, x, y_bot_grid, lw=0.55, color="white", zorder=4)
    for j in range(31):
        y = y_top_grid + j * row_h
        line(ax, x_left, y, x_right, y, lw=0.55, color="white", zorder=4)

    ax.text(460.9, 42.6, "Genes", ha="left", va="center", fontsize=6)
    for j, row in enumerate(top30):
        ax.text(460.9, y_top_grid + (j + 0.5) * row_h, row["Gene"], ha="left",
                va="center", fontsize=5, style="italic")
    for x0, x1, colour, label in ((197.69, 327.37, RIBBON_RED, "COVID-19"),
                                  (327.37, 457.06, RIBBON_TEAL, "HC")):
        ax.add_patch(Rectangle((x0, 190.12), x1 - x0, 2.17, facecolor=colour,
                               edgecolor="none", clip_on=False, zorder=3))
        ax.text((x0 + x1) / 2, 196.3, label, ha="center", va="top", fontsize=7)
    save(fig, outdir, "Figure_3b_v2")
    report["figure_3b"] = {
        "max_genes_per_sample": int(counts.max()),
        "covid_samples_with_any": int((counts[:39] > 0).sum()),
        "hc_samples_with_any": int((counts[39:] > 0).sum()),
        "column_order": "COVID-19 block then HC block, each sorted by descending "
                        "number of the 30 genes detected",
    }


# ================================================================== main =====
# Regions of each submitted page that a rebuilt panel owns.  Used only to make
# the in-context previews: the old panel is painted out, the new one placed on
# top.  The delivered panel PDFs themselves are transparent.
PANEL_BOX = {
    "Figure_1b_v2": (1, (0.0, 252.0, 155.0, 412.0)),
    "Figure_1c_v2": (1, (155.0, 252.0, 280.0, 412.0)),
    "Figure_1e_v2": (1, (0.0, 410.0, 518.74, 545.0)),
    "Figure_2b_v2": (2, (0.0, 295.0, 356.0, 448.0)),
    "Figure_2d_v2": (2, (356.0, 295.0, 518.74, 448.0)),
    "Figure_3a_v2": (3, (0.0, 0.0, 170.0, 206.0)),
    "Figure_3b_v2": (3, (170.0, 0.0, 518.74, 206.0)),
}


def previews(outdir: Path) -> None:
    """Drop each rebuilt panel into a copy of the submitted figure, in place."""
    try:
        import fitz  # PyMuPDF
    except ImportError:
        print("  PyMuPDF not available -- skipping in-context previews")
        return
    by_figure: dict[int, list[str]] = {}
    for stem, (fignum, _) in PANEL_BOX.items():
        by_figure.setdefault(fignum, []).append(stem)
    for fignum, stems in sorted(by_figure.items()):
        source = PROJECT / "Submit" / f"Figure_{fignum}.pdf"
        if not source.is_file():
            print(f"  {source.name} not found -- skipping preview")
            continue
        doc = fitz.open(source)
        page = doc[0]
        for stem in stems:
            page.draw_rect(fitz.Rect(*PANEL_BOX[stem][1]), color=None,
                           fill=(1, 1, 1), overlay=True)
        for stem in stems:
            panel = fitz.open(outdir / f"{stem}.pdf")
            page.show_pdf_page(page.rect, panel, 0, overlay=True)
        out = outdir / f"preview_Figure_{fignum}_with_v2_panels.pdf"
        doc.save(out)
        fitz.open(out)[0].get_pixmap(dpi=200).save(out.with_suffix(".png"))
        print(f"  wrote {out.name} (+ .png)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path,
                        default=RUN / "figures" / "illustrator_panels_v2")
    parser.add_argument("--matrix-dir", type=Path,
                        default=TABLES / "eccGene" / "matrices",
                        help="directory holding junction_EA.tsv and "
                             "junction_detected.tsv from the v2 eccGene run")
    parser.add_argument("--no-preview", action="store_true",
                        help="skip the in-context previews built off Submit/")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("v1 statistics module")
    stats = load_v1_stats_module()

    report: dict = {}

    # Figure 1b/1c group tests, same method as the v1 sample_burden family.
    rows = read_tsv(PANEL_DATA / "Figure_1bce_per_sample.tsv")
    covid = [r for r in rows if r["group"] == "COVID-19"]
    hc = [r for r in rows if r["group"] == "HC"]
    for panel, column in (("figure_1b", "count_v2"), ("figure_1c", "epm_v2")):
        x = [float(r[column]) for r in covid]
        y = [float(r[column]) for r in hc]
        p = stats.mann_whitney_two_sided(x, y)[1]
        stars = "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
        report[panel] = {
            "p": p, "stars": STAR * len(stars) if stars != "ns" else "ns",
            "cliffs_delta": stats.cliffs_delta(x, y),
            "covid_hinges": hinges(x), "hc_hinges": hinges(y),
            "covid_range": [min(x), max(x)], "hc_range": [min(y), max(y)],
        }

    print("Figure 1")
    figure1(args.outdir, report)

    # Figure 2 families.
    chrom_rows = read_tsv(PANEL_DATA / "Figure_2b_chromosome_per_sample.tsv")
    chrom_values = {}
    for chrom in CHROMS:
        chrom_values[chrom] = (
            [float(r[chrom]) * 100 for r in chrom_rows if r["group"] == "COVID-19"],
            [float(r[chrom]) * 100 for r in chrom_rows if r["group"] == "HC"])
    report["figure_2b"] = group_test(stats, chrom_values)

    elem_rows = read_tsv(PANEL_DATA / "Figure_2d_element_oe_per_sample.tsv")
    elem_values = {}
    for element in ELEMENTS:
        sub = [r for r in elem_rows if r["element"] == element]
        elem_values[element] = (
            [float(r["oe_ratio"]) for r in sub if r["group"] == "COVID-19"],
            [float(r["oe_ratio"]) for r in sub if r["group"] == "HC"])
    report["figure_2d"] = group_test(stats, elem_values)

    print("Figure 2")
    figure2(args.outdir, report)

    print("Figure 3")
    figure3(args.outdir, report, args.matrix_dir)

    # ---- traceability -------------------------------------------------------
    stats_path = args.outdir / "panel_statistics_v2.tsv"
    with stats_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["panel", "family", "feature", "covid_hinge_low",
                         "covid_median", "covid_hinge_high", "hc_hinge_low",
                         "hc_median", "hc_hinge_high", "cliffs_delta",
                         "p_value", "q_value_BH", "marked_on_panel"])
        for panel, column in (("Figure_1b", "count_v2"), ("Figure_1c", "epm_v2")):
            key = "figure_1b" if panel == "Figure_1b" else "figure_1c"
            entry = report[key]
            writer.writerow([panel, "sample_burden", column,
                             *[f"{v:.6g}" for v in entry["covid_hinges"]],
                             *[f"{v:.6g}" for v in entry["hc_hinges"]],
                             f"{entry['cliffs_delta']:.6g}",
                             f"{entry['p']:.6g}", "NA",
                             "****" if entry["p"] < 1e-4 else "see p"])
        for chrom in CHROMS:
            e = report["figure_2b"][chrom]
            c_lo, c_med, c_hi = hinges(chrom_values[chrom][0])
            h_lo, h_med, h_hi = hinges(chrom_values[chrom][1])
            writer.writerow(["Figure_2b", "chromosome_distribution", chrom,
                             f"{c_lo:.6g}", f"{c_med:.6g}", f"{c_hi:.6g}",
                             f"{h_lo:.6g}", f"{h_med:.6g}", f"{h_hi:.6g}",
                             f"{e['cliffs_delta']:.6g}", f"{e['p']:.6g}",
                             f"{e['q']:.6g}", "*" if e["q"] < 0.05 else ""])
        for element in ELEMENTS:
            e = report["figure_2d"][element]
            c_lo, c_med, c_hi = hinges(elem_values[element][0])
            h_lo, h_med, h_hi = hinges(elem_values[element][1])
            writer.writerow(["Figure_2d", "gene_feature_normalized_ratio", element,
                             f"{c_lo:.6g}", f"{c_med:.6g}", f"{c_hi:.6g}",
                             f"{h_lo:.6g}", f"{h_med:.6g}", f"{h_hi:.6g}",
                             f"{e['cliffs_delta']:.6g}", f"{e['p']:.6g}",
                             f"{e['q']:.6g}", "*" if e["q"] < 0.05 else ""])
    print(f"  wrote {stats_path.name}")

    def jsonable(obj):
        if isinstance(obj, dict):
            return {k: jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [jsonable(v) for v in obj]
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        return obj

    (args.outdir / "panel_rebuild_report.json").write_text(
        json.dumps(jsonable(report), indent=2, ensure_ascii=False) + "\n")
    print("  wrote panel_rebuild_report.json")

    if not args.no_preview:
        print("In-context previews")
        previews(args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
