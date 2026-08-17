#!/usr/bin/env python3
"""Rebuild main-text Figure 4 and Figure 5 (panels a-c) from the hg38 analysis.

The published panels were built before the coordinate harmonization and before
the enrichment framework was corrected, so they no longer match the Methods or
the Results:

  - Figure 4a and 5a plotted ChIP-seq signal around eccDNA midpoints in hg19, and
    showed HC above COVID-19, which is the relative-enrichment claim that the
    burden analysis retracted.
  - Figure 4b and 5b used a peak-set-internal EPM denominator, giving values 8-85
    times larger than the definition now stated in Methods.
  - Figure 5a and 5b included a SARS-CoV-2 H3K27me3 comparison for which
    GSE179184 provides no processed peak set.
  - Figure 5c and 5d were IGV screenshots labelled with hg19 coordinates, and
    their tracks were labelled "COVID-19 H3K27ac/H3K27me3" although they are
    A549-ACE2 cell-line data, not patient ChIP-seq.

These replacements are drawn directly from the analysed hg38 call set. Panel d of
Figure 5 (PROCR) is dropped: the window contains no CD14+ H3K27ac or H3K27me3
peak, so the locus does not illustrate what the original legend claimed.

Formatted to the Nature figure specification.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402

ALLOWED_BP = 2831289217
CD14_MARKS = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K9ac", "H3K27me3", "H3K9me3"]
A549_SETS = [
    ("A549_GSE179184_H3K27ac_official", "H3K27ac"),
    ("A549_GSE179184_H3K4me3_official", "H3K4me3"),
]

HC_COLOR, COVID_COLOR = "#0272B2", "#EC6F00"
EXPECT_COLOR = "#6F7B91"
ZERO_LINE = "#49566D"
GENE_COLOR = "#253247"
BIN_BP = 25
FLANK_BP = 1500

# Figure 4a deliberately reuses the exact Figure 5a vector palette so cohort
# identity is encoded consistently across the two main figures.  The Figure 5
# source and output are not changed by this Figure-4-only revision.
FIG4_LINE_COLORS = {"COVID": "#C84F50", "HC": "#36617B"}
FIG4_BOX_COLORS = {"COVID": "#FF8080", "HC": "#8BABD3"}
FIG4_POINT_COLORS = {"COVID": "#D70000", "HC": "#034E61"}
FIG4_REFERENCE_GREY = "#B2B2B2"
FIG4_GROUP_OFFSET = 0.18

PEAK_TRACKS = [
    ("A549_GSE179184_H3K27ac_official", "A549-ACE2 H3K27ac", "#C93E3F"),
    ("A549_GSE179184_H3K4me3_official", "A549-ACE2 H3K4me3", "#EC6F00"),
    ("CD14_H3K27ac", "CD14$^+$ H3K27ac", "#459434"),
    ("CD14_H3K27me3", "CD14$^+$ H3K27me3", "#6F7B91"),
]


def sig_tier(q) -> str:
    """Manuscript convention: * q<0.05, ** q<0.01, *** q<0.001, **** q<1e-4."""
    if q is None or pd.isna(q):
        return ""
    q = float(q)
    return "****" if q < 1e-4 else "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "NS"


def setup(mpl) -> None:
    base.setup_matplotlib()
    body = mpl.rcParams["font.family"]
    body = body[0] if isinstance(body, (list, tuple)) else body
    mpl.rcParams.update(
        {
            "mathtext.fontset": "custom",
            "mathtext.rm": body,
            "mathtext.it": f"{body}:italic",
            "mathtext.bf": f"{body}:bold",
            "mathtext.cal": f"{body}:italic",
            "mathtext.tt": body,
            "mathtext.default": "regular",
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "hatch.linewidth": 0.4,
            "legend.handlelength": 1.6,
            "legend.handletextpad": 0.5,
            "legend.labelspacing": 0.35,
            "legend.borderpad": 0.2,
        }
    )


def density_fold(density: pd.DataFrame, peak_set_id: str) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Cohort mean fold-enrichment of eccDNA breakpoints over the genome average.

    density_per_bp_per_centre is breakpoints per bp per peak centre. Dividing by
    each sample's genome-wide breakpoint density (total units / allowed bp) gives
    a fold value that is comparable between cohorts despite the sevenfold
    difference in eccDNA burden, and whose null value is exactly 1.
    """
    sub = density[density["peak_set_id"] == peak_set_id].copy()
    sub["fold"] = sub["density_per_bp_per_centre"] * ALLOWED_BP / sub["total_breakpoint_units"]
    out = {}
    for group in ("COVID", "HC"):
        g = sub[sub["group"] == group]
        piv = g.pivot_table(index="offset_bp", columns="sample_id", values="fold")
        offsets = piv.index.to_numpy(float)
        values = piv.to_numpy(float)
        mean = np.nanmean(values, axis=1)
        sem = np.nanstd(values, axis=1, ddof=1) / np.sqrt(values.shape[1])
        out[group] = (offsets, mean, sem)
    return out


def profile_panel(
    ax,
    folds,
    title,
    label_size,
    tick_size,
    legend_size,
    show_legend=False,
    group_colors=None,
    group_order=("HC", "COVID"),
    expectation_color=EXPECT_COLOR,
    band_alpha=0.25,
) -> None:
    colors = group_colors or {"HC": HC_COLOR, "COVID": COVID_COLOR}
    labels = {"COVID": "COVID-19", "HC": "HC"}
    for group in group_order:
        label = labels[group]
        color = colors[group]
        offsets, mean, sem = folds[group]
        ax.plot(offsets / 1000.0, mean, "-", color=color, linewidth=0.9, label=label)
        ax.fill_between(
            offsets / 1000.0,
            mean - sem,
            mean + sem,
            color=color,
            alpha=band_alpha,
            linewidth=0,
        )
    ax.axhline(1.0, color=expectation_color, linewidth=0.7, linestyle="--", zorder=0,
               label="Genome average" if show_legend else None)
    ax.axvline(0.0, color=ZERO_LINE, linewidth=0.5, zorder=0)
    ax.set_title(title, fontsize=label_size, pad=2.5)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xticklabels(["-2", "", "peak\ncentre", "", "2"])
    ax.tick_params(labelsize=tick_size, length=2.0, pad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_legend:
        ax.legend(frameon=False, fontsize=legend_size, loc="lower left",
                  bbox_to_anchor=(0.0, 0.02))


def box_by_group(
    ax,
    frame,
    value_col,
    marks,
    mark_col="mark",
    log=False,
    widths=0.3,
    group_order=("HC", "COVID"),
    fill_colors=None,
    point_colors=None,
    show_points=False,
    use_hatch=True,
    reference_style=False,
):
    x = np.arange(len(marks))
    fills = fill_colors or {"HC": HC_COLOR, "COVID": COVID_COLOR}
    points = point_colors or fills
    offsets = np.linspace(-FIG4_GROUP_OFFSET, FIG4_GROUP_OFFSET, len(group_order))
    for group_index, (offset, group) in enumerate(zip(offsets, group_order)):
        color = fills[group]
        hatch = "///" if use_hatch and group == "COVID" else ""
        data = []
        for m in marks:
            vals = (
                frame.loc[(frame[mark_col] == m) & (frame["group"] == group), value_col]
                .replace([np.inf, -np.inf], np.nan)
                .dropna()
                .to_numpy(float)
            )
            data.append(vals)
        bp = ax.boxplot(data, positions=x + offset, widths=widths, patch_artist=True,
                        showfliers=False,
                        medianprops={
                            "color": "black",
                            "linewidth": 0.75 if reference_style else 0.7,
                        },
                        whiskerprops={
                            "color": "black",
                            "linewidth": 0.55 if reference_style else 0.5,
                        },
                        capprops={
                            "color": "black",
                            "linewidth": 0.55 if reference_style else 0.5,
                        })
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(1.0 if reference_style else 0.85)
            patch.set_hatch(hatch)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.55 if reference_style else 0.4)
        if show_points:
            for mark_index, (position, vals) in enumerate(zip(x + offset, data)):
                # Fixed seeds make the visual jitter exactly reproducible while
                # preserving every source value unchanged.
                rng = np.random.default_rng(4400 + group_index * 100 + mark_index)
                jitter = rng.uniform(-0.075, 0.075, size=len(vals))
                ax.scatter(
                    np.full(len(vals), position) + jitter,
                    vals,
                    s=3.0,
                    color=points[group],
                    edgecolors="none",
                    alpha=0.88,
                    zorder=3,
                    rasterized=False,
                )
    if log:
        ax.set_yscale("log")
    return x


def locus_panel(fig, gs_cell, locus, coverage, grid, peaks, gene, sizes) -> None:
    """BCL3 locus: cohort eccDNA coverage, peak tracks and gene model, all hg38."""
    import matplotlib.pyplot as plt

    LAB, TICK, TRACK = sizes
    heights = [1.35, 1.35] + [0.30] * len(PEAK_TRACKS) + [0.45]
    inner = gs_cell.subgridspec(len(heights), 1, height_ratios=heights, hspace=0.2)
    centres = (grid[:-1] + grid[1:]) / 2.0 / 1000.0
    vmax = max(float(np.max(coverage["COVID"])), float(np.max(coverage["HC"])), 1e-9)

    for row, (group, color, label) in enumerate(
        (("COVID", COVID_COLOR, "COVID-19"), ("HC", HC_COLOR, "HC"))
    ):
        ax = fig.add_subplot(inner[row, 0])
        ax.fill_between(centres, 0, coverage[group], color=color, linewidth=0, alpha=0.9, step="mid")
        ax.set_ylim(0, vmax * 1.12)
        ax.set_xlim(locus[1] / 1000.0, locus[2] / 1000.0)
        ax.set_yticks([0, vmax])
        ax.set_yticklabels(["0", f"{vmax:.2g}"])
        ax.set_xticks([])
        ax.tick_params(labelsize=TICK, length=1.8, pad=1.2)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.text(0.988, 0.86, label, transform=ax.transAxes, fontsize=TRACK, ha="right",
                va="top", color=color, fontweight="bold")
        if row == 0:
            ax.set_ylabel("eccDNA per 10$^6$\nmapped alignments", fontsize=TRACK, labelpad=2)
            ax.yaxis.set_label_coords(-0.085, -0.05)
            ax.set_title(f"{gene['gene']}   {locus[0]}:{locus[1]:,}-{locus[2]:,} (hg38)",
                         fontsize=LAB, pad=3)

    for offset, (peak_id, peak_label, peak_color) in enumerate(PEAK_TRACKS):
        ax = fig.add_subplot(inner[2 + offset, 0])
        for start, end in peaks[peak_id]:
            ax.axvspan(start / 1000.0, end / 1000.0, color=peak_color, linewidth=0)
        ax.set_xlim(locus[1] / 1000.0, locus[2] / 1000.0)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.text(-0.008, 0.5, peak_label, transform=ax.transAxes, fontsize=TRACK,
                ha="right", va="center")

    ax = fig.add_subplot(inner[len(heights) - 1, 0])
    lo = max(gene["start"], locus[1]) / 1000.0
    hi = min(gene["end"], locus[2]) / 1000.0
    ax.plot([lo, hi], [0.5, 0.5], "-", color=GENE_COLOR, linewidth=0.6)
    for exon_start, exon_end in gene["exons"]:
        if exon_end <= locus[1] or exon_start >= locus[2]:
            continue
        left = max(exon_start, locus[1]) / 1000.0
        width = (min(exon_end, locus[2]) - max(exon_start, locus[1])) / 1000.0
        ax.add_patch(plt.Rectangle((left, 0.18), width, 0.64, facecolor=GENE_COLOR, edgecolor="none"))
    span = hi - lo
    for frac in np.linspace(0.08, 0.92, 8):
        xpos = lo + frac * span
        ax.annotate("", xy=(xpos + (0.012 if gene["strand"] == "+" else -0.012) * span, 0.5),
                    xytext=(xpos, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="white", linewidth=0.4,
                                    mutation_scale=3, shrinkA=0, shrinkB=0))
    ax.set_xlim(locus[1] / 1000.0, locus[2] / 1000.0)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(labelsize=TICK, length=1.8, pad=1.2)
    import matplotlib as mpl
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=6, steps=[1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _p: f"{v:,.0f}"))
    ax.set_xlabel(f"{locus[0]} position (kb, hg38)", fontsize=LAB)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.text(-0.008, 0.5, f"{gene['gene']} ({gene['strand']})", transform=ax.transAxes,
            fontsize=TRACK, ha="right", va="center", style="italic")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument("--refgene", default="/gpfs/data/gao/Covid-eccDNA/eccDNA_abundance3/refGene.txt")
    parser.add_argument(
        "--figure4-only",
        action="store_true",
        help="Stop after rendering Figure 4; the default remains Figures 4 and 5.",
    )
    args = parser.parse_args(argv)

    dirs = base.ensure_dirs(Path(args.outdir))
    t = dirs.tables

    density = pd.read_csv(t / "ext_peak_centred_density.tsv", sep="\t")
    absolute_all = pd.read_csv(t / "sample_absolute_burden.tsv", sep="\t")
    absolute_all = absolute_all[absolute_all["mode"] == "full_interval"]
    relative_all = pd.read_csv(t / "sample_relative_enrichment.tsv", sep="\t")
    relative_all = relative_all[relative_all["mode"] == "full_interval"]
    abs_stats = pd.read_csv(t / "absolute_burden_group_stats.tsv", sep="\t")
    abs_stats = abs_stats[abs_stats["mode"] == "full_interval"]
    relative_stats = pd.read_csv(t / "relative_enrichment_group_stats.tsv", sep="\t")
    relative_stats = relative_stats[relative_stats["mode"] == "full_interval"]

    # Figure 4 is the CD14+ reference figure. Several peak sets share a mark name
    # (CD14, A549 official and A549 legacy all contain H3K27ac), so these must be
    # selected on the dataset, not on the mark.
    absolute = absolute_all[absolute_all["dataset"] == "CD14_reference"]
    relative = relative_all[relative_all["dataset"] == "CD14_reference"]

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch

    setup(mpl)
    LAB, TICK, LEG, PAN = 6.0, 5.5, 5.0, 8.0
    sizes = (LAB, TICK, LEG)

    def panel_label(ax, letter, dx=0.0, dy=0.03):
        box = ax.get_position()
        fig = ax.get_figure()
        fig.text(max(0.004, box.x0 + dx), box.y1 + dy, letter, fontsize=PAN,
                 fontweight="bold", va="bottom", ha="left")

    def tidy(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=TICK, length=2.0, pad=1.5)
        ax.yaxis.label.set_size(LAB)
        ax.xaxis.label.set_size(LAB)

    # ---------------------------------------------------------------- Figure 4
    fig = plt.figure(figsize=(183 / 25.4, 155 / 25.4))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.25, 1.0], hspace=0.42,
                  left=0.075, right=0.985, top=0.94, bottom=0.075)
    top = gs[0, 0].subgridspec(2, 3, hspace=0.62, wspace=0.30)
    first_ax = None
    for k, mark in enumerate(CD14_MARKS):
        ax = fig.add_subplot(top[k // 3, k % 3])
        if first_ax is None:
            first_ax = ax
        folds = density_fold(density, f"CD14_{mark}")
        profile_panel(
            ax,
            folds,
            mark,
            LAB,
            TICK,
            LEG,
            show_legend=(k == 0),
            group_colors=FIG4_LINE_COLORS,
            group_order=("COVID", "HC"),
            expectation_color=FIG4_REFERENCE_GREY,
            band_alpha=0.18,
        )
        if k % 3 == 0:
            ax.set_ylabel("eccDNA breakpoint density\n(fold over genome average)", fontsize=LAB)
        if k // 3 == 1:
            ax.set_xlabel("Distance from peak centre (kb)", fontsize=LAB)
    panel_label(first_ax, "a", dx=-0.062)

    bottom = gs[1, 0].subgridspec(1, 2, wspace=0.24)
    ax = fig.add_subplot(bottom[0, 0])
    x = box_by_group(
        ax,
        absolute,
        "absolute_epm_per_total_mapped_alignments",
        CD14_MARKS,
        log=True,
        group_order=("COVID", "HC"),
        fill_colors=FIG4_BOX_COLORS,
        point_colors=FIG4_POINT_COLORS,
        show_points=True,
        use_hatch=False,
        reference_style=True,
    )
    stats_by_mark = abs_stats[abs_stats["dataset"] == "CD14_reference"].set_index("mark")
    top_y = ax.get_ylim()[1]
    for xi, m in zip(x, CD14_MARKS):
        bracket_y = top_y * 1.12
        bracket_low = top_y * 1.04
        ax.plot(
            [xi - FIG4_GROUP_OFFSET, xi - FIG4_GROUP_OFFSET, xi + FIG4_GROUP_OFFSET, xi + FIG4_GROUP_OFFSET],
            [bracket_low, bracket_y, bracket_y, bracket_low],
            color="black",
            linewidth=0.55,
            clip_on=False,
        )
        ax.text(xi, top_y * 1.19, sig_tier(stats_by_mark.loc[m, "bh_fdr_all_tests"]),
                ha="center", va="bottom", fontsize=5.5, fontweight="bold")
    ax.set_ylim(ax.get_ylim()[0], top_y * 2.6)
    ax.set_xticks(x)
    ax.set_xticklabels(CD14_MARKS, rotation=40, ha="right")
    ax.set_xlim(-0.6, len(CD14_MARKS) - 0.4)
    ax.set_ylabel("eccDNAs overlapping peak set\nper 10$^6$ mapped alignments")
    tidy(ax)
    panel_label(ax, "b", dx=-0.075)

    ax = fig.add_subplot(bottom[0, 1])
    x = box_by_group(
        ax,
        relative,
        "log2_enrichment_ratio",
        CD14_MARKS,
        group_order=("COVID", "HC"),
        fill_colors=FIG4_BOX_COLORS,
        point_colors=FIG4_POINT_COLORS,
        show_points=True,
        use_hatch=False,
        reference_style=True,
    )
    ax.axhline(0.0, color=ZERO_LINE, linewidth=0.7, zorder=0)
    relative_stats_by_mark = (
        relative_stats[relative_stats["dataset"] == "CD14_reference"].set_index("mark")
    )
    low_y, top_y = ax.get_ylim()
    span_y = top_y - low_y
    bracket_low = top_y + 0.035 * span_y
    bracket_high = top_y + 0.070 * span_y
    for xi, mark in zip(x, CD14_MARKS):
        tier = sig_tier(relative_stats_by_mark.loc[mark, "bh_fdr_all_tests"])
        if tier == "NS":
            continue
        ax.plot(
            [
                xi - FIG4_GROUP_OFFSET,
                xi - FIG4_GROUP_OFFSET,
                xi + FIG4_GROUP_OFFSET,
                xi + FIG4_GROUP_OFFSET,
            ],
            [bracket_low, bracket_high, bracket_high, bracket_low],
            color="black",
            linewidth=0.55,
            clip_on=False,
        )
        ax.text(
            xi,
            top_y + 0.090 * span_y,
            tier,
            ha="center",
            va="bottom",
            fontsize=5.5,
            fontweight="bold",
        )
    ax.set_ylim(low_y, top_y + 0.17 * span_y)
    ax.set_xticks(x)
    ax.set_xticklabels(CD14_MARKS, rotation=40, ha="right")
    ax.set_xlim(-0.6, len(CD14_MARKS) - 0.4)
    ax.set_ylabel("log$_2$(observed / matched expected)")
    tidy(ax)
    ax.legend(
        handles=[
            Patch(
                facecolor=FIG4_BOX_COLORS["COVID"],
                edgecolor="black",
                linewidth=0.45,
                label="COVID-19",
            ),
            Patch(
                facecolor=FIG4_BOX_COLORS["HC"],
                edgecolor="black",
                linewidth=0.45,
                label="HC",
            ),
        ],
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        ncol=2,
        frameon=False,
        fontsize=5.5,
        handlelength=1.2,
        handleheight=0.8,
        columnspacing=1.1,
        borderaxespad=0.2,
    )
    panel_label(ax, "c", dx=-0.075)

    for ext in ("pdf", "svg", "png"):
        fig.savefig(dirs.figures / f"figure4_hg38.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    base.log("wrote figures/figure4_hg38.{pdf,svg,png}")

    if args.figure4_only:
        return 0

    # ---------------------------------------------------------------- Figure 5
    genes = {}
    with open(args.refgene, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16 or parts[12] != "BCL3" or parts[2] not in base.CANON_SET:
                continue
            record = {
                "gene": parts[12], "chrom": parts[2], "strand": parts[3],
                "start": int(parts[4]), "end": int(parts[5]),
                "exons": list(zip(
                    [int(v) for v in parts[9].rstrip(",").split(",") if v],
                    [int(v) for v in parts[10].rstrip(",").split(",") if v],
                )),
            }
            if "BCL3" not in genes or (record["end"] - record["start"]) > (
                genes["BCL3"]["end"] - genes["BCL3"]["start"]):
                genes["BCL3"] = record
    gene = genes["BCL3"]
    locus = (gene["chrom"], gene["start"] - FLANK_BP, gene["end"] + FLANK_BP)

    cov = pd.read_csv(t / "figure5_cd_locus_coverage.tsv", sep="\t")
    cov = cov[cov["gene"] == "BCL3"]
    grid = np.arange(locus[1], locus[2] + BIN_BP, BIN_BP)
    coverage = {
        g: cov[cov["group"] == g].sort_values("bin_start")["mean_eccdna_per_million_mapped"].to_numpy(float)
        for g in ("COVID", "HC")
    }
    peaks = {}
    for peak_id, _, _ in PEAK_TRACKS:
        bed = dirs.peaks / f"{base.sanitize_id(peak_id)}.liftover.hg38.canonical.filtered.merged.bed"
        found = []
        with open(bed, "rt", encoding="utf-8") as handle:
            for line in handle:
                parts = line.split("\t")
                if parts[0] != locus[0]:
                    continue
                s, e = int(parts[1]), int(parts[2])
                if e > locus[1] and s < locus[2]:
                    found.append((max(s, locus[1]), min(e, locus[2])))
        peaks[peak_id] = found
        base.log(f"BCL3 / {peak_id}: {len(found)} peaks in window")

    fig = plt.figure(figsize=(183 / 25.4, 140 / 25.4))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[0.85, 1.15], hspace=0.30,
                  left=0.135, right=0.985, top=0.93, bottom=0.075)
    top = gs[0, 0].subgridspec(1, 3, wspace=0.34)
    first_ax = None
    for k, (peak_id, mark) in enumerate(A549_SETS):
        ax = fig.add_subplot(top[0, k])
        if first_ax is None:
            first_ax = ax
        profile_panel(ax, density_fold(density, peak_id), f"A549-ACE2 {mark}", LAB, TICK, LEG,
                      show_legend=(k == 0))
        ax.set_xlabel("Distance from peak centre (kb)", fontsize=LAB)
        if k == 0:
            ax.set_ylabel("eccDNA breakpoint density\n(fold over genome average)", fontsize=LAB)
    panel_label(first_ax, "a", dx=-0.078)

    ax = fig.add_subplot(top[0, 2])
    a549_abs = absolute_all[absolute_all["peak_set_id"].isin([p for p, _ in A549_SETS])].copy()
    a549_abs["label"] = a549_abs["mark"]
    labels = [m for _, m in A549_SETS]
    x = box_by_group(ax, a549_abs, "absolute_epm_per_total_mapped_alignments", labels,
                     mark_col="label", log=True, widths=0.26)
    a549_stats = abs_stats[abs_stats["peak_set_id"].isin([p for p, _ in A549_SETS])].set_index("mark")
    top_y = ax.get_ylim()[1]
    for xi, m in zip(x, labels):
        ax.text(xi, top_y * 1.15, sig_tier(a549_stats.loc[m, "bh_fdr_all_tests"]),
                ha="center", va="bottom", fontsize=5.5, fontweight="bold")
    ax.set_ylim(ax.get_ylim()[0], top_y * 3.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylabel("eccDNAs overlapping peak set\nper 10$^6$ mapped alignments")
    tidy(ax)
    panel_label(ax, "b", dx=-0.085)

    locus_ax_holder = fig.add_subplot(gs[1, 0])
    locus_ax_holder.axis("off")
    locus_panel(fig, gs[1, 0], locus, coverage, grid, peaks, gene, sizes)
    panel_label(locus_ax_holder, "c", dx=-0.082, dy=0.005)

    for ext in ("pdf", "svg", "png"):
        fig.savefig(dirs.figures / f"figure5_hg38.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    base.log("wrote figures/figure5_hg38.{pdf,svg,png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
