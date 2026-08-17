#!/usr/bin/env python3
"""Rebuild Figure 5 with the supplied reference layout and vector colours.

The quantitative values, transformations, group definitions, statistics and
genomic coordinates are read from the existing hg38 analysis outputs.  This
script changes presentation only:

  - panel a stacks the two peak-centred profiles;
  - the two absolute-burden comparisons become panels b and c;
  - the validated hg38 BCL3 and PROCR locus views become panels d and e.

All colours below were extracted from vector objects in the supplied reference
PDF rather than estimated from a raster screenshot.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import make_main_figures as shared  # noqa: E402
import run_chromatin_hg38_analysis as base  # noqa: E402
from make_locus_panels import read_refgene  # noqa: E402


A549_SETS = [
    ("A549_GSE179184_H3K27ac_official", "H3K27ac"),
    ("A549_GSE179184_H3K4me3_official", "H3K4me3"),
]

# Exact reference-PDF palette.
LINE_COLORS = {"COVID": "#C84F50", "HC": "#36617B"}
REFERENCE_GREY = "#B2B2B2"
BOX_COLORS = {"COVID": "#FF8080", "HC": "#8BABD3"}
POINT_COLORS = {"COVID": "#D70000", "HC": "#034E61"}
LOCUS_COLORS = {"COVID": "#C84F50", "HC": "#36617B"}
GENE_COLOR = "#0000B2"
PEAK_TRACKS = [
    ("A549_GSE179184_H3K27ac_official", "A549-ACE2 H3K27ac", "#009999"),
    ("A549_GSE179184_H3K4me3_official", "A549-ACE2 H3K4me3", "#006666"),
    ("CD14_H3K27ac", "CD14$^+$ H3K27ac", "#BFBFBF"),
    ("CD14_H3K27me3", "CD14$^+$ H3K27me3", "#3F3F3F"),
]

FLANK_BP = 1500
GROUP_POSITIONS = {"COVID": 0.0, "HC": 1.0}


def peaks_in_window(path: Path, chrom: str, start: int, end: int) -> List[Tuple[int, int]]:
    """Return merged peak intervals clipped to one locus."""
    found: List[Tuple[int, int]] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            parts = line.split("\t")
            if parts[0] != chrom:
                continue
            left, right = int(parts[1]), int(parts[2])
            if right > start and left < end:
                found.append((max(left, start), min(right, end)))
    return found


def style_axis(ax, tick_size: float, label_size: float) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=tick_size, length=2.0, pad=1.5)
    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)


def add_panel_label(ax, label: str, x: float = -0.18, y: float = 1.10) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=8.0,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def absolute_box_panel(
    ax,
    frame: pd.DataFrame,
    stats: pd.DataFrame,
    mark: str,
    show_ylabel: bool,
) -> None:
    """Reference-panel-B box style with all 39 source values per group."""
    data = []
    for group in ("COVID", "HC"):
        values = (
            frame.loc[(frame["mark"] == mark) & (frame["group"] == group),
                      "absolute_epm_per_total_mapped_alignments"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
        )
        data.append(values)

    bp = ax.boxplot(
        data,
        positions=[GROUP_POSITIONS["COVID"], GROUP_POSITIONS["HC"]],
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 0.75},
        whiskerprops={"color": "black", "linewidth": 0.55},
        capprops={"color": "black", "linewidth": 0.55},
    )
    for patch, group in zip(bp["boxes"], ("COVID", "HC")):
        patch.set_facecolor(BOX_COLORS[group])
        patch.set_edgecolor("black")
        patch.set_linewidth(0.55)

    for group_index, (group, values) in enumerate(zip(("COVID", "HC"), data)):
        rng = np.random.default_rng(5500 + group_index + (0 if mark == "H3K27ac" else 10))
        jitter = rng.uniform(-0.13, 0.13, size=len(values))
        ax.scatter(
            np.full(len(values), GROUP_POSITIONS[group]) + jitter,
            values,
            s=3.0,
            color=POINT_COLORS[group],
            edgecolors="none",
            alpha=0.90,
            zorder=3,
            rasterized=False,
        )

    ax.set_yscale("log")
    low, high = ax.get_ylim()
    bracket_low, bracket_high = high * 1.04, high * 1.16
    ax.plot(
        [0, 0, 1, 1],
        [bracket_low, bracket_high, bracket_high, bracket_low],
        color="black",
        linewidth=0.55,
        clip_on=False,
    )
    q = float(stats.loc[mark, "bh_fdr_all_tests"])
    ax.text(
        0.5,
        high * 1.25,
        shared.sig_tier(q),
        ha="center",
        va="bottom",
        fontsize=5.5,
        fontweight="bold",
    )
    ax.set_ylim(low, high * 2.25)
    ax.set_xlim(-0.55, 1.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["COVID-19", "HC"], rotation=45, ha="right")
    ax.set_title(mark, fontsize=6.0, pad=2.5)
    if show_ylabel:
        ax.set_ylabel("eccDNAs overlapping peak set\nper 10$^6$ mapped alignments")
    style_axis(ax, 5.5, 6.0)


def locus_panel(
    fig,
    cell,
    panel_label: str,
    gene: str,
    record: dict,
    coverage: pd.DataFrame,
    peaks: Dict[str, List[Tuple[int, int]]],
) -> None:
    """Compact hg38 locus view, following the reference D/E track hierarchy."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    start = int(record["start"]) - FLANK_BP
    end = int(record["end"]) + FLANK_BP
    sub = coverage[coverage["gene"] == gene]
    vmax = max(float(sub["mean_eccdna_per_million_mapped"].max()), 1e-9)
    heights = [1.05, 1.05] + [0.20] * len(PEAK_TRACKS) + [0.45]
    inner = cell.subgridspec(len(heights), 1, height_ratios=heights, hspace=0.12)

    first_ax = None
    for row, group in enumerate(("COVID", "HC")):
        ax = fig.add_subplot(inner[row, 0])
        if first_ax is None:
            first_ax = ax
        values = sub[sub["group"] == group].sort_values("bin_start")
        centres = (values["bin_start"].to_numpy(float) + values["bin_end"].to_numpy(float)) / 2000.0
        signal = values["mean_eccdna_per_million_mapped"].to_numpy(float)
        ax.fill_between(
            centres,
            0,
            signal,
            color=LOCUS_COLORS[group],
            linewidth=0,
            alpha=0.96,
            step="mid",
        )
        ax.set_xlim(start / 1000.0, end / 1000.0)
        ax.set_ylim(0, vmax * 1.12)
        ax.set_yticks([0, vmax])
        ax.set_yticklabels(["0", f"{vmax:.4f}"])
        ax.set_xticks([])
        ax.tick_params(labelsize=4.7, length=1.5, pad=1.0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(
            0.992,
            0.86,
            "COVID-19" if group == "COVID" else "HC",
            transform=ax.transAxes,
            fontsize=4.8,
            color=LOCUS_COLORS[group],
            fontweight="bold",
            ha="right",
            va="top",
        )
        if row == 0:
            ax.set_ylabel("eccDNA per 10$^6$\nmapped alignments", fontsize=4.7, labelpad=2)
            ax.yaxis.set_label_coords(-0.075, -0.06)
            ax.set_title(
                rf"$\it{{{gene}}}$   {record['chrom']}:{start:,}-{end:,} (hg38)",
                fontsize=5.8,
                pad=2.5,
            )

    assert first_ax is not None
    add_panel_label(first_ax, panel_label, x=-0.09, y=1.15)

    for track_index, (peak_id, track_label, color) in enumerate(PEAK_TRACKS):
        ax = fig.add_subplot(inner[2 + track_index, 0])
        for left, right in peaks[peak_id]:
            ax.axvspan(left / 1000.0, right / 1000.0, color=color, linewidth=0)
        ax.set_xlim(start / 1000.0, end / 1000.0)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.text(
            0.992,
            0.50,
            track_label,
            transform=ax.transAxes,
            fontsize=4.4,
            ha="right",
            va="center",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.4},
        )

    ax = fig.add_subplot(inner[-1, 0])
    tx_start = max(int(record["start"]), start)
    tx_end = min(int(record["end"]), end)
    ax.plot([tx_start / 1000.0, tx_end / 1000.0], [0.5, 0.5], color=GENE_COLOR, linewidth=0.6)
    for exon_start, exon_end in record["exons"]:
        if exon_end <= start or exon_start >= end:
            continue
        left = max(exon_start, start) / 1000.0
        width = (min(exon_end, end) - max(exon_start, start)) / 1000.0
        ax.add_patch(
            plt.Rectangle((left, 0.18), width, 0.64, facecolor=GENE_COLOR, edgecolor="none")
        )
    span = (tx_end - tx_start) / 1000.0
    for frac in np.linspace(0.08, 0.92, 7):
        xpos = tx_start / 1000.0 + frac * span
        ax.annotate(
            "",
            xy=(xpos + 0.012 * span, 0.5),
            xytext=(xpos, 0.5),
            arrowprops={
                "arrowstyle": "-|>",
                "color": "white",
                "linewidth": 0.35,
                "mutation_scale": 2.8,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )
    ax.set_xlim(start / 1000.0, end / 1000.0)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(labelsize=4.7, length=1.5, pad=1.0)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda value, _pos: f"{value:,.0f}"))
    ax.set_xlabel(f"{record['chrom']} position (kb, hg38)", fontsize=4.8, labelpad=1.2)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.text(
        -0.008,
        0.5,
        rf"$\it{{{gene}}}$ (+)",
        transform=ax.transAxes,
        fontsize=4.6,
        ha="right",
        va="center",
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(SCRIPT_DIR.parent))
    parser.add_argument(
        "--refgene",
        default=str(SCRIPT_DIR.parent / "inputs" / "figure5_refGene_hg38.txt"),
    )
    parser.add_argument("--stem", default="figure5_hg38")
    args = parser.parse_args(argv)

    dirs = base.ensure_dirs(Path(args.outdir))
    tables = dirs.tables
    density = pd.read_csv(tables / "ext_peak_centred_density.tsv", sep="\t")
    absolute = pd.read_csv(tables / "sample_absolute_burden.tsv", sep="\t")
    absolute = absolute[
        (absolute["mode"] == "full_interval")
        & (absolute["dataset"] == "A549_GSE179184_official")
    ].copy()
    stats = pd.read_csv(tables / "absolute_burden_group_stats.tsv", sep="\t")
    stats = stats[
        (stats["mode"] == "full_interval")
        & (stats["dataset"] == "A549_GSE179184_official")
    ].set_index("mark")
    coverage = pd.read_csv(tables / "figure5_cd_locus_coverage.tsv", sep="\t")
    genes = read_refgene(Path(args.refgene), ["BCL3", "PROCR"])

    locus_peaks: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    for gene in ("BCL3", "PROCR"):
        record = genes[gene]
        start, end = int(record["start"]) - FLANK_BP, int(record["end"]) + FLANK_BP
        locus_peaks[gene] = {}
        for peak_id, _, _ in PEAK_TRACKS:
            path = dirs.peaks / (
                f"{base.sanitize_id(peak_id)}.liftover.hg38.canonical.filtered.merged.bed"
            )
            locus_peaks[gene][peak_id] = peaks_in_window(
                path, record["chrom"], start, end
            )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    shared.setup(mpl)
    mpl.rcParams["hatch.linewidth"] = 0.4
    label_size, tick_size, legend_size = 6.0, 5.5, 5.0

    fig = plt.figure(figsize=(183 / 25.4, 162.75 / 25.4))
    outer = GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[0.36, 0.64],
        height_ratios=[1.0, 1.0],
        wspace=0.20,
        hspace=0.28,
        left=0.075,
        right=0.99,
        top=0.96,
        bottom=0.075,
    )

    profile_grid = outer[0, 0].subgridspec(2, 1, hspace=0.46)
    first_profile = None
    for index, (peak_id, mark) in enumerate(A549_SETS):
        ax = fig.add_subplot(profile_grid[index, 0])
        if first_profile is None:
            first_profile = ax
        shared.profile_panel(
            ax,
            shared.density_fold(density, peak_id),
            mark,
            label_size,
            tick_size,
            legend_size,
            show_legend=(index == 0),
            group_colors=LINE_COLORS,
            group_order=("COVID", "HC"),
            expectation_color=REFERENCE_GREY,
            band_alpha=0.18,
        )
        if index == 0:
            legend = ax.get_legend()
            if legend is not None:
                legend.set_loc("upper right")
                legend.set_bbox_to_anchor((1.0, 1.0))
        ax.set_ylabel("eccDNA breakpoint density\n(fold over genome average)")
        if index == 1:
            ax.set_xlabel("Distance from peak centre (kb)")
        style_axis(ax, tick_size, label_size)
    assert first_profile is not None
    add_panel_label(first_profile, "a", x=-0.24, y=1.14)

    box_grid = outer[1, 0].subgridspec(1, 2, wspace=0.62)
    for index, mark in enumerate(("H3K27ac", "H3K4me3")):
        ax = fig.add_subplot(box_grid[0, index])
        absolute_box_panel(ax, absolute, stats, mark, show_ylabel=(index == 0))
        add_panel_label(ax, "b" if index == 0 else "c", x=-0.36, y=1.08)

    locus_panel(
        fig,
        outer[0, 1],
        "d",
        "BCL3",
        genes["BCL3"],
        coverage,
        locus_peaks["BCL3"],
    )
    locus_panel(
        fig,
        outer[1, 1],
        "e",
        "PROCR",
        genes["PROCR"],
        coverage,
        locus_peaks["PROCR"],
    )

    for extension in ("pdf", "svg", "png"):
        fig.savefig(
            dirs.figures / f"{args.stem}.{extension}",
            dpi=600 if extension == "png" else None,
        )
    plt.close(fig)
    base.log(f"wrote figures/{args.stem}.{{pdf,svg,png}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
