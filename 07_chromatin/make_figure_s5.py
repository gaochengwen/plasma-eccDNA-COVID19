#!/usr/bin/env python3
"""Supplementary Figure S5: cross-cell-type and mappability checks.

Answers Reviewer 3 major comment 4, second point (is the chromatin association
specific to the CD14+ monocyte reference?) and tests whether the H3K9me3
depletion is an alignment-accessibility artefact.

Formatted to the Nature figure specification: 183 mm double column, Arial or its
metric-compatible substitute, 5-7 pt labels, 8 pt bold lowercase panel labels,
palette Hex values from the Nature illustration guide, and non-colour encoding
wherever a group is distinguished.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402

ALLOWED_BP = 2831289217
MARKS = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K9ac", "H3K27me3", "H3K9me3"]
CELL_ORDER = [
    ("CD14_positive_monocyte", "CD14$^+$ monocyte"),
    ("neutrophil", "Neutrophil"),
    ("CD4_positive_alpha_beta_T_cell", "CD4$^+$ T cell"),
    ("B_cell", "B cell"),
    ("upper_lobe_of_left_lung", "Lung"),
    ("endothelial_cell_of_umbilical_vein", "HUVEC"),
]

BAR_PRIMARY = "#6F7B91"        # Nature Grey 4 (medium slate grey)
BAR_SECONDARY = "#C8CEDA"      # Nature Grey 2 (soft light slate grey, non-white)
ZERO_LINE = "#49566D"


def sig_tier(q) -> str:
    if q is None or pd.isna(q):
        return ""
    q = float(q)
    return "****" if q < 1e-4 else "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "NS"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument("--stem", default="figure_s5_crosscell_mappability")
    args = parser.parse_args(argv)

    dirs = base.ensure_dirs(Path(args.outdir))
    t = dirs.tables

    mc = pd.read_csv(t / "ext_multicell_group_stats.tsv", sep="\t")
    mc = mc[mc["mode"] == "full_interval"].copy()
    mc["cell"] = mc["dataset"].str.replace("ENCODE_hg38_", "", regex=False)
    mcq = pd.read_csv(t / "ext_multicell_peak_qc.tsv", sep="\t")
    mcq["cell"] = mcq["dataset"].str.replace("ENCODE_hg38_", "", regex=False)
    lift = pd.read_csv(t / "chipseq_liftover_qc.tsv", sep="\t")
    prim = pd.read_csv(t / "relative_enrichment_group_stats.tsv", sep="\t")
    prim = prim[(prim["analysis_tier"] == "primary") & (prim["mode"] == "full_interval")]
    mapq = pd.read_csv(t / "ext_peak_set_mappability.tsv", sep="\t")
    mapr = pd.read_csv(t / "ext_mappability_group_stats.tsv", sep="\t")
    mapr = mapr[(mapr["analysis_tier"] == "primary") & (mapr["mode"] == "full_interval")]

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt

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
            "hatch.linewidth": 0.4,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "legend.handlelength": 1.6,
            "legend.handletextpad": 0.5,
            "legend.labelspacing": 0.35,
            "legend.borderpad": 0.2,
        }
    )
    LAB, TICK, LEG, PAN = 6.0, 5.5, 5.0, 8.0

    def tidy(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=TICK, length=2.0, pad=1.5)
        ax.yaxis.label.set_size(LAB)
        ax.xaxis.label.set_size(LAB)

    panel_slots: list[tuple] = []

    def panel(ax, letter):
        # Only records the panel; the letter is drawn by place_panel_labels() once the
        # final layout is known. Reading ax.get_position() here would use the default
        # subplot geometry, not the one fig.subplots_adjust() sets below.
        panel_slots.append((ax, letter))

    def place_panel_labels():
        """Draw the panel letters in the margin left of each panel's y-axis.

        Anchored to the y-axis tight bbox (tick labels plus axis label) so the letter
        never lands on the axis text, then shared column-wise. The column takes the
        rightmost of its candidate positions, which keeps the two letters of a column
        aligned without either drifting over the panel to its left.
        """
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        to_fig = fig.transFigure.inverted()
        column_x: dict[int, float] = {}
        for ax, _ in panel_slots:
            box = ax.get_position()
            bb = ax.yaxis.get_tightbbox(renderer)
            axis_left = to_fig.transform((bb.x0, 0))[0] if bb is not None else box.x0
            key = round(box.x0, 3)
            column_x[key] = max(column_x.get(key, -np.inf), axis_left - 0.016)
        for ax, letter in panel_slots:
            box = ax.get_position()
            fig.text(max(column_x[round(box.x0, 3)], 0.004), box.y1 + 0.014, letter,
                     fontsize=PAN, fontweight="bold", va="bottom", ha="left")

    fig, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 122 / 25.4))
    x = np.arange(len(MARKS))

    # (a) within-cohort enrichment across six GRCh38-native cell types
    ax = axes[0][0]
    grid = np.full((len(CELL_ORDER), len(MARKS)), np.nan)
    for i, (cell, _) in enumerate(CELL_ORDER):
        for j, mark in enumerate(MARKS):
            sub = mc[(mc["cell"] == cell) & (mc["mark"] == mark)]
            if len(sub):
                grid[i, j] = float(sub["median_ratio_covid"].iloc[0])
    # Drawn as rectangles rather than imshow: imshow embeds a raster image, and the
    # specification requires vector art with no rasterized elements.
    vmax = float(np.nanmax(np.abs(np.log2(grid))))
    cmap = mpl.colormaps["RdBu_r"]
    norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
    for i in range(len(CELL_ORDER)):
        for j in range(len(MARKS)):
            value = grid[i, j]
            face = "#FFFFFF" if np.isnan(value) else cmap(norm(np.log2(value)))
            ax.add_patch(
                plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=face,
                              edgecolor="white", linewidth=0.4)
            )
            ax.text(j, i, "n.a." if np.isnan(value) else f"{value:.2f}",
                    ha="center", va="center", fontsize=4.4,
                    color="0.35" if np.isnan(value) else "black")
    ax.set_xlim(-0.5, len(MARKS) - 0.5)
    ax.set_ylim(len(CELL_ORDER) - 0.5, -0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(MARKS, rotation=40, ha="right")
    ax.set_yticks(range(len(CELL_ORDER)))
    ax.set_yticklabels([label for _, label in CELL_ORDER])
    ax.tick_params(labelsize=TICK, length=0, pad=1.5)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)

    # Vector colour bar, also drawn as rectangles.
    cax = ax.inset_axes([1.03, 0.0, 0.045, 1.0])
    steps = 64
    for k in range(steps):
        frac = k / (steps - 1)
        cax.add_patch(
            plt.Rectangle((0, frac), 1, 1.0 / steps,
                          facecolor=cmap(frac), edgecolor="none")
        )
    cax.set_xlim(0, 1)
    cax.set_ylim(0, 1)
    cax.set_xticks([])
    cax.set_yticks([0, 0.5, 1])
    cax.set_yticklabels([f"{-vmax:.2f}", "0", f"{vmax:.2f}"])
    cax.yaxis.tick_right()
    cax.tick_params(labelsize=4.6, length=1.5, pad=1.0)
    for side in ("top", "right", "left", "bottom"):
        cax.spines[side].set_visible(False)
    cax.set_ylabel("log$_2$ enrichment", fontsize=LEG, labelpad=1)
    cax.yaxis.set_label_position("right")
    panel(ax, "a")

    # (b) the same six marks from two independent CD14+ peak-set sources
    ax = axes[0][1]
    lifted = [2 ** float(prim[prim["mark"] == m]["median_covid"].iloc[0]) for m in MARKS]
    native = [
        float(mc[(mc["cell"] == "CD14_positive_monocyte") & (mc["mark"] == m)]["median_ratio_covid"].iloc[0])
        for m in MARKS
    ]
    w = 0.36
    ax.bar(x - w / 2, lifted, width=w, color=BAR_PRIMARY, edgecolor="black", linewidth=0.3,
           label="Lifted hg19 (this study)")
    ax.bar(x + w / 2, native, width=w, color=BAR_SECONDARY, edgecolor="black",
           linewidth=0.3, label="ENCODE GRCh38-native")
    ax.axhline(1.0, color=ZERO_LINE, linewidth=0.6, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(MARKS, rotation=40, ha="right")
    ax.set_ylabel("Enrichment ratio, COVID-19")
    ax.set_ylim(0, 2.15)
    ax.legend(frameon=False, fontsize=LEG, loc="upper right")
    tidy(ax)
    panel(ax, "b")

    # (c) how much of the genome each peak set covers
    ax = axes[0][2]
    lift_pct = [
        100 * int(lift[lift["peak_set_id"] == f"CD14_{m}"]["analysis_merged_coverage_bp"].iloc[0]) / ALLOWED_BP
        for m in MARKS
    ]
    nat_pct = [
        100 * int(mcq[(mcq["cell"] == "CD14_positive_monocyte") & (mcq["mark"] == m)]["merged_coverage_bp"].iloc[0])
        / ALLOWED_BP
        for m in MARKS
    ]
    ax.bar(x - w / 2, lift_pct, width=w, color=BAR_PRIMARY, edgecolor="black", linewidth=0.3,
           label="Lifted hg19 (this study)")
    ax.bar(x + w / 2, nat_pct, width=w, color=BAR_SECONDARY, edgecolor="black",
           linewidth=0.3, label="ENCODE GRCh38-native")
    ax.set_xticks(x)
    ax.set_xticklabels(MARKS, rotation=40, ha="right")
    ax.set_ylabel("Genome covered by peak set (%)")
    ax.set_ylim(0, 34)
    ax.legend(frameon=False, fontsize=LEG, loc="upper left")
    tidy(ax)
    panel(ax, "c")

    # (d) is the H3K9me3 depletion an alignment-accessibility artefact?
    ax = axes[1][0]
    order = [f"CD14_{m}" for m in MARKS] + [
        "A549_GSE179184_H3K27ac_official",
        "A549_GSE179184_H3K4me3_official",
        "A549_legacy_H3K27me3_thresholded",
    ]
    labels = MARKS + ["A549 H3K27ac", "A549 H3K4me3", "A549 legacy H3K27me3"]
    vals = []
    for pid in order:
        row = mapq[mapq["peak_set_id"] == pid]
        vals.append(float(row["mappable_fraction"].iloc[0]) if len(row) else np.nan)
    bg = float(mapq["genome_background_mappable_fraction"].iloc[0])
    xs = np.arange(len(order))
    ax.bar(xs, vals, color=BAR_PRIMARY, edgecolor="black", linewidth=0.3)
    ax.axhline(bg, color=ZERO_LINE, linewidth=0.7, linestyle="--")
    # Held above the tallest bar so the label sits in clear space, not across the bars
    # and the dashed reference line it annotates.
    ax.text(len(order) - 0.4, 1.04, f"genome background {bg:.2f}", fontsize=LEG,
            ha="right", va="bottom", color=ZERO_LINE)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=40, ha="right")
    ax.set_ylabel("Fraction uniquely mappable (Umap k36)")
    ax.set_ylim(0, 1.12)
    tidy(ax)
    panel(ax, "d")

    # (e) enrichment before and after restricting to uniquely mappable sequence
    ax = axes[1][1]
    unres = [2 ** float(prim[prim["mark"] == m]["median_covid"].iloc[0]) for m in MARKS]
    res = [float(mapr[mapr["mark"] == m]["median_ratio_covid"].iloc[0]) for m in MARKS]
    ax.bar(x - w / 2, unres, width=w, color=BAR_PRIMARY, edgecolor="black", linewidth=0.3,
           label="All allowed sequence")
    ax.bar(x + w / 2, res, width=w, color=BAR_SECONDARY, edgecolor="black",
           linewidth=0.3, label="Uniquely mappable only")
    ax.axhline(1.0, color=ZERO_LINE, linewidth=0.6, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(MARKS, rotation=40, ha="right")
    ax.set_ylabel("Enrichment ratio, COVID-19")
    ax.set_ylim(0, 2.15)
    ax.legend(frameon=False, fontsize=LEG, loc="upper right")
    tidy(ax)
    panel(ax, "e")

    # (f) between-cohort effect size is not consistent across cell types
    ax = axes[1][2]
    ys = np.arange(len(CELL_ORDER))
    mark_colors = {
        "H3K27ac": "#0272B2",   # Bright Blue
        "H3K4me1": "#459434",   # Vibrant Green
        "H3K4me3": "#D99B00",   # Warm Gold / Amber
        "H3K9ac": "#019AA3",    # Bright Teal
        "H3K27me3": "#A84E94",  # Vibrant Purple
        "H3K9me3": "#D92121",   # Vibrant Red
    }
    markers = {
        "H3K27ac": "o", "H3K4me1": "s", "H3K4me3": "^",
        "H3K9ac": "D", "H3K27me3": "v", "H3K9me3": "p"
    }
    for j, mark in enumerate(MARKS):
        deltas, pos = [], []
        for i, (cell, _) in enumerate(CELL_ORDER):
            sub = mc[(mc["cell"] == cell) & (mc["mark"] == mark)]
            if len(sub):
                deltas.append(float(sub["cliffs_delta_covid_vs_hc"].iloc[0]))
                pos.append(i + (j - 2.5) * 0.13)
        ax.plot(deltas, pos, markers[mark], color=mark_colors[mark], markersize=3.2,
                linestyle="none", markeredgecolor="black", markeredgewidth=0.20, label=mark)
    ax.axvline(0.0, color=ZERO_LINE, linewidth=0.6, zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([label for _, label in CELL_ORDER])
    ax.set_xlabel("Cliff's $\\it{\\delta}$ (COVID-19 vs HC)")
    ax.set_ylim(-1.9, len(CELL_ORDER) - 0.4)
    ax.legend(frameon=False, fontsize=4.6, loc="lower center", ncol=3, columnspacing=0.8)
    tidy(ax)
    panel(ax, "f")

    fig.subplots_adjust(left=0.135, right=0.985, top=0.925, bottom=0.125, wspace=0.58, hspace=0.52)
    place_panel_labels()
    for ext in ("pdf", "svg", "png"):
        fig.savefig(dirs.figures / f"{args.stem}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    base.log(f"wrote figures/{args.stem}.{{pdf,svg,png}}")

    import shutil
    revise_dir = Path("/home/gao/eccDNA/revise")
    if revise_dir.exists():
        for ext in ("pdf", "svg", "png"):
            shutil.copy2(dirs.figures / f"{args.stem}.{ext}", revise_dir / f"Figure_S5.{ext}")
        base.log("copied to revise/Figure_S5.{pdf,svg,png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
