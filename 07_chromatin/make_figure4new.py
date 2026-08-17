#!/usr/bin/env python3
"""Build the main reference-chromatin figure from the corrected hg38 analysis.

The figure follows one evidence sequence:

  a. CD14+ peak-centred breakpoint profiles;
  b. CD14+ relative enrichment and absolute overlap abundance;
  c. reproduction across GRCh38-native reference biosamples;
  d. high-mappability sensitivity analysis;
  e. SARS-CoV-2-infected A549-ACE2 reference results;
  f. cohort coefficients before and after eccDNA-burden adjustment; and
  g. enrichment in relation to detected eccDNA burden.

All displayed values are read from the deposited analysis tables. The script
does not refit models, redefine localization, or alter any statistical result.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


ALLOWED_BP = 2_831_289_217
MARKS = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K9ac", "H3K27me3", "H3K9me3"]
A549_PEAKS = {
    "H3K27ac": "A549_GSE179184_H3K27ac_official",
    "H3K4me3": "A549_GSE179184_H3K4me3_official",
}
CELL_ORDER = [
    ("CD14_positive_monocyte", "CD14$^+$ monocyte"),
    ("neutrophil", "Neutrophil"),
    ("CD4_positive_alpha_beta_T_cell", "CD4$^+$ T cell"),
    ("B_cell", "B cell"),
    ("upper_lobe_of_left_lung", "Lung"),
    ("endothelial_cell_of_umbilical_vein", "HUVEC"),
]
CATEGORIES = ("Humoral immune response", "Endothelial/coagulation")
SELECTED_GENES = ("HLA-E", "MCAM")

# Exact Nature scientific-illustration palette values.
COVID = "#C93E3F"
COVID_LIGHT = "#FAD0CE"
HC = "#0272B2"
HC_LIGHT = "#C8E7FB"
TEAL = "#019AA3"
TEAL_LIGHT = "#CCE7EE"
PURPLE = "#A84E94"
PURPLE_DARK = "#792C74"
GREEN = "#459434"
ORANGE = "#EC6F00"
YELLOW = "#CCA02C"
GREY_1 = "#E6E6ED"
GREY_2 = "#C8CEDA"
GREY_3 = "#99A3B4"
GREY_4 = "#6F7B91"
GREY_5 = "#49566D"
GREY_6 = "#253247"
STONE_1 = "#F7F5EF"


def read_table(path: Path, required: Iterable[str]) -> pd.DataFrame:
    """Read a TSV and fail early if its schema no longer matches."""
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    missing = set(required) - set(frame.columns)
    if missing:
        raise RuntimeError(f"{path.name} lacks columns: {', '.join(sorted(missing))}")
    return frame


def require_rows(frame: pd.DataFrame, expected: int, label: str) -> None:
    if len(frame) != expected:
        raise RuntimeError(f"{label}: expected {expected} rows, found {len(frame)}")


def density_fold(
    density: pd.DataFrame, peak_set_id: str
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Return cohort mean ± SEM breakpoint density relative to genome average."""
    subset = density[density["peak_set_id"].eq(peak_set_id)].copy()
    subset["fold"] = (
        subset["density_per_bp_per_centre"]
        * ALLOWED_BP
        / subset["total_breakpoint_units"]
    )
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for group in ("COVID", "HC"):
        group_frame = subset[subset["group"].eq(group)]
        pivot = group_frame.pivot_table(
            index="offset_bp", columns="sample_id", values="fold"
        ).sort_index()
        if pivot.shape[1] != 39:
            raise RuntimeError(
                f"{peak_set_id}/{group}: expected 39 samples, found {pivot.shape[1]}"
            )
        values = pivot.to_numpy(float)
        out[group] = (
            pivot.index.to_numpy(float),
            np.nanmean(values, axis=1),
            np.nanstd(values, axis=1, ddof=1) / np.sqrt(values.shape[1]),
        )
    return out


def quantiles(values: pd.Series) -> tuple[float, float, float]:
    clean = pd.to_numeric(values, errors="raise").dropna().to_numpy(float)
    if len(clean) != 39:
        raise RuntimeError(f"expected 39 sample values, found {len(clean)}")
    q1, median, q3 = np.quantile(clean, [0.25, 0.50, 0.75])
    return float(q1), float(median), float(q3)


def setup_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 5.8,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.cal": "Arial:italic",
            "mathtext.tt": "Arial",
            "mathtext.default": "regular",
            "axes.linewidth": 0.5,
            "axes.labelsize": 5.7,
            "axes.titlesize": 6.0,
            "xtick.labelsize": 5.2,
            "ytick.labelsize": 5.2,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "legend.fontsize": 5.1,
            "legend.frameon": False,
            "legend.handlelength": 1.5,
            "legend.handletextpad": 0.45,
            "legend.labelspacing": 0.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "none",
        }
    )


def tidy(ax, left: bool = True, bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(length=2.0, width=0.5, pad=1.4)


def panel_label(fig, ax, label: str, dx: float = -0.030, dy: float = 0.008) -> None:
    box = ax.get_position()
    fig.text(
        box.x0 + dx,
        box.y1 + dy,
        label,
        fontsize=8.0,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="black",
    )


def q_to_stars(value: float) -> str:
    """Return the manuscript-wide BH-q significance tier."""
    if value < 0.0001:
        return "****"
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return ""


def summary_errorbar(
    ax,
    y: float,
    median: float,
    low: float,
    high: float,
    color: str,
    marker: str,
    label: str | None = None,
    zorder: int = 3,
) -> None:
    ax.errorbar(
        median,
        y,
        xerr=np.array([[median - low], [high - median]]),
        fmt=marker,
        markersize=3.7,
        markerfacecolor=color,
        markeredgecolor="white",
        markeredgewidth=0.35,
        color=color,
        ecolor=color,
        elinewidth=0.75,
        capsize=1.7,
        capthick=0.65,
        label=label,
        zorder=zorder,
    )


def draw_profile_panel(fig, cell, density: pd.DataFrame) -> None:
    sub = cell.subgridspec(2, 3, hspace=0.58, wspace=0.31)
    first = None
    for index, mark in enumerate(MARKS):
        ax = fig.add_subplot(sub[index // 3, index % 3])
        if first is None:
            first = ax
        folds = density_fold(density, f"CD14_{mark}")
        for group, color, label in (
            ("COVID", COVID, "COVID-19"),
            ("HC", HC, "HC"),
        ):
            offsets, mean, sem = folds[group]
            x = offsets / 1_000.0
            ax.plot(x, mean, color=color, linewidth=0.85, label=label)
            ax.fill_between(
                x, mean - sem, mean + sem, color=color, alpha=0.18, linewidth=0
            )
        ax.axhline(
            1.0,
            color=GREY_4,
            linewidth=0.55,
            linestyle=(0, (3, 2)),
            label="Genome average" if index == 0 else None,
            zorder=0,
        )
        ax.axvline(0.0, color=GREY_3, linewidth=0.45, zorder=0)
        ax.set_title(mark, pad=1.8)
        ax.set_xlim(-2.0, 2.0)
        ax.set_ylim(0.72, 1.45)
        ax.set_xticks([-2, 0, 2])
        ax.set_xticklabels(["−2", "peak\ncentre", "2"])
        if index % 3 == 0:
            ax.set_ylabel("Breakpoint density\n(fold over genome average)")
        if index // 3 == 1:
            ax.set_xlabel("Distance (kb)", labelpad=1.1)
        tidy(ax)
        if index == 0:
            ax.legend(
                loc="lower left",
                bbox_to_anchor=(0.0, 0.015),
                fontsize=5.0,
                borderaxespad=0,
            )
    assert first is not None
    panel_label(fig, first, "a", dx=-0.044, dy=0.006)


def draw_cd14_summary(
    fig,
    cell,
    within: pd.DataFrame,
    absolute_samples: pd.DataFrame,
    relative_stats: pd.DataFrame,
    absolute_stats: pd.DataFrame,
) -> None:
    sub = cell.subgridspec(2, 1, hspace=0.42, height_ratios=[1.0, 1.0])
    y = np.arange(len(MARKS))
    offsets = {"COVID": -0.12, "HC": 0.12}

    ax_rel = fig.add_subplot(sub[0, 0])
    for group, color, marker, legend_label in (
        ("COVID", COVID, "o", "COVID-19"),
        ("HC", HC, "s", "HC"),
    ):
        for index, mark in enumerate(MARKS):
            row = within[
                within["peak_set_id"].eq(f"CD14_{mark}")
                & within["mode"].eq("full_interval")
                & within["group"].eq(group)
            ]
            require_rows(row, 1, f"CD14 {mark}/{group} within-group summary")
            row = row.iloc[0]
            summary_errorbar(
                ax_rel,
                index + offsets[group],
                float(row["median_enrichment_ratio"]),
                float(row["ratio_ci_low"]),
                float(row["ratio_ci_high"]),
                color,
                marker,
                legend_label if index == 0 else None,
            )
    ax_rel.axvline(1.0, color=GREY_4, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_rel.set_yticks(y)
    ax_rel.set_yticklabels(MARKS)
    ax_rel.set_ylim(len(MARKS) - 0.5, -0.5)
    ax_rel.set_xlim(0.70, 1.36)
    ax_rel.set_xlabel("Observed / matched expected")
    ax_rel.set_title("Relative enrichment", loc="left", pad=1.8)
    ax_rel.legend(loc="upper left", ncol=1, borderaxespad=0.2)
    rel_tests = relative_stats[
        relative_stats["dataset"].eq("CD14_reference")
        & relative_stats["mode"].eq("full_interval")
    ].set_index("mark")
    for index, mark in enumerate(MARKS):
        stars = q_to_stars(float(rel_tests.loc[mark, "bh_fdr_all_tests"]))
        if stars:
            ax_rel.text(
                1.345,
                index,
                stars,
                fontsize=5.2,
                fontweight="bold",
                ha="right",
                va="center",
            )
    tidy(ax_rel)

    ax_abs = fig.add_subplot(sub[1, 0])
    absolute = absolute_samples[
        absolute_samples["dataset"].eq("CD14_reference")
        & absolute_samples["mode"].eq("full_interval")
    ]
    for group, color, marker in (("COVID", COVID, "o"), ("HC", HC, "s")):
        for index, mark in enumerate(MARKS):
            q1, median, q3 = quantiles(
                absolute[
                    absolute["mark"].eq(mark) & absolute["group"].eq(group)
                ]["absolute_epm_per_total_mapped_alignments"]
            )
            summary_errorbar(
                ax_abs,
                index + offsets[group],
                median,
                q1,
                q3,
                color,
                marker,
            )
    ax_abs.set_xscale("log")
    ax_abs.set_yticks(y)
    ax_abs.set_yticklabels(MARKS)
    ax_abs.set_ylim(len(MARKS) - 0.5, -0.5)
    ax_abs.set_xlim(25, 1_500)
    ax_abs.set_xlabel("Peak-overlapping eccDNAs per 10$^6$ alignments")
    ax_abs.set_title("Absolute overlap abundance", loc="left", pad=1.8)
    abs_tests = absolute_stats[
        absolute_stats["dataset"].eq("CD14_reference")
        & absolute_stats["mode"].eq("full_interval")
    ].set_index("mark")
    for index, mark in enumerate(MARKS):
        stars = q_to_stars(float(abs_tests.loc[mark, "bh_fdr_all_tests"]))
        if stars:
            ax_abs.text(
                1_420,
                index,
                stars,
                fontsize=5.2,
                fontweight="bold",
                ha="right",
                va="center",
            )
    tidy(ax_abs)

    panel_label(fig, ax_rel, "b", dx=-0.055, dy=0.006)


def draw_crosscell_heatmap(fig, cell, multicell: pd.DataFrame) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    ax = fig.add_subplot(cell)
    position = ax.get_position()
    ax.set_position(
        [
            position.x0 + 0.018,
            position.y0,
            position.width - 0.018,
            position.height,
        ]
    )
    frame = multicell[multicell["mode"].eq("full_interval")].copy()
    frame["cell"] = frame["dataset"].str.replace("ENCODE_hg38_", "", regex=False)
    grid = np.full((len(CELL_ORDER), len(MARKS)), np.nan)
    for row_index, (cell_id, _) in enumerate(CELL_ORDER):
        for col_index, mark in enumerate(MARKS):
            row = frame[frame["cell"].eq(cell_id) & frame["mark"].eq(mark)]
            if len(row) == 1:
                grid[row_index, col_index] = float(row.iloc[0]["median_ratio_covid"])
            elif len(row) > 1:
                raise RuntimeError(f"duplicate cross-cell summary for {cell_id}/{mark}")

    # Keep the heat-map scale distinct from the red/blue cohort mapping.
    colours = ["#00394E", "#4BBCBD", GREY_1, "#EBD6E9", "#BB7CB4", "#792C74"]
    bounds = [-0.50, -0.25, -0.08, 0.08, 0.20, 0.35, 0.50]
    cmap = mpl.colors.ListedColormap(colours)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    for i in range(len(CELL_ORDER)):
        for j in range(len(MARKS)):
            ratio = grid[i, j]
            face = "white" if np.isnan(ratio) else cmap(norm(np.log2(ratio)))
            patch = plt.Rectangle(
                (j - 0.5, i - 0.5),
                1,
                1,
                facecolor=face,
                edgecolor="white",
                linewidth=0.45,
                hatch="///" if np.isnan(ratio) else None,
            )
            ax.add_patch(patch)
            label = "n.a." if np.isnan(ratio) else f"{ratio:.2f}"
            text_color = (
                "white"
                if not np.isnan(ratio)
                and (np.log2(ratio) <= -0.25 or np.log2(ratio) >= 0.35)
                else GREY_6
            )
            ax.text(j, i, label, ha="center", va="center", fontsize=5.0, color=text_color)
    ax.set_xlim(-0.5, len(MARKS) - 0.5)
    ax.set_ylim(len(CELL_ORDER) - 0.5, -0.5)
    ax.set_xticks(range(len(MARKS)))
    ax.set_xticklabels(MARKS, rotation=38, ha="right")
    ax.set_yticks(range(len(CELL_ORDER)))
    ax.set_yticklabels([label for _, label in CELL_ORDER])
    ax.tick_params(length=0, pad=1.3)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "GRCh38-native biosamples",
        loc="left",
        pad=2.0,
        fontweight="bold",
        color=GREY_6,
    )
    cax = ax.inset_axes([1.025, 0.02, 0.032, 0.96])
    cb = mpl.colorbar.ColorbarBase(
        cax,
        cmap=cmap,
        norm=norm,
        boundaries=bounds,
        ticks=[-0.4, 0.0, 0.4],
        orientation="vertical",
    )
    cb.set_ticklabels(["−0.4", "0", "0.4"])
    cb.set_label("log$_2$ ratio", fontsize=5.0, labelpad=1.2)
    cb.ax.tick_params(labelsize=5.0, length=1.5, pad=1.0)
    cb.outline.set_linewidth(0.4)
    panel_label(fig, ax, "c", dx=-0.052, dy=0.006)


def draw_mappability(fig, cell, relative_stats: pd.DataFrame, mappable: pd.DataFrame) -> None:
    from matplotlib.lines import Line2D

    ax = fig.add_subplot(cell)
    primary = relative_stats[
        relative_stats["dataset"].eq("CD14_reference")
        & relative_stats["mode"].eq("full_interval")
    ].set_index("mark")
    restricted = mappable[
        mappable["dataset"].eq("CD14_reference")
        & mappable["mode"].eq("full_interval")
    ].set_index("mark")
    if set(primary.index) != set(MARKS) or set(restricted.index) != set(MARKS):
        raise RuntimeError("mappability comparison does not contain all six CD14 marks")

    y = np.arange(len(MARKS))
    group_styles = [
        ("covid", COVID, -0.11),
        ("hc", HC, 0.11),
    ]
    for index, mark in enumerate(MARKS):
        for group, color, offset in group_styles:
            unrestricted = 2 ** float(primary.loc[mark, f"median_{group}"])
            masked = float(restricted.loc[mark, f"median_ratio_{group}"])
            ax.plot(
                [unrestricted, masked],
                [index + offset, index + offset],
                color=color,
                alpha=0.55,
                linewidth=0.65,
                zorder=1,
            )
            ax.plot(
                unrestricted,
                index + offset,
                "o",
                markerfacecolor="white",
                markeredgecolor=color,
                markeredgewidth=0.75,
                markersize=3.7,
                zorder=3,
            )
            ax.plot(
                masked,
                index + offset,
                "D",
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=0.35,
                markersize=3.5,
                zorder=3,
            )
    ax.axvline(1.0, color=GREY_4, linewidth=0.6, linestyle=(0, (3, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels(MARKS)
    ax.set_ylim(len(MARKS) - 0.5, -0.5)
    ax.set_xlim(0.70, 1.36)
    ax.set_xlabel("Median enrichment ratio")
    ax.set_title(
        "High-mappability mask",
        loc="left",
        pad=2.0,
        fontweight="bold",
        color=GREY_6,
    )
    handles = [
        Line2D([0], [0], color=COVID, linewidth=1.0, label="COVID-19"),
        Line2D([0], [0], color=HC, linewidth=1.0, label="HC"),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor=GREY_5,
            markersize=3.8,
            label="Allowed",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            linestyle="none",
            markerfacecolor=GREY_5,
            markeredgecolor="white",
            markersize=3.8,
            label="High-mappability",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.23),
        ncol=2,
        fontsize=5.0,
        columnspacing=0.9,
        borderaxespad=0,
    )
    tidy(ax)
    panel_label(fig, ax, "d", dx=-0.054, dy=0.006)


def draw_adjusted_effects(fig, cell, models: pd.DataFrame) -> None:
    ax = fig.add_subplot(cell)
    frame = models[
        models["dataset"].eq("CD14_reference")
        & models["mode"].eq("full_interval")
        & models["term"].eq("COVID_vs_HC")
        & models["model"].isin(
            ["m1_group_only", "m2_group_rca_age_sex", "m3_plus_eccdna_burden"]
        )
    ].copy()
    model_styles = [
        ("m1_group_only", "Group only", GREY_4, "o", -0.18),
        ("m2_group_rca_age_sex", "+ RCA, age, sex", PURPLE, "s", 0.00),
        ("m3_plus_eccdna_burden", "+ eccDNA burden", GREEN, "D", 0.18),
    ]
    y = np.arange(len(MARKS))
    for model, label, color, marker, offset in model_styles:
        sub = frame[frame["model"].eq(model)].set_index("mark")
        if set(sub.index) != set(MARKS):
            raise RuntimeError(f"{model} does not contain all six CD14 marks")
        for index, mark in enumerate(MARKS):
            row = sub.loc[mark]
            q_value = float(row["bh_fdr_within_family"])
            ax.errorbar(
                float(row["beta"]),
                index + offset,
                xerr=np.array(
                    [
                        [float(row["beta"]) - float(row["ci_low"])],
                        [float(row["ci_high"]) - float(row["beta"])],
                    ]
                ),
                fmt=marker,
                markersize=3.5,
                markerfacecolor=color,
                markeredgecolor=color,
                markeredgewidth=0.75,
                color=color,
                ecolor=color,
                elinewidth=0.65,
                capsize=1.5,
                capthick=0.6,
                label=label if index == 0 else None,
                zorder=3,
            )
            stars = q_to_stars(q_value)
            if stars:
                beta = float(row["beta"])
                if beta >= 0:
                    star_x = float(row["ci_high"]) + 0.006
                    horizontal_alignment = "left"
                else:
                    star_x = float(row["ci_low"]) - 0.006
                    horizontal_alignment = "right"
                ax.text(
                    star_x,
                    index + offset,
                    stars,
                    fontsize=5.0,
                    fontweight="bold",
                    color=color,
                    ha=horizontal_alignment,
                    va="center",
                )
    ax.axvline(0.0, color=GREY_5, linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(MARKS)
    ax.set_ylim(len(MARKS) - 0.55, -0.55)
    ax.set_xlim(-0.20, 0.18)
    ax.set_xlabel("COVID-19 − HC coefficient\nin log$_2$ enrichment (95% CI)")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.58, 1.015),
        ncol=3,
        fontsize=5.0,
        columnspacing=0.9,
        borderaxespad=0,
    )
    tidy(ax)
    panel_label(fig, ax, "f", dx=-0.050, dy=0.006)


def draw_a549_summary(
    fig,
    cell,
    within: pd.DataFrame,
    absolute_samples: pd.DataFrame,
    absolute_stats: pd.DataFrame,
) -> None:
    outer = cell.subgridspec(1, 2, wspace=0.58)
    y = np.arange(2)
    marks = ["H3K27ac", "H3K4me3"]
    offsets = {"COVID": -0.12, "HC": 0.12}

    ax_rel = fig.add_subplot(outer[0, 0])
    for group, color, marker, legend_label in (
        ("COVID", COVID, "o", "COVID-19"),
        ("HC", HC, "s", "HC"),
    ):
        for index, mark in enumerate(marks):
            row = within[
                within["peak_set_id"].eq(A549_PEAKS[mark])
                & within["mode"].eq("full_interval")
                & within["group"].eq(group)
            ]
            require_rows(row, 1, f"A549 {mark}/{group} within-group summary")
            row = row.iloc[0]
            summary_errorbar(
                ax_rel,
                index + offsets[group],
                float(row["median_enrichment_ratio"]),
                float(row["ratio_ci_low"]),
                float(row["ratio_ci_high"]),
                color,
                marker,
                legend_label if index == 0 else None,
            )
    ax_rel.axvline(1.0, color=GREY_4, linewidth=0.6, linestyle=(0, (3, 2)))
    ax_rel.set_yticks(y)
    ax_rel.set_yticklabels(marks)
    ax_rel.set_ylim(1.55, -0.55)
    ax_rel.set_xlim(0.82, 1.25)
    ax_rel.set_xlabel("Observed / expected")
    ax_rel.set_title("A549-ACE2", loc="left", pad=1.8, fontweight="bold", color=GREY_6)
    ax_rel.legend(
        loc="upper center",
        bbox_to_anchor=(1.37, -0.23),
        ncol=2,
        fontsize=5.0,
        columnspacing=0.9,
        borderaxespad=0,
    )
    tidy(ax_rel)

    ax_abs = fig.add_subplot(outer[0, 1])
    absolute = absolute_samples[
        absolute_samples["dataset"].eq("A549_GSE179184_official")
        & absolute_samples["mode"].eq("full_interval")
    ]
    for group, color, marker in (("COVID", COVID, "o"), ("HC", HC, "s")):
        for index, mark in enumerate(marks):
            q1, median, q3 = quantiles(
                absolute[
                    absolute["mark"].eq(mark) & absolute["group"].eq(group)
                ]["absolute_epm_per_total_mapped_alignments"]
            )
            summary_errorbar(
                ax_abs,
                index + offsets[group],
                median,
                q1,
                q3,
                color,
                marker,
            )
    ax_abs.set_xscale("log")
    ax_abs.set_yticks(y)
    ax_abs.set_yticklabels([])
    ax_abs.set_ylim(1.55, -0.55)
    ax_abs.set_xlim(5, 220)
    ax_abs.set_xlabel("eccDNAs per 10$^6$\nalignments")
    ax_abs.set_title("")
    abs_frame = absolute_stats[
        absolute_stats["dataset"].eq("A549_GSE179184_official")
        & absolute_stats["mode"].eq("full_interval")
    ].set_index("mark")
    for index, mark in enumerate(marks):
        stars = q_to_stars(float(abs_frame.loc[mark, "bh_fdr_all_tests"]))
        if stars:
            ax_abs.text(
                210,
                index,
                stars,
                fontsize=5.2,
                fontweight="bold",
                ha="right",
                va="center",
            )
    tidy(ax_abs)

    panel_label(fig, ax_rel, "e", dx=-0.055, dy=0.006)


def draw_burden_relationship(
    fig,
    cell,
    sample_enrichment: pd.DataFrame,
    associations: pd.DataFrame,
    matched_stats: pd.DataFrame,
) -> None:
    sub = cell.subgridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.10)
    marks = ["H3K27me3", "H3K9me3"]
    axes = []
    for index, mark in enumerate(marks):
        ax = fig.add_subplot(sub[index, 0])
        axes.append(ax)
        frame = sample_enrichment[
            sample_enrichment["dataset"].eq("CD14_reference")
            & sample_enrichment["mark"].eq(mark)
            & sample_enrichment["mode"].eq("full_interval")
        ].copy()
        require_rows(frame, 78, f"{mark} burden-association samples")
        for group, color, label in (
            ("COVID", COVID, "COVID-19"),
            ("HC", HC, "HC"),
        ):
            group_frame = frame[frame["group"].eq(group)]
            ax.scatter(
                np.log10(group_frame["eligible_eccdna_count"].to_numpy(float)),
                group_frame["log2_enrichment_ratio"].to_numpy(float),
                s=6.5,
                facecolor=color,
                edgecolor="white",
                linewidth=0.25,
                alpha=0.88,
                label=label if index == 0 else None,
                zorder=3,
            )

        match = matched_stats[
            matched_stats["dataset"].eq("CD14_reference")
            & matched_stats["mark"].eq(mark)
            & matched_stats["mode"].eq("full_interval")
        ]
        require_rows(match, 1, f"{mark} burden-matched comparison")
        match = match.iloc[0]
        window_low = np.log10(float(match["matched_window_low"]))
        window_high = np.log10(float(match["matched_window_high"]))
        ax.axvspan(window_low, window_high, color=GREY_1, alpha=0.85, zorder=0)
        ax.axhline(0.0, color=GREY_5, linewidth=0.6, zorder=1)

        mark_stats = associations[
            associations["dataset"].eq("CD14_reference")
            & associations["mark"].eq(mark)
            & associations["mode"].eq("full_interval")
        ].set_index("stratum")
        if not {"COVID", "HC"}.issubset(mark_stats.index):
            raise RuntimeError(f"{mark}: missing within-cohort burden association")
        covid_row = mark_stats.loc["COVID"]
        hc_row = mark_stats.loc["HC"]
        annotation_y = (0.20, 0.07) if index == 0 else (0.88, 0.74)
        for y_position, row, color, label in (
            (annotation_y[0], covid_row, COVID, "COVID-19"),
            (annotation_y[1], hc_row, HC, "HC"),
        ):
            ax.text(
                0.985,
                y_position,
                (
                    f"{label}: ρ={float(row['spearman_rho']):+.2f}, "
                    f"P={float(row['spearman_p']):.3f}, "
                    f"q={float(row['bh_fdr_within_family']):.3f}"
                ),
                transform=ax.transAxes,
                fontsize=4.8,
                color=color,
                ha="right",
                va="center",
            )
        ax.text(
            0.015,
            0.96,
            mark,
            transform=ax.transAxes,
            fontsize=5.7,
            color=GREY_5,
            ha="left",
            va="top",
        )
        ax.set_xlim(3.5, 6.25)
        ax.set_ylim(-0.52, 0.47)
        ax.set_yticks([-0.4, 0.0, 0.4])
        tidy(ax)
        if index == 0:
            ax.tick_params(labelbottom=False)
            ax.set_ylabel("log$_2$ enrichment")
        else:
            ax.set_xlabel("log$_{10}$ detected eccDNAs per sample")

    panel_label(fig, axes[0], "g", dx=-0.048, dy=0.006)


def candidate_point_area(values: pd.Series | np.ndarray) -> np.ndarray:
    return 5.0 + 0.45 * np.asarray(values, dtype=float)


def draw_candidate_landscape(fig, cell, candidates: pd.DataFrame) -> None:
    ax = fig.add_subplot(cell)
    frame = candidates.copy()
    frame["minus_log10_fdr"] = -np.log10(
        pd.to_numeric(frame["FDR"], errors="raise")
    )
    markers = {
        "Humoral immune response": "o",
        "Endothelial/coagulation": "^",
    }
    states = [
        ("Not strict-eligible", ~frame["Strict_selection_eligible"].eq("Yes")),
        (
            "Strict-eligible",
            frame["Strict_selection_eligible"].eq("Yes")
            & ~frame["Dual_reference_H3K27ac_eligible"].eq("Yes"),
        ),
        (
            "Dual-reference H3K27ac",
            frame["Dual_reference_H3K27ac_eligible"].eq("Yes")
            & ~frame["Selected_for_Figure_5_display"].eq("Yes"),
        ),
        ("Selected locus", frame["Selected_for_Figure_5_display"].eq("Yes")),
    ]
    style = {
        "Not strict-eligible": dict(
            facecolors="white", edgecolors=GREY_2, linewidths=0.45, alpha=0.78
        ),
        "Strict-eligible": dict(
            facecolors=GREY_3, edgecolors="white", linewidths=0.35, alpha=0.72
        ),
        "Dual-reference H3K27ac": dict(
            facecolors=TEAL_LIGHT, edgecolors=TEAL, linewidths=0.65, alpha=0.95
        ),
        "Selected locus": dict(
            facecolors=PURPLE, edgecolors=PURPLE_DARK, linewidths=0.75, alpha=1.0
        ),
    }
    for state, state_mask in states:
        for category, marker in markers.items():
            subset = frame[state_mask & frame["Category_for_row"].eq(category)]
            if subset.empty:
                continue
            ax.scatter(
                subset["Effect_size"],
                subset["minus_log10_fdr"],
                s=candidate_point_area(subset["COVID_detected_n"]),
                marker=marker,
                zorder=5 if state == "Selected locus" else 2,
                **style[state],
            )
    label_offsets = {"HLA-E": (-30, 13), "MCAM": (8, 9)}
    for gene in SELECTED_GENES:
        row = frame[
            frame["Gene"].eq(gene)
            & frame["Selected_for_Figure_5_display"].eq("Yes")
        ]
        require_rows(row, 1, f"selected candidate {gene}")
        row = row.iloc[0]
        ax.annotate(
            gene,
            xy=(float(row["Effect_size"]), float(row["minus_log10_fdr"])),
            xytext=label_offsets[gene],
            textcoords="offset points",
            fontsize=5.2,
            fontstyle="italic",
            fontweight="bold",
            color=PURPLE_DARK,
            ha="left",
            va="center",
            arrowprops=dict(
                arrowstyle="-",
                color=PURPLE_DARK,
                linewidth=0.55,
                shrinkA=2,
                shrinkB=2,
            ),
            bbox=dict(
                boxstyle="round,pad=0.12",
                facecolor="white",
                edgecolor="none",
                alpha=0.92,
            ),
            zorder=8,
        )
    ax.set_xlim(0.10, 0.66)
    ax.set_ylim(1.25, 4.48)
    ax.set_xlabel("Cliff's δ (COVID-19 vs HC)")
    ax.set_ylabel("−log$_{10}$(BH FDR)")
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    ax.set_yticks([1.3, 2.0, 3.0, 4.0])
    ax.grid(axis="y", color=GREY_1, linewidth=0.45, zorder=0)
    ax.set_title(
        "Candidate-locus landscape",
        loc="left",
        pad=2.0,
        fontweight="bold",
        color=GREY_6,
    )
    tidy(ax)

    # Compact, non-colour-reliant key.
    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=GREY_3,
            markeredgecolor="white",
            markersize=4.0,
            label="Humoral immune",
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            linestyle="none",
            markerfacecolor=GREY_3,
            markeredgecolor="white",
            markersize=4.2,
            label="Endothelial/coagulation",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor=GREY_2,
            markersize=4.0,
            label="Not strict-eligible",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=GREY_3,
            markeredgecolor="white",
            markersize=4.0,
            label="Strict-eligible",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=TEAL_LIGHT,
            markeredgecolor=TEAL,
            markersize=4.0,
            label="Dual-reference H3K27ac",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=PURPLE,
            markeredgecolor=PURPLE_DARK,
            markersize=4.0,
            label="Selected locus",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.19),
        ncol=3,
        fontsize=4.4,
        columnspacing=0.85,
        handletextpad=0.35,
        borderaxespad=0,
    )
    ax.text(
        0.99,
        0.02,
        "n=281; strict=149; dual-reference=54",
        transform=ax.transAxes,
        fontsize=4.5,
        color=GREY_5,
        ha="right",
        va="bottom",
    )
    panel_label(fig, ax, "g", dx=-0.048, dy=0.006)


def format_q(value: float) -> str:
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    return f"{mantissa:.2g}×10$^{{{exponent}}}$"


def draw_locus_axis(
    ax,
    samples: pd.DataFrame,
    candidates: pd.DataFrame,
    gene: str,
    show_ylabel: bool,
) -> None:
    row = candidates[
        candidates["Gene"].eq(gene)
        & candidates["Selected_for_Figure_5_display"].eq("Yes")
    ]
    require_rows(row, 1, f"selected candidate {gene}")
    row = row.iloc[0]
    groups = ("COVID-19", "HC")
    colours = (COVID, HC)
    transformed: list[np.ndarray] = []
    for group in groups:
        values = samples[
            samples["Gene"].eq(gene) & samples["group"].eq(group)
        ]["primary_junction_EA"].to_numpy(float)
        if len(values) != 39:
            raise RuntimeError(f"{gene}/{group}: expected 39 samples, found {len(values)}")
        transformed.append(np.log10(values + 1.0))

    box = ax.boxplot(
        transformed,
        positions=[0, 1],
        widths=0.42,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 0.75},
        boxprops={"color": "black", "linewidth": 0.50},
        whiskerprops={"color": "black", "linewidth": 0.50},
        capprops={"color": "black", "linewidth": 0.50},
    )
    for patch, color in zip(box["boxes"], colours):
        patch.set_facecolor(color)
        patch.set_alpha(0.24)

    for index, (values, color) in enumerate(zip(transformed, colours)):
        rng = np.random.default_rng(6200 + index + (0 if gene == "HLA-E" else 20))
        ax.scatter(
            np.full(len(values), index) + rng.uniform(-0.16, 0.16, len(values)),
            values,
            s=5.0,
            facecolors=color,
            edgecolors="white",
            linewidths=0.20,
            alpha=0.88,
            zorder=3,
        )

    tick_values = np.array([0, 10, 100, 1_000], dtype=float)
    ax.set_yticks(np.log10(tick_values + 1.0))
    ax.set_yticklabels(["0", "10", "100", "1,000"])
    ax.set_ylim(-0.08, 3.27)
    ax.set_xlim(-0.52, 1.52)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(groups)
    ax.grid(axis="y", color=GREY_1, linewidth=0.45, zorder=0)
    if show_ylabel:
        ax.set_ylabel("Junction-level EA\n(log$_{10}$[EA + 1])")
    else:
        ax.tick_params(labelleft=False)
        ax.spines["left"].set_visible(False)
    tidy(ax, left=show_ylabel)
    ax.set_title(
        gene,
        fontsize=6.7,
        fontstyle="italic",
        fontweight="bold",
        color=GREY_6,
        pad=24,
    )
    ax.text(
        0.5,
        1.16,
        f"BH q={format_q(float(row['FDR']))}; δ={float(row['Effect_size']):.2f}",
        transform=ax.transAxes,
        fontsize=4.8,
        color="black",
        ha="center",
        va="bottom",
    )
    ax.text(
        0.5,
        1.09,
        (
            f"Detected: {int(row['COVID_detected_n'])}/39 vs "
            f"{int(row['HC_detected_n'])}/39   |   "
            f"masked: {int(row['quality_COVID_detected_n'])}/39 vs "
            f"{int(row['quality_HC_detected_n'])}/39"
        ),
        transform=ax.transAxes,
        fontsize=4.7,
        color=GREY_5,
        ha="center",
        va="bottom",
    )
    ax.text(
        0.5,
        1.025,
        "dual-reference H3K27ac overlap",
        transform=ax.transAxes,
        fontsize=4.6,
        color=TEAL,
        ha="center",
        va="bottom",
    )


def draw_selected_loci(
    fig, cell, samples: pd.DataFrame, candidates: pd.DataFrame
) -> None:
    sub = cell.subgridspec(1, 2, wspace=0.22)
    axes = []
    for index, gene in enumerate(SELECTED_GENES):
        ax = fig.add_subplot(sub[0, index])
        draw_locus_axis(ax, samples, candidates, gene, show_ylabel=index == 0)
        axes.append(ax)
    panel_label(fig, axes[0], "h", dx=-0.040, dy=0.050)
    box = axes[0].get_position()
    fig.text(
        box.x0,
        box.y1 + 0.054,
        "Selected loci: sample-level primary abundance",
        fontsize=6.2,
        fontweight="bold",
        color=GREY_6,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.52,
        0.018,
        (
            "Points are individual samples; intervals are IQRs or 95% CIs as labelled. "
            "Masked denotes the prespecified breakpoint-quality sensitivity analysis."
        ),
        fontsize=4.7,
        color=GREY_5,
        ha="center",
        va="bottom",
    )


def validate_inputs(
    density: pd.DataFrame,
    within: pd.DataFrame,
    relative_stats: pd.DataFrame,
    absolute_stats: pd.DataFrame,
    absolute_samples: pd.DataFrame,
    multicell: pd.DataFrame,
    mappable: pd.DataFrame,
    models: pd.DataFrame,
    sample_enrichment: pd.DataFrame,
    associations: pd.DataFrame,
    matched_stats: pd.DataFrame,
) -> None:
    for mark in MARKS:
        if not density["peak_set_id"].eq(f"CD14_{mark}").any():
            raise RuntimeError(f"peak-centred density lacks CD14_{mark}")
    if not {
        "CD14_reference",
        "A549_GSE179184_official",
    }.issubset(set(within["dataset"])):
        raise RuntimeError("within-group table lacks a required reference dataset")
    for frame, label in (
        (relative_stats, "relative group statistics"),
        (absolute_stats, "absolute group statistics"),
        (absolute_samples, "absolute sample table"),
        (multicell, "cross-cell statistics"),
        (mappable, "mappability statistics"),
        (models, "adjusted model table"),
        (sample_enrichment, "sample enrichment table"),
        (associations, "burden-association statistics"),
        (matched_stats, "burden-matched statistics"),
    ):
        if frame.empty:
            raise RuntimeError(f"{label} is empty")
    for mark in ("H3K27me3", "H3K9me3"):
        require_rows(
            sample_enrichment[
                sample_enrichment["dataset"].eq("CD14_reference")
                & sample_enrichment["mark"].eq(mark)
                & sample_enrichment["mode"].eq("full_interval")
            ],
            78,
            f"{mark} full-interval sample enrichment",
        )


def main(argv: Sequence[str] | None = None) -> int:
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tables",
        type=Path,
        default=repo / "analyse" / "chromatin" / "tables",
    )
    parser.add_argument("--output-dir", type=Path, default=repo / "Submit")
    parser.add_argument("--stem", default="Figure4new")
    args = parser.parse_args(argv)

    tables = args.tables
    density = read_table(
        tables / "ext_peak_centred_density.tsv",
        [
            "sample_id",
            "group",
            "peak_set_id",
            "offset_bp",
            "density_per_bp_per_centre",
            "total_breakpoint_units",
        ],
    )
    within = read_table(
        tables / "p0_within_group_enrichment.tsv",
        [
            "peak_set_id",
            "dataset",
            "mark",
            "mode",
            "group",
            "median_enrichment_ratio",
            "ratio_ci_low",
            "ratio_ci_high",
        ],
    )
    relative_stats = read_table(
        tables / "relative_enrichment_group_stats.tsv",
        [
            "dataset",
            "mark",
            "mode",
            "median_covid",
            "median_hc",
            "bh_fdr_all_tests",
        ],
    )
    absolute_stats = read_table(
        tables / "absolute_burden_group_stats.tsv",
        ["dataset", "mark", "mode", "bh_fdr_all_tests"],
    )
    absolute_samples = read_table(
        tables / "sample_absolute_burden.tsv",
        [
            "sample_id",
            "group",
            "dataset",
            "mark",
            "mode",
            "absolute_epm_per_total_mapped_alignments",
        ],
    )
    multicell = read_table(
        tables / "ext_multicell_group_stats.tsv",
        ["dataset", "mark", "mode", "median_ratio_covid"],
    )
    mappable = read_table(
        tables / "ext_mappability_group_stats.tsv",
        [
            "dataset",
            "mark",
            "mode",
            "median_ratio_covid",
            "median_ratio_hc",
        ],
    )
    models = read_table(
        tables / "p0_covariate_adjusted_hc3.tsv",
        [
            "dataset",
            "mark",
            "mode",
            "model",
            "term",
            "beta",
            "ci_low",
            "ci_high",
            "bh_fdr_within_family",
        ],
    )
    sample_enrichment = read_table(
        tables / "sample_relative_enrichment.tsv",
        [
            "sample_id",
            "group",
            "dataset",
            "mark",
            "mode",
            "eligible_eccdna_count",
            "log2_enrichment_ratio",
        ],
    )
    associations = read_table(
        tables / "p0_burden_association.tsv",
        [
            "dataset",
            "mark",
            "mode",
            "stratum",
            "spearman_rho",
            "spearman_p",
            "bh_fdr_within_family",
        ],
    )
    matched_stats = read_table(
        tables / "p0_burden_matched_subset.tsv",
        [
            "dataset",
            "mark",
            "mode",
            "matched_window_low",
            "matched_window_high",
            "n_covid",
            "n_hc",
            "cliffs_delta_covid_vs_hc",
            "bh_fdr_within_family",
        ],
    )
    validate_inputs(
        density,
        within,
        relative_stats,
        absolute_stats,
        absolute_samples,
        multicell,
        mappable,
        models,
        sample_enrichment,
        associations,
        matched_stats,
    )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    setup_style()
    fig = plt.figure(figsize=(183 / 25.4, 190 / 25.4))
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.38, 0.92, 0.94],
        hspace=0.52,
        left=0.078,
        right=0.985,
        top=0.975,
        bottom=0.075,
    )

    top = outer[0, 0].subgridspec(1, 2, width_ratios=[1.78, 1.0], wspace=0.30)
    draw_profile_panel(fig, top[0, 0], density)
    draw_cd14_summary(
        fig,
        top[0, 1],
        within,
        absolute_samples,
        relative_stats,
        absolute_stats,
    )

    controls = outer[1, 0].subgridspec(
        1, 3, width_ratios=[1.32, 0.91, 1.18], wspace=0.57
    )
    draw_crosscell_heatmap(fig, controls[0, 0], multicell)
    draw_mappability(fig, controls[0, 1], relative_stats, mappable)
    draw_a549_summary(
        fig,
        controls[0, 2],
        within,
        absolute_samples,
        absolute_stats,
    )

    bottom = outer[2, 0].subgridspec(
        1, 2, width_ratios=[1.02, 1.38], wspace=0.38
    )
    draw_adjusted_effects(fig, bottom[0, 0], models)
    draw_burden_relationship(
        fig,
        bottom[0, 1],
        sample_enrichment,
        associations,
        matched_stats,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "svg", "png"):
        path = args.output_dir / f"{args.stem}.{extension}"
        fig.savefig(
            path,
            dpi=600 if extension == "png" else None,
            metadata={
                "Title": "Plasma eccDNA overlap with reference chromatin annotations",
                "Creator": "make_figure4new.py",
            },
        )
    plt.close(fig)
    print(
        f"Wrote {args.output_dir / (args.stem + '.pdf')}, "
        f"{args.stem}.svg and {args.stem}.png"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
