#!/usr/bin/env python3
"""Merged chromatin main figure: replaces both Figure 4 and Figure 5.

Figure contract
---------------
Core conclusion
    Plasma eccDNA is non-randomly distributed with respect to histone-marked
    chromatin - enriched at active marks and depleted at H3K9me3, by the same
    amount in both cohorts - and the uniform Cliff's delta of about +0.75 in
    COVID-19 across every peak set is an eccDNA burden effect that disappears
    once burden is matched.

Archetype
    Asymmetric composite: a small-multiple block on top, then a three-column
    forest layout in which every column shares the same eight peak-set rows.
    Panel d is the hero.

Panel map and evidence hierarchy
    a  establish the system   peak-centred breakpoint density, 8 peak sets,
                              one shared y axis so the mark contrast is visible
    b  what the values are    absolute abundance, per-sample strip + median/IQR
    c  the positive claim     enrichment over the matched placement expectation,
                              per-sample strip + median with the locked 95% CI
    d  HERO, the contrast     Cliff's delta for COVID-19 vs HC across all eight
                              prespecified analysis variants, per peak set

Design rationale
    The first merged draft used sixteen box plots per metric and a 64-cell
    numeric heat map; it read as a spreadsheet.  This version keeps the marker
    idiom that the rest of the figure already uses: dots, medians, reference
    rules and shaded bands.  It also separates description from inference.
    Panels b and c say what the values are and carry no significance stars at
    all; every inferential statement lives in panel d, where a filled marker
    means q < 0.05.  In c the locked 95% interval does that work directly - an
    interval clear of zero is the significant within-cohort enrichment - so no
    asterisk column is needed.

Why the two figures were merged
    Old Figure 4 (CD14+ reference chromatin) and old Figure 5 (A549-ACE2
    infected-cell chromatin) ran the same three analyses on two peak-set
    families and reached the same conclusion.  Merged, the uniformity of the
    absolute effect across all eight peak sets becomes the argument rather than
    a coincidence: a genuine chromatin-specific redistribution could not give
    active marks and heterochromatin the same effect size.

Statistical integrity
    Every plotted number is read from the locked hg38 analysis tables.  No test,
    transformation, group definition, threshold or significance tier is
    recomputed anywhere in this figure.

Backend
    Python only (matplotlib).  183 mm double column, Arial, editable vector text.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Locked analysis constants

ALLOWED_BP = 2831289217

# Row order: CD14+ reference marks, active first, then repressive; then the two
# infected-cell sets.  Panels b, c and d all reuse this order.
PEAK_SETS: List[Tuple[str, str, str]] = [
    ("CD14_H3K27ac", "H3K27ac", "CD14"),
    ("CD14_H3K4me1", "H3K4me1", "CD14"),
    ("CD14_H3K4me3", "H3K4me3", "CD14"),
    ("CD14_H3K9ac", "H3K9ac", "CD14"),
    ("CD14_H3K27me3", "H3K27me3", "CD14"),
    ("CD14_H3K9me3", "H3K9me3", "CD14"),
    ("A549_GSE179184_H3K27ac_official", "H3K27ac", "A549"),
    ("A549_GSE179184_H3K4me3_official", "H3K4me3", "A549"),
]
SET_IDS = [pid for pid, _, _ in PEAK_SETS]
N_CD14 = sum(1 for _, _, family in PEAK_SETS if family == "CD14")

# ---------------------------------------------------------------------------
# Palette: one signal family.  Red always means COVID-19 and blue always means
# HC, in panel d too, where the sign of Cliff's delta is therefore encoded by
# the cohort it favours.
COVID_LINE, HC_LINE = "#C84F50", "#36617B"
COVID_POINT, HC_POINT = "#D70000", "#034E61"
COHORT_LINE = {"COVID": COVID_LINE, "HC": HC_LINE}
COHORT_POINT = {"COVID": COVID_POINT, "HC": HC_POINT}
REFERENCE_GREY = "#B2B2B2"
ZERO_LINE = "#49566D"
NOTE_GREY = "#5A6673"
LABEL_SLATE = "#3C4C5A"
BAND_FILL = "#F2F4F6"
DARK_SLATE = "#22303C"

ROW_OFFSET = 0.19          # cohort bands within one peak-set row
JITTER = 0.105

# Typography: three sizes.
LAB, TICK, PAN = 6.0, 5.5, 8.0

# Panel d variants.  The first is the absolute metric; the rest are relative
# enrichment under the prespecified sensitivity analyses.
VARIANTS: List[Tuple[str, str, str, str, Optional[Tuple[str, int]]]] = [
    ("Absolute abundance", "absolute_burden_group_stats.tsv", "full_interval", "bh_fdr_all_tests", None),
    ("Full interval", "relative_enrichment_group_stats.tsv", "full_interval", "bh_fdr_all_tests", None),
    ("Junction breakpoints", "relative_enrichment_group_stats.tsv", "junction_start_end", "bh_fdr_all_tests", None),
    ("Midpoint", "relative_enrichment_group_stats.tsv", "midpoint", "bh_fdr_all_tests", None),
    ("Downsampled to 4,000", "p0_downsampled_group_stats.tsv", "full_interval", "bh_fdr_within_family",
     ("downsample_target", 4000)),
    ("Downsampled to 20,000", "p0_downsampled_group_stats.tsv", "full_interval", "bh_fdr_within_family",
     ("downsample_target", 20000)),
    ("Burden-matched subset", "p0_burden_matched_subset.tsv", "full_interval", "bh_fdr_within_family", None),
    ("Uniquely mappable only", "ext_mappability_group_stats.tsv", "full_interval",
     "bh_fdr_group_within_family", None),
]
ABSOLUTE_ROW = 0
BURDEN_MATCHED_ROW = 6


def setup_matplotlib(mpl) -> None:
    """Locked project rcParams, unchanged from the published Figure 4."""
    from matplotlib import font_manager

    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next(
        (name for name in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans") if name in available),
        "DejaVu Sans",
    )
    mpl.rcParams.update(
        {
            "font.family": family,
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "font.size": LAB,
            "axes.labelsize": LAB,
            "xtick.labelsize": TICK,
            "ytick.labelsize": TICK,
            "legend.fontsize": TICK,
            "legend.frameon": False,
            "mathtext.fontset": "custom",
            "mathtext.rm": family,
            "mathtext.it": f"{family}:italic",
            "mathtext.bf": f"{family}:bold",
            "mathtext.cal": f"{family}:italic",
            "mathtext.tt": family,
            "mathtext.sf": family,
            "mathtext.default": "regular",
            "legend.handlelength": 1.2,
            "legend.handletextpad": 0.4,
            "legend.labelspacing": 0.25,
            "legend.borderpad": 0.15,
        }
    )


def tidy(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=TICK, length=2.0, pad=1.5)
    ax.xaxis.label.set_size(LAB)
    ax.yaxis.label.set_size(LAB)


LABEL_X = 0.010


def panel_label(fig, ax, letter: str, dy: float = 0.010, x: Optional[float] = None) -> None:
    box = ax.get_position()
    fig.text(LABEL_X if x is None else x, min(0.999, box.y1 + dy), letter,
             fontsize=PAN, fontweight="bold", va="bottom", ha="left")


def row_bands(ax) -> None:
    """Alternating background bands, so eight rows stay readable across columns."""
    for index in range(len(PEAK_SETS)):
        if index % 2 == 0:
            ax.axhspan(index - 0.5, index + 0.5, color=BAND_FILL, linewidth=0, zorder=0)
    ax.axhline(N_CD14 - 0.5, color=LABEL_SLATE, linewidth=0.7, zorder=2)
    ax.set_ylim(len(PEAK_SETS) - 0.5, -0.5)


def strip_row(ax, index: int, values: np.ndarray, group: str, seed: int,
              interval: Optional[Tuple[float, float]] = None) -> None:
    """One cohort inside one peak-set row: dots, median tick and a range bar."""
    centre = index + (-ROW_OFFSET if group == "COVID" else ROW_OFFSET)
    rng = np.random.default_rng(seed)
    ax.scatter(values, centre + rng.uniform(-JITTER, JITTER, size=len(values)),
               s=1.6, color=COHORT_POINT[group], edgecolors="none", alpha=0.7, zorder=3)
    low, high = interval if interval is not None else (np.percentile(values, 25), np.percentile(values, 75))
    ax.plot([low, high], [centre, centre], color=COHORT_LINE[group], linewidth=0.9,
            solid_capstyle="butt", zorder=4)
    median = float(np.median(values))
    ax.plot([median, median], [centre - 0.155, centre + 0.155], color=COHORT_LINE[group],
            linewidth=1.5, solid_capstyle="butt", zorder=5)


# ---------------------------------------------------------------------------
# Panel a


def density_fold(density: pd.DataFrame, peak_set_id: str) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Cohort mean fold-enrichment of eccDNA breakpoints over the genome average.

    ``density_per_bp_per_centre`` is breakpoints per bp per peak centre.  Dividing
    by each sample's genome-wide breakpoint density (total units / allowed bp)
    yields a fold value whose null is exactly 1 and which is comparable between
    cohorts despite the sevenfold difference in eccDNA burden.
    """
    sub = density[density["peak_set_id"] == peak_set_id].copy()
    sub["fold"] = sub["density_per_bp_per_centre"] * ALLOWED_BP / sub["total_breakpoint_units"]
    out: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for group in ("COVID", "HC"):
        pivot = sub[sub["group"] == group].pivot_table(index="offset_bp", columns="sample_id", values="fold")
        values = pivot.to_numpy(float)
        out[group] = (
            pivot.index.to_numpy(float),
            np.nanmean(values, axis=1),
            np.nanstd(values, axis=1, ddof=1) / np.sqrt(values.shape[1]),
        )
    return out


# ---------------------------------------------------------------------------
# Panel d


def robustness_matrix(tables: Path) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Locked Cliff's delta and BH-FDR for every variant x peak set."""
    deltas = np.full((len(VARIANTS), len(PEAK_SETS)), np.nan)
    qvalues = np.full((len(VARIANTS), len(PEAK_SETS)), np.nan)
    labels = []
    for row, (label, filename, mode, q_column, extra) in enumerate(VARIANTS):
        labels.append(label)
        frame = pd.read_csv(tables / filename, sep="\t")
        frame = frame[frame["peak_set_id"].isin(SET_IDS) & (frame["mode"] == mode)]
        if extra is not None:
            frame = frame[frame[extra[0]] == extra[1]]
        frame = frame.set_index("peak_set_id")
        for column, pid in enumerate(SET_IDS):
            if pid not in frame.index:
                raise SystemExit(f"panel d: {label} is missing {pid}")
            deltas[row, column] = float(frame.loc[pid, "cliffs_delta_covid_vs_hc"])
            qvalues[row, column] = float(frame.loc[pid, q_column])
    return deltas, qvalues, labels


# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--outdir", default=None, help="default: <analysis-dir>/figures")
    parser.add_argument("--stem", default="figure_chromatin_merged")
    args = parser.parse_args(argv)

    analysis = Path(args.analysis_dir)
    tables = analysis / "tables"
    outdir = Path(args.outdir) if args.outdir else analysis / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    density = pd.read_csv(tables / "ext_peak_centred_density.tsv", sep="\t")
    per_sample = pd.read_csv(tables / "sample_absolute_burden.tsv", sep="\t")
    per_sample = per_sample[
        per_sample["peak_set_id"].isin(SET_IDS) & (per_sample["mode"] == "full_interval")
    ].copy()
    within = pd.read_csv(tables / "p0_within_group_enrichment.tsv", sep="\t")
    within = within[within["mode"] == "full_interval"].set_index(["peak_set_id", "group"])

    deltas, qvalues, variant_labels = robustness_matrix(tables)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.lines import Line2D

    setup_matplotlib(mpl)

    fig = plt.figure(figsize=(183 / 25.4, 190 / 25.4))
    LEFT, RIGHT = 0.118, 0.988
    a_grid = GridSpec(2, 4, figure=fig, left=LEFT, right=RIGHT, top=0.962, bottom=0.716,
                      wspace=0.16, hspace=0.40)
    bcd = GridSpec(1, 3, figure=fig, left=LEFT, right=RIGHT, top=0.606, bottom=0.178,
                   width_ratios=[0.235, 0.285, 0.480], wspace=0.10)

    # ------------------------------------------------------------------ a
    a_axes = []
    for index, (peak_set_id, label, family) in enumerate(PEAK_SETS):
        ax = fig.add_subplot(a_grid[index // 4, index % 4])
        a_axes.append(ax)
        folds = density_fold(density, peak_set_id)
        for group in ("COVID", "HC"):
            offsets, mean, sem = folds[group]
            ax.plot(offsets / 1000.0, mean, "-", color=COHORT_LINE[group], linewidth=0.8,
                    label="COVID-19" if group == "COVID" else "HC", solid_capstyle="round")
            ax.fill_between(offsets / 1000.0, mean - sem, mean + sem,
                            color=COHORT_LINE[group], alpha=0.18, linewidth=0)
        ax.axhline(1.0, color=REFERENCE_GREY, linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
        ax.axvline(0.0, color=ZERO_LINE, linewidth=0.4, zorder=0)
        prefix = "CD14$^+$ " if family == "CD14" else "A549-ACE2 "
        ax.set_title(prefix + label, fontsize=LAB, pad=2.0)
        ax.set_xticks([-2, 0, 2])
        tidy(ax)
    # One shared y range across all eight peak sets, so the depletion at
    # H3K9me3 and the flat A549 H3K4me3 profile are visual facts rather than
    # artefacts of eight independently scaled axes.
    low = min(ax.get_ylim()[0] for ax in a_axes)
    high = max(ax.get_ylim()[1] for ax in a_axes)
    for index, ax in enumerate(a_axes):
        ax.set_ylim(low, high)
        ax.set_yticks([0.8, 1.0, 1.2, 1.4])
        if index % 4 != 0:
            ax.set_yticklabels([])
        if index // 4 == 0:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(["-2", "0", "2"])
    fig.text((a_axes[4].get_position().x0 + a_axes[7].get_position().x1) / 2.0,
             a_axes[4].get_position().y0 - 0.030,
             "Distance from peak centre (kb)", ha="center", va="top", fontsize=LAB)
    fig.text(a_axes[0].get_position().x0 - 0.048,
             (a_axes[0].get_position().y1 + a_axes[4].get_position().y0) / 2.0,
             "eccDNA breakpoint density\n(fold over genome average)",
             rotation=90, ha="center", va="center", fontsize=LAB)
    a_axes[0].legend(loc="lower left", bbox_to_anchor=(-0.03, -0.04), fontsize=TICK,
                     handlelength=0.9, handletextpad=0.3, labelspacing=0.18)
    panel_label(fig, a_axes[0], "a", dy=0.012)

    # ------------------------------------------------------------------ b
    ax_b = fig.add_subplot(bcd[0, 0])
    row_bands(ax_b)
    for index, pid in enumerate(SET_IDS):
        for group in ("COVID", "HC"):
            values = per_sample.loc[
                (per_sample["peak_set_id"] == pid) & (per_sample["group"] == group),
                "absolute_epm_per_total_mapped_alignments"].to_numpy(float)
            strip_row(ax_b, index, values, group, seed=9100 + index * 10 + (group == "HC"))
    ax_b.set_xscale("log")
    ax_b.set_xlim(1.6, 6000)
    ax_b.set_xticks([10, 100, 1000])
    ax_b.set_yticks(np.arange(len(PEAK_SETS)))
    ax_b.set_yticklabels([("CD14$^+$ " if family == "CD14" else "A549-ACE2 ") + label
                          for _, label, family in PEAK_SETS])
    ax_b.tick_params(axis="y", length=0, pad=2.0)
    ax_b.set_xlabel("eccDNAs per 10$^6$\nmapped alignments")
    ax_b.set_title("Absolute abundance", fontsize=LAB, pad=3.0)
    tidy(ax_b)
    panel_label(fig, ax_b, "b", dy=0.026)

    # ------------------------------------------------------------------ c
    ax_c = fig.add_subplot(bcd[0, 1])
    row_bands(ax_c)
    ax_c.axvline(0.0, color=ZERO_LINE, linewidth=0.8, zorder=2)
    for index, pid in enumerate(SET_IDS):
        for group in ("COVID", "HC"):
            values = per_sample.loc[
                (per_sample["peak_set_id"] == pid) & (per_sample["group"] == group),
                "log2_enrichment_ratio"].to_numpy(float)
            record = within.loc[(pid, group)]
            # The locked 95% interval of the median does the inferential work:
            # an interval clear of 0 is the significant within-cohort enrichment.
            interval = (float(np.log2(record["ratio_ci_low"])), float(np.log2(record["ratio_ci_high"])))
            strip_row(ax_c, index, values, group, seed=9300 + index * 10 + (group == "HC"),
                      interval=interval)
    ax_c.set_xlim(-0.62, 0.72)
    ax_c.set_xticks([-0.4, 0.0, 0.4])
    ax_c.set_yticks(np.arange(len(PEAK_SETS)))
    ax_c.set_yticklabels([])
    ax_c.tick_params(axis="y", length=0)
    ax_c.set_xlabel("log$_2$(observed /\nmatched expected)")
    ax_c.set_title("Relative enrichment", fontsize=LAB, pad=3.0)
    tidy(ax_c)
    panel_label(fig, ax_c, "c", dy=0.026, x=ax_c.get_position().x0 - 0.014)

    # ------------------------------------------------------------------ d, hero
    ax_d = fig.add_subplot(bcd[0, 2])
    row_bands(ax_d)
    ax_d.axvline(0.0, color=ZERO_LINE, linewidth=0.8, zorder=2)
    for column in range(len(PEAK_SETS)):
        relative = deltas[1:, column]
        # A thin bar spans where the seven relative-enrichment variants land, so
        # the spread of the sensitivity analyses is visible before any marker.
        ax_d.plot([relative.min(), relative.max()], [column, column], color="#AAB4BE",
                  linewidth=0.8, solid_capstyle="round", zorder=3)
        for row in range(len(VARIANTS)):
            value = deltas[row, column]
            significant = qvalues[row, column] < 0.05
            colour = COVID_LINE if value >= 0 else HC_LINE
            if row == ABSOLUTE_ROW:
                ax_d.scatter([value], [column], s=20, marker="D", zorder=6,
                             facecolors=colour, edgecolors=DARK_SLATE, linewidths=0.45)
            elif row == BURDEN_MATCHED_ROW:
                ax_d.scatter([value], [column], s=17, marker="o", zorder=6,
                             facecolors=colour if significant else "white",
                             edgecolors=DARK_SLATE, linewidths=0.85)
            else:
                ax_d.scatter([value], [column], s=7, marker="o", zorder=5,
                             facecolors=colour if significant else "white",
                             edgecolors=colour, linewidths=0.45)
    ax_d.set_xlim(-0.88, 0.94)
    ax_d.set_xticks([-0.8, -0.4, 0.0, 0.4, 0.8])
    ax_d.set_yticks(np.arange(len(PEAK_SETS)))
    ax_d.set_yticklabels([])
    ax_d.tick_params(axis="y", length=0)
    ax_d.set_xlabel("Cliff's $\\delta$, COVID-19 vs HC\n(red, COVID-19 higher; blue, HC higher)")
    ax_d.set_title("Contrast between cohorts, across all eight prespecified analyses",
                   fontsize=LAB, pad=3.0)
    tidy(ax_d)
    panel_label(fig, ax_d, "d", dy=0.026, x=ax_d.get_position().x0 - 0.014)

    # ------------------------------------------------------------------ keys and notes
    box_b, box_d = ax_b.get_position(), ax_d.get_position()
    # Panel d's key sits below the axes: inside it would cover the data in the
    # bottom two rows.
    fig.legend(handles=[
        Line2D([], [], marker="D", linestyle="none", markersize=3.0, markerfacecolor=COVID_LINE,
               markeredgecolor=DARK_SLATE, markeredgewidth=0.45, label="Absolute abundance"),
        Line2D([], [], marker="o", linestyle="none", markersize=2.0, markerfacecolor="white",
               markeredgecolor=HC_LINE, markeredgewidth=0.45, label="Relative enrichment, 6 variants"),
        Line2D([], [], marker="o", linestyle="none", markersize=2.7, markerfacecolor="white",
               markeredgecolor=DARK_SLATE, markeredgewidth=0.85, label="Burden-matched subset"),
        Line2D([], [], marker="o", linestyle="none", markersize=2.4, markerfacecolor=COVID_LINE,
               markeredgecolor=COVID_LINE, markeredgewidth=0.45, label="filled, $q$ < 0.05")],
        loc="upper center", bbox_to_anchor=((box_b.x0 + box_d.x1) / 2.0, 0.092),
        ncol=4, columnspacing=1.1, handletextpad=0.35, fontsize=TICK)

    fig.text((box_b.x0 + box_d.x1) / 2.0, 0.038,
             "b and c, every sample is drawn as a dot; the vertical tick is the median, "
             "the bar is the interquartile range in b and the locked 95% confidence "
             "interval of the median in c.\n"
             "An interval clear of 0 in c is a significant within-cohort enrichment or "
             "depletion.  All $q$ are Benjamini-Hochberg corrected within each metric family.  "
             "$n$ = 39 COVID-19 and 39 HC throughout.",
             ha="center", va="top", fontsize=TICK, color=NOTE_GREY, linespacing=1.35)

    for extension in ("pdf", "svg", "png", "tiff"):
        kwargs = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        fig.savefig(outdir / f"{args.stem}.{extension}",
                    dpi=600 if extension in ("png", "tiff") else None, **kwargs)
    plt.close(fig)
    print(f"wrote {outdir}/{args.stem}.{{pdf,svg,png,tiff}}")

    export = []
    for row, label in enumerate(variant_labels):
        for column, (pid, set_label, family) in enumerate(PEAK_SETS):
            export.append({
                "variant": label,
                "peak_set_id": pid,
                "peak_set": ("CD14+ " if family == "CD14" else "A549-ACE2 ") + set_label,
                "cliffs_delta_covid_vs_hc": deltas[row, column],
                "bh_fdr": qvalues[row, column],
            })
    pd.DataFrame(export).to_csv(outdir / f"{args.stem}_panel_d_values.tsv", sep="\t", index=False)
    print(f"wrote {outdir}/{args.stem}_panel_d_values.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
