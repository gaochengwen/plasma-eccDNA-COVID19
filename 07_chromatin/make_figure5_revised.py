#!/usr/bin/env python3
"""Main-text Figure 5: the COVID-19 eccDNA change is quantitative, not positional.

Figure contract
---------------
Core conclusion
    In COVID-19, plasma eccDNA output from SARS-CoV-2-infected-cell regulatory
    chromatin rises about fourfold, but the chromatin compartment it comes from
    does not change: the disease signal is quantitative, not positional.

Archetype
    Quantitative grid with one hero panel (asymmetric).  Panel d spans two
    columns and carries the decisive evidence; e and f are deliberately
    quieter controls.

Panel map and evidence hierarchy
    a  establish the system   peak-centred breakpoint density, both marks
    b  main effect            absolute abundance inside the peak sets
    c  localization           enrichment over the matched placement expectation
    d  HERO                   exact waterfall decomposition of the difference
    e  design limitation      enrichment against detection burden (quiet)
    f  robustness             ten prespecified analysis variants (quiet)

Why the previous version was replaced
    The published five-panel figure spent its whole right column on single-locus
    views of BCL3 and PROCR.  The objective, prespecified locus-selection audit
    in ``analyse/BCL3`` finds neither locus eligible: BCL3 fails the primary
    direction and recurrence rule, PROCR fails the artifact-masked sensitivity
    analysis and keeps a qualifying call in only one COVID-19 sample.  The old
    layout also omitted the relative-enrichment panel that Figure 4 carries, so
    a reader could not tell that the large cohort difference in absolute
    abundance is a detection-burden effect.

Statistical integrity
    Every plotted number is read from the locked hg38 analysis tables.  No test,
    transformation, group definition, threshold or significance tier is
    recomputed, with one exception that is pure algebra: panel d splits each
    sample's log2 absolute abundance into three additive terms using

        absolute_epm = EPM_total * matched_expected_fraction * enrichment_ratio

    which holds exactly because the full-interval denominator equals the
    eligible eccDNA count.  Both the identity and the additivity of the three
    component differences are asserted at run time.  The third term equals the
    locked ``m1_group_only`` HC3 group effect to all reported digits, and the
    interval drawn on it is the locked HC3 95% interval.

Reviewer risk handled in the layout
    The A549-ACE2 peak sets are a secondary tier covering about 1% of the
    allowed genome, so the figure never invites a magnitude comparison against
    the Figure 4 CD14+ primary sets, and no panel implies infection specificity.

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

# Allowed hg38 placement space after blacklist and assembly-gap removal, the
# same denominator the analysis used for genome-wide breakpoint density.
ALLOWED_BP = 2831289217

A549_DATASET = "A549_GSE179184_official"
A549_SETS: List[Tuple[str, str]] = [
    ("A549_GSE179184_H3K27ac_official", "H3K27ac"),
    ("A549_GSE179184_H3K4me3_official", "H3K4me3"),
]
MARKS = [mark for _, mark in A549_SETS]

# ---------------------------------------------------------------------------
# Palette: one signal family (cohort red/blue, inherited unchanged from
# Figure 4 so a colour never means two things across the manuscript) and one
# neutral slate family for the accounting and robustness panels.  Every
# encoding is also carried by position, lightness or marker fill, so the figure
# survives greyscale printing.
LINE_COLORS = {"COVID": "#C84F50", "HC": "#36617B"}
BOX_COLORS = {"COVID": "#FF8080", "HC": "#8BABD3"}
POINT_COLORS = {"COVID": "#D70000", "HC": "#034E61"}
REFERENCE_GREY = "#B2B2B2"
ZERO_LINE = "#49566D"

SLATE = {
    "total": "#22303C",
    "burden": "#4A6070",
    "expectation": "#8FA0AC",
    "positional": "#C6D0D7",
}
MARK_MARKERS = {"H3K27ac": "o", "H3K4me3": "s"}
ACCENT = "#1F5B7A"
BAND_FILL = "#F0F2F4"
NOTE_GREY = "#5A6673"
LABEL_SLATE = "#3C4C5A"
CONNECTOR = "#8A96A2"

GROUP_OFFSET = 0.19
BOX_WIDTH = 0.30

ABSOLUTE_BLOCK = "Absolute abundance"
RELATIVE_BLOCK = "Relative enrichment (observed / matched expected)"

# Typography: three sizes only, plus one smaller size reserved for the dense
# category labels of panel f.
LAB, TICK, PAN, DENSE = 6.0, 5.5, 8.0, 5.0


def sig_tier(q) -> str:
    """Manuscript convention: * q<0.05, ** q<0.01, *** q<0.001, **** q<1e-4."""
    if q is None or pd.isna(q):
        return ""
    q = float(q)
    return "****" if q < 1e-4 else "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "NS"


def setup_matplotlib(mpl) -> None:
    """Locked project rcParams, so this figure matches Figure 4 exactly."""
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
            "pdf.fonttype": 42,          # editable TrueType text in PDF
            "ps.fonttype": 42,
            "svg.fonttype": "none",      # editable text nodes in SVG
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
            "legend.handlelength": 1.4,
            "legend.handletextpad": 0.4,
            "legend.labelspacing": 0.28,
            "legend.borderpad": 0.15,
        }
    )


def tidy(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=TICK, length=2.0, pad=1.5)
    ax.xaxis.label.set_size(LAB)
    ax.yaxis.label.set_size(LAB)


def panel_label(fig, ax, letter: str, dx: float = -0.055, dy: float = 0.018) -> None:
    box = ax.get_position()
    fig.text(
        max(0.003, box.x0 + dx),
        min(0.999, box.y1 + dy),
        letter,
        fontsize=PAN,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def note_under(fig, axes, text: str, drop: float = 0.042) -> None:
    """Reading note on a baseline shared by every panel in the same block."""
    x0 = min(ax.get_position().x0 for ax in axes)
    x1 = max(ax.get_position().x1 for ax in axes)
    y0 = min(ax.get_position().y0 for ax in axes)
    fig.text(
        (x0 + x1) / 2.0,
        y0 - drop,
        text,
        ha="center",
        va="top",
        fontsize=TICK,
        color=NOTE_GREY,
        linespacing=1.3,
    )


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
        pivot = sub[sub["group"] == group].pivot_table(
            index="offset_bp", columns="sample_id", values="fold"
        )
        values = pivot.to_numpy(float)
        out[group] = (
            pivot.index.to_numpy(float),
            np.nanmean(values, axis=1),
            np.nanstd(values, axis=1, ddof=1) / np.sqrt(values.shape[1]),
        )
    return out


def profile_panel(ax, folds, title: str, show_legend: bool, show_ylabel: bool) -> None:
    for group in ("COVID", "HC"):
        offsets, mean, sem = folds[group]
        ax.plot(offsets / 1000.0, mean, "-", color=LINE_COLORS[group], linewidth=0.85,
                label="COVID-19" if group == "COVID" else "HC", solid_capstyle="round")
        ax.fill_between(offsets / 1000.0, mean - sem, mean + sem,
                        color=LINE_COLORS[group], alpha=0.18, linewidth=0)
    ax.axhline(1.0, color=REFERENCE_GREY, linewidth=0.7, linestyle=(0, (3, 2)), zorder=0,
               label="Genome average" if show_legend else None)
    ax.axvline(0.0, color=ZERO_LINE, linewidth=0.45, zorder=0)
    ax.set_title(title, fontsize=LAB, pad=2.5)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xticklabels(["-2", "", "peak\ncentre", "", "2"])
    ax.set_xlabel("Distance from peak centre (kb)")
    if show_ylabel:
        ax.set_ylabel("eccDNA breakpoint density\n(fold over genome average)")
    if show_legend:
        ax.legend(loc="lower left", bbox_to_anchor=(-0.01, -0.02), fontsize=TICK)
    tidy(ax)


# ---------------------------------------------------------------------------
# Panels b and c


def grouped_boxes(ax, frame: pd.DataFrame, value_col: str, seed_base: int) -> np.ndarray:
    """Two cohorts per mark, box plus every source value as a jittered point."""
    x = np.arange(len(MARKS), dtype=float)
    for group_index, group in enumerate(("COVID", "HC")):
        offset = (-GROUP_OFFSET) if group == "COVID" else GROUP_OFFSET
        data = [
            frame.loc[(frame["mark"] == mark) & (frame["group"] == group), value_col]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
            for mark in MARKS
        ]
        bp = ax.boxplot(
            data,
            positions=x + offset,
            widths=BOX_WIDTH,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 0.75},
            whiskerprops={"color": "black", "linewidth": 0.55},
            capprops={"color": "black", "linewidth": 0.55},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(BOX_COLORS[group])
            patch.set_edgecolor("black")
            patch.set_linewidth(0.55)
        for mark_index, (position, values) in enumerate(zip(x + offset, data)):
            # Fixed seeds keep the visual jitter reproducible while every source
            # value stays exactly as analysed.
            rng = np.random.default_rng(seed_base + group_index * 100 + mark_index)
            ax.scatter(
                np.full(len(values), position) + rng.uniform(-0.078, 0.078, size=len(values)),
                values,
                s=2.8,
                color=POINT_COLORS[group],
                edgecolors="none",
                alpha=0.88,
                zorder=3,
            )
    ax.set_xticks(x)
    ax.set_xticklabels(MARKS)
    ax.set_xlim(-0.55, len(MARKS) - 0.45)
    return x


def bracket(ax, xi: float, low: float, high: float, tier: str, text_y: float) -> None:
    ax.plot(
        [xi - GROUP_OFFSET, xi - GROUP_OFFSET, xi + GROUP_OFFSET, xi + GROUP_OFFSET],
        [low, high, high, low],
        color="black",
        linewidth=0.55,
        clip_on=False,
    )
    ax.text(xi, text_y, tier, ha="center", va="bottom", fontsize=TICK, fontweight="bold")


# ---------------------------------------------------------------------------
# Panel d, the hero


DECOMP_STEPS = ["burden", "expectation", "positional"]
DECOMP_ORDER = DECOMP_STEPS + ["total"]
DECOMP_LABELS = {
    "burden": "eccDNA\nburden",
    "expectation": "Matched\nexpectation",
    "positional": "Positional\npreference",
    "total": "Observed\ntotal",
}


def decompose(absolute: pd.DataFrame) -> pd.DataFrame:
    """Exact additive split of the cohort difference in log2 absolute abundance.

    absolute_epm = EPM_total * expected_fraction * enrichment_ratio, so the
    difference of cohort means of log2 absolute abundance is the sum of the
    differences of cohort means of the three log2 terms.  Nothing is estimated
    here; the assertions below fail if the identity does not hold.
    """
    frame = absolute.copy()
    frame["epm_total"] = frame["eligible_eccdna_count"] / frame["mapped_alignments_idxstats"] * 1e6
    residual = np.abs(
        frame["observed_fraction"] * frame["epm_total"]
        - frame["absolute_epm_per_total_mapped_alignments"]
    ) / frame["absolute_epm_per_total_mapped_alignments"]
    assert float(residual.max()) < 1e-10, "absolute-abundance identity does not hold"

    terms = {
        "total": "absolute_epm_per_total_mapped_alignments",
        "burden": "epm_total",
        "expectation": "mean_shuffled_overlap_fraction",
        "positional": "enrichment_ratio",
    }
    rows = []
    for mark in MARKS:
        sub = frame[frame["mark"] == mark]
        values = {}
        for name, column in terms.items():
            covid = np.log2(sub.loc[sub["group"] == "COVID", column].to_numpy(float))
            hc = np.log2(sub.loc[sub["group"] == "HC", column].to_numpy(float))
            values[name] = float(covid.mean() - hc.mean())
        assert abs(values["burden"] + values["expectation"] + values["positional"] - values["total"]) < 1e-9
        for name, value in values.items():
            rows.append({"mark": mark, "term": name, "delta_log2": value})
    return pd.DataFrame(rows)


def waterfall(ax, values: Dict[str, float], hc3_row, mark: str, show_ylabel: bool) -> None:
    """One mark: three additive steps landing exactly on the observed total."""
    positions = np.arange(len(DECOMP_ORDER), dtype=float)

    running = 0.0
    tops: List[float] = []
    for index, term in enumerate(DECOMP_STEPS):
        delta = values[term]
        bottom, top = (running, running + delta) if delta >= 0 else (running + delta, running)
        ax.bar(positions[index], top - bottom, bottom=bottom, width=0.62,
               color=SLATE[term], edgecolor="black", linewidth=0.4, zorder=2)
        running += delta
        tops.append(running)
        # Connector to the next bar makes the additivity visible at a glance.
        ax.plot([positions[index] - 0.31, positions[index + 1] + 0.31], [running, running],
                linestyle=(0, (2, 1.6)), linewidth=0.45, color=CONNECTOR, zorder=1)

    ax.bar(positions[-1], values["total"], width=0.62, color=SLATE["total"],
           edgecolor="black", linewidth=0.4, zorder=2)
    tops.append(values["total"])

    # The positional term is the panel's punch line, so it alone carries the
    # locked HC3 interval and its significance tier.
    beta = float(hc3_row["beta"])
    ci_low, ci_high = float(hc3_row["ci_low"]), float(hc3_row["ci_high"])
    positional_index = DECOMP_STEPS.index("positional")
    ax.errorbar(
        positions[positional_index], beta,
        yerr=[[beta - ci_low], [ci_high - beta]],
        fmt="none", ecolor="black", elinewidth=0.55, capsize=1.4, capthick=0.55, zorder=4,
    )
    ax.text(positions[positional_index], ci_low - 0.09,
            sig_tier(hc3_row["bh_fdr_within_family"]), ha="center", va="top",
            fontsize=TICK, color=LABEL_SLATE)

    for index, term in enumerate(DECOMP_ORDER):
        anchor = ci_high if term == "positional" else tops[index]
        ax.text(positions[index], anchor + 0.09, f"{values[term]:+.2f}",
                ha="center", va="bottom", fontsize=TICK, color=LABEL_SLATE)
    # Fold change only where it is interpretable: the driver and the result.
    for term in ("burden", "total"):
        ax.text(positions[DECOMP_ORDER.index(term)], -0.36, f"{2 ** values[term]:.1f}×",
                ha="center", va="top", fontsize=TICK, color=SLATE["total"], fontweight="bold")

    ax.axhline(0.0, color=ZERO_LINE, linewidth=0.7, zorder=1)
    ax.set_xticks(positions)
    ax.set_xticklabels([DECOMP_LABELS[term] for term in DECOMP_ORDER], linespacing=1.15)
    ax.set_xlim(-0.62, len(DECOMP_ORDER) - 0.38)
    ax.set_ylim(-0.70, 2.72)
    ax.set_title(mark, fontsize=LAB, pad=2.5)
    if show_ylabel:
        ax.set_ylabel("COVID-19 − HC difference in\nmean log$_2$ abundance")
    else:
        ax.set_yticklabels([])
    tidy(ax)


# ---------------------------------------------------------------------------
# Panel f


def robustness_rows(tables: Path) -> pd.DataFrame:
    """Collect every locked COVID-19 vs HC Cliff's delta for the A549 peak sets."""

    def read(name: str) -> pd.DataFrame:
        frame = pd.read_csv(tables / name, sep="\t")
        return frame[frame["peak_set_id"].str.contains("A549_GSE179184", na=False)]

    absolute = read("absolute_burden_group_stats.tsv")
    relative = read("relative_enrichment_group_stats.tsv")
    downsampled = read("p0_downsampled_group_stats.tsv")
    matched = read("p0_burden_matched_subset.tsv")
    mappability = read("ext_mappability_group_stats.tsv")

    specs = [
        (ABSOLUTE_BLOCK, "Full\ninterval", absolute, {"mode": "full_interval"}, "bh_fdr_all_tests"),
        (ABSOLUTE_BLOCK, "Junction\nbreakpoints", absolute, {"mode": "junction_start_end"}, "bh_fdr_all_tests"),
        (ABSOLUTE_BLOCK, "Midpoint", absolute, {"mode": "midpoint"}, "bh_fdr_all_tests"),
        (RELATIVE_BLOCK, "Full\ninterval", relative, {"mode": "full_interval"}, "bh_fdr_all_tests"),
        (RELATIVE_BLOCK, "Junction\nbreakpoints", relative, {"mode": "junction_start_end"}, "bh_fdr_all_tests"),
        (RELATIVE_BLOCK, "Midpoint", relative, {"mode": "midpoint"}, "bh_fdr_all_tests"),
        (RELATIVE_BLOCK, "Down-\nsampled\n4,000", downsampled,
         {"mode": "full_interval", "downsample_target": 4000}, "bh_fdr_within_family"),
        (RELATIVE_BLOCK, "Down-\nsampled\n20,000", downsampled,
         {"mode": "full_interval", "downsample_target": 20000}, "bh_fdr_within_family"),
        (RELATIVE_BLOCK, "Burden-\nmatched", matched, {"mode": "full_interval"}, "bh_fdr_within_family"),
        (RELATIVE_BLOCK, "Uniquely\nmappable", mappability, {"mode": "full_interval"},
         "bh_fdr_group_within_family"),
    ]

    rows = []
    for order, (block, label, frame, selector, q_column) in enumerate(specs):
        for mark in MARKS:
            sub = frame[frame["mark"] == mark]
            for key, value in selector.items():
                sub = sub[sub[key] == value]
            if len(sub) != 1:
                raise SystemExit(f"panel f: {block} / {label} / {mark} matched {len(sub)} rows")
            record = sub.iloc[0]
            rows.append({
                "order": order,
                "block": block,
                "label": label,
                "mark": mark,
                "delta": float(record["cliffs_delta_covid_vs_hc"]),
                "q": float(record[q_column]),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--outdir", default=None, help="default: <analysis-dir>/figures")
    parser.add_argument("--stem", default="figure5_revised")
    args = parser.parse_args(argv)

    analysis = Path(args.analysis_dir)
    tables = analysis / "tables"
    outdir = Path(args.outdir) if args.outdir else analysis / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    density = pd.read_csv(tables / "ext_peak_centred_density.tsv", sep="\t")
    absolute = pd.read_csv(tables / "sample_absolute_burden.tsv", sep="\t")
    absolute = absolute[
        (absolute["dataset"] == A549_DATASET) & (absolute["mode"] == "full_interval")
    ].copy()

    abs_stats = pd.read_csv(tables / "absolute_burden_group_stats.tsv", sep="\t")
    abs_stats = abs_stats[
        (abs_stats["dataset"] == A549_DATASET) & (abs_stats["mode"] == "full_interval")
    ].set_index("mark")
    rel_stats = pd.read_csv(tables / "relative_enrichment_group_stats.tsv", sep="\t")
    rel_stats = rel_stats[
        (rel_stats["dataset"] == A549_DATASET) & (rel_stats["mode"] == "full_interval")
    ].set_index("mark")

    within = pd.read_csv(tables / "p0_within_group_enrichment.tsv", sep="\t")
    within = within[
        (within["dataset"] == A549_DATASET) & (within["mode"] == "full_interval")
    ].set_index(["mark", "group"])

    hc3 = pd.read_csv(tables / "p0_covariate_adjusted_hc3.tsv", sep="\t")
    hc3 = hc3[
        (hc3["dataset"] == A549_DATASET)
        & (hc3["mode"] == "full_interval")
        & (hc3["model"] == "m1_group_only")
        & (hc3["term"] == "COVID_vs_HC")
    ].set_index("mark")

    burden_assoc = pd.read_csv(tables / "p0_burden_association.tsv", sep="\t")
    burden_assoc = burden_assoc[
        (burden_assoc["dataset"] == A549_DATASET) & (burden_assoc["mode"] == "full_interval")
    ].set_index(["mark", "stratum"])

    matched_window = pd.read_csv(tables / "p0_burden_matched_subset.tsv", sep="\t")
    matched_window = matched_window[
        (matched_window["dataset"] == A549_DATASET) & (matched_window["mode"] == "full_interval")
    ].iloc[0]

    decomposition = decompose(absolute)
    robust = robustness_rows(tables)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    setup_matplotlib(mpl)

    fig = plt.figure(figsize=(183 / 25.4, 184 / 25.4))
    # Asymmetric hierarchy: the hero row is the tallest and the control row the
    # shortest, so evidence weight is visible before any label is read.
    outer = GridSpec(3, 1, figure=fig, height_ratios=[1.00, 1.22, 0.96], hspace=0.64,
                     left=0.062, right=0.988, top=0.952, bottom=0.086)
    row1 = outer[0, 0].subgridspec(1, 3, wspace=0.42)
    row2 = outer[1, 0].subgridspec(1, 3, wspace=0.42)
    row3 = outer[2, 0].subgridspec(1, 3, wspace=0.34)

    # ------------------------------------------------------------------ a
    profile_axes = []
    for index, (peak_set_id, mark) in enumerate(A549_SETS):
        ax = fig.add_subplot(row1[0, index])
        profile_axes.append(ax)
        profile_panel(ax, density_fold(density, peak_set_id), f"A549-ACE2 {mark}",
                      show_legend=(index == 0), show_ylabel=(index == 0))
    # One shared y range so the mark-specific contrast is a visual fact rather
    # than an artefact of two independently scaled axes.
    low = min(ax.get_ylim()[0] for ax in profile_axes)
    high = max(ax.get_ylim()[1] for ax in profile_axes)
    for ax in profile_axes:
        ax.set_ylim(low, high)
    profile_axes[1].set_yticklabels([])
    panel_label(fig, profile_axes[0], "a", dx=-0.052)

    # ------------------------------------------------------------------ b
    ax_b = fig.add_subplot(row1[0, 2])
    x = grouped_boxes(ax_b, absolute, "absolute_epm_per_total_mapped_alignments", seed_base=7100)
    ax_b.set_yscale("log")
    low, high = ax_b.get_ylim()
    for xi, mark in zip(x, MARKS):
        bracket(ax_b, xi, high * 1.06, high * 1.30,
                sig_tier(abs_stats.loc[mark, "bh_fdr_all_tests"]), high * 1.42)
    ax_b.set_ylim(low, high * 5.5)
    ax_b.set_ylabel("eccDNAs overlapping peak set\nper 10$^6$ mapped alignments")
    ax_b.set_xlabel("A549-ACE2 peak set")
    ax_b.set_title("Absolute abundance", fontsize=LAB, pad=2.5)
    tidy(ax_b)
    # The only cohort key for the box panels: panel a already defines the line
    # colours and panel c inherits this one.
    ax_b.legend(
        handles=[
            Patch(facecolor=BOX_COLORS["COVID"], edgecolor="black", linewidth=0.45, label="COVID-19"),
            Patch(facecolor=BOX_COLORS["HC"], edgecolor="black", linewidth=0.45, label="HC"),
        ],
        loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=2,
        handlelength=1.1, handleheight=0.8, columnspacing=1.2, borderaxespad=0.1,
    )
    panel_label(fig, ax_b, "b", dx=-0.058)

    # ------------------------------------------------------------------ c
    ax_c = fig.add_subplot(row2[0, 0])
    x = grouped_boxes(ax_c, absolute, "log2_enrichment_ratio", seed_base=7300)
    ax_c.axhline(0.0, color=ZERO_LINE, linewidth=0.7, zorder=0)
    low, high = ax_c.get_ylim()
    span = high - low
    # Tier above each individual box: one-sample signed-rank test of that cohort
    # against its own matched placement expectation.
    for xi, mark in zip(x, MARKS):
        for group in ("COVID", "HC"):
            offset = (-GROUP_OFFSET) if group == "COVID" else GROUP_OFFSET
            ax_c.text(xi + offset, high + 0.02 * span,
                      sig_tier(within.loc[(mark, group), "bh_fdr_within_family"]),
                      ha="center", va="bottom", fontsize=TICK, color=POINT_COLORS[group])
        bracket(ax_c, xi, high + 0.16 * span, high + 0.24 * span,
                sig_tier(rel_stats.loc[mark, "bh_fdr_all_tests"]), high + 0.26 * span)
    ax_c.set_ylim(low, high + 0.46 * span)
    ax_c.set_ylabel("log$_2$(observed / matched expected)")
    ax_c.set_xlabel("A549-ACE2 peak set")
    ax_c.set_title("Relative enrichment", fontsize=LAB, pad=2.5)
    tidy(ax_c)
    panel_label(fig, ax_c, "c", dx=-0.052)

    # ------------------------------------------------------------------ d, hero
    d_grid = row2[0, 1:3].subgridspec(1, 2, wspace=0.09)
    d_axes = []
    for index, mark in enumerate(MARKS):
        ax = fig.add_subplot(d_grid[0, index])
        d_axes.append(ax)
        values = {
            row.term: float(row.delta_log2)
            for row in decomposition[decomposition["mark"] == mark].itertuples()
        }
        waterfall(ax, values, hc3.loc[mark], mark, show_ylabel=(index == 0))
    covid_median = absolute.loc[(absolute["mark"] == "H3K27ac") & (absolute["group"] == "COVID"),
                                "eligible_eccdna_count"].median()
    hc_median = absolute.loc[(absolute["mark"] == "H3K27ac") & (absolute["group"] == "HC"),
                             "eligible_eccdna_count"].median()
    d_box_left = d_axes[0].get_position()
    d_box_right = d_axes[1].get_position()
    fig.text(
        (d_box_left.x0 + d_box_right.x1) / 2.0,
        max(d_box_left.y1, d_box_right.y1) + 0.030,
        "Decomposition of the cohort difference in absolute abundance",
        ha="center", va="bottom", fontsize=LAB, color=LABEL_SLATE,
    )
    panel_label(fig, d_axes[0], "d", dx=-0.058, dy=0.046)

    # ------------------------------------------------------------------ e
    window_low = float(matched_window["matched_window_low"])
    window_high = float(matched_window["matched_window_high"])
    e_grid = row3[0, 0].subgridspec(2, 1, hspace=0.18)
    e_axes = []
    for index, mark in enumerate(MARKS):
        ax = fig.add_subplot(e_grid[index, 0])
        e_axes.append(ax)
        ax.axvspan(np.log10(window_low), np.log10(window_high), color=BAND_FILL,
                   linewidth=0, zorder=0)
        ax.axhline(0.0, color=ZERO_LINE, linewidth=0.6, zorder=1)
        for group in ("COVID", "HC"):
            sub = absolute[(absolute["mark"] == mark) & (absolute["group"] == group)]
            ax.scatter(np.log10(sub["eligible_eccdna_count"].to_numpy(float)),
                       sub["log2_enrichment_ratio"].to_numpy(float),
                       s=4.0, marker="o", facecolors=POINT_COLORS[group], edgecolors="none",
                       alpha=0.9, zorder=3)
        ax.set_xlim(3.55, 6.42)
        ax.set_ylim(-0.58, 1.06)
        ax.set_yticks([-0.4, 0.0, 0.4])
        ax.text(0.015, 0.97, mark, transform=ax.transAxes, ha="left", va="top",
                fontsize=TICK, color=LABEL_SLATE)
        covid = burden_assoc.loc[(mark, "COVID")]
        hc = burden_assoc.loc[(mark, "HC")]
        ax.text(
            0.015, 0.022,
            f"$\\rho_{{COVID}}$ = {float(covid['spearman_rho']):+.2f}"
            f"{sig_tier(covid['bh_fdr_within_family']).replace('NS', '')}"
            f"    $\\rho_{{HC}}$ = {float(hc['spearman_rho']):+.2f}"
            f"{sig_tier(hc['bh_fdr_within_family']).replace('NS', '')}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=TICK, color=LABEL_SLATE,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.8},
        )
        if index == 0:
            ax.set_xticklabels([])
            ax.set_title("Burden separates the cohorts", fontsize=LAB, pad=2.5)
            # Direct labels: the cohorts are spatially fixed on this axis, so a
            # legend would only add eye travel.
            ax.text(4.42, 0.72, "HC", color=POINT_COLORS["HC"], fontsize=TICK,
                    ha="center", va="center", fontweight="bold")
            ax.text(5.90, 0.72, "COVID-19", color=POINT_COLORS["COVID"], fontsize=TICK,
                    ha="center", va="center", fontweight="bold")
        else:
            ax.set_xlabel("log$_{10}$ eccDNAs detected per sample")
        tidy(ax)
    fig.text(
        e_axes[0].get_position().x0 - 0.036,
        (e_axes[0].get_position().y1 + e_axes[1].get_position().y0) / 2.0,
        "log$_2$(observed / matched expected)",
        rotation=90, ha="center", va="center", fontsize=LAB,
    )
    panel_label(fig, e_axes[0], "e", dx=-0.052)

    # ------------------------------------------------------------------ f
    ax_f = fig.add_subplot(row3[0, 1:3])
    labels = robust.drop_duplicates("order").sort_values("order")
    n_absolute = int((labels["block"] == ABSOLUTE_BLOCK).sum())
    xs = np.arange(len(labels), dtype=float)
    ax_f.axvspan(n_absolute - 0.5, len(labels) - 0.45, color=BAND_FILL, linewidth=0, zorder=0)
    ax_f.axhline(0.0, color=ZERO_LINE, linewidth=0.7, zorder=1)
    for mark_index, mark in enumerate(MARKS):
        offset = (mark_index - 0.5) * 0.30
        sub = robust[robust["mark"] == mark].sort_values("order")
        significant = sub["q"].to_numpy(float) < 0.05
        deltas = sub["delta"].to_numpy(float)
        for xi, delta, sig in zip(xs + offset, deltas, significant):
            ax_f.plot([xi, xi], [0.0, delta], color=ACCENT, linewidth=0.5,
                      alpha=1.0 if sig else 0.45, zorder=3)
        ax_f.scatter(xs[significant] + offset, deltas[significant], s=12.0,
                     marker=MARK_MARKERS[mark], facecolors=ACCENT, edgecolors=ACCENT,
                     linewidths=0.5, zorder=4)
        ax_f.scatter(xs[~significant] + offset, deltas[~significant], s=12.0,
                     marker=MARK_MARKERS[mark], facecolors="white", edgecolors=ACCENT,
                     linewidths=0.5, zorder=4)
    ax_f.set_xticks(xs)
    ax_f.set_xticklabels(labels["label"].tolist(), linespacing=1.15, fontsize=DENSE)
    ax_f.set_xlim(-0.55, len(labels) - 0.45)
    ax_f.set_ylim(-0.68, 0.94)
    ax_f.set_yticks([-0.5, -0.25, 0.0, 0.25, 0.5, 0.75])
    ax_f.set_ylabel("Cliff's $\\delta$, COVID-19 vs HC")
    # Direct block labels instead of another legend entry for the metric family.
    ax_f.text((n_absolute - 1) / 2.0, 0.88, ABSOLUTE_BLOCK, ha="center", va="center",
              fontsize=LAB, color=LABEL_SLATE)
    ax_f.text((n_absolute + len(labels) - 1) / 2.0, 0.88, RELATIVE_BLOCK, ha="center",
              va="center", fontsize=LAB, color=LABEL_SLATE)
    ax_f.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", markersize=2.8, markerfacecolor=ACCENT,
                   markeredgecolor=ACCENT, markeredgewidth=0.5, label="H3K27ac"),
            Line2D([], [], marker="s", linestyle="none", markersize=2.6, markerfacecolor=ACCENT,
                   markeredgecolor=ACCENT, markeredgewidth=0.5, label="H3K4me3"),
            Line2D([], [], marker="o", linestyle="none", markersize=2.8, markerfacecolor=ACCENT,
                   markeredgecolor=ACCENT, markeredgewidth=0.5, label="filled, $q$ < 0.05"),
            Line2D([], [], marker="o", linestyle="none", markersize=2.8, markerfacecolor="white",
                   markeredgecolor=ACCENT, markeredgewidth=0.5, label="open, $q$ ≥ 0.05"),
        ],
        loc="upper right", bbox_to_anchor=(1.005, 0.79), ncol=2,
        handlelength=0.9, handletextpad=0.35, columnspacing=1.2, borderaxespad=0.1,
    )
    ax_f.set_title("Robustness of the cohort contrast", fontsize=LAB, pad=2.5)
    tidy(ax_f)
    panel_label(fig, ax_f, "f", dx=-0.052)

    # ------------------------------------------------------------------ notes
    note_under(fig, [ax_c], "coloured tiers, each cohort vs its matched\n"
                            "expectation; bracket, COVID-19 vs HC")
    note_under(fig, d_axes,
               "burden + expectation + positional = observed total (exact); "
               "error bar, locked HC3 95% CI of the positional term\n"
               f"median detected eccDNAs {covid_median:,.0f} (COVID-19) vs {hc_median:,.0f} (HC)")
    note_under(fig, e_axes[1:], "shaded band, burden-matched window used in f;\n"
                                "$\\rho$, within-cohort Spearman correlation")
    note_under(fig, [ax_f], "shaded block, relative-enrichment variants")

    for extension in ("pdf", "svg", "png", "tiff"):
        kwargs = {}
        if extension == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(outdir / f"{args.stem}.{extension}",
                    dpi=600 if extension in ("png", "tiff") else None, **kwargs)
    plt.close(fig)
    print(f"wrote {outdir}/{args.stem}.{{pdf,svg,png,tiff}}")

    # Machine-readable record of every derived value, so the figure can be
    # audited without rerunning it.
    decomposition.to_csv(outdir / f"{args.stem}_panel_d_values.tsv", sep="\t", index=False)
    robust_export = robust.copy()
    robust_export["label"] = robust_export["label"].str.replace("\n", " ", regex=False)
    robust_export.to_csv(outdir / f"{args.stem}_panel_f_values.tsv", sep="\t", index=False)
    print(f"wrote {outdir}/{args.stem}_panel_d_values.tsv and _panel_f_values.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
