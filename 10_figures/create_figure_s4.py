#!/usr/bin/env python3
"""Create Supplementary Figure S4 (Submit numbering): eccDNA-burden robustness
across callers and call sets.

Every plotted number is read from the locked Circle-Map/Circle_finder robustness
tables under ``analyse/Circle_finder/tables``; nothing is recomputed. The only
derived quantities are order statistics (median, quartiles) of the locked
per-sample columns, and the script asserts that they reproduce the values quoted
in ``analyse/Circle_finder/RESULTS_SUMMARY.md`` and in the manuscript.

Panel a follows the format and inks of Submit/Figure_S2 panel c: open spines,
a frameless patch legend above the axes, filled boxes with black outlines and
jittered cohort-coloured points.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
from matplotlib.transforms import blended_transform_factory


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "analyse" / "Circle_finder" / "tables"
OUTPUT_DIR = ROOT / "revise"
OUTPUT_STEM = OUTPUT_DIR / "Figure_S4_callset_robustness"

MM = 1.0 / 25.4
FIG_W_MM = 183.0
# Taller canvas at unchanged width. Type sizes are absolute (points), so the extra
# height is spent entirely on the panel boxes; the margins and the inter-panel gaps,
# which are sized by the text that sits in them, keep their previous values.
FIG_H_MM = 200.0

# Cohort inks are the locked project palette (Figure 1b/1c vector colours, also
# used by Figure 4, Figure 5 and Submit/Figure_S2 panel c).
POINT = {"COVID-19": "#D70000", "HC": "#034E61"}
BOX_FILL = {"COVID-19": "#FF8080", "HC": "#8BABD3"}
BOX_ORDER = ("COVID-19", "HC")
# Neutral family reserved by the project for quantities that are not cohorts.
NEUTRAL_DARK = "#3C4C5A"
NEUTRAL_MID = "#6E7A87"
NEUTRAL_LIGHT = "#AFBAC4"
# Row band for the panel b table: a pale tint of the neutral family, light enough that
# the black row labels and the numerals keep full contrast on it.
BAND = "#F1F3F6"
ACCENT = "#1F5B7A"
RULE = "#49566D"
SEED = 20260803

# Call sets in the order used throughout the figure and Table S5.
CALLSETS = (
    "circlemap_current",
    "circlemap_methods",
    "circlemap_strict",
    "circlemap_very_strict",
    "circlemap_artifact_masked",
    "circlemap_high_support_masked",
    "consensus_t10",
)
SHORT_LABEL = {
    "circlemap_current": "Current\nCircle-Map",
    "circlemap_methods": "Methods\nfilter",
    "circlemap_strict": "split ≥10\nscore >1,000",
    "circlemap_very_strict": "split ≥20\nscore >2,000",
    "circlemap_artifact_masked": "Artifact\nmasked",
    "circlemap_high_support_masked": "High-support\n+ masked",
    "consensus_t10": "Two-caller\nconsensus",
}
LONG_LABEL = {
    "circlemap_current": "Current Circle-Map",
    "circlemap_methods": "Methods filter",
    "circlemap_strict": "split ≥10, score >1,000",
    "circlemap_very_strict": "split ≥20, score >2,000",
    "circlemap_artifact_masked": "Artifact masked",
    "circlemap_high_support_masked": "High-support + masked",
    "consensus_t10": "Two-caller consensus",
}

LENGTH_BINS = ("lt200", "200_399", "400_599", "600_999", "1000_1999", "ge2000")
BIN_LABEL = {
    "lt200": "<200",
    "200_399": "200–\n399",
    "400_599": "400–\n599",
    "600_999": "600–\n999",
    "1000_1999": "1,000–\n1,999",
    "ge2000": "≥2,000",
}

TOLERANCES = (2, 5, 10, 20)
PRIMARY_TOLERANCE = 10

# Values quoted in RESULTS_SUMMARY.md and in the Results paragraph on call-set
# robustness. The script refuses to emit a figure that disagrees with them.
EXPECTED_MEDIAN_EPM = {
    ("circlemap_current", "COVID-19"): 3143.4,
    ("circlemap_current", "HC"): 814.2,
    ("circlemap_methods", "COVID-19"): 2052.6,
    ("circlemap_methods", "HC"): 551.5,
    ("consensus_t10", "COVID-19"): 1882.1,
    ("consensus_t10", "HC"): 527.0,
}
EXPECTED_CONCORDANCE = {
    2: (0.927, 0.703, 0.660, 9_076_736),
    5: (0.932, 0.705, 0.663, 9_109_452),
    10: (0.935, 0.705, 0.665, 9_130_043),
    20: (0.939, 0.706, 0.667, 9_152_395),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def configure_plotting() -> None:
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        family = "Arial"
    except ValueError:
        family = "Liberation Sans"
    mpl.rcParams.update(
        {
            "font.family": family,
            "font.size": 5.5,
            "axes.labelsize": 6.5,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.3,
            "ytick.major.size": 2.3,
            "xtick.minor.size": 1.2,
            "ytick.minor.size": 1.2,
            "lines.solid_capstyle": "round",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "custom",
            "mathtext.rm": family,
            "mathtext.it": f"{family}:italic",
            "mathtext.bf": f"{family}:bold",
        }
    )


def minus(text: str) -> str:
    """Replace ASCII hyphen-minus with the typographic minus sign."""
    return text.replace("-", "−")


def scientific(value: float, digits: int = 1) -> str:
    """Mantissa x 10^exponent, laid out by mathtext so the spacing is even."""
    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def asterisks(q_value: float) -> str:
    for threshold, mark in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")):
        if q_value < threshold:
            return mark
    return "n.s."


def panel_label(figure: plt.Figure, x_mm: float, y_mm: float, label: str) -> None:
    figure.text(
        x_mm / FIG_W_MM,
        y_mm / FIG_H_MM,
        label,
        fontsize=8,
        fontweight="bold",
        va="baseline",
        ha="left",
    )


def add_axes_mm(
    figure: plt.Figure, left: float, bottom: float, width: float, height: float
) -> plt.Axes:
    return figure.add_axes(
        [left / FIG_W_MM, bottom / FIG_H_MM, width / FIG_W_MM, height / FIG_H_MM]
    )


# --------------------------------------------------------------------------- #
# Panel a: per-sample EPM for the seven call sets
# --------------------------------------------------------------------------- #
def panel_burden(
    axis: plt.Axes, metrics: list[dict[str, str]], tests: list[dict[str, str]]
) -> None:
    rng = np.random.default_rng(SEED)
    offset = 0.185
    values: dict[tuple[str, str], np.ndarray] = {}
    for callset in CALLSETS:
        for group in BOX_ORDER:
            values[(callset, group)] = np.asarray(
                [
                    float(row["epm"])
                    for row in metrics
                    if row["callset"] == callset and row["group"] == group
                ]
            )

    for (callset, group), sample_values in values.items():
        assert sample_values.size == 39, f"{callset}/{group}: n = {sample_values.size}"
        expected = EXPECTED_MEDIAN_EPM.get((callset, group))
        if expected is not None:
            assert abs(float(np.median(sample_values)) - expected) < 0.1, (
                f"{callset}/{group} median EPM disagrees with RESULTS_SUMMARY.md"
            )

    positions = []
    ordered_values = []
    ordered_groups = []
    for index, callset in enumerate(CALLSETS):
        for sign, group in zip((-1, 1), BOX_ORDER):
            positions.append(index + sign * offset)
            ordered_values.append(values[(callset, group)])
            ordered_groups.append(group)

    boxplot = axis.boxplot(
        ordered_values,
        positions=positions,
        widths=0.26,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 0.7},
        whiskerprops={"color": "black", "linewidth": 0.6},
        capprops={"color": "black", "linewidth": 0.6},
        boxprops={"edgecolor": "black", "linewidth": 0.6},
    )
    for patch, group in zip(boxplot["boxes"], ordered_groups):
        patch.set_facecolor(BOX_FILL[group])
        patch.set_zorder(2)

    for position, group, sample_values in zip(positions, ordered_groups, ordered_values):
        jitter = rng.uniform(-0.095, 0.095, sample_values.size)
        axis.scatter(
            position + jitter,
            sample_values,
            s=4.0,
            marker="o",
            color=POINT[group],
            linewidth=0,
            alpha=0.9,
            zorder=3,
        )

    # Pairwise BH q values are the locked two-sided Mann-Whitney results.
    for index, callset in enumerate(CALLSETS):
        test = next(
            row for row in tests if row["callset"] == callset and row["metric"] == "epm"
        )
        top = max(
            float(values[(callset, group)].max()) for group in BOX_ORDER
        )
        bracket_y = top * 1.55
        axis.plot(
            [index - offset, index - offset, index + offset, index + offset],
            [bracket_y / 1.22, bracket_y, bracket_y, bracket_y / 1.22],
            color="black",
            linewidth=0.6,
            solid_joinstyle="miter",
        )
        axis.text(
            index,
            bracket_y * 1.1,
            asterisks(float(test["q_value_bh_within_metric"])),
            ha="center",
            va="bottom",
            fontsize=5.5,
        )

    axis.set_yscale("log")
    axis.set_ylim(4, 40_000)
    axis.set_yticks([10, 100, 1_000, 10_000])
    axis.set_yticklabels(["10", "100", "1,000", "10,000"])
    axis.set_ylabel("eccDNA per million mapped reads (EPM)")
    axis.set_xlim(-0.55, len(CALLSETS) - 0.45)
    axis.set_xticks(range(len(CALLSETS)))
    axis.set_xticklabels([SHORT_LABEL[c] for c in CALLSETS], linespacing=1.25)
    axis.tick_params(axis="x", length=0, pad=2.5)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(
        handles=[
            Patch(facecolor=BOX_FILL[group], edgecolor="black", linewidth=0.6, label=group)
            for group in BOX_ORDER
        ],
        loc="lower left",
        bbox_to_anchor=(0.0, 1.005),
        ncol=2,
        frameon=False,
        handlelength=1.5,
        handleheight=0.85,
        columnspacing=1.1,
        handletextpad=0.5,
    )


# --------------------------------------------------------------------------- #
# Panel b: covariate-adjusted COVID-19/HC fold change
# --------------------------------------------------------------------------- #
def panel_adjusted_effect(axis: plt.Axes, regressions: list[dict[str, str]]) -> None:
    rows = []
    for callset in CALLSETS:
        rows.append(
            next(
                row
                for row in regressions
                if row["callset"] == callset and row["term"] == "COVID-19_vs_HC"
            )
        )

    y_positions = np.arange(len(CALLSETS))[::-1]
    effects = np.asarray([float(row["multiplicative_effect"]) for row in rows])
    low = np.asarray([float(row["multiplicative_ci_low"]) for row in rows])
    high = np.asarray([float(row["multiplicative_ci_high"]) for row in rows])

    assert 3.30 < effects.min() < 3.32 and 5.07 < effects.max() < 5.09, (
        "adjusted fold changes disagree with the 3.31-5.08 range quoted in Results"
    )

    # The row label, the interval and the two numeric columns are ~150 mm apart end to
    # end, so each row is banded across that whole span. Alternating fill carries the
    # eye from a call set to its numbers without seven hairlines of extra ink; the rule
    # under the headers is what makes the block read as a table rather than as a plot
    # with captions floating beside it.
    transform = blended_transform_factory(axis.transAxes, axis.transData)
    band_left, band_right = -0.355, 1.68
    for banded_row in y_positions[::2]:
        axis.add_patch(
            Rectangle(
                (band_left, banded_row - 0.5),
                band_right - band_left,
                1.0,
                transform=transform,
                facecolor=BAND,
                edgecolor="none",
                zorder=0,
                clip_on=False,
            )
        )
    axis.plot(
        [band_left, band_right],
        [len(CALLSETS) - 0.38] * 2,
        transform=transform,
        color=NEUTRAL_LIGHT,
        linewidth=0.5,
        clip_on=False,
        zorder=1,
    )

    axis.axvline(1.0, color=RULE, linewidth=0.6, linestyle=(0, (2.5, 2.0)), zorder=1)
    axis.hlines(
        y_positions,
        low,
        high,
        color=NEUTRAL_DARK,
        linewidth=0.7,
        zorder=2,
    )
    for y_position, lower, upper in zip(y_positions, low, high):
        for bound in (lower, upper):
            axis.plot(
                [bound, bound],
                [y_position - 0.16, y_position + 0.16],
                color=NEUTRAL_DARK,
                linewidth=0.7,
                zorder=2,
            )
    axis.scatter(
        effects,
        y_positions,
        s=13,
        marker="o",
        facecolor=POINT["COVID-19"],
        edgecolor="black",
        linewidth=0.4,
        zorder=3,
    )

    axis.set_xscale("log")
    axis.set_xlim(0.85, 12.0)
    axis.set_xticks([1, 2, 3, 5, 10])
    axis.set_xticklabels(["1", "2", "3", "5", "10"])
    axis.minorticks_off()
    axis.set_ylim(-0.65, len(CALLSETS) - 0.35)
    axis.set_yticks(y_positions)
    axis.set_yticklabels([LONG_LABEL[c] for c in CALLSETS])
    axis.tick_params(axis="y", length=0, pad=2.0)
    axis.set_xlabel("Adjusted COVID-19/HC fold change in EPM (95% CI)")
    axis.spines[["top", "right", "left"]].set_visible(False)

    # Numeric read-out columns to the right of the plotting area.
    columns = ((1.06, "Fold change (95% CI)"), (1.50, "Adjusted $P$"))
    for x_position, header in columns:
        axis.text(
            x_position,
            len(CALLSETS) - 0.22,
            header,
            transform=transform,
            fontsize=5.5,
            ha="left",
            va="bottom",
            clip_on=False,
        )
    for y_position, row, effect, lower, upper in zip(y_positions, rows, effects, low, high):
        axis.text(
            columns[0][0],
            y_position,
            f"{effect:.2f} ({lower:.2f}–{upper:.2f})",
            transform=transform,
            fontsize=5.5,
            ha="left",
            va="center",
            clip_on=False,
        )
        axis.text(
            columns[1][0],
            y_position,
            scientific(float(row["p_value_two_sided"])),
            transform=transform,
            fontsize=5.5,
            ha="left",
            va="center",
            clip_on=False,
        )


# --------------------------------------------------------------------------- #
# Panel c: two-caller concordance as a function of breakpoint tolerance
# --------------------------------------------------------------------------- #
def panel_concordance(axis: plt.Axes, concordance: list[dict[str, str]]) -> None:
    # Each series is dodged slightly in x so that the interquartile bars of the
    # two lower series, which overlap in y, stay separable.
    # Each entry carries its own direct-label anchor. The Jaccard label is set inside
    # the panel rather than at the left margin: at the margin it is crossed by the
    # tolerance-2 interquartile bar, whose lower quartile (0.603) reaches below it.
    series = (
        ("circlemap_retention_fraction", "Circle-Map calls retained", NEUTRAL_DARK, "o", -0.07, -0.20, 0.968),
        ("circlefinder_retention_fraction", "Circle_finder calls retained", ACCENT, "s", 0.0, -0.20, 0.762),
        ("jaccard_index", "Jaccard index", NEUTRAL_MID, "^", 0.07, 0.30, 0.639),
    )
    x = np.arange(len(TOLERANCES), dtype=float)

    by_tolerance = {
        tolerance: [row for row in concordance if int(row["tolerance_bp"]) == tolerance]
        for tolerance in TOLERANCES
    }
    for tolerance, rows in by_tolerance.items():
        assert len(rows) == 78, f"tolerance {tolerance}: {len(rows)} samples"
        medians = [
            float(np.median([float(row[entry[0]]) for row in rows])) for entry in series
        ]
        consensus_calls = sum(int(row["consensus_count"]) for row in rows)
        expected = EXPECTED_CONCORDANCE[tolerance]
        for observed, target in zip(medians, expected[:3]):
            assert abs(round(observed, 3) - target) < 1e-9, (
                f"tolerance {tolerance}: median {observed} disagrees with RESULTS_SUMMARY.md"
            )
        assert consensus_calls == expected[3]

    primary = TOLERANCES.index(PRIMARY_TOLERANCE)
    axis.axvline(primary, color=NEUTRAL_LIGHT, linewidth=0.6, linestyle=(0, (2.0, 1.8)), zorder=1)

    for column, label, color, marker, dodge, label_x, label_y in series:
        median = np.asarray(
            [
                float(np.median([float(row[column]) for row in by_tolerance[t]]))
                for t in TOLERANCES
            ]
        )
        q1 = np.asarray(
            [
                float(np.percentile([float(row[column]) for row in by_tolerance[t]], 25))
                for t in TOLERANCES
            ]
        )
        q3 = np.asarray(
            [
                float(np.percentile([float(row[column]) for row in by_tolerance[t]], 75))
                for t in TOLERANCES
            ]
        )
        axis.vlines(x + dodge, q1, q3, color=color, linewidth=0.6, alpha=0.75, zorder=2)
        axis.plot(
            x + dodge,
            median,
            color=color,
            linewidth=0.9,
            marker=marker,
            markersize=2.8,
            markeredgecolor="white",
            markeredgewidth=0.35,
            zorder=3,
        )
        # Direct labels instead of a legend: the three series never cross.
        axis.text(
            label_x,
            label_y,
            label,
            fontsize=5,
            color=color,
            ha="left",
            va="bottom" if marker != "^" else "top",
        )

    axis.set_xlim(-0.28, len(TOLERANCES) - 0.72)
    axis.set_xticks(x)
    axis.set_xticklabels([str(t) for t in TOLERANCES])
    axis.set_xlabel("Breakpoint tolerance (bp)")
    axis.set_ylim(0.56, 1.02)
    axis.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    axis.set_yticklabels(["0.6", "0.7", "0.8", "0.9", "1.0"])
    axis.set_ylabel("Per-sample fraction (median, IQR)")
    axis.spines[["top", "right"]].set_visible(False)
    axis.annotate(
        "prespecified primary tolerance\n9,130,043 consensus calls",
        xy=(primary, 0.828),
        xytext=(primary - 0.12, 0.828),
        fontsize=4.8,
        color=NEUTRAL_MID,
        ha="right",
        va="center",
        linespacing=1.25,
    )


# --------------------------------------------------------------------------- #
# Panel d: fragment-size composition difference by call set
# --------------------------------------------------------------------------- #
def panel_fragment_size(
    axis: plt.Axes, colorbar_axis: plt.Axes, fragment_tests: list[dict[str, str]]
) -> None:
    callsets = [c for c in CALLSETS if c != "circlemap_current"]
    difference = np.full((len(callsets), len(LENGTH_BINS)), np.nan)
    significant = np.zeros_like(difference, dtype=bool)
    for row_index, callset in enumerate(callsets):
        for column_index, length_bin in enumerate(LENGTH_BINS):
            row = next(
                r
                for r in fragment_tests
                if r["callset"] == callset and r["length_bin"] == length_bin
            )
            # Locked proportions are fractions; the panel reports percentage points.
            difference[row_index, column_index] = (
                100.0 * float(row["median_difference_COVID_minus_HC"])
            )
            significant[row_index, column_index] = (
                float(row["q_value_bh_across_callsets_and_bins"]) < 0.05
            )

    assert significant[:, LENGTH_BINS.index("ge2000")].all(), (
        "the >=2 kb depletion is expected to be BH-significant in every call set"
    )
    assert not significant[:, : LENGTH_BINS.index("ge2000")].any(), (
        "no other length bin is expected to be BH-significant"
    )

    # Diverging map built from the locked cohort inks: HC-ward blue, COVID-ward red.
    cmap = LinearSegmentedColormap.from_list(
        "cohort_diverging",
        ["#034E61", "#8BABD3", "#FFFFFF", "#FF8080", "#D70000"],
    )
    limit = 5.5
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    mesh = axis.pcolormesh(
        np.arange(len(LENGTH_BINS) + 1),
        np.arange(len(callsets) + 1),
        difference,
        cmap=cmap,
        norm=norm,
        edgecolors="white",
        linewidth=0.7,
    )

    for row_index in range(len(callsets)):
        for column_index in range(len(LENGTH_BINS)):
            value = difference[row_index, column_index]
            text = minus(f"{value:+.1f}")
            if significant[row_index, column_index]:
                text += "*"
            axis.text(
                column_index + 0.5,
                row_index + 0.5,
                text,
                ha="center",
                va="center",
                fontsize=4.6,
                color="white" if abs(value) > 3.6 else "black",
            )

    axis.set_xlim(0, len(LENGTH_BINS))
    axis.set_ylim(len(callsets), 0)
    axis.set_xticks(np.arange(len(LENGTH_BINS)) + 0.5)
    axis.set_xticklabels([BIN_LABEL[b] for b in LENGTH_BINS], linespacing=1.2)
    axis.set_yticks(np.arange(len(callsets)) + 0.5)
    axis.set_yticklabels([LONG_LABEL[c] for c in callsets])
    axis.set_xlabel("eccDNA length bin (bp)")
    axis.tick_params(length=0, pad=2.0)
    axis.spines[:].set_visible(False)

    colorbar = axis.figure.colorbar(mesh, cax=colorbar_axis)
    # Matplotlib rasterizes colorbar solids by default; the journal requires
    # fully vector art, so the gradient is emitted as filled paths instead.
    colorbar.solids.set_rasterized(False)
    colorbar.solids.set_edgecolor("face")
    colorbar.outline.set_linewidth(0.5)
    colorbar.set_ticks([-5, -2.5, 0, 2.5, 5])
    colorbar.set_ticklabels([minus("-5"), minus("-2.5"), "0", "+2.5", "+5"])
    colorbar_axis.tick_params(labelsize=5, length=1.6, width=0.5, pad=1.5)
    colorbar_axis.set_ylabel(
        "Median difference, COVID-19 − HC\n(percentage points)",
        fontsize=5.5,
        labelpad=3.0,
    )


def main() -> None:
    configure_plotting()
    metrics = read_tsv(TABLE_DIR / "robustness_callset_sample_metrics.tsv")
    tests = read_tsv(TABLE_DIR / "robustness_group_tests.tsv")
    regressions = read_tsv(TABLE_DIR / "robustness_regression_HC3.tsv")
    concordance = read_tsv(TABLE_DIR / "caller_concordance_all.tsv")
    fragment_tests = read_tsv(TABLE_DIR / "robustness_fragment_size_group_tests.tsv")

    figure = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))

    axis_a = add_axes_mm(figure, 15.5, 131.0, 164.0, 59.0)
    axis_b = add_axes_mm(figure, 46.0, 74.0, 74.0, 37.0)
    axis_c = add_axes_mm(figure, 15.5, 14.0, 74.0, 40.0)
    axis_d = add_axes_mm(figure, 122.0, 14.0, 43.0, 40.0)
    axis_cbar = add_axes_mm(figure, 168.0, 14.0, 2.6, 40.0)

    panel_burden(axis_a, metrics, tests)
    panel_adjusted_effect(axis_b, regressions)
    panel_concordance(axis_c, concordance)
    panel_fragment_size(axis_d, axis_cbar, fragment_tests)

    panel_label(figure, 3.0, 195.5, "a")
    panel_label(figure, 3.0, 115.0, "b")
    panel_label(figure, 3.0, 57.0, "c")
    panel_label(figure, 93.0, 57.0, "d")

    OUTPUT_DIR.mkdir(exist_ok=True)
    figure.savefig(OUTPUT_STEM.with_suffix(".pdf"))
    figure.savefig(OUTPUT_STEM.with_suffix(".svg"))
    figure.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=600)
    plt.close(figure)
    print(f"Created {OUTPUT_STEM}.pdf/.svg/.png")


if __name__ == "__main__":
    main()
