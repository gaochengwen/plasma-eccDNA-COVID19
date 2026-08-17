#!/usr/bin/env python3
"""Create Supplementary Figure S2 for the RCA technical-bias revision."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "analyse" / "RCA" / "results"
OUTPUT_DIR = ROOT / "revise"

GROUP_ORDER = ("HC", "COVID-19")
# Panel a mirrors Figure 1b, which plots the patient cohort first.
BOX_ORDER = ("COVID-19", "HC")
# Cohort colours sampled from the vector content of Submit/Figure_1.pdf so that
# this supplementary figure and Figure 1b use exactly the same inks.
COLORS = {"HC": "#034E61", "COVID-19": "#D70000"}
BOX_FILL = {"HC": "#8BABD3", "COVID-19": "#FF8080"}
MARKERS = {"HC": "o", "COVID-19": "^"}
GREY = "#6F7B91"
SEED = 20260723

LABELS = {
    "eccdna_count": "Detected eccDNA count",
    "total_support_events": "Circle-Map support events",
    "mapped_reads": "Mapped reads (millions)",
    "epm": "eccDNA per million mapped reads (EPM)",
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
            "axes.labelsize": 5.5,
            "xtick.labelsize": 5,
            "ytick.labelsize": 5,
            "legend.fontsize": 5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.3,
            "ytick.major.size": 2.3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            # Italic statistical symbols in the same face as the rest of the figure.
            "mathtext.fontset": "custom",
            "mathtext.rm": family,
            "mathtext.it": f"{family}:italic",
            "mathtext.bf": f"{family}:bold",
        }
    )


def stat_label(value: float, symbol: str = "P") -> str:
    """Journal-style label with an italic symbol and a superscript exponent."""
    if value >= 0.01:
        return f"${symbol}$ = {value:.3f}"
    if value >= 0.001:
        # four decimals so the panel matches the precision used in the main text
        return f"${symbol}$ = {value:.4f}"
    mantissa, exponent = f"{value:.2e}".split("e")
    # one mathtext expression, so the multiplication sign and the exponent are spaced
    # by the math layout rather than by the surrounding text
    return rf"${symbol} = {mantissa}\times10^{{{int(exponent)}}}$"


def signed(value: float, decimals: int = 2) -> str:
    """Format with a true minus sign (U+2212) rather than a hyphen."""
    return f"{value:.{decimals}f}".replace("-", "−")


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.18,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
    )


def rca_group_panel(
    axis: plt.Axes,
    samples: list[dict[str, str]],
    comparisons: list[dict[str, str]],
) -> None:
    rng = np.random.default_rng(SEED)
    values = [
        np.asarray(
            [float(row["rca_concentration_ng_ul"]) for row in samples if row["group"] == group]
        )
        for group in BOX_ORDER
    ]
    # Whiskers span the full data range and the box carries a black outline, as in
    # Figure 1b; no axis break is used because the concentrations fit one linear axis.
    boxplot = axis.boxplot(
        values,
        positions=np.arange(2),
        widths=0.38,
        patch_artist=True,
        showfliers=False,
        whis=(0, 100),
        medianprops={"color": "black", "linewidth": 0.7},
        whiskerprops={"color": "black", "linewidth": 0.7},
        capprops={"color": "black", "linewidth": 0.7},
        boxprops={"edgecolor": "black", "linewidth": 0.7},
    )
    for patch, group in zip(boxplot["boxes"], BOX_ORDER):
        patch.set_facecolor(BOX_FILL[group])
    for position, (group, group_values) in enumerate(zip(BOX_ORDER, values)):
        jitter = rng.uniform(-0.13, 0.13, len(group_values))
        axis.scatter(
            position + jitter,
            group_values,
            s=9,
            marker=MARKERS[group],
            color=COLORS[group],
            linewidth=0,
            zorder=3,
        )
    test = next(row for row in comparisons if row["metric"] == "rca_concentration_ng_ul")
    bracket_y = 236.0
    tick = 6.0
    axis.plot(
        [0, 0, 1, 1],
        [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
        color="black",
        linewidth=0.7,
        clip_on=False,
    )
    axis.text(
        0.5,
        bracket_y + 2.0,
        stat_label(float(test["p_value_two_sided"])),
        ha="center",
        va="bottom",
        fontsize=5,
    )
    axis.set_xticks([0, 1], BOX_ORDER)
    axis.set_ylabel("RCA concentration (ng/µL)")
    axis.set_ylim(0, 260)
    axis.spines[["top", "right"]].set_visible(False)


def correlation_result(
    correlations: list[dict[str, str]], scope: str, metric: str
) -> dict[str, str]:
    return next(
        row
        for row in correlations
        if row["scope"] == scope and row["y_metric"] == metric
    )


def correlation_panel(
    axis: plt.Axes,
    samples: list[dict[str, str]],
    correlations: list[dict[str, str]],
    metric: str,
) -> None:
    log_y = metric in {"eccdna_count", "total_support_events", "epm"}
    for group in GROUP_ORDER:
        subset = [row for row in samples if row["group"] == group]
        x = np.asarray([float(row["rca_concentration_ng_ul"]) for row in subset])
        y = np.asarray([float(row[metric]) for row in subset])
        plot_y = y / 1_000_000.0 if metric == "mapped_reads" else y
        axis.scatter(
            x,
            plot_y,
            s=11,
            marker=MARKERS[group],
            color=COLORS[group],
            edgecolor="white",
            linewidth=0.25,
            alpha=0.82,
            zorder=3,
        )
        fit_y = np.log10(plot_y) if log_y else plot_y
        slope, intercept = np.polyfit(x, fit_y, 1)
        fit_x = np.linspace(np.min(x), np.max(x), 100)
        fitted = 10 ** (intercept + slope * fit_x) if log_y else intercept + slope * fit_x
        axis.plot(fit_x, fitted, color=COLORS[group], linewidth=0.8)
    if log_y:
        axis.set_yscale("log")
    axis.set_xlabel("RCA concentration (ng/µL)")
    axis.set_ylabel(LABELS[metric])
    axis.spines[["top", "right"]].set_visible(False)
    lines = []
    for scope in ("All", "HC", "COVID-19"):
        result = correlation_result(correlations, scope, metric)
        lines.append(
            f"{scope}: $\\rho$ = {signed(float(result['spearman_rho']))}, "
            f"{stat_label(float(result['p_value_two_sided']))}"
        )
    axis.text(
        0.03,
        0.97,
        "\n".join(lines),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=4.3,
        linespacing=1.18,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.4},
    )


def fit_hc3(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xtx_inverse = np.linalg.pinv(x.T @ x)
    beta = xtx_inverse @ x.T @ y
    residuals = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, xtx_inverse, x)
    scaled_residual_sq = (residuals / np.maximum(1.0 - leverage, np.finfo(float).eps)) ** 2
    meat = x.T @ (x * scaled_residual_sq[:, None])
    covariance = xtx_inverse @ meat @ xtx_inverse
    return beta, covariance


def adjusted_prediction_panel(
    axis: plt.Axes,
    samples: list[dict[str, str]],
    regressions: list[dict[str, str]],
    outcome: str,
    y_label: str,
) -> None:
    group = np.asarray([row["group"] == "COVID-19" for row in samples], dtype=float)
    rca = np.asarray([float(row["rca_concentration_ng_ul"]) for row in samples])
    log2_rca = np.log2(rca)
    mean_log2_rca = np.mean(log2_rca)
    centered_log2_rca = log2_rca - mean_log2_rca
    age = (np.asarray([float(row["age_years"]) for row in samples]) - 70.0) / 10.0
    male = np.asarray([row["sex"] == "Male" for row in samples], dtype=float)
    outcome_values = np.asarray([float(row[outcome]) for row in samples])
    design = np.column_stack(
        [
            np.ones(len(samples)),
            group,
            centered_log2_rca,
            age,
            male,
            group * centered_log2_rca,
        ]
    )
    beta, covariance = fit_hc3(np.log2(outcome_values + 1.0), design)

    shared_min = max(np.min(rca[group == 0]), np.min(rca[group == 1]))
    shared_max = min(np.max(rca[group == 0]), np.max(rca[group == 1]))
    grid = np.linspace(shared_min, shared_max, 200)
    mean_age = float(np.mean(age))
    male_proportion = float(np.mean(male))
    critical = stats.t.ppf(0.975, len(samples) - design.shape[1])

    for group_name, group_value in (("HC", 0.0), ("COVID-19", 1.0)):
        centered_grid = np.log2(grid) - mean_log2_rca
        new_design = np.column_stack(
            [
                np.ones(len(grid)),
                np.full(len(grid), group_value),
                centered_grid,
                np.full(len(grid), mean_age),
                np.full(len(grid), male_proportion),
                group_value * centered_grid,
            ]
        )
        predicted = new_design @ beta
        standard_error = np.sqrt(np.einsum("ij,jk,ik->i", new_design, covariance, new_design))
        lower = predicted - critical * standard_error
        upper = predicted + critical * standard_error
        axis.plot(
            grid,
            np.maximum(2**predicted - 1.0, np.finfo(float).tiny),
            color=COLORS[group_name],
            linewidth=1.1,
            label=group_name,
        )
        axis.fill_between(
            grid,
            np.maximum(2**lower - 1.0, np.finfo(float).tiny),
            np.maximum(2**upper - 1.0, np.finfo(float).tiny),
            color=COLORS[group_name],
            alpha=0.16,
            linewidth=0,
        )

    interaction = next(
        row
        for row in regressions
        if row["outcome"] == outcome
        and row["model"] == "full_plus_group_RCA_interaction"
        and row["term"] == "COVID-19_x_log2_RCA"
    )
    axis.set_yscale("log")
    axis.set_xlabel("RCA concentration (ng/µL)")
    axis.set_ylabel(y_label)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="upper left", frameon=False, handlelength=1.4)
    axis.text(
        0.03,
        0.72,
        "95% CI for adjusted mean\n"
        f"Group × RCA: {stat_label(float(interaction['p_value_two_sided']))}, "
        f"{stat_label(float(interaction['q_value_bh_within_term']), 'q')}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=4.5,
        linespacing=1.22,
    )


def main() -> None:
    configure_plotting()
    samples = read_tsv(INPUT_DIR / "RCA_sample_metrics.tsv")
    correlations = read_tsv(INPUT_DIR / "RCA_correlations.tsv")
    comparisons = read_tsv(INPUT_DIR / "RCA_group_comparisons.tsv")
    regressions = read_tsv(INPUT_DIR / "RCA_regression_HC3.tsv")

    figure, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 126 / 25.4))
    flat = axes.ravel()
    rca_group_panel(flat[0], samples, comparisons)
    correlation_panel(flat[1], samples, correlations, "eccdna_count")
    correlation_panel(flat[2], samples, correlations, "total_support_events")
    correlation_panel(flat[3], samples, correlations, "epm")
    adjusted_prediction_panel(
        flat[4],
        samples,
        regressions,
        "eccdna_count",
        "Age- and sex-adjusted eccDNA count",
    )
    adjusted_prediction_panel(
        flat[5],
        samples,
        regressions,
        "epm",
        "Age- and sex-adjusted EPM",
    )
    for axis, label in zip(flat, "abcdef"):
        panel_label(axis, label)
    figure.subplots_adjust(
        left=0.085,
        right=0.985,
        bottom=0.105,
        top=0.96,
        wspace=0.46,
        hspace=0.47,
    )

    output_stem = OUTPUT_DIR / "Figure_S2"
    figure.savefig(output_stem.with_suffix(".pdf"))
    figure.savefig(output_stem.with_suffix(".svg"))
    figure.savefig(output_stem.with_suffix(".png"), dpi=300)
    plt.close(figure)

    source_dir = OUTPUT_DIR / "source_data"
    source_dir.mkdir(exist_ok=True)
    for filename in (
        "RCA_sample_metrics.tsv",
        "RCA_correlations.tsv",
        "RCA_group_comparisons.tsv",
        "RCA_regression_HC3.tsv",
    ):
        destination = source_dir / filename
        destination.write_bytes((INPUT_DIR / filename).read_bytes())
    print(f"Created {output_stem}.pdf/.svg/.png")


if __name__ == "__main__":
    main()
