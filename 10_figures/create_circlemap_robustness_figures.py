#!/usr/bin/env python3
"""Create supplementary figures for Circle-Map robustness analyses.

The script reads only finalized tables from analyse/Circle_finder and writes
publication-ready Figure S6/S7 outputs to revise/.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


REVISE = Path(__file__).resolve().parent.parent
ROOT = REVISE.parent / "analyse" / "Circle_finder"
TABLES = ROOT / "tables"

# Figure S7d,e report only the three targets selected from the computational
# results for experimental validation. The v1 fidelity run rebinds this list
# to the historical four-target set; the v2 output uses the three-target set.
VALIDATION_TARGETS = ["CTNNA2", "SDK1", "TCF7L1"]

CALLSET_ORDER = [
    "circlemap_current",
    "circlemap_methods",
    "circlemap_strict",
    "circlemap_very_strict",
    "circlemap_artifact_masked",
    "circlemap_high_support_masked",
    "consensus_t10",
]

CALLSET_LABELS = {
    "circlemap_current": "Current\nCM",
    "circlemap_methods": "Methods\nCM",
    "circlemap_strict": "CM split>=10\nscore>1,000",
    "circlemap_very_strict": "CM split>=20\nscore>2,000",
    "circlemap_artifact_masked": "CM artifact\nmasked",
    "circlemap_high_support_masked": "CM high-support\n+ masked",
    "consensus_t10": "CM/CF\nconsensus",
}

CALLSET_FULL_LABELS = {
    "circlemap_current": "Current Circle-Map",
    "circlemap_methods": "Methods-consistent Circle-Map",
    "circlemap_strict": "Circle-Map split>=10, score>1,000",
    "circlemap_very_strict": "Circle-Map split>=20, score>2,000",
    "circlemap_artifact_masked": "Circle-Map artifact-masked",
    "circlemap_high_support_masked": "Circle-Map high-support + artifact-masked",
    "consensus_t10": "Circle-Map/Circle_finder consensus (10 bp)",
}

LENGTH_ORDER = [
    "lt200",
    "200_399",
    "400_599",
    "600_999",
    "1000_1999",
    "ge2000",
]

LENGTH_LABELS = {
    "lt200": "<200",
    "200_399": "200-399",
    "400_599": "400-599",
    "600_999": "600-999",
    "1000_1999": "1,000-1,999",
    "ge2000": ">=2,000",
}

PALETTE = {
    "hc": "#0272B2",
    "covid": "#C93E3F",
    "orange": "#EC6F00",
    "green": "#459434",
    "purple": "#A84E94",
    "teal": "#019AA3",
    "yellow": "#CCA02C",
    "grey1": "#E6E6ED",
    "grey2": "#C8CEDA",
    "grey3": "#99A3B4",
    "grey4": "#6F7B91",
    "grey5": "#49566D",
    "grey6": "#253247",
}

FIGURE4B_GROUP_STYLE = {
    "COVID-19": {"fill": "#FF8080", "point": "#D70000"},
    "HC": {"fill": "#8BABD3", "point": "#034E61"},
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "xtick.labelsize": 5,
            "ytick.labelsize": 5,
            "legend.fontsize": 5,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
        }
    )


def read_table(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def fmt_p(p: float) -> str:
    if not math.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.4f}".rstrip("0").rstrip(".")


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
    )


def strip_boxplot(ax: plt.Axes, sample_metrics: pd.DataFrame) -> None:
    rng = np.random.default_rng(1729)
    width = 0.28
    offsets = {"COVID-19": -width / 1.5, "HC": width / 1.5}

    for i, callset in enumerate(CALLSET_ORDER):
        for group in ["COVID-19", "HC"]:
            values = sample_metrics.loc[
                (sample_metrics["callset"] == callset)
                & (sample_metrics["group"] == group),
                "epm",
            ].to_numpy(float)
            xpos = i + offsets[group]
            style = FIGURE4B_GROUP_STYLE[group]
            ax.boxplot(
                np.log10(values + 1),
                positions=[xpos],
                widths=width,
                patch_artist=True,
                manage_ticks=False,
                showfliers=False,
                medianprops={"color": "#000000", "linewidth": 0.60},
                boxprops={
                    "facecolor": style["fill"],
                    "edgecolor": PALETTE["grey6"],
                    "linewidth": 0.45,
                },
                whiskerprops={"color": PALETTE["grey6"], "linewidth": 0.45},
                capprops={"color": PALETTE["grey6"], "linewidth": 0.45},
            )
            jitter = rng.normal(0, 0.025, len(values))
            ax.scatter(
                np.full_like(values, xpos, dtype=float) + jitter,
                np.log10(values + 1),
                s=5,
                facecolors=style["point"],
                edgecolors=style["point"],
                linewidths=0.20,
                alpha=0.75,
                zorder=3,
            )

    ax.set_xlim(-0.55, len(CALLSET_ORDER) - 0.45)
    ax.set_xticks(range(len(CALLSET_ORDER)))
    ax.set_xticklabels([CALLSET_LABELS[c] for c in CALLSET_ORDER], rotation=35, ha="right")
    ax.set_ylabel("log10(EPM + 1)")
    ax.set_title("Sample-level eccDNA burden")
    ax.grid(axis="y", color=PALETTE["grey2"], linewidth=0.35, alpha=0.7)
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=4,
                markerfacecolor=FIGURE4B_GROUP_STYLE["COVID-19"]["point"],
                markeredgecolor=PALETTE["grey6"],
                color="none",
                label="COVID-19",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markersize=4,
                markerfacecolor=FIGURE4B_GROUP_STYLE["HC"]["point"],
                markeredgecolor=PALETTE["grey6"],
                color="none",
                label="HC",
            ),
        ],
        frameon=False,
        loc="upper right",
        handletextpad=0.3,
    )
    panel_label(ax, "a")


def forest_plot(ax: plt.Axes, regressions: pd.DataFrame) -> None:
    rows = (
        regressions.loc[
            (regressions["term"] == "COVID-19_vs_HC")
            & (regressions["callset"].isin(CALLSET_ORDER))
        ]
        .set_index("callset")
        .loc[CALLSET_ORDER]
        .reset_index()
    )
    y = np.arange(len(rows))[::-1]
    lo = rows["multiplicative_ci_low"].to_numpy(float)
    hi = rows["multiplicative_ci_high"].to_numpy(float)
    eff = rows["multiplicative_effect"].to_numpy(float)
    ax.hlines(y, lo, hi, color=PALETTE["grey5"], linewidth=0.8)
    ax.scatter(eff, y, s=26, color=PALETTE["covid"], edgecolor="black", linewidth=0.35, zorder=3)
    ax.axvline(1, color="black", linestyle=":", linewidth=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.8, 12)
    ax.set_yticks(y)
    ax.set_yticklabels([CALLSET_FULL_LABELS[c] for c in rows["callset"]])
    ax.set_xlabel("Adjusted COVID-19/HC fold change in EPM")
    ax.set_title("Adjusted burden effect")
    ax.grid(axis="x", color=PALETTE["grey2"], linewidth=0.35, alpha=0.7)
    for yi, p in zip(y, rows["p_value_two_sided"].to_numpy(float)):
        ax.text(11.7, yi, f"P={fmt_p(p)}", va="center", ha="right", fontsize=5)
    panel_label(ax, "b", x=-0.22)


def concordance_plot(ax: plt.Axes, concordance: pd.DataFrame) -> None:
    summary = (
        concordance.groupby("tolerance_bp", as_index=False)
        .agg(
            circlemap_retention=("circlemap_retention_fraction", "median"),
            circlefinder_retention=("circlefinder_retention_fraction", "median"),
            jaccard=("jaccard_index", "median"),
            consensus_calls=("consensus_count", "sum"),
        )
        .sort_values("tolerance_bp")
    )
    series = [
        ("Circle-Map retained", "circlemap_retention", PALETTE["covid"], "o"),
        ("Circle_finder retained", "circlefinder_retention", PALETTE["orange"], "s"),
        ("Jaccard", "jaccard", PALETTE["hc"], "^"),
    ]
    for label, column, color, marker in series:
        ax.plot(
            summary["tolerance_bp"],
            summary[column],
            marker=marker,
            markersize=3.5,
            linewidth=0.9,
            color=color,
            label=label,
        )
    ax.set_ylim(0.55, 1.00)
    ax.set_xticks(summary["tolerance_bp"])
    ax.set_xlabel("Breakpoint tolerance (bp)")
    ax.set_ylabel("Median fraction")
    ax.set_title("Two-caller concordance")
    ax.grid(axis="y", color=PALETTE["grey2"], linewidth=0.35, alpha=0.7)
    ax.legend(frameon=False, loc="lower right", handlelength=1.5)
    primary = summary.loc[summary["tolerance_bp"] == 10].iloc[0]
    ax.text(
        0.03,
        0.08,
        f"10-bp consensus: {primary['consensus_calls'] / 1e6:.2f}M calls",
        transform=ax.transAxes,
        fontsize=5,
        ha="left",
        va="bottom",
    )
    panel_label(ax, "c")


def fragment_heatmap(ax: plt.Axes, fragment_tests: pd.DataFrame) -> None:
    callsets = [callset for callset in CALLSET_ORDER if callset in set(fragment_tests["callset"])]
    table = (
        fragment_tests.pivot(index="callset", columns="length_bin", values="median_difference_COVID_minus_HC")
        .loc[callsets, LENGTH_ORDER]
        .astype(float)
    )
    q = (
        fragment_tests.pivot(index="callset", columns="length_bin", values="q_value_bh_across_callsets_and_bins")
        .loc[callsets, LENGTH_ORDER]
        .astype(float)
    )
    vmax = float(np.nanmax(np.abs(table.to_numpy())))
    cmap = LinearSegmentedColormap.from_list(
        "blue_white_red", [PALETTE["hc"], "#FFFFFF", PALETTE["covid"]]
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for r in range(table.shape[0]):
        for c in range(table.shape[1]):
            value = table.iat[r, c]
            ax.add_patch(
                Rectangle(
                    (c - 0.5, r - 0.5),
                    1.0,
                    1.0,
                    facecolor=cmap(norm(value)),
                    edgecolor="white",
                    linewidth=0.45,
                )
            )
    ax.set_xlim(-0.5, table.shape[1] - 0.5)
    ax.set_ylim(table.shape[0] - 0.5, -0.5)
    ax.set_yticks(range(len(callsets)))
    ax.set_yticklabels([CALLSET_LABELS[c].replace("\n", " ") for c in callsets])
    ax.set_xticks(range(len(LENGTH_ORDER)))
    ax.set_xticklabels([LENGTH_LABELS[b] for b in LENGTH_ORDER], rotation=35, ha="right")
    ax.set_xlabel("eccDNA length bin (bp)")
    ax.set_title("Fragment-size composition")
    for r in range(table.shape[0]):
        for c in range(table.shape[1]):
            value = table.iat[r, c]
            star = "*" if q.iat[r, c] < 0.05 else ""
            text_color = "white" if abs(value) > vmax * 0.55 else "black"
            ax.text(c, r, f"{value:+.2f}{star}", ha="center", va="center", fontsize=4.5, color=text_color)
    panel_label(ax, "d", x=-0.18)


def save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix in ("pdf", "svg", "png"):
        path = REVISE / f"{stem}.{suffix}"
        if suffix == "png":
            fig.savefig(path, dpi=600)
        else:
            fig.savefig(path)


def create_figure_s6() -> None:
    sample_metrics = read_table("robustness_callset_sample_metrics.tsv")
    regressions = read_table("robustness_regression_HC3.tsv")
    concordance = read_table("caller_concordance_all.tsv")
    fragment_tests = read_table("robustness_fragment_size_group_tests.tsv")

    fig = plt.figure(figsize=(183 / 25.4, 168 / 25.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.28, 1.0], height_ratios=[1.06, 1.0])
    strip_boxplot(fig.add_subplot(gs[0, 0]), sample_metrics)
    forest_plot(fig.add_subplot(gs[0, 1]), regressions)
    concordance_plot(fig.add_subplot(gs[1, 0]), concordance)
    fragment_heatmap(fig.add_subplot(gs[1, 1]), fragment_tests)
    save_figure(fig, "Figure_S6")
    plt.close(fig)


def clip_for_display(values: pd.Series, limit: float = 6.0) -> pd.Series:
    return values.astype(float).clip(-limit, limit)


def eccgene_scatter(ax: plt.Axes) -> None:
    methods = pd.read_csv(
        ROOT / "downstream" / "circlemap_methods" / "eccGene" / "results" / "junction_abundance_wilcoxon.tsv",
        sep="\t",
    )
    consensus = pd.read_csv(
        ROOT / "downstream" / "consensus_t10" / "eccGene" / "results" / "junction_abundance_wilcoxon.tsv",
        sep="\t",
    )
    merged = methods.merge(consensus, on="Gene", suffixes=("_methods", "_consensus"))
    both_sig = (merged["wilcoxon_q_BH_methods"] < 0.05) & (merged["wilcoxon_q_BH_consensus"] < 0.05)
    x = clip_for_display(merged["log2FC_mean_EA_COVID_vs_HC_methods"])
    y = clip_for_display(merged["log2FC_mean_EA_COVID_vs_HC_consensus"])
    ax.scatter(x[~both_sig], y[~both_sig], s=3, color=PALETTE["grey3"], alpha=0.28, linewidths=0)
    ax.scatter(x[both_sig], y[both_sig], s=4, color=PALETTE["covid"], alpha=0.32, linewidths=0)
    ax.axhline(0, color=PALETTE["grey4"], linewidth=0.45)
    ax.axvline(0, color=PALETTE["grey4"], linewidth=0.45)
    ax.plot([-6, 6], [-6, 6], color="black", linewidth=0.6, linestyle=":")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_xlabel("Methods Circle-Map eccGene log2FC")
    ax.set_ylabel("Consensus eccGene log2FC")
    ax.set_title("Genome-wide eccGene effect concordance")
    rho = read_table("eccgene_callset_concordance.tsv")
    primary = rho[(rho["comparison_callset"] == "consensus_t10") & (rho["definition"] == "junction")].iloc[0]
    ax.text(
        0.03,
        0.95,
        f"junction definition\nrho={primary['log2FC_spearman_rho']:.3f}; Jaccard={primary['significant_gene_jaccard']:.3f}\naxes clipped at +/-6",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5,
    )
    panel_label(ax, "a")


def eccgene_concordance_summary(ax: plt.Axes) -> None:
    data = read_table("eccgene_callset_concordance.tsv")
    data = data[data["comparison_callset"] == "consensus_t10"].copy()
    data["definition"] = pd.Categorical(data["definition"], ["junction", "interval", "midpoint"], ordered=True)
    data = data.sort_values("definition")
    x = np.arange(len(data))
    width = 0.28
    ax.bar(
        x - width / 2,
        data["log2FC_spearman_rho"],
        width=width,
        color=PALETTE["hc"],
        edgecolor="black",
        linewidth=0.35,
        label="log2FC rho",
    )
    ax.bar(
        x + width / 2,
        data["significant_gene_jaccard"],
        width=width,
        color=PALETTE["orange"],
        edgecolor="black",
        linewidth=0.35,
        label="significant-gene Jaccard",
    )
    ax.set_ylim(0, 1.20)
    ax.set_xticks(x)
    ax.set_xticklabels(["Junction", "Interval", "Midpoint"], rotation=25, ha="right")
    ax.set_ylabel("Concordance")
    ax.set_title("eccGene robustness by definition")
    ax.legend(frameon=False, loc="upper left", handletextpad=0.4)
    ax.grid(axis="y", color=PALETTE["grey2"], linewidth=0.35, alpha=0.7)
    panel_label(ax, "b")


def chromatin_scatter(ax: plt.Axes) -> None:
    data = read_table("chromatin_callset_concordance.tsv")
    data = data[data["comparison_callset"] == "consensus_t10"].copy()
    color_map = {"absolute": PALETTE["covid"], "relative": PALETTE["hc"]}
    marker_map = {"absolute": "o", "relative": "^"}
    for analysis, sub in data.groupby("analysis"):
        ax.scatter(
            sub["reference_median_difference"],
            sub["comparison_median_difference"],
            s=22,
            color=color_map[analysis],
            marker=marker_map[analysis],
            alpha=0.75,
            edgecolor="black",
            linewidth=0.25,
            label=analysis.capitalize(),
        )
    lim = max(
        abs(data["reference_median_difference"]).max(),
        abs(data["comparison_median_difference"]).max(),
    )
    lim = max(0.5, float(lim) * 1.08)
    ax.plot([-lim, lim], [-lim, lim], color="black", linewidth=0.6, linestyle=":")
    ax.axhline(0, color=PALETTE["grey4"], linewidth=0.45)
    ax.axvline(0, color=PALETTE["grey4"], linewidth=0.45)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Methods Circle-Map median difference")
    ax.set_ylabel("Consensus median difference")
    ax.set_title("Chromatin-effect concordance")
    direction = data["direction_concordant"].mean() * 100
    significance = data["significance_concordant"].mean() * 100
    ax.text(
        0.03,
        0.95,
        f"n={len(data)} tests\ndirection {direction:.1f}%\nsignificance {significance:.1f}%",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5,
    )
    ax.legend(frameon=False, loc="lower right", handletextpad=0.3)
    panel_label(ax, "c", y=1.13)


def target_support(ax: plt.Axes) -> None:
    data = read_table("validation_target_support_summary.tsv")
    selected_callsets = [
        "circlemap_methods",
        "circlemap_strict",
        "circlemap_very_strict",
        "circlemap_high_support_masked",
        "consensus_t10",
    ]
    selected_targets = list(VALIDATION_TARGETS)
    covid = data[
        (data["group"] == "COVID-19")
        & (data["tolerance_bp"] == 10)
        & (data["callset"].isin(selected_callsets))
        & (data["target"].isin(selected_targets))
    ].copy()
    matrix = (
        covid.pivot(index="target", columns="callset", values="supported_samples")
        .loc[selected_targets, selected_callsets]
        .fillna(0)
    )
    ax.set_xlim(-0.5, len(selected_callsets) - 0.5)
    ax.set_ylim(-0.5, len(selected_targets) - 0.5)
    for y, target in enumerate(selected_targets):
        for x, callset in enumerate(selected_callsets):
            value = int(matrix.loc[target, callset])
            size = 20 + (value / 39) * 260
            ax.scatter(x, y, s=size, color=PALETTE["teal"], alpha=0.82, edgecolor="black", linewidth=0.35)
            ax.text(x, y, str(value), ha="center", va="center", fontsize=5, color="black")
    ax.set_xticks(range(len(selected_callsets)))
    ax.set_xticklabels([CALLSET_LABELS[c].replace("\n", " ") for c in selected_callsets], rotation=35, ha="right")
    ax.set_yticks(range(len(selected_targets)))
    ax.set_yticklabels(["chr22 target" if t == "chr22_target" else t for t in selected_targets])
    ax.invert_yaxis()
    ax.set_title("Validation-target support", pad=8)
    ax.set_xlabel("Callset; count shown out of 39 samples")
    ax.set_facecolor("#FFFFFF")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    for x in range(len(selected_callsets) + 1):
        ax.axvline(x - 0.5, color=PALETTE["grey1"], linewidth=0.5, zorder=0)
    for y in range(len(selected_targets) + 1):
        ax.axhline(y - 0.5, color=PALETTE["grey1"], linewidth=0.5, zorder=0)
    panel_label(ax, "d", x=-0.24, y=1.17)


def target_mask_audit(ax: plt.Axes) -> None:
    data = read_table("validation_target_mask_audit.tsv")
    targets = list(VALIDATION_TARGETS)
    columns = ["Umap >=0.90", "No blacklist/gap", "No segmental dup.", "No RepeatMasker", "Strict target"]
    status = {}
    for target, sub in data.groupby("target"):
        status[target] = [
            bool((sub["umap_k100_unique_overlap_fraction"].astype(float) >= 0.90).all()),
            bool((sub["blacklist_gap_overlap_bp"].astype(float) == 0).all()),
            bool((sub["segmental_duplication_overlap_bp"].astype(float) == 0).all()),
            bool((sub["repeatmasker_overlap_bp"].astype(float) == 0).all()),
        ]
        status[target].append(all(status[target]))

    ax.set_xlim(-0.5, len(columns) - 0.5)
    ax.set_ylim(-0.5, len(targets) - 0.5)
    for y, target in enumerate(targets):
        for x, ok in enumerate(status[target]):
            color = PALETTE["green"] if ok else PALETTE["covid"]
            ax.add_patch(Rectangle((x - 0.45, y - 0.36), 0.9, 0.72, facecolor=color, edgecolor="black", linewidth=0.35))
            ax.text(x, y, "Pass" if ok else "Fail", ha="center", va="center", fontsize=5, color="white")
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(columns, rotation=35, ha="right")
    ax.set_yticks(range(len(targets)))
    ax.set_yticklabels(["chr22 target" if t == "chr22_target" else t for t in targets])
    ax.invert_yaxis()
    ax.set_title("Breakpoint artefact-mask audit", pad=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "e", x=-0.22, y=1.17)


def create_figure_s7() -> None:
    fig = plt.figure(figsize=(183 / 25.4, 170 / 25.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.0, 1.22], height_ratios=[1.0, 1.0])
    eccgene_scatter(fig.add_subplot(gs[0, :2]))
    eccgene_concordance_summary(fig.add_subplot(gs[0, 2]))
    chromatin_scatter(fig.add_subplot(gs[1, 0]))
    target_support(fig.add_subplot(gs[1, 1]))
    target_mask_audit(fig.add_subplot(gs[1, 2]))
    save_figure(fig, "Figure_S7")
    plt.close(fig)


def write_captions() -> None:
    (REVISE / "Figure_S6_caption.md").write_text(
        """Supplementary Figure 6. Robustness of the COVID-19 versus HC eccDNA-burden conclusion to caller and call-set definitions. (a) Sample-level EPM distributions for seven call sets: the current Circle-Map call set, the manuscript methods-consistent Circle-Map call set, two progressively stricter Circle-Map split-read/score thresholds, a stringent artifact-masked Circle-Map call set, a high-support artifact-masked call set, and a Circle-Map/Circle_finder two-caller consensus call set requiring both breakpoints to match within 10 bp. Points are samples (n = 39 per group), boxes show the median and interquartile range and whiskers extend to 1.5 times the interquartile range. (b) Multiplicative COVID-19 versus HC effects from robust linear models of log2(EPM + 1), adjusted for log2 RCA product concentration, age and sex, with HC3 standard errors; points show adjusted fold changes and bars show 95% confidence intervals. (c) One-to-one Circle-Map/Circle_finder concordance across breakpoint tolerances. Lines show the sample-level median Circle-Map retention fraction, Circle_finder retention fraction and Jaccard index; the 10-bp consensus was used for downstream analyses. (d) Fragment-size sensitivity analysis. Heat-map values are COVID-19 minus HC differences in median per-sample length-bin proportions; asterisks denote BH q < 0.05 across call set and length-bin tests.""",
        encoding="utf-8",
    )
    (REVISE / "Figure_S7_caption.md").write_text(
        """Supplementary Figure 7. Robustness of downstream eccGene, chromatin and validation-target analyses to high-support filtering and two-caller consensus calling. (a) Genome-wide concordance of junction-defined eccGene log2 fold changes between the methods-consistent Circle-Map call set and the 10-bp Circle-Map/Circle_finder consensus call set. Points are genes tested in both call sets; red points were BH significant in both call sets. Axes are clipped at +/-6 for display only. (b) Methods-Circle-Map versus consensus concordance across the three eccGene assignment definitions, summarized by Spearman correlation of log2 fold changes and Jaccard overlap of BH-significant gene sets. (c) Concordance of chromatin group effects between methods Circle-Map and the 10-bp consensus call set across 60 matched absolute and relative enrichment tests. (d) COVID-19 support for the three candidates selected for experimental validation across call sets at a 10-bp breakpoint tolerance; numbers are supported samples out of 39. No HC sample supported any of the three candidates. (e) Strict breakpoint artifact-mask audit for the three candidates. CTNNA2 and TCF7L1 passed all mask components, retaining support in 14/39 and 8/39 COVID-19 samples after masking. SDK1 had one repeat-overlapping breakpoint window but full Umap k = 100 mappability at both breakpoints.""",
        encoding="utf-8",
    )


def write_palette_and_qa(stems: list[str] | None = None) -> None:
    s6_palette_rows = [
        ("S6a COVID-19 box fill", FIGURE4B_GROUP_STYLE["COVID-19"]["fill"], "Matched to Figure 4b"),
        ("S6a HC box fill", FIGURE4B_GROUP_STYLE["HC"]["fill"], "Matched to Figure 4b"),
        ("S6a COVID-19 raw points", FIGURE4B_GROUP_STYLE["COVID-19"]["point"], "Matched to Figure 4b"),
        ("S6a HC raw points", FIGURE4B_GROUP_STYLE["HC"]["point"], "Matched to Figure 4b"),
        ("S6b/S6d COVID-19 or positive effect", PALETTE["covid"], "Robustness effect and diverging heat-map positive side"),
        ("Circle_finder retained", PALETTE["orange"], "Caller concordance"),
        ("Neutral grey", PALETTE["grey2"], "Grid/context"),
    ]
    s7_palette_rows = [
        ("HC", PALETTE["hc"], "Group colour"),
        ("COVID-19", PALETTE["covid"], "Group colour and positive/failed mask"),
        ("Circle_finder retained", PALETTE["orange"], "Caller concordance"),
        ("Pass", PALETTE["green"], "Artifact-mask pass"),
        ("Target support", PALETTE["teal"], "Target-support points"),
        ("Neutral grey", PALETTE["grey2"], "Grid/context"),
    ]
    entries = {
        "Figure_S6": (
            s6_palette_rows,
            "Colours: panel a uses the exact group-colour format extracted from Figure 4b; the remaining panels use the Nature palette mappings recorded in the palette file.",
        ),
        "Figure_S7": (
            s7_palette_rows,
            "Colours: exact Hex values from the bundled Nature palette; colour is paired with labels, shapes, or direct annotations.",
        ),
    }
    selected = stems if stems is not None else list(entries)
    for stem in selected:
        palette_rows, colour_note = entries[stem]
        (REVISE / f"{stem}_palette.tsv").write_text(
            "meaning\thex\trole\n"
            + "\n".join("\t".join(row) for row in palette_rows)
            + "\n",
            encoding="utf-8",
        )
        (REVISE / f"{stem}_style_QA.md").write_text(
            f"""# {stem} style QA

- Width: 183 mm double-column.
- Font: Arial requested through matplotlib; PDF export uses TrueType font embedding.
- Text: live text in PDF/SVG; no labels converted to paths.
- Line widths: axes and plotted lines >=0.35 pt.
- {colour_note}
- Outputs: PDF and SVG for vector editing; 600 dpi PNG for quick review.
""",
            encoding="utf-8",
        )


def main() -> None:
    configure_matplotlib()
    create_figure_s6()
    create_figure_s7()
    write_captions()
    write_palette_and_qa()


if __name__ == "__main__":
    main()
