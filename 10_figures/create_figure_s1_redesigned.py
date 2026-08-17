#!/usr/bin/env python3
"""Create Supplementary Figure S1, the age-matched sensitivity analysis.

The pairing is not derived here.  It is read from age_matched_pairs.json, the
single canonical assignment produced by build_age_matched_pairs.py under the
1-year age caliper described in the Methods, so that this figure and
Supplementary Table S19 can never diverge.  This script only visualises the
sample-level eccDNA burden and length densities of those pairs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np
import pandas as pd
import json

from scipy.stats import mannwhitneyu, wilcoxon


ROOT = Path(__file__).resolve().parents[2]
METADATA_PATH = ROOT / "analyse" / "EPM" / "input" / "sample_metadata.tsv"
BURDEN_PATH = (
    ROOT
    / "analyse"
    / "EPM"
    / "results"
    / "sample_burden_recomputed_from_qduh_bed.tsv"
)
LENGTH_CURVES_PATH = (
    ROOT
    / "analyse"
    / "fragment-length"
    / "source_data"
    / "sample_level_primary_curves.tsv"
)
PEAK_PATH = (
    ROOT
    / "analyse"
    / "fragment-length"
    / "tables"
    / "objective_peak_identification_and_stability.tsv"
)
OUTPUT_BASE = ROOT / "revise" / "Figure_S1_redesigned"
SOURCE_DATA_DIR = ROOT / "revise" / "source_data"

GROUP_COVID = "COVID-19"
GROUP_HC = "HC"
PAIRS_PATH = ROOT / "revise" / "source_data" / "age_matched_pairs.json"
AGE_CALIPER = 1
EXPECT_STABLE_PEAKS = [196, 365, 571]
BOOTSTRAP_ITERATIONS = 5000
RANDOM_SEED = 20260731

# Nature scientific-illustration palette.
COVID_COLOR = "#C93E3F"
COVID_LIGHT = "#FAD0CE"
HC_COLOR = "#0272B2"
HC_LIGHT = "#C8E7FB"
GREY_DARK = "#6F7B91"
GREY_MID = "#99A3B4"
GREY_LIGHT = "#E6E6ED"
BLACK = "#000000"
PEAK_LABEL_COLOR = "#253247"   # matches the peak labels in Figure 1d


def configure_matplotlib() -> None:
    """Set final-size, editable publication defaults."""
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
            "font.size": 6.0,
            "axes.labelsize": 6.5,
            "axes.linewidth": 0.6,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.unicode_minus": False,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.minor.size": 1.5,
            "ytick.minor.size": 1.5,
            "legend.fontsize": 5.5,
            "legend.frameon": False,
            "lines.linewidth": 0.9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_sample_data() -> pd.DataFrame:
    """Load and validate sample metadata and recomputed eccDNA burden."""
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    burden = pd.read_csv(BURDEN_PATH, sep="\t")

    required_metadata = {"sample_id", "gender", "age", "group"}
    required_burden = {
        "sample_id",
        "group",
        "eccdna_count_qduh_bed",
        "epm_qduh_bed",
    }
    if not required_metadata.issubset(metadata.columns):
        raise ValueError(
            f"Missing metadata columns: {required_metadata - set(metadata.columns)}"
        )
    if not required_burden.issubset(burden.columns):
        raise ValueError(
            f"Missing burden columns: {required_burden - set(burden.columns)}"
        )

    data = metadata[list(required_metadata)].merge(
        burden[list(required_burden)],
        on=["sample_id", "group"],
        how="inner",
        validate="one_to_one",
    )
    data = data.rename(
        columns={
            "eccdna_count_qduh_bed": "eccdna_count",
            "epm_qduh_bed": "epm",
            "gender": "sex",
        }
    )
    if len(data) != 78 or data["sample_id"].nunique() != 78:
        raise ValueError(f"Expected 78 unique samples, found {len(data)}")
    if set(data["group"]) != {GROUP_COVID, GROUP_HC}:
        raise ValueError(f"Unexpected groups: {sorted(data['group'].unique())}")
    if (data[["age", "eccdna_count", "epm"]] <= 0).any().any():
        raise ValueError("Age, eccDNA count and EPM must all be positive")
    return data


def match_samples(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the canonical age-matched pairing built by build_age_matched_pairs.py.

    That script matches each COVID-19 patient 1:1 and without replacement to a
    healthy control whose age differs by at most one year, maximising the number
    of pairs, then minimising the total absolute age difference, then maximising
    sex concordance.  Reading it here rather than re-deriving it keeps the figure
    and Supplementary Table S19 on one source of truth.
    """
    records = json.loads(PAIRS_PATH.read_text())
    pairs = pd.DataFrame(
        {
            "pair_id": np.arange(1, len(records) + 1),
            "covid_sample_id": [r[0] for r in records],
            "covid_sex": [r[2] for r in records],
            "covid_age": [r[1] for r in records],
            "hc_sample_id": [r[5] for r in records],
            "hc_sex": [r[7] for r in records],
            "hc_age": [r[6] for r in records],
        }
    )
    pairs["age_difference_covid_minus_hc"] = pairs["covid_age"] - pairs["hc_age"]
    pairs["absolute_age_difference"] = pairs["age_difference_covid_minus_hc"].abs()
    if pairs["absolute_age_difference"].max() > AGE_CALIPER:
        raise ValueError("A pair exceeds the 1-year age caliper")

    idx = data.set_index("sample_id")
    covid = idx.loc[pairs["covid_sample_id"]].reset_index()
    hc = idx.loc[pairs["hc_sample_id"]].reset_index()
    if (covid["group"] != GROUP_COVID).any() or (hc["group"] != GROUP_HC).any():
        raise ValueError("Pair file does not respect the cohort labels")
    if (covid["age"].to_numpy() != pairs["covid_age"].to_numpy()).any() or (
        hc["age"].to_numpy() != pairs["hc_age"].to_numpy()
    ).any():
        raise ValueError("Ages in the pair file disagree with Supplementary Table S1")

    matched = pd.concat(
        [
            covid.assign(pair_id=pairs["pair_id"].to_numpy(), matched_group_order=0),
            hc.assign(pair_id=pairs["pair_id"].to_numpy(), matched_group_order=1),
        ],
        ignore_index=True,
    ).sort_values(["pair_id", "matched_group_order"])
    return pairs, matched


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Return Cliff's delta, positive when x tends to exceed y."""
    u_value = mannwhitneyu(x, y, alternative="two-sided").statistic
    return float(2.0 * u_value / (len(x) * len(y)) - 1.0)


def build_statistics(
    data: pd.DataFrame, pairs: pd.DataFrame, matched: pd.DataFrame
) -> pd.DataFrame:
    """Create the auditable statistics table used by panels a--c."""
    rows: list[dict[str, object]] = []

    for label, subset in (
        ("full_cohort", data),
        ("matched_subset", matched),
    ):
        covid = subset.loc[subset["group"] == GROUP_COVID]
        hc = subset.loc[subset["group"] == GROUP_HC]
        age_test = mannwhitneyu(
            covid["age"].to_numpy(),
            hc["age"].to_numpy(),
            alternative="two-sided",
        )
        rows.append(
            {
                "analysis_set": label,
                "metric": "age",
                "test": "two-sided Mann-Whitney U",
                "n_covid": len(covid),
                "n_hc": len(hc),
                "covid_median": float(covid["age"].median()),
                "hc_median": float(hc["age"].median()),
                "effect": cliffs_delta(
                    covid["age"].to_numpy(), hc["age"].to_numpy()
                ),
                "effect_definition": "Cliff's delta (COVID-19 minus HC)",
                "p_value": float(age_test.pvalue),
            }
        )

    covid_by_pair = (
        matched.loc[matched["group"] == GROUP_COVID]
        .set_index("pair_id")
        .sort_index()
    )
    hc_by_pair = (
        matched.loc[matched["group"] == GROUP_HC]
        .set_index("pair_id")
        .sort_index()
    )
    for metric in ("eccdna_count", "epm"):
        covid_values = covid_by_pair[metric].to_numpy(dtype=float)
        hc_values = hc_by_pair[metric].to_numpy(dtype=float)
        test = wilcoxon(
            covid_values,
            hc_values,
            alternative="two-sided",
            method="auto",
        )
        rows.append(
            {
                "analysis_set": "matched_subset",
                "metric": metric,
                "test": "two-sided paired Wilcoxon signed-rank",
                "n_covid": len(covid_values),
                "n_hc": len(hc_values),
                "covid_median": float(np.median(covid_values)),
                "hc_median": float(np.median(hc_values)),
                "effect": float(np.median(covid_values / hc_values)),
                "effect_definition": "median within-pair fold change",
                "p_value": float(test.pvalue),
            }
        )

    rows.append(
        {
            "analysis_set": "matched_subset",
            "metric": "matching_quality",
            "test": "1:1 matching without replacement, 1-year age caliper",
            "n_covid": len(pairs),
            "n_hc": len(pairs),
            "covid_median": float(pairs["covid_age"].median()),
            "hc_median": float(pairs["hc_age"].median()),
            "effect": float(pairs["absolute_age_difference"].median()),
            "effect_definition": (
                "median absolute age difference; maximum = "
                f"{int(pairs['absolute_age_difference'].max())} year"
            ),
            "p_value": np.nan,
        }
    )
    return pd.DataFrame(rows)


def bootstrap_matched_length_curves(
    pairs: pd.DataFrame,
) -> tuple[pd.DataFrame, list[int]]:
    """Compute equal-sample-weight mean curves and pair-bootstrap intervals."""
    curves = pd.read_csv(LENGTH_CURVES_PATH, sep="\t")
    required = {"sample_id", "group", "length_bp", "density_sigma15"}
    if not required.issubset(curves.columns):
        raise ValueError(
            f"Missing length-curve columns: {required - set(curves.columns)}"
        )

    covid_ids = pairs["covid_sample_id"].tolist()
    hc_ids = pairs["hc_sample_id"].tolist()
    selected_ids = covid_ids + hc_ids
    selected = curves.loc[curves["sample_id"].isin(selected_ids)].copy()
    found_ids = set(selected["sample_id"])
    if found_ids != set(selected_ids):
        raise ValueError(
            f"Missing sample-level length curves for: {sorted(set(selected_ids) - found_ids)}"
        )

    pivot = selected.pivot(
        index="sample_id", columns="length_bp", values="density_sigma15"
    )
    lengths = pivot.columns.to_numpy(dtype=int)
    if lengths[0] != 50 or lengths[-1] != 1000:
        raise ValueError(
            f"Expected length grid 50--1000 bp, found {lengths[0]}--{lengths[-1]}"
        )
    covid_matrix = pivot.loc[covid_ids].to_numpy(dtype=float)
    hc_matrix = pivot.loc[hc_ids].to_numpy(dtype=float)
    if covid_matrix.shape != hc_matrix.shape:
        raise ValueError("Matched groups have inconsistent length-curve dimensions")

    rng = np.random.default_rng(RANDOM_SEED)
    resampled_pairs = rng.integers(
        0, len(pairs), size=(BOOTSTRAP_ITERATIONS, len(pairs))
    )
    covid_bootstrap = covid_matrix[resampled_pairs].mean(axis=1)
    hc_bootstrap = hc_matrix[resampled_pairs].mean(axis=1)

    summary = pd.DataFrame(
        {
            "length_bp": lengths,
            "covid_mean_density": covid_matrix.mean(axis=0),
            "covid_pair_bootstrap_ci_low": np.quantile(
                covid_bootstrap, 0.025, axis=0
            ),
            "covid_pair_bootstrap_ci_high": np.quantile(
                covid_bootstrap, 0.975, axis=0
            ),
            "hc_mean_density": hc_matrix.mean(axis=0),
            "hc_pair_bootstrap_ci_low": np.quantile(
                hc_bootstrap, 0.025, axis=0
            ),
            "hc_pair_bootstrap_ci_high": np.quantile(
                hc_bootstrap, 0.975, axis=0
            ),
        }
    )

    peaks = pd.read_csv(PEAK_PATH, sep="\t")
    stable = peaks.loc[
        peaks["stable_consensus_peak"].astype(str).isin({"1", "True", "true"})
    ]
    stable_positions = stable["primary_peak_position_bp"].astype(int).tolist()
    # v1 pinned the peak positions it was built with. Under the methods call
    # set the stable consensus peaks move by at most 1 bp (196/366/571 vs
    # 196/365/571), so the pin is expressed as a module constant that the v2
    # driver rebinds, rather than a literal that aborts the run.
    if EXPECT_STABLE_PEAKS is not None and stable_positions != EXPECT_STABLE_PEAKS:
        raise ValueError(f"Unexpected stable consensus peaks: {stable_positions}")
    print(f"stable consensus peaks: {stable_positions}", flush=True)
    return summary, stable_positions


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.18) -> None:
    ax.text(
        x,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="top",
    )


def draw_matching_panel(
    ax: plt.Axes, data: pd.DataFrame, pairs: pd.DataFrame
) -> None:
    """Panel a: full age distributions and selected matched pairs."""
    group_x = {GROUP_COVID: 0.0, GROUP_HC: 1.0}
    rng = np.random.default_rng(RANDOM_SEED + 1)

    ax.axhspan(
        min(pairs["covid_age"].min(), pairs["hc_age"].min()),
        max(pairs["covid_age"].max(), pairs["hc_age"].max()),
        color=GREY_LIGHT, alpha=0.48, zorder=0,
    )
    for group in (GROUP_COVID, GROUP_HC):
        values = (
            data.loc[data["group"] == group]
            .sort_values(["age", "sample_id"])
            .reset_index(drop=True)
        )
        jitter = rng.uniform(-0.065, 0.065, len(values))
        ax.scatter(
            group_x[group] + jitter,
            values["age"],
            s=6,
            facecolor="white",
            edgecolor=GREY_MID,
            linewidth=0.4,
            alpha=0.88,
            zorder=1,
        )

    ordered = pairs.sort_values(["covid_age", "covid_sample_id"]).reset_index(drop=True)
    pair_jitter = np.linspace(-0.045, 0.045, len(ordered))
    for index, row in ordered.iterrows():
        xs = np.array([0.0, 1.0]) + pair_jitter[index]
        ys = np.array([row["covid_age"], row["hc_age"]], dtype=float)
        ax.plot(xs, ys, color=GREY_DARK, linewidth=0.45, alpha=0.72, zorder=2)
        ax.scatter(
            xs[0],
            ys[0],
            s=12,
            color=COVID_COLOR,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
        ax.scatter(
            xs[1],
            ys[1],
            s=12,
            color=HC_COLOR,
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )

    ax.set_xlim(-0.30, 1.30)
    ax.set_ylim(46, 96)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [
            f"COVID-19\n39 total; {len(pairs)} matched",
            f"HC\n39 total; {len(pairs)} matched",
        ]
    )
    ax.set_ylabel("Age (years)")
    ax.text(
        0.52,
        0.98,
        "Open grey, full cohort; colour, matched pairs",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=5.0,
        linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 1.0},
    )
    add_panel_label(ax, "a", x=-0.20)


def draw_paired_burden_panel(
    ax: plt.Axes,
    matched: pd.DataFrame,
    metric: str,
    ylabel: str,
    p_value: float,
    median_fold_change: float,
    panel_label: str,
) -> None:
    """Panels b/c: paired sample trajectories on a logarithmic y axis."""
    n_pairs = int(matched["pair_id"].nunique())
    covid = (
        matched.loc[matched["group"] == GROUP_COVID]
        .set_index("pair_id")
        .sort_index()
    )
    hc = (
        matched.loc[matched["group"] == GROUP_HC]
        .set_index("pair_id")
        .sort_index()
    )
    jitter = np.linspace(-0.035, 0.035, len(covid))
    covid_values = covid[metric].to_numpy(dtype=float)
    hc_values = hc[metric].to_numpy(dtype=float)

    for index in range(len(covid)):
        xs = np.array([0.0, 1.0]) + jitter[index]
        ys = np.array([covid_values[index], hc_values[index]])
        ax.plot(xs, ys, color=GREY_MID, linewidth=0.48, alpha=0.78, zorder=1)
        ax.scatter(
            xs[0],
            ys[0],
            s=13,
            color=COVID_COLOR,
            edgecolor="white",
            linewidth=0.35,
            zorder=2,
        )
        ax.scatter(
            xs[1],
            ys[1],
            s=13,
            color=HC_COLOR,
            edgecolor="white",
            linewidth=0.35,
            zorder=2,
        )

    medians = [np.median(covid_values), np.median(hc_values)]
    for x_position, median in enumerate(medians):
        ax.plot(
            [x_position - 0.12, x_position + 0.12],
            [median, median],
            color=BLACK,
            linewidth=1.15,
            solid_capstyle="butt",
            zorder=4,
        )

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(-0.28, 1.28)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"COVID-19\n(n = {n_pairs})", f"HC\n(n = {n_pairs})"]
    )
    ax.set_ylabel(ylabel)
    ax.text(
        0.5,
        0.985,
        f"Paired P = {p_value:.4f}\nmedian paired FC = {median_fold_change:.2g}×",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=5.2,
        linespacing=1.2,
    )
    add_panel_label(ax, panel_label, x=-0.22)


def draw_length_panel(
    ax: plt.Axes,
    curve_summary: pd.DataFrame,
    stable_positions: list[int],
    n_pairs: int,
) -> None:
    """Panel d: sample-weighted length density in the matched subset."""
    x = curve_summary["length_bp"].to_numpy(dtype=float)
    scale = 100.0
    covid_mean = curve_summary["covid_mean_density"].to_numpy() * scale
    covid_low = curve_summary["covid_pair_bootstrap_ci_low"].to_numpy() * scale
    covid_high = curve_summary["covid_pair_bootstrap_ci_high"].to_numpy() * scale
    hc_mean = curve_summary["hc_mean_density"].to_numpy() * scale
    hc_low = curve_summary["hc_pair_bootstrap_ci_low"].to_numpy() * scale
    hc_high = curve_summary["hc_pair_bootstrap_ci_high"].to_numpy() * scale

    ax.fill_between(
        x, covid_low, covid_high, color=COVID_LIGHT, alpha=0.60, linewidth=0
    )
    ax.fill_between(x, hc_low, hc_high, color=HC_LIGHT, alpha=0.60, linewidth=0)
    ax.plot(x, covid_mean, color=COVID_COLOR, label=f"COVID-19 (n = {n_pairs})")
    ax.plot(x, hc_mean, color=HC_COLOR, label=f"HC (n = {n_pairs})")

    y_max = max(float(covid_high.max()), float(hc_high.max()))
    for position in stable_positions:
        i = int(np.argmin(np.abs(x - float(position))))
        apex = max(float(covid_high[i]), float(hc_high[i]))
        ax.text(
            position,
            apex + y_max * 0.030,
            f"{position} bp",
            ha="center",
            va="bottom",
            fontsize=5.0,
            color=PEAK_LABEL_COLOR,
        )

    ax.set_xlim(50, 1000)
    ax.set_ylim(0, y_max * 1.12)
    ax.set_xticks([100, 300, 500, 700, 900])
    ax.set_xlabel("eccDNA length (bp)")
    ax.set_ylabel("Probability density (% per bp)")
    ax.legend(
        loc="upper right",
        ncol=2,
        handlelength=2.2,
        columnspacing=1.0,
        borderaxespad=0.2,
    )
    ax.text(
        0.995,
        0.73,
        "Mean ± pair-bootstrap 95% CI",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.0,
        color=GREY_DARK,
    )
    add_panel_label(ax, "d", x=-0.055)


def make_figure(
    data: pd.DataFrame,
    pairs: pd.DataFrame,
    matched: pd.DataFrame,
    statistics: pd.DataFrame,
    curve_summary: pd.DataFrame,
    stable_positions: list[int],
) -> plt.Figure:
    """Compose the final 183-mm multi-panel figure."""
    configure_matplotlib()
    fig = plt.figure(figsize=(183 / 25.4, 130 / 25.4))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=(1.03, 0.92),
        left=0.065,
        right=0.985,
        bottom=0.105,
        top=0.965,
        hspace=0.50,
    )
    top = outer[0].subgridspec(
        1, 3, width_ratios=(1.13, 1.0, 1.0), wspace=0.44
    )

    ax_a = fig.add_subplot(top[0, 0])
    ax_b = fig.add_subplot(top[0, 1])
    ax_c = fig.add_subplot(top[0, 2])
    ax_d = fig.add_subplot(outer[1, 0])

    draw_matching_panel(ax_a, data, pairs)

    burden_stats = statistics.set_index(["analysis_set", "metric"])
    count_row = burden_stats.loc[("matched_subset", "eccdna_count")]
    epm_row = burden_stats.loc[("matched_subset", "epm")]
    draw_paired_burden_panel(
        ax_b,
        matched,
        metric="eccdna_count",
        ylabel="Detected eccDNA count",
        p_value=float(count_row["p_value"]),
        median_fold_change=float(count_row["effect"]),
        panel_label="b",
    )
    draw_paired_burden_panel(
        ax_c,
        matched,
        metric="epm",
        ylabel="EPM",
        p_value=float(epm_row["p_value"]),
        median_fold_change=float(epm_row["effect"]),
        panel_label="c",
    )
    draw_length_panel(
        ax_d, curve_summary, stable_positions, int(pairs["pair_id"].nunique())
    )
    return fig


def save_outputs(fig: plt.Figure) -> None:
    """Save editable and submission-ready versions at exact physical size."""
    OUTPUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_BASE.with_suffix(".pdf"), facecolor="white")
    fig.savefig(OUTPUT_BASE.with_suffix(".svg"), facecolor="white")
    fig.savefig(OUTPUT_BASE.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(
        OUTPUT_BASE.with_suffix(".tiff"),
        dpi=600,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def write_caption(pairs: pd.DataFrame) -> None:
    n = len(pairs)
    caption = (
        "Supplementary Figure S1. Age-matched sensitivity analysis of plasma "
        f"eccDNA. (a) Age distributions in the full cohort (open grey points) and "
        f"in the {n} matched pairs (coloured points and connecting lines); shading "
        "marks the age span covered by the matched pairs. Each patient was matched "
        "1:1 and without replacement to a control whose age differed by at most one "
        "year, maximising the number of pairs, then minimising the total absolute "
        "age difference, then maximising sex concordance (Methods). "
        "(b, c) Paired comparisons of detected eccDNA count and eccDNA per million "
        "mapped reads (EPM), displayed on logarithmic axes. Horizontal bars denote "
        "medians; P values are from two-sided paired Wilcoxon signed-rank tests. FC "
        "denotes the median within-pair fold change. (d) Equally weighted "
        "sample-level eccDNA length-density curves in the matched subset. Lines show "
        "group means and shaded bands show 95% confidence intervals from 5,000 "
        "matched-pair bootstrap resamples. Labels mark the three stable consensus "
        "length peaks identified in the prespecified full-cohort peak analysis."
    )
    OUTPUT_BASE.with_name(f"{OUTPUT_BASE.name}_caption.txt").write_text(
        caption + "\n", encoding="utf-8"
    )


def main() -> None:
    data = load_sample_data()
    pairs, matched = match_samples(data)
    statistics = build_statistics(data, pairs, matched)
    curve_summary, stable_positions = bootstrap_matched_length_curves(pairs)

    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(
        SOURCE_DATA_DIR / "Figure_S1_redesigned_matched_pairs.tsv",
        sep="\t",
        index=False,
    )
    statistics.to_csv(
        SOURCE_DATA_DIR / "Figure_S1_redesigned_statistics.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    curve_summary.to_csv(
        SOURCE_DATA_DIR / "Figure_S1_redesigned_length_curves.tsv",
        sep="\t",
        index=False,
    )

    figure = make_figure(
        data, pairs, matched, statistics, curve_summary, stable_positions
    )
    save_outputs(figure)
    write_caption(pairs)

    print(f"Wrote {OUTPUT_BASE.with_suffix('.pdf')}")
    print(f"Matched pairs: {len(pairs)}")
    print(
        "Maximum absolute age difference: "
        f"{pairs['absolute_age_difference'].max():.0f} year"
    )


if __name__ == "__main__":
    main()
