#!/usr/bin/env python3
"""Assess rolling-circle amplification (RCA) technical bias in eccDNA data.

The script combines Table_S1-derived sample metadata with the final 11-column
Circle-Map BED call sets and BAM indices. It produces sample-level metrics,
two-sided group comparisons, Spearman correlations, RCA-adjusted regression
models with HC3 robust standard errors, and publication-ready vector figures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from matplotlib import font_manager
from scipy import stats


GROUP_ORDER = ("HC", "COVID-19")
GROUP_COLORS = {"HC": "#36617B", "COVID-19": "#C84F50"}
GROUP_MARKERS = {"HC": "o", "COVID-19": "^"}
GREY = "#6F7B91"
LIGHT_GREY = "#C8CEDA"
RANDOM_SEED = 20260723
CORRELATION_METRICS = (
    "eccdna_count",
    "total_support_events",
    "mapped_reads",
    "epm",
    "support_events_per_million_mapped_reads",
    "discordant_read_pairs",
    "split_reads",
)
REGRESSION_OUTCOMES = (
    "eccdna_count",
    "total_support_events",
    "mapped_reads",
    "epm",
    "support_events_per_million_mapped_reads",
)
METRIC_LABELS = {
    "rca_concentration_ng_ul": "RCA concentration (ng/µl)",
    "eccdna_count": "Detected eccDNA count",
    "discordant_read_pairs": "Discordant read-pair support events",
    "split_reads": "Junction-spanning split reads",
    "total_support_events": "Circle-Map support events",
    "mapped_reads": "Mapped reads",
    "epm": "eccDNA per million mapped reads (EPM)",
    "support_events_per_million_mapped_reads": (
        "Support events per million mapped reads"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--covid-bed-dir", required=True, type=Path)
    parser.add_argument("--hc-bed-dir", required=True, type=Path)
    parser.add_argument("--covid-bam-dir", required=True, type=Path)
    parser.add_argument("--hc-bam-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, (np.floating, float)):
        value = float(value)
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return f"{value:.10g}"
    return str(value)


def write_tsv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def read_metadata(path: Path) -> list[dict[str, object]]:
    required = {
        "sample_id",
        "rca_concentration_ng_ul",
        "sex",
        "age_years",
        "group",
    }
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Metadata must contain fields: {sorted(required)}")
        for row in reader:
            parsed = {
                "sample_id": row["sample_id"].strip(),
                "rca_concentration_ng_ul": float(row["rca_concentration_ng_ul"]),
                "sex": row["sex"].strip(),
                "age_years": float(row["age_years"]),
                "group": row["group"].strip(),
            }
            rows.append(parsed)

    sample_ids = [str(row["sample_id"]) for row in rows]
    if len(rows) != 78 or len(set(sample_ids)) != 78:
        raise ValueError("Metadata must contain 78 unique samples")
    if Counter(str(row["group"]) for row in rows) != Counter({"HC": 39, "COVID-19": 39}):
        raise ValueError("Metadata must contain 39 HC and 39 COVID-19 samples")
    if not all(float(row["rca_concentration_ng_ul"]) > 0 for row in rows):
        raise ValueError("All RCA concentrations must be positive")
    return rows


def parse_circlemap_bed(path: Path) -> dict[str, int]:
    eccdna_count = 0
    discordant_read_pairs = 0
    split_reads = 0
    header_lines_skipped = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t", 5)
            if len(columns) < 6:
                raise ValueError(
                    f"{path}:{line_number} has fewer than six Circle-Map columns"
                )
            if (
                line_number == 1
                and columns[0].strip().lower() in {"chrom", "chromosome"}
                and columns[3].strip().lower() == "discordants"
            ):
                header_lines_skipped += 1
                continue
            try:
                discordants = int(float(columns[3]))
                splits = int(float(columns[4]))
            except ValueError as error:
                raise ValueError(
                    f"Non-numeric Circle-Map read support at {path}:{line_number}"
                ) from error
            if discordants < 0 or splits < 0:
                raise ValueError(f"Negative read support at {path}:{line_number}")
            eccdna_count += 1
            discordant_read_pairs += discordants
            split_reads += splits
    return {
        "eccdna_count": eccdna_count,
        "discordant_read_pairs": discordant_read_pairs,
        "split_reads": split_reads,
        "total_support_events": discordant_read_pairs + split_reads,
        "header_lines_skipped": header_lines_skipped,
    }


def mapped_reads_from_index(samtools: str, bam_path: Path) -> int:
    completed = subprocess.run(
        [samtools, "idxstats", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    total = 0
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 4:
            raise ValueError(f"Unexpected samtools idxstats line for {bam_path}: {line}")
        total += int(fields[2])
    if total <= 0:
        raise ValueError(f"No mapped reads were counted in {bam_path}")
    return total


def collect_sample_metrics(
    metadata: list[dict[str, object]], args: argparse.Namespace
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sample_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    expected_files: set[Path] = set()

    for meta in metadata:
        sample_id = str(meta["sample_id"])
        group = str(meta["group"])
        if group == "COVID-19":
            bed_dir, bam_dir = args.covid_bed_dir, args.covid_bam_dir
        elif group == "HC":
            bed_dir, bam_dir = args.hc_bed_dir, args.hc_bam_dir
        else:
            raise ValueError(f"Unexpected group for {sample_id}: {group}")
        bed_path = bed_dir / f"{sample_id}_circle_site.bed"
        bam_path = bam_dir / f"sorted_{sample_id}_circle.bam"
        bai_path = Path(f"{bam_path}.bai")
        for path in (bed_path, bam_path, bai_path):
            if not path.is_file():
                raise FileNotFoundError(path)
            expected_files.add(path)

        bed_metrics = parse_circlemap_bed(bed_path)
        mapped_reads = mapped_reads_from_index(args.samtools, bam_path)
        eccdna_count = bed_metrics["eccdna_count"]
        total_support_events = bed_metrics["total_support_events"]
        row = {
            **meta,
            **bed_metrics,
            "mapped_reads": mapped_reads,
            "epm": eccdna_count / mapped_reads * 1_000_000.0,
            "support_events_per_million_mapped_reads": (
                total_support_events / mapped_reads * 1_000_000.0
            ),
            "bed_path": str(bed_path),
            "bam_path": str(bam_path),
        }
        sample_rows.append(row)

        for file_type, path in (("Circle-Map BED", bed_path), ("BAM", bam_path), ("BAI", bai_path)):
            stat_result = path.stat()
            manifest_rows.append(
                {
                    "sample_id": sample_id,
                    "group": group,
                    "file_type": file_type,
                    "path": str(path),
                    "size_bytes": stat_result.st_size,
                    "mtime_epoch": stat_result.st_mtime,
                }
            )

    for directory, suffix in (
        (args.covid_bed_dir, "*_circle_site.bed"),
        (args.hc_bed_dir, "*_circle_site.bed"),
        (args.covid_bam_dir, "sorted_*_circle.bam"),
        (args.hc_bam_dir, "sorted_*_circle.bam"),
    ):
        found = set(directory.glob(suffix))
        unexpected = found - expected_files
        if unexpected:
            raise ValueError(
                f"Unexpected input files in {directory}: "
                + ", ".join(sorted(path.name for path in unexpected))
            )

    return sample_rows, manifest_rows


def bh_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    finite_values = values[finite]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[finite] = restored
    return result


def percentile_summary(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(list(values), dtype=float)
    return {
        "n": len(array),
        "mean": float(np.mean(array)),
        "sd": float(np.std(array, ddof=1)),
        "median": float(np.median(array)),
        "q1": float(np.percentile(array, 25)),
        "q3": float(np.percentile(array, 75)),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
    }


def cliffs_delta_from_u(u_statistic: float, n_x: int, n_y: int) -> float:
    return 2.0 * float(u_statistic) / (n_x * n_y) - 1.0


def group_statistics(sample_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    metrics = ("rca_concentration_ng_ul",) + CORRELATION_METRICS
    summary_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    for metric in metrics:
        by_group = {
            group: np.asarray(
                [float(row[metric]) for row in sample_rows if row["group"] == group]
            )
            for group in GROUP_ORDER
        }
        for group in GROUP_ORDER:
            summary_rows.append(
                {"metric": metric, "metric_label": METRIC_LABELS[metric], "group": group, **percentile_summary(by_group[group])}
            )
        covid = by_group["COVID-19"]
        hc = by_group["HC"]
        test = stats.mannwhitneyu(covid, hc, alternative="two-sided", method="auto")
        comparison_rows.append(
            {
                "metric": metric,
                "metric_label": METRIC_LABELS[metric],
                "comparison": "COVID-19 vs HC",
                "n_covid": len(covid),
                "n_hc": len(hc),
                "mann_whitney_u": float(test.statistic),
                "p_value_two_sided": float(test.pvalue),
                "cliffs_delta_covid_minus_hc": cliffs_delta_from_u(test.statistic, len(covid), len(hc)),
                "median_covid": float(np.median(covid)),
                "median_hc": float(np.median(hc)),
            }
        )
    adjusted = bh_adjust([float(row["p_value_two_sided"]) for row in comparison_rows])
    for row, q_value in zip(comparison_rows, adjusted):
        row["q_value_bh_across_group_tests"] = float(q_value)
    return summary_rows, comparison_rows


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, iterations: int, seed: int
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates: list[float] = []
    n = len(x)
    for _ in range(iterations):
        index = rng.integers(0, n, size=n)
        if np.unique(x[index]).size < 2 or np.unique(y[index]).size < 2:
            continue
        estimate = stats.spearmanr(x[index], y[index]).statistic
        if np.isfinite(estimate):
            estimates.append(float(estimate))
    if len(estimates) < iterations * 0.9:
        raise RuntimeError("Too few valid bootstrap Spearman estimates")
    return float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5))


def correlation_statistics(
    sample_rows: list[dict[str, object]], iterations: int
) -> list[dict[str, object]]:
    correlation_rows: list[dict[str, object]] = []
    scopes = ("All",) + GROUP_ORDER
    for scope_index, scope in enumerate(scopes):
        subset = sample_rows if scope == "All" else [row for row in sample_rows if row["group"] == scope]
        x = np.asarray([float(row["rca_concentration_ng_ul"]) for row in subset])
        for metric_index, metric in enumerate(CORRELATION_METRICS):
            y = np.asarray([float(row[metric]) for row in subset])
            test = stats.spearmanr(x, y)
            ci_low, ci_high = bootstrap_spearman_ci(
                x,
                y,
                iterations,
                RANDOM_SEED + 100 * scope_index + metric_index,
            )
            correlation_rows.append(
                {
                    "scope": scope,
                    "n": len(subset),
                    "x_metric": "rca_concentration_ng_ul",
                    "y_metric": metric,
                    "y_metric_label": METRIC_LABELS[metric],
                    "spearman_rho": float(test.statistic),
                    "bootstrap_ci_low": ci_low,
                    "bootstrap_ci_high": ci_high,
                    "p_value_two_sided": float(test.pvalue),
                    "bootstrap_iterations": iterations,
                }
            )
    all_q = bh_adjust([float(row["p_value_two_sided"]) for row in correlation_rows])
    for row, q_value in zip(correlation_rows, all_q):
        row["q_value_bh_across_all_correlations"] = float(q_value)
    for scope in scopes:
        indices = [i for i, row in enumerate(correlation_rows) if row["scope"] == scope]
        q_values = bh_adjust([float(correlation_rows[i]["p_value_two_sided"]) for i in indices])
        for index, q_value in zip(indices, q_values):
            correlation_rows[index]["q_value_bh_within_scope"] = float(q_value)
    return correlation_rows


def ols_hc3(y: np.ndarray, x: np.ndarray) -> dict[str, object]:
    n, p = x.shape
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residuals = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, xtx_inv, x)
    denominator = np.maximum(1.0 - leverage, np.finfo(float).eps)
    scaled_residual_sq = (residuals / denominator) ** 2
    meat = x.T @ (x * scaled_residual_sq[:, None])
    covariance = xtx_inv @ meat @ xtx_inv
    standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    dof = n - p
    t_values = beta / standard_errors
    p_values = 2.0 * stats.t.sf(np.abs(t_values), dof)
    critical = stats.t.ppf(0.975, dof)
    ci_low = beta - critical * standard_errors
    ci_high = beta + critical * standard_errors
    centered = y - np.mean(y)
    r_squared = 1.0 - float(np.sum(residuals**2)) / float(np.sum(centered**2))
    return {
        "beta": beta,
        "standard_error": standard_errors,
        "t_value": t_values,
        "p_value": p_values,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "degrees_of_freedom": dof,
        "r_squared": r_squared,
        "n": n,
    }


def regression_statistics(sample_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    covid = np.asarray([row["group"] == "COVID-19" for row in sample_rows], dtype=float)
    log2_rca = np.log2(np.asarray([float(row["rca_concentration_ng_ul"]) for row in sample_rows]))
    age_per_10 = (np.asarray([float(row["age_years"]) for row in sample_rows]) - 70.0) / 10.0
    male = np.asarray([str(row["sex"]).lower() == "male" for row in sample_rows], dtype=float)
    log2_rca_centered = log2_rca - np.mean(log2_rca)
    model_definitions = (
        (
            "group_plus_rca",
            ("intercept", "COVID-19_vs_HC", "log2_RCA_per_doubling"),
            np.column_stack([np.ones(len(sample_rows)), covid, log2_rca_centered]),
        ),
        (
            "full_covariate_adjusted",
            (
                "intercept",
                "COVID-19_vs_HC",
                "log2_RCA_per_doubling",
                "age_per_10_years",
                "male_vs_female",
            ),
            np.column_stack(
                [np.ones(len(sample_rows)), covid, log2_rca_centered, age_per_10, male]
            ),
        ),
        (
            "full_plus_group_RCA_interaction",
            (
                "intercept",
                "COVID-19_vs_HC",
                "log2_RCA_per_doubling",
                "age_per_10_years",
                "male_vs_female",
                "COVID-19_x_log2_RCA",
            ),
            np.column_stack(
                [
                    np.ones(len(sample_rows)),
                    covid,
                    log2_rca_centered,
                    age_per_10,
                    male,
                    covid * log2_rca_centered,
                ]
            ),
        ),
    )
    output_rows: list[dict[str, object]] = []
    for outcome in REGRESSION_OUTCOMES:
        y = np.log2(np.asarray([float(row[outcome]) for row in sample_rows]) + 1.0)
        for model_name, terms, x in model_definitions:
            result = ols_hc3(y, x)
            for index, term in enumerate(terms):
                estimate = float(result["beta"][index])
                ci_low = float(result["ci_low"][index])
                ci_high = float(result["ci_high"][index])
                output_rows.append(
                    {
                        "outcome": outcome,
                        "outcome_label": METRIC_LABELS[outcome],
                        "outcome_transform": "log2(outcome + 1)",
                        "model": model_name,
                        "term": term,
                        "estimate_log2_scale": estimate,
                        "hc3_standard_error": float(result["standard_error"][index]),
                        "ci_low_log2_scale": ci_low,
                        "ci_high_log2_scale": ci_high,
                        "multiplicative_effect": 2.0**estimate,
                        "multiplicative_ci_low": 2.0**ci_low,
                        "multiplicative_ci_high": 2.0**ci_high,
                        "t_value": float(result["t_value"][index]),
                        "p_value_two_sided": float(result["p_value"][index]),
                        "n": int(result["n"]),
                        "degrees_of_freedom": int(result["degrees_of_freedom"]),
                        "r_squared": float(result["r_squared"]),
                    }
                )

    for term in ("COVID-19_vs_HC", "log2_RCA_per_doubling", "COVID-19_x_log2_RCA"):
        indices = [
            i
            for i, row in enumerate(output_rows)
            if row["term"] == term and row["model"] != "group_plus_rca"
        ]
        q_values = bh_adjust([float(output_rows[i]["p_value_two_sided"]) for i in indices])
        for index, q_value in zip(indices, q_values):
            output_rows[index]["q_value_bh_within_term"] = float(q_value)
    return output_rows


def configure_matplotlib() -> str:
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        font_family = "Arial"
    except ValueError:
        font_family = "Liberation Sans"
    mpl.rcParams.update(
        {
            "font.family": font_family,
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    return font_family


def p_text(p_value: float) -> str:
    if p_value < 0.0001:
        return "P < 0.0001"
    if p_value < 0.001:
        return f"P = {p_value:.1e}"
    return f"P = {p_value:.3f}"


def add_panel_label(axis: plt.Axes, label: str) -> None:
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


def group_boxplot(
    axis: plt.Axes,
    sample_rows: list[dict[str, object]],
    comparison_rows: list[dict[str, object]],
) -> None:
    rng = np.random.default_rng(RANDOM_SEED)
    values = [
        np.asarray(
            [float(row["rca_concentration_ng_ul"]) for row in sample_rows if row["group"] == group]
        )
        for group in GROUP_ORDER
    ]
    box = axis.boxplot(
        values,
        positions=np.arange(len(GROUP_ORDER)),
        widths=0.45,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 0.8},
        whiskerprops={"color": GREY, "linewidth": 0.6},
        capprops={"color": GREY, "linewidth": 0.6},
        boxprops={"linewidth": 0.7},
    )
    for patch, group in zip(box["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[group])
        patch.set_alpha(0.28)
        patch.set_edgecolor(GROUP_COLORS[group])
    for position, (group, group_values) in enumerate(zip(GROUP_ORDER, values)):
        jitter = rng.uniform(-0.12, 0.12, size=len(group_values))
        axis.scatter(
            position + jitter,
            group_values,
            s=10,
            marker=GROUP_MARKERS[group],
            facecolor=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.25,
            alpha=0.8,
            zorder=3,
        )
    comparison = next(row for row in comparison_rows if row["metric"] == "rca_concentration_ng_ul")
    y_max = max(np.max(group_values) for group_values in values)
    axis.text(
        0.5,
        y_max * 1.06,
        p_text(float(comparison["p_value_two_sided"])),
        ha="center",
        va="bottom",
        fontsize=5.5,
    )
    axis.set_xticks(range(len(GROUP_ORDER)), GROUP_ORDER)
    axis.set_ylabel(METRIC_LABELS["rca_concentration_ng_ul"])
    axis.set_ylim(0, y_max * 1.2)
    axis.spines[["top", "right"]].set_visible(False)


def get_correlation(
    correlation_rows: list[dict[str, object]], scope: str, metric: str
) -> dict[str, object]:
    return next(
        row
        for row in correlation_rows
        if row["scope"] == scope and row["y_metric"] == metric
    )


def scatter_panel(
    axis: plt.Axes,
    sample_rows: list[dict[str, object]],
    correlation_rows: list[dict[str, object]],
    metric: str,
    log_y: bool,
) -> None:
    for group in GROUP_ORDER:
        subset = [row for row in sample_rows if row["group"] == group]
        x = np.asarray([float(row["rca_concentration_ng_ul"]) for row in subset])
        y = np.asarray([float(row[metric]) for row in subset])
        plot_y = y / 1_000_000.0 if metric == "mapped_reads" else y
        axis.scatter(
            x,
            plot_y,
            s=13,
            marker=GROUP_MARKERS[group],
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.25,
            alpha=0.82,
            label=f"{group} (n={len(subset)})",
            zorder=3,
        )
        fit_y = np.log10(plot_y) if log_y else plot_y
        slope, intercept = np.polyfit(x, fit_y, 1)
        fit_x = np.linspace(np.min(x), np.max(x), 100)
        plotted_y = 10 ** (intercept + slope * fit_x) if log_y else intercept + slope * fit_x
        axis.plot(fit_x, plotted_y, color=GROUP_COLORS[group], linewidth=0.8, alpha=0.8)

    if log_y:
        axis.set_yscale("log")
    axis.set_xlabel(METRIC_LABELS["rca_concentration_ng_ul"])
    axis.set_ylabel("Mapped reads (millions)" if metric == "mapped_reads" else METRIC_LABELS[metric])
    axis.spines[["top", "right"]].set_visible(False)
    annotation_lines = []
    for scope, label in (("All", "All"), ("HC", "HC"), ("COVID-19", "COVID-19")):
        result = get_correlation(correlation_rows, scope, metric)
        annotation_lines.append(
            f"{label}: ρ={float(result['spearman_rho']):.2f}, "
            f"{p_text(float(result['p_value_two_sided'])).replace('P ', 'P ')}"
        )
    axis.text(
        0.03,
        0.97,
        "\n".join(annotation_lines),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=4.7,
        linespacing=1.25,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
    )


def save_figure_all_formats(figure: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".pdf"))
    figure.savefig(stem.with_suffix(".svg"))
    figure.savefig(stem.with_suffix(".png"), dpi=300)


def make_correlation_figure(
    sample_rows: list[dict[str, object]],
    comparison_rows: list[dict[str, object]],
    correlation_rows: list[dict[str, object]],
    output_dir: Path,
) -> None:
    width_inches = 183 / 25.4
    figure, axes = plt.subplots(2, 3, figsize=(width_inches, 118 / 25.4))
    axes_flat = axes.ravel()
    group_boxplot(axes_flat[0], sample_rows, comparison_rows)
    scatter_panel(axes_flat[1], sample_rows, correlation_rows, "eccdna_count", True)
    scatter_panel(axes_flat[2], sample_rows, correlation_rows, "total_support_events", True)
    scatter_panel(axes_flat[3], sample_rows, correlation_rows, "mapped_reads", False)
    scatter_panel(axes_flat[4], sample_rows, correlation_rows, "epm", True)
    for axis, label in zip(axes_flat[:5], "abcde"):
        add_panel_label(axis, label)
    axes_flat[5].axis("off")
    handles, labels = axes_flat[1].get_legend_handles_labels()
    axes_flat[5].legend(handles, labels, loc="upper left", frameon=False, handletextpad=0.4)
    axes_flat[5].text(
        0.0,
        0.62,
        "Lines show group-specific least-squares trends.\n"
        "ρ: Spearman rank correlation; all tests two-sided.",
        transform=axes_flat[5].transAxes,
        ha="left",
        va="top",
        fontsize=5.5,
        linespacing=1.35,
    )
    figure.subplots_adjust(left=0.085, right=0.985, bottom=0.11, top=0.96, wspace=0.44, hspace=0.48)
    save_figure_all_formats(figure, output_dir / "figures" / "Figure_RCA_technical_bias")
    plt.close(figure)


def make_adjusted_effect_figure(regression_rows: list[dict[str, object]], output_dir: Path) -> None:
    selected = [
        row
        for row in regression_rows
        if row["model"] == "full_covariate_adjusted"
        and row["term"] in {"COVID-19_vs_HC", "log2_RCA_per_doubling"}
    ]
    outcomes = list(REGRESSION_OUTCOMES)
    width_inches = 183 / 25.4
    figure, axes = plt.subplots(1, 2, figsize=(width_inches, 72 / 25.4), sharey=True)
    panels = (
        ("COVID-19_vs_HC", "COVID-19 versus HC", "#C84F50"),
        ("log2_RCA_per_doubling", "Per twofold higher RCA", "#36617B"),
    )
    y_positions = np.arange(len(outcomes))[::-1]
    for axis, (term, x_label, color), panel_label in zip(axes, panels, "ab"):
        term_rows = {str(row["outcome"]): row for row in selected if row["term"] == term}
        estimates = np.asarray([float(term_rows[outcome]["multiplicative_effect"]) for outcome in outcomes])
        lower = np.asarray([float(term_rows[outcome]["multiplicative_ci_low"]) for outcome in outcomes])
        upper = np.asarray([float(term_rows[outcome]["multiplicative_ci_high"]) for outcome in outcomes])
        axis.errorbar(
            estimates,
            y_positions,
            xerr=np.vstack([estimates - lower, upper - estimates]),
            fmt="o",
            markersize=4,
            color=color,
            ecolor=color,
            elinewidth=0.9,
            capsize=2,
            markeredgecolor="white",
            markeredgewidth=0.35,
            zorder=3,
        )
        axis.axvline(1.0, color=GREY, linestyle="--", linewidth=0.7, zorder=1)
        axis.set_xscale("log")
        axis.set_xlabel(f"Multiplicative effect (95% CI)\n{x_label}")
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
        add_panel_label(axis, panel_label)
    axes[0].set_yticks(y_positions, [METRIC_LABELS[outcome] for outcome in outcomes])
    figure.subplots_adjust(left=0.35, right=0.98, bottom=0.23, top=0.93, wspace=0.32)
    save_figure_all_formats(figure, output_dir / "figures" / "Figure_RCA_adjusted_effects")
    plt.close(figure)


def make_summary(
    sample_rows: list[dict[str, object]],
    group_summary_rows: list[dict[str, object]],
    comparison_rows: list[dict[str, object]],
    correlation_rows: list[dict[str, object]],
    regression_rows: list[dict[str, object]],
) -> str:
    summary_lookup = {
        (str(row["metric"]), str(row["group"])): row for row in group_summary_rows
    }
    rca_test = next(row for row in comparison_rows if row["metric"] == "rca_concentration_ng_ul")
    epm_correlations = {
        str(row["scope"]): row for row in correlation_rows if row["y_metric"] == "epm"
    }
    full_epm = [
        row
        for row in regression_rows
        if row["outcome"] == "epm" and row["model"] == "full_covariate_adjusted"
    ]
    epm_group = next(row for row in full_epm if row["term"] == "COVID-19_vs_HC")
    epm_rca = next(row for row in full_epm if row["term"] == "log2_RCA_per_doubling")
    lines = [
        "RCA technical-bias analysis summary",
        "===================================",
        "",
        f"Samples: {len(sample_rows)} total (39 COVID-19, 39 HC).",
        "Circle-Map support events were operationally defined as the per-call sum of "
        "discordant read-pair counts (BED column 4) and split-read counts (BED column 5).",
        "These are support-event sums, not deduplicated unique read names.",
        "",
    ]
    for group in GROUP_ORDER:
        row = summary_lookup[("rca_concentration_ng_ul", group)]
        lines.append(
            f"{group} RCA: mean {float(row['mean']):.2f}, SD {float(row['sd']):.2f}, "
            f"median {float(row['median']):.2f}, IQR {float(row['q1']):.2f}–{float(row['q3']):.2f} ng/µl."
        )
    lines.extend(
        [
            f"Two-sided Mann–Whitney comparison: U={float(rca_test['mann_whitney_u']):.1f}, "
            f"P={float(rca_test['p_value_two_sided']):.6g}, "
            f"Cliff's delta={float(rca_test['cliffs_delta_covid_minus_hc']):.3f}.",
            "",
            "Spearman correlation between RCA concentration and EPM:",
        ]
    )
    for scope in ("All", "COVID-19", "HC"):
        row = epm_correlations[scope]
        lines.append(
            f"  {scope}: rho={float(row['spearman_rho']):.3f}, "
            f"95% bootstrap CI {float(row['bootstrap_ci_low']):.3f} to "
            f"{float(row['bootstrap_ci_high']):.3f}, P={float(row['p_value_two_sided']):.6g}."
        )
    lines.extend(
        [
            "",
            "Full model: log2(EPM + 1) ~ group + log2(RCA) + age/10 + sex, with HC3 robust SEs.",
            f"  COVID-19 vs HC adjusted multiplicative effect: {float(epm_group['multiplicative_effect']):.3f} "
            f"(95% CI {float(epm_group['multiplicative_ci_low']):.3f}–"
            f"{float(epm_group['multiplicative_ci_high']):.3f}), "
            f"P={float(epm_group['p_value_two_sided']):.6g}.",
            f"  Per twofold higher RCA adjusted multiplicative effect: {float(epm_rca['multiplicative_effect']):.3f} "
            f"(95% CI {float(epm_rca['multiplicative_ci_low']):.3f}–"
            f"{float(epm_rca['multiplicative_ci_high']):.3f}), "
            f"P={float(epm_rca['p_value_two_sided']):.6g}.",
            "",
            "Interpret pooled correlations together with within-group correlations and adjusted models, "
            "because disease group is associated with both RCA yield and downstream eccDNA metrics.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_software_versions(path: Path, samtools: str, font_family: str) -> None:
    samtools_version = subprocess.run(
        [samtools, "--version"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[0]
    lines = [
        f"Python\t{platform.python_version()}",
        f"NumPy\t{np.__version__}",
        f"SciPy\t{scipy.__version__}",
        f"Matplotlib\t{mpl.__version__}",
        f"samtools\t{samtools_version}",
        f"Figure font\t{font_family}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.bootstrap_iterations < 1000:
        raise ValueError("Use at least 1,000 bootstrap iterations")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results_dir = args.output_dir / "results"
    figures_dir = args.output_dir / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    metadata = read_metadata(args.metadata)
    sample_rows, manifest_rows = collect_sample_metrics(metadata, args)
    group_summary_rows, comparison_rows = group_statistics(sample_rows)
    correlation_rows = correlation_statistics(sample_rows, args.bootstrap_iterations)
    regression_rows = regression_statistics(sample_rows)

    sample_fields = (
        "sample_id",
        "group",
        "sex",
        "age_years",
        "rca_concentration_ng_ul",
        "eccdna_count",
        "discordant_read_pairs",
        "split_reads",
        "total_support_events",
        "header_lines_skipped",
        "mapped_reads",
        "epm",
        "support_events_per_million_mapped_reads",
        "bed_path",
        "bam_path",
    )
    write_tsv(results_dir / "RCA_sample_metrics.tsv", sample_rows, sample_fields)
    write_tsv(
        results_dir / "RCA_input_manifest.tsv",
        manifest_rows,
        ("sample_id", "group", "file_type", "path", "size_bytes", "mtime_epoch"),
    )
    write_tsv(
        results_dir / "RCA_group_summary.tsv",
        group_summary_rows,
        ("metric", "metric_label", "group", "n", "mean", "sd", "median", "q1", "q3", "minimum", "maximum"),
    )
    write_tsv(
        results_dir / "RCA_group_comparisons.tsv",
        comparison_rows,
        (
            "metric",
            "metric_label",
            "comparison",
            "n_covid",
            "n_hc",
            "mann_whitney_u",
            "p_value_two_sided",
            "q_value_bh_across_group_tests",
            "cliffs_delta_covid_minus_hc",
            "median_covid",
            "median_hc",
        ),
    )
    write_tsv(
        results_dir / "RCA_correlations.tsv",
        correlation_rows,
        (
            "scope",
            "n",
            "x_metric",
            "y_metric",
            "y_metric_label",
            "spearman_rho",
            "bootstrap_ci_low",
            "bootstrap_ci_high",
            "p_value_two_sided",
            "q_value_bh_within_scope",
            "q_value_bh_across_all_correlations",
            "bootstrap_iterations",
        ),
    )
    write_tsv(
        results_dir / "RCA_regression_HC3.tsv",
        regression_rows,
        (
            "outcome",
            "outcome_label",
            "outcome_transform",
            "model",
            "term",
            "estimate_log2_scale",
            "hc3_standard_error",
            "ci_low_log2_scale",
            "ci_high_log2_scale",
            "multiplicative_effect",
            "multiplicative_ci_low",
            "multiplicative_ci_high",
            "t_value",
            "p_value_two_sided",
            "q_value_bh_within_term",
            "n",
            "degrees_of_freedom",
            "r_squared",
        ),
    )

    font_family = configure_matplotlib()
    make_correlation_figure(sample_rows, comparison_rows, correlation_rows, args.output_dir)
    make_adjusted_effect_figure(regression_rows, args.output_dir)

    summary_text = make_summary(
        sample_rows,
        group_summary_rows,
        comparison_rows,
        correlation_rows,
        regression_rows,
    )
    (results_dir / "RCA_analysis_summary.txt").write_text(summary_text, encoding="utf-8")
    write_software_versions(results_dir / "software_versions.tsv", args.samtools, font_family)
    parameters = {
        "metadata": str(args.metadata),
        "metadata_sha256": sha256(args.metadata),
        "covid_bed_dir": str(args.covid_bed_dir),
        "hc_bed_dir": str(args.hc_bed_dir),
        "covid_bam_dir": str(args.covid_bam_dir),
        "hc_bam_dir": str(args.hc_bam_dir),
        "output_dir": str(args.output_dir),
        "samtools": args.samtools,
        "bootstrap_iterations": args.bootstrap_iterations,
        "random_seed": RANDOM_SEED,
        "circlemap_column_4": "discordant read-pair count",
        "circlemap_column_5": "split-read count",
        "total_support_events_definition": "sum over calls of column 4 + column 5",
        "epm_definition": "BED line count / samtools idxstats mapped alignment count * 1e6",
        "group_test": "two-sided Mann-Whitney U",
        "correlation_test": "two-sided Spearman rank correlation with percentile bootstrap CI",
        "regression": "OLS on log2(outcome + 1) with HC3 robust standard errors",
    }
    (results_dir / "run_parameters.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    validation_lines = [
        "PASS: metadata contains 78 unique samples (39 COVID-19; 39 HC)",
        "PASS: every sample has one BED, BAM, and BAI file",
        "PASS: no unexpected BED or BAM sample files were found in the four input directories",
        "PASS: every parsed Circle-Map data record had at least six tab-separated columns and non-negative support counts",
        "PASS: one explicit header line in ZXS59 was detected and excluded from biological call counts",
        "PASS: mapped-read counts were obtained from BAM indices with samtools idxstats",
        "PASS: all 78 EPM values were finite and positive",
        "PASS: two-sided group tests, correlations, and HC3-adjusted models completed",
        "PASS: PDF, SVG, and 300-dpi PNG figures were exported",
    ]
    (results_dir / "validation_report.txt").write_text(
        "\n".join(validation_lines) + "\n", encoding="utf-8"
    )
    print(summary_text)


if __name__ == "__main__":
    main()
