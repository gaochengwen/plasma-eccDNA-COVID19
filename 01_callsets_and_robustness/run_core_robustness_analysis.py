#!/usr/bin/env python3
"""Statistical synthesis and figure for Circle-Map robustness callsets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy import stats


SEED = 20260728
BOOTSTRAPS = 10_000
LENGTH_BINS = ("lt200", "200_399", "400_599", "600_999", "1000_1999", "ge2000")
LENGTH_LABELS = ("<200", "200–399", "400–599", "600–999", "1,000–1,999", "≥2,000")
CALLSET_ORDER = (
    "circlemap_current",
    "circlemap_methods",
    "circlemap_strict",
    "circlemap_very_strict",
    "circlemap_artifact_masked",
    "circlemap_high_support_masked",
    "consensus_t10",
)
CALLSET_LABELS = {
    "circlemap_current": "Current\nCircle-Map",
    "circlemap_methods": "Methods\nfilter",
    "circlemap_strict": "split≥10,\nscore>1,000",
    "circlemap_very_strict": "split≥20,\nscore>2,000",
    "circlemap_artifact_masked": "Artifact-\nmasked",
    "circlemap_high_support_masked": "High-support\n+ masked",
    "consensus_t10": "Two-caller\nconsensus",
}
TARGETS = {
    "CTNNA2": ("chr2", 79_885_063, 79_885_444),
    "CAB39": ("chr2", 230_727_472, 230_727_831),
    "chr22_target": ("chr22", 42_359_179, 42_359_558),
    "SDK1": ("chr7", 3_619_742, 3_620_122),
    "TCF7L1": ("chr2", 85_213_855, 85_214_235),
}
BLUE = "#0272B2"
ORANGE = "#EC6F00"
RED = "#C93E3F"
PURPLE = "#A84E94"
GREEN = "#459434"
TEAL = "#019AA3"
YELLOW = "#CCA02C"
GREY1 = "#E6E6ED"
GREY2 = "#C8CEDA"
GREY4 = "#6F7B91"
GREY6 = "#253247"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_group(value: str) -> str:
    return "COVID-19" if value.strip().lower() in {"covid", "covid-19", "covid19"} else "HC"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    fieldnames: list[str] = []
    for row in rows:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(values: list[float]) -> list[float]:
    output = [math.nan] * len(values)
    valid = [index for index, value in enumerate(values) if math.isfinite(value)]
    if not valid:
        return output
    ordered = sorted(valid, key=lambda index: values[index])
    running = 1.0
    for reverse_rank, index in enumerate(reversed(ordered), start=1):
        rank = len(ordered) - reverse_rank + 1
        running = min(running, values[index] * len(ordered) / rank)
        output[index] = min(1.0, running)
    return output


def cliffs_delta(covid: np.ndarray, hc: np.ndarray) -> float:
    difference = covid[:, None] - hc[None, :]
    return float(
        (np.count_nonzero(difference > 0) - np.count_nonzero(difference < 0))
        / difference.size
    )


def bootstrap_median_ratio(covid: np.ndarray, hc: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    ratios = np.empty(BOOTSTRAPS)
    for index in range(BOOTSTRAPS):
        covid_median = np.median(rng.choice(covid, size=len(covid), replace=True))
        hc_median = np.median(rng.choice(hc, size=len(hc), replace=True))
        ratios[index] = (covid_median + 1e-12) / (hc_median + 1e-12)
    low, high = np.percentile(ratios, [2.5, 97.5])
    return float(low), float(high)


def summarize_bed(path: Path) -> tuple[int, dict[str, int]]:
    count = 0
    bins = defaultdict(int)
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                length = int(fields[2]) - int(fields[1])
            except ValueError:
                continue
            if length <= 0:
                continue
            count += 1
            if length < 200:
                name = "lt200"
            elif length < 400:
                name = "200_399"
            elif length < 600:
                name = "400_599"
            elif length < 1000:
                name = "600_999"
            elif length < 2000:
                name = "1000_1999"
            else:
                name = "ge2000"
            bins[name] += 1
    return count, dict(bins)


def load_metrics(root: Path, metadata: list[dict[str, str]]) -> list[dict]:
    mapped = {row["sample_id"]: int(row["mapped_alignments_idxstats"]) for row in metadata}
    groups = {row["sample_id"]: normalize_group(row["group"]) for row in metadata}
    output = []
    standard = root / "tables" / "circlemap_output_callset_sample_metrics.tsv"
    artifact = root / "tables" / "artifact_masked_sample_metrics.tsv"
    for path in (standard, artifact):
        if not path.is_file():
            continue
        for row in read_tsv(path):
            sample = row["sample_id"]
            output.append(
                {
                    "sample_id": sample,
                    "group": normalize_group(row["group"]),
                    "callset": row["callset"],
                    "call_count": int(row["call_count"]),
                    "mapped_alignments": int(row["mapped_alignments"]),
                    "epm": float(row["epm"]),
                    "retention_fraction": float(
                        row.get("retention_fraction_of_current")
                        or row.get("retention_fraction_of_methods")
                        or "nan"
                    ),
                    **{
                        f"length_{name}_count": int(row[f"length_{name}_count"])
                        for name in LENGTH_BINS
                    },
                }
            )

    rca_path = root / "metadata" / "RCA_sample_metrics.tsv"
    if rca_path.is_file():
        for row in read_tsv(rca_path):
            output.append(
                {
                    "sample_id": row["sample_id"],
                    "group": normalize_group(row["group"]),
                    "callset": "circlemap_current",
                    "call_count": int(row["eccdna_count"]),
                    "mapped_alignments": int(row["mapped_reads"]),
                    "epm": float(row["epm"]),
                    "retention_fraction": 1.0,
                    **{f"length_{name}_count": math.nan for name in LENGTH_BINS},
                }
            )

    consensus_summary = root / "tables" / "caller_concordance_all.tsv"
    if consensus_summary.is_file():
        summary = {
            row["sample_id"]: row
            for row in read_tsv(consensus_summary)
            if int(row["tolerance_bp"]) == 10
        }
        for sample in sorted(summary):
            bed = root / "callsets" / "consensus_t10" / f"{sample}_circle_site.bed"
            count, bins = summarize_bed(bed)
            if count != int(summary[sample]["consensus_count"]):
                raise RuntimeError(f"Consensus count mismatch: {sample}")
            output.append(
                {
                    "sample_id": sample,
                    "group": groups[sample],
                    "callset": "consensus_t10",
                    "call_count": count,
                    "mapped_alignments": mapped[sample],
                    "epm": count / mapped[sample] * 1_000_000,
                    "retention_fraction": float(
                        summary[sample]["circlemap_retention_fraction"]
                    ),
                    **{f"length_{name}_count": bins.get(name, 0) for name in LENGTH_BINS},
                }
            )
    return output


def group_statistics(metrics: list[dict]) -> tuple[list[dict], list[dict]]:
    summary_rows = []
    test_rows = []
    callsets = [name for name in CALLSET_ORDER if any(row["callset"] == name for row in metrics)]
    for callset in callsets:
        subset = [row for row in metrics if row["callset"] == callset]
        for metric in ("call_count", "epm", "retention_fraction"):
            if not all(math.isfinite(float(row[metric])) for row in subset):
                continue
            by_group = {
                group: np.asarray(
                    [float(row[metric]) for row in subset if row["group"] == group]
                )
                for group in ("COVID-19", "HC")
            }
            for group, values in by_group.items():
                summary_rows.append(
                    {
                        "callset": callset,
                        "metric": metric,
                        "group": group,
                        "n": len(values),
                        "mean": float(np.mean(values)),
                        "sd": float(np.std(values, ddof=1)),
                        "median": float(np.median(values)),
                        "q1": float(np.percentile(values, 25)),
                        "q3": float(np.percentile(values, 75)),
                        "minimum": float(np.min(values)),
                        "maximum": float(np.max(values)),
                    }
                )
            covid, hc = by_group["COVID-19"], by_group["HC"]
            result = stats.mannwhitneyu(covid, hc, alternative="two-sided", method="asymptotic")
            ci_low, ci_high = bootstrap_median_ratio(
                covid, hc, SEED + 100 * callsets.index(callset) + len(test_rows)
            )
            test_rows.append(
                {
                    "callset": callset,
                    "metric": metric,
                    "covid_n": len(covid),
                    "hc_n": len(hc),
                    "covid_median": float(np.median(covid)),
                    "hc_median": float(np.median(hc)),
                    "median_ratio_COVID_over_HC": float(
                        (np.median(covid) + 1e-12) / (np.median(hc) + 1e-12)
                    ),
                    "median_ratio_bootstrap_ci_low": ci_low,
                    "median_ratio_bootstrap_ci_high": ci_high,
                    "mannwhitney_u": float(result.statistic),
                    "p_value_two_sided": float(result.pvalue),
                    "cliffs_delta_COVID_vs_HC": cliffs_delta(covid, hc),
                    "bootstrap_iterations": BOOTSTRAPS,
                }
            )
    for metric in ("call_count", "epm", "retention_fraction"):
        indices = [i for i, row in enumerate(test_rows) if row["metric"] == metric]
        adjusted = bh_adjust([test_rows[index]["p_value_two_sided"] for index in indices])
        for index, value in zip(indices, adjusted):
            test_rows[index]["q_value_bh_within_metric"] = value
    return summary_rows, test_rows


def fragment_statistics(metrics: list[dict]) -> tuple[list[dict], list[dict]]:
    sample_rows = []
    test_rows = []
    for row in metrics:
        counts = [row[f"length_{name}_count"] for name in LENGTH_BINS]
        if not all(math.isfinite(float(value)) for value in counts):
            continue
        denominator = sum(int(value) for value in counts)
        for name, count in zip(LENGTH_BINS, counts):
            sample_rows.append(
                {
                    "sample_id": row["sample_id"],
                    "group": row["group"],
                    "callset": row["callset"],
                    "length_bin": name,
                    "count": int(count),
                    "proportion": int(count) / denominator if denominator else math.nan,
                }
            )
    for callset in [name for name in CALLSET_ORDER if any(row["callset"] == name for row in sample_rows)]:
        for name in LENGTH_BINS:
            subset = [
                row
                for row in sample_rows
                if row["callset"] == callset
                and row["length_bin"] == name
                and math.isfinite(float(row["proportion"]))
            ]
            covid = np.asarray(
                [row["proportion"] for row in subset if row["group"] == "COVID-19"]
            )
            hc = np.asarray([row["proportion"] for row in subset if row["group"] == "HC"])
            result = stats.mannwhitneyu(covid, hc, alternative="two-sided", method="asymptotic")
            test_rows.append(
                {
                    "callset": callset,
                    "length_bin": name,
                    "covid_n": len(covid),
                    "hc_n": len(hc),
                    "covid_median_proportion": float(np.median(covid)),
                    "hc_median_proportion": float(np.median(hc)),
                    "median_difference_COVID_minus_HC": float(
                        np.median(covid) - np.median(hc)
                    ),
                    "mannwhitney_u": float(result.statistic),
                    "p_value_two_sided": float(result.pvalue),
                    "cliffs_delta_COVID_vs_HC": cliffs_delta(covid, hc),
                }
            )
    adjusted = bh_adjust([row["p_value_two_sided"] for row in test_rows])
    for row, value in zip(test_rows, adjusted):
        row["q_value_bh_across_callsets_and_bins"] = value
    return sample_rows, test_rows


def ols_hc3(y: np.ndarray, x: np.ndarray) -> dict:
    n, p = x.shape
    inverse = np.linalg.pinv(x.T @ x)
    beta = inverse @ x.T @ y
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, inverse, x)
    scaled = (residual / np.maximum(1 - leverage, np.finfo(float).eps)) ** 2
    covariance = inverse @ (x.T @ (x * scaled[:, None])) @ inverse
    standard_error = np.sqrt(np.maximum(np.diag(covariance), 0))
    dof = n - p
    t_value = beta / standard_error
    p_value = 2 * stats.t.sf(np.abs(t_value), dof)
    critical = stats.t.ppf(0.975, dof)
    centered = y - np.mean(y)
    return {
        "beta": beta,
        "se": standard_error,
        "t": t_value,
        "p": p_value,
        "low": beta - critical * standard_error,
        "high": beta + critical * standard_error,
        "n": n,
        "dof": dof,
        "r2": 1 - np.sum(residual**2) / np.sum(centered**2),
    }


def adjusted_regressions(metrics: list[dict], rca_path: Path) -> list[dict]:
    covariates = {row["sample_id"]: row for row in read_tsv(rca_path)}
    output = []
    for callset in [name for name in CALLSET_ORDER if any(row["callset"] == name for row in metrics)]:
        rows = sorted(
            [row for row in metrics if row["callset"] == callset],
            key=lambda row: row["sample_id"],
        )
        if len(rows) != 78 or any(row["sample_id"] not in covariates for row in rows):
            continue
        covid = np.asarray([row["group"] == "COVID-19" for row in rows], dtype=float)
        log_rca = np.log2(
            np.asarray(
                [float(covariates[row["sample_id"]]["rca_concentration_ng_ul"]) for row in rows]
            )
        )
        age = np.asarray(
            [float(covariates[row["sample_id"]]["age_years"]) for row in rows]
        )
        male = np.asarray(
            [covariates[row["sample_id"]]["sex"].lower() == "male" for row in rows],
            dtype=float,
        )
        x = np.column_stack(
            [
                np.ones(len(rows)),
                covid,
                log_rca - np.mean(log_rca),
                (age - np.mean(age)) / 10,
                male,
            ]
        )
        y = np.log2(np.asarray([float(row["epm"]) for row in rows]) + 1)
        result = ols_hc3(y, x)
        for index, term in enumerate(
            ("intercept", "COVID-19_vs_HC", "log2_RCA_per_doubling", "age_per_10_years", "male_vs_female")
        ):
            estimate = float(result["beta"][index])
            low = float(result["low"][index])
            high = float(result["high"][index])
            output.append(
                {
                    "callset": callset,
                    "outcome": "log2(EPM + 1)",
                    "model": "group_plus_log2RCA_age_sex_HC3",
                    "term": term,
                    "estimate_log2_scale": estimate,
                    "hc3_standard_error": float(result["se"][index]),
                    "ci_low_log2_scale": low,
                    "ci_high_log2_scale": high,
                    "multiplicative_effect": 2**estimate,
                    "multiplicative_ci_low": 2**low,
                    "multiplicative_ci_high": 2**high,
                    "t_value": float(result["t"][index]),
                    "p_value_two_sided": float(result["p"][index]),
                    "n": int(result["n"]),
                    "degrees_of_freedom": int(result["dof"]),
                    "r_squared": float(result["r2"]),
                }
            )
    indices = [i for i, row in enumerate(output) if row["term"] == "COVID-19_vs_HC"]
    adjusted = bh_adjust([output[index]["p_value_two_sided"] for index in indices])
    for index, value in zip(indices, adjusted):
        output[index]["q_value_bh_across_callsets_for_group_term"] = value
    return output


def threshold_statistics(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    rows = read_tsv(path)
    output = []
    combinations = sorted(
        {(int(row["split_reads_min"]), float(row["score_strictly_greater_than"])) for row in rows}
    )
    for split_min, score_gt in combinations:
        subset = [
            row
            for row in rows
            if int(row["split_reads_min"]) == split_min
            and float(row["score_strictly_greater_than"]) == score_gt
        ]
        covid = np.asarray([float(row["epm"]) for row in subset if normalize_group(row["group"]) == "COVID-19"])
        hc = np.asarray([float(row["epm"]) for row in subset if normalize_group(row["group"]) == "HC"])
        result = stats.mannwhitneyu(covid, hc, alternative="two-sided", method="asymptotic")
        output.append(
            {
                "split_reads_min": split_min,
                "score_strictly_greater_than": score_gt,
                "covid_median_epm": float(np.median(covid)),
                "hc_median_epm": float(np.median(hc)),
                "median_ratio_COVID_over_HC": float(np.median(covid) / np.median(hc)),
                "mannwhitney_u": float(result.statistic),
                "p_value_two_sided": float(result.pvalue),
                "cliffs_delta_COVID_vs_HC": cliffs_delta(covid, hc),
            }
        )
    adjusted = bh_adjust([row["p_value_two_sided"] for row in output])
    for row, value in zip(output, adjusted):
        row["q_value_bh_across_threshold_grid"] = value
    return output


def target_support(root: Path, metadata: list[dict[str, str]], callsets: list[str]) -> tuple[list[dict], list[dict]]:
    sample_rows = []
    for callset in callsets:
        directory = root / "callsets" / callset
        if not directory.is_dir():
            continue
        for metadata_row in metadata:
            sample = metadata_row["sample_id"]
            path = directory / f"{sample}_circle_site.bed"
            if not path.is_file():
                continue
            calls = []
            with path.open() as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    try:
                        calls.append((fields[0], int(fields[1]), int(fields[2]), fields))
                    except ValueError:
                        continue
            for tolerance in (10, 20):
                for target, (chrom, target_start, target_end) in TARGETS.items():
                    matches = [
                        call
                        for call in calls
                        if call[0] == chrom
                        and abs(call[1] - target_start) <= tolerance
                        and abs(call[2] - target_end) <= tolerance
                    ]
                    best = min(
                        matches,
                        key=lambda call: (
                            abs(call[1] - target_start) + abs(call[2] - target_end),
                            call[1],
                            call[2],
                        ),
                        default=None,
                    )
                    sample_rows.append(
                        {
                            "sample_id": sample,
                            "group": normalize_group(metadata_row["group"]),
                            "callset": callset,
                            "tolerance_bp": tolerance,
                            "target": target,
                            "target_chrom": chrom,
                            "target_start": target_start,
                            "target_end": target_end,
                            "supported": int(best is not None),
                            "matched_call_count": len(matches),
                            "best_match_chrom": best[0] if best else "",
                            "best_match_start": best[1] if best else "",
                            "best_match_end": best[2] if best else "",
                            "best_match_total_breakpoint_distance": (
                                abs(best[1] - target_start) + abs(best[2] - target_end)
                                if best
                                else ""
                            ),
                            "best_match_split_reads": (
                                best[3][4] if best and len(best[3]) > 5 else ""
                            ),
                            "best_match_score": (
                                best[3][5] if best and len(best[3]) > 5 else ""
                            ),
                        }
                    )
    summary_rows = []
    keys = sorted(
        {(row["callset"], row["tolerance_bp"], row["target"], row["group"]) for row in sample_rows}
    )
    for callset, tolerance, target, group in keys:
        subset = [
            row
            for row in sample_rows
            if row["callset"] == callset
            and row["tolerance_bp"] == tolerance
            and row["target"] == target
            and row["group"] == group
        ]
        summary_rows.append(
            {
                "callset": callset,
                "tolerance_bp": tolerance,
                "target": target,
                "group": group,
                "n_samples": len(subset),
                "supported_samples": sum(row["supported"] for row in subset),
                "support_fraction": sum(row["supported"] for row in subset) / len(subset),
                "matched_call_count": sum(row["matched_call_count"] for row in subset),
            }
        )
    return sample_rows, summary_rows


def configure_matplotlib() -> str:
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        family = "Arial"
    except ValueError:
        family = "Liberation Sans"
    mpl.rcParams.update(
        {
            "font.family": family,
            "font.size": 6,
            "axes.labelsize": 6,
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
    return family


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.16, 1.07, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")


def make_figure(
    root: Path,
    metrics: list[dict],
    tests: list[dict],
    fragment_tests: list[dict],
    regressions: list[dict],
    target_summary: list[dict],
) -> str:
    family = configure_matplotlib()
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    callsets = [name for name in CALLSET_ORDER if any(row["callset"] == name for row in metrics)]
    fig, axes = plt.subplots(2, 2, figsize=(183 / 25.4, 137 / 25.4))
    ax_a, ax_b, ax_c, ax_d = axes.flat

    positions = np.arange(len(callsets))
    for group, color, marker, offset in (
        ("COVID-19", ORANGE, "o", -0.08),
        ("HC", BLUE, "s", 0.08),
    ):
        medians = [
            np.median(
                [row["epm"] for row in metrics if row["callset"] == callset and row["group"] == group]
            )
            for callset in callsets
        ]
        ax_a.plot(
            positions + offset,
            medians,
            color=color,
            marker=marker,
            markersize=3.5,
            linewidth=1,
            label=group,
        )
    ax_a.set_yscale("log")
    ax_a.set_ylabel("Median circles per million mapped reads")
    ax_a.set_xticks(positions, [CALLSET_LABELS[name] for name in callsets])
    ax_a.legend(frameon=False, ncol=2, loc="upper right")
    ax_a.grid(axis="y", color=GREY1, linewidth=0.5)
    panel_label(ax_a, "a")

    group_effects = [
        row
        for row in regressions
        if row["term"] == "COVID-19_vs_HC" and row["callset"] in callsets
    ]
    effect_by_callset = {row["callset"]: row for row in group_effects}
    effect_callsets = [name for name in callsets if name in effect_by_callset]
    y = np.arange(len(effect_callsets))
    effects = np.asarray([effect_by_callset[name]["multiplicative_effect"] for name in effect_callsets])
    lows = np.asarray([effect_by_callset[name]["multiplicative_ci_low"] for name in effect_callsets])
    highs = np.asarray([effect_by_callset[name]["multiplicative_ci_high"] for name in effect_callsets])
    ax_b.errorbar(
        effects,
        y,
        xerr=np.vstack([effects - lows, highs - effects]),
        fmt="o",
        color=RED,
        ecolor=GREY4,
        elinewidth=0.8,
        capsize=2,
        markersize=3.5,
    )
    ax_b.axvline(1, color=GREY6, linewidth=0.6, linestyle="--")
    ax_b.set_xscale("log")
    ax_b.set_yticks(y, [CALLSET_LABELS[name].replace("\n", " ") for name in effect_callsets])
    ax_b.invert_yaxis()
    ax_b.set_xlabel("Adjusted COVID-19/HC multiplicative effect (95% CI)")
    ax_b.grid(axis="x", color=GREY1, linewidth=0.5)
    panel_label(ax_b, "b")

    fragment_callsets = [
        name for name in callsets if any(row["callset"] == name for row in fragment_tests)
    ]
    matrix = np.full((len(fragment_callsets), len(LENGTH_BINS)), np.nan)
    for i, callset in enumerate(fragment_callsets):
        for j, name in enumerate(LENGTH_BINS):
            candidates = [
                row
                for row in fragment_tests
                if row["callset"] == callset and row["length_bin"] == name
            ]
            if candidates:
                matrix[i, j] = candidates[0]["median_difference_COVID_minus_HC"] * 100
    maximum = max(1.0, float(np.nanmax(np.abs(matrix))))
    cmap = mpl.colors.LinearSegmentedColormap.from_list("nature_diverging", [BLUE, "#FFFFFF", RED])
    image = ax_c.imshow(matrix, aspect="auto", cmap=cmap, vmin=-maximum, vmax=maximum)
    ax_c.set_xticks(np.arange(len(LENGTH_BINS)), LENGTH_LABELS, rotation=35, ha="right")
    ax_c.set_yticks(
        np.arange(len(fragment_callsets)),
        [CALLSET_LABELS[name].replace("\n", " ") for name in fragment_callsets],
    )
    ax_c.set_xlabel("Circle length (bp)")
    colorbar = fig.colorbar(image, ax=ax_c, fraction=0.045, pad=0.03)
    colorbar.set_label("Median proportion difference (COVID-19 − HC, pp)")
    panel_label(ax_c, "c")

    target_callsets = [
        name
        for name in ("circlemap_methods", "circlemap_high_support_masked", "consensus_t10")
        if any(row["callset"] == name for row in target_summary)
    ]
    x = np.arange(len(TARGETS))
    total_width = 0.72
    series = [(callset, group) for callset in target_callsets for group in ("COVID-19", "HC")]
    colors = {
        "circlemap_methods": GREY4,
        "circlemap_high_support_masked": GREEN,
        "consensus_t10": PURPLE,
    }
    for index, (callset, group) in enumerate(series):
        values = []
        for target in TARGETS:
            candidates = [
                row
                for row in target_summary
                if row["callset"] == callset
                and row["group"] == group
                and row["target"] == target
                and row["tolerance_bp"] == 10
            ]
            values.append(candidates[0]["support_fraction"] * 100 if candidates else 0)
        width = total_width / max(1, len(series))
        offset = (index - (len(series) - 1) / 2) * width
        ax_d.bar(
            x + offset,
            values,
            width=width * 0.92,
            color=colors[callset],
            edgecolor="black",
            linewidth=0.35,
            hatch="" if group == "COVID-19" else "////",
            label=f"{CALLSET_LABELS[callset].replace(chr(10), ' ')}; {group}",
        )
    ax_d.set_xticks(x, list(TARGETS), rotation=20, ha="right")
    ax_d.set_ylabel("Samples supporting target (%)")
    ax_d.set_ylim(0, 100)
    ax_d.legend(frameon=False, fontsize=4.8, ncol=2, loc="upper right")
    ax_d.grid(axis="y", color=GREY1, linewidth=0.5)
    panel_label(ax_d, "d")

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for label in ax_a.get_xticklabels():
        label.set_rotation(18)
        label.set_ha("right")
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.15, top=0.98, wspace=0.50, hspace=0.50)
    stem = figures / "Figure_CircleMap_robustness"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    return family


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--rca-metadata", required=True, type=Path)
    args = parser.parse_args()
    tables = args.root / "tables"
    metadata_dir = args.root / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    rca_copy = metadata_dir / "RCA_sample_metrics.tsv"
    if args.rca_metadata.resolve() != rca_copy.resolve():
        rca_copy.write_bytes(args.rca_metadata.read_bytes())

    metadata = read_tsv(args.manifest)
    metrics = load_metrics(args.root, metadata)
    callsets = [name for name in CALLSET_ORDER if any(row["callset"] == name for row in metrics)]
    for callset in callsets:
        observed = [row for row in metrics if row["callset"] == callset]
        if len(observed) != 78:
            raise RuntimeError(f"Expected 78 rows for {callset}; observed {len(observed)}")
    write_tsv(tables / "robustness_callset_sample_metrics.tsv", metrics)
    summaries, tests = group_statistics(metrics)
    write_tsv(tables / "robustness_group_summary.tsv", summaries)
    write_tsv(tables / "robustness_group_tests.tsv", tests)
    fragment_samples, fragment_tests = fragment_statistics(metrics)
    write_tsv(tables / "robustness_fragment_size_sample_proportions.tsv", fragment_samples)
    write_tsv(tables / "robustness_fragment_size_group_tests.tsv", fragment_tests)
    regressions = adjusted_regressions(metrics, rca_copy)
    write_tsv(tables / "robustness_regression_HC3.tsv", regressions)
    threshold_rows = threshold_statistics(
        tables / "circlemap_threshold_grid_sample_metrics.tsv"
    )
    write_tsv(tables / "robustness_threshold_grid_tests.tsv", threshold_rows)
    target_samples, target_summary = target_support(args.root, metadata, callsets)
    write_tsv(tables / "validation_target_sample_support.tsv", target_samples)
    write_tsv(tables / "validation_target_support_summary.tsv", target_summary)
    family = make_figure(
        args.root, metrics, tests, fragment_tests, regressions, target_summary
    )

    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(metadata),
        "callsets": callsets,
        "primary_group_metric": "EPM (calls per 1,000,000 mapped alignments)",
        "unadjusted_test": "two-sided Mann-Whitney U with asymptotic tie correction",
        "effect_size": "Cliff's delta, COVID-19 minus HC",
        "median_ratio_ci": f"{BOOTSTRAPS}-replicate percentile bootstrap",
        "regression": "OLS on log2(EPM+1), adjusted for log2 RCA, age, sex; HC3 standard errors",
        "multiple_testing": "Benjamini-Hochberg within prespecified result families",
        "target_match": "same chromosome and both breakpoints within tolerance; primary 10 bp, sensitivity 20 bp",
        "random_seed": SEED,
        "font_family": family,
        "figure_width_mm": 183,
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.root / "core_robustness_run_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")
    with (args.root / "core_robustness_validation.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write("samples=78\n")
        handle.write(f"callsets={len(callsets)}\n")
        handle.write(f"sample_metric_rows={len(metrics)}\n")
        handle.write(f"group_test_rows={len(tests)}\n")
        handle.write(f"target_sample_rows={len(target_samples)}\n")


if __name__ == "__main__":
    main()
