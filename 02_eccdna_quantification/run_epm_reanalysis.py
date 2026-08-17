#!/usr/bin/env python3
"""Redo COVID-19 versus HC basic statistics for the eccDNA revision.

The script uses only the Python standard library so it can run on QDUH without
adding dependencies.  It consumes TSV matrices exported from the manuscript
supplementary workbook and, when BED directories are supplied, recalculates
total eccDNA counts and fragment-length bin proportions directly from the
final Circle-Map BED files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Optional


GROUP_COVID = "COVID-19"
GROUP_HC = "HC"

FRAGMENT_BINS = (
    ("lt_1kb", "<1 kb", lambda length: length < 1000),
    ("1_2kb", "1-2 kb", lambda length: 1000 <= length <= 2000),
    ("gt_2kb", ">2 kb", lambda length: length > 2000),
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key, "")) for key in fieldnames})


def format_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return f"{value:.10g}"
    return value


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def normalize_group(value: object) -> str:
    text = "" if value is None else str(value).strip()
    lowered = text.lower()
    if lowered in {"covid", "covid-19", "covid_19", "case", "cases"}:
        return GROUP_COVID
    if lowered in {"hc", "healthy", "healthy control", "healthy controls", "normal", "control", "controls"}:
        return GROUP_HC
    return text


def canonical_feature_label(value: object) -> str:
    text = "" if value is None else str(value).strip()
    lowered = text.lower()
    if lowered.startswith("chr"):
        suffix = text[3:]
        if suffix.lower() in {"x", "y", "m", "mt"}:
            return "chr" + suffix.upper()
        return "chr" + suffix
    if lowered.endswith(".th3.bed"):
        text = text[: -len(".th3.bed")]
    return text


def to_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.upper() in {"NA", "NAN", "NONE"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if math.isnan(number):
        return None
    return number


def quantile(values: list[float], q: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[int(position)]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def group_summary(values: list[float]) -> dict[str, float]:
    clean = [value for value in values if value is not None and not math.isnan(value)]
    if not clean:
        return {
            "n": 0,
            "median": math.nan,
            "q1": math.nan,
            "q3": math.nan,
            "iqr": math.nan,
            "mean": math.nan,
            "min": math.nan,
            "max": math.nan,
        }
    q1 = quantile(clean, 0.25)
    q3 = quantile(clean, 0.75)
    return {
        "n": len(clean),
        "median": statistics.median(clean),
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "mean": statistics.fmean(clean),
        "min": min(clean),
        "max": max(clean),
    }


def average_ranks(values: list[float]) -> tuple[list[float], list[int]]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    tie_sizes: list[int] = []
    idx = 0
    while idx < len(indexed):
        end = idx + 1
        while end < len(indexed) and indexed[end][1] == indexed[idx][1]:
            end += 1
        average_rank = (idx + 1 + end) / 2.0
        for pos in range(idx, end):
            ranks[indexed[pos][0]] = average_rank
        tie_sizes.append(end - idx)
        idx = end
    return ranks, tie_sizes


def mann_whitney_two_sided(x: list[float], y: list[float]) -> tuple[float, float, float]:
    """Return U for x, two-sided P value, and z score.

    The P value uses the normal approximation with tie correction.  The cohort
    sizes are 39 versus 39 for the primary manuscript comparisons, where this is
    the standard Wilcoxon rank-sum large-sample approach.
    """

    n1 = len(x)
    n2 = len(y)
    if n1 == 0 or n2 == 0:
        return math.nan, math.nan, math.nan
    combined = x + y
    ranks, tie_sizes = average_ranks(combined)
    rank_sum_x = sum(ranks[:n1])
    u_x = rank_sum_x - n1 * (n1 + 1) / 2.0
    total_n = n1 + n2
    mean_u = n1 * n2 / 2.0
    if total_n <= 1:
        return u_x, math.nan, math.nan
    tie_sum = sum(tie**3 - tie for tie in tie_sizes if tie > 1)
    variance = n1 * n2 / 12.0 * ((total_n + 1) - tie_sum / (total_n * (total_n - 1)))
    if variance <= 0:
        return u_x, 1.0, 0.0
    z = (u_x - mean_u) / math.sqrt(variance)
    p_value = math.erfc(abs(z) / math.sqrt(2.0))
    return u_x, min(max(p_value, 0.0), 1.0), z


def cliffs_delta(x: list[float], y: list[float]) -> float:
    if not x or not y:
        return math.nan
    greater = 0
    less = 0
    for x_value in x:
        for y_value in y:
            if x_value > y_value:
                greater += 1
            elif x_value < y_value:
                less += 1
    return (greater - less) / (len(x) * len(y))


def bh_adjust(p_values: list[Optional[float]]) -> list[float]:
    indexed = [(idx, p) for idx, p in enumerate(p_values) if p is not None and not math.isnan(p)]
    adjusted = [math.nan] * len(p_values)
    if not indexed:
        return adjusted
    indexed.sort(key=lambda item: item[1])
    m = len(indexed)
    running = 1.0
    for rank_from_end, (idx, p_value) in enumerate(reversed(indexed), start=1):
        rank = m - rank_from_end + 1
        candidate = p_value * m / rank
        running = min(running, candidate)
        adjusted[idx] = min(running, 1.0)
    return adjusted


def clean_numeric_pairs(rows: Iterable[dict[str, object]], value_key: str) -> dict[str, list[float]]:
    values = {GROUP_COVID: [], GROUP_HC: []}
    for row in rows:
        group = normalize_group(row.get("group"))
        value = to_float(row.get(value_key))
        if group in values and value is not None:
            values[group].append(value)
    return values


def build_comparison(
    *,
    family: str,
    fdr_family: str,
    source: str,
    metric: str,
    feature: str,
    unit: str,
    covid_values: list[float],
    hc_values: list[float],
    transform: str = "none",
    notes: str = "",
) -> dict[str, object]:
    covid_summary = group_summary(covid_values)
    hc_summary = group_summary(hc_values)
    u_value, p_value, z_score = mann_whitney_two_sided(covid_values, hc_values)
    delta = cliffs_delta(covid_values, hc_values)
    median_difference = covid_summary["median"] - hc_summary["median"]
    if math.isnan(delta):
        direction = "NA"
    elif delta > 0:
        direction = "COVID-19 higher"
    elif delta < 0:
        direction = "COVID-19 lower"
    else:
        direction = "no stochastic shift"
    return {
        "family": family,
        "fdr_family": fdr_family,
        "source": source,
        "metric": metric,
        "feature": feature,
        "unit": unit,
        "transform": transform,
        "test": "two-sided Wilcoxon rank-sum (Mann-Whitney U, tie-corrected normal approximation)",
        "n_covid": covid_summary["n"],
        "n_hc": hc_summary["n"],
        "covid_median": covid_summary["median"],
        "covid_q1": covid_summary["q1"],
        "covid_q3": covid_summary["q3"],
        "covid_iqr": covid_summary["iqr"],
        "hc_median": hc_summary["median"],
        "hc_q1": hc_summary["q1"],
        "hc_q3": hc_summary["q3"],
        "hc_iqr": hc_summary["iqr"],
        "median_difference_covid_minus_hc": median_difference,
        "cliffs_delta": delta,
        "rank_biserial_correlation": delta,
        "mann_whitney_u_covid": u_value,
        "z_score": z_score,
        "p_value": p_value,
        "p_adj_bh": math.nan,
        "direction": direction,
        "notes": notes,
    }


def apply_fdr(comparisons: list[dict[str, object]]) -> None:
    by_family: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(comparisons):
        by_family[str(row["fdr_family"])].append(idx)
    for _, indices in by_family.items():
        adjusted = bh_adjust([comparisons[idx]["p_value"] for idx in indices])
        for idx, p_adj in zip(indices, adjusted):
            comparisons[idx]["p_adj_bh"] = p_adj


def bed_path_for_sample(sample_id: str, group: str, covid_bed_dir: Optional[Path], hc_bed_dir: Optional[Path]) -> Optional[Path]:
    if group == GROUP_COVID and covid_bed_dir is not None:
        return covid_bed_dir / f"{sample_id}_circle_site.bed"
    if group == GROUP_HC and hc_bed_dir is not None:
        return hc_bed_dir / f"{sample_id}_circle_site.bed"
    return None


def iter_bed_lengths(path: Path) -> tuple[list[int], int, int]:
    lengths: list[int] = []
    skipped_header_or_invalid = 0
    empty_lines = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip():
                empty_lines += 1
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                skipped_header_or_invalid += 1
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                skipped_header_or_invalid += 1
                continue
            length = end - start
            if length < 0:
                skipped_header_or_invalid += 1
                continue
            lengths.append(length)
    return lengths, skipped_header_or_invalid, empty_lines


def recalculate_bed_counts_and_fragments(
    burden_rows: list[dict[str, str]],
    covid_bed_dir: Optional[Path],
    hc_bed_dir: Optional[Path],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    sample_rows: list[dict[str, object]] = []
    fragment_rows: list[dict[str, object]] = []
    warnings: list[dict[str, object]] = []
    for row in burden_rows:
        sample_id = str(row.get("sample_id", "")).strip()
        group = normalize_group(row.get("group"))
        mapped_reads = to_float(row.get("mapped_reads"))
        input_count = to_float(row.get("eccdna_count"))
        input_epm = to_float(row.get("epm"))
        bed_path = bed_path_for_sample(sample_id, group, covid_bed_dir, hc_bed_dir)
        computed_count: Optional[int] = None
        skipped = 0
        empty = 0
        if bed_path is not None and bed_path.exists():
            lengths, skipped, empty = iter_bed_lengths(bed_path)
            computed_count = len(lengths)
            bin_counts = {bin_id: 0 for bin_id, _, _ in FRAGMENT_BINS}
            for length in lengths:
                for bin_id, _, predicate in FRAGMENT_BINS:
                    if predicate(length):
                        bin_counts[bin_id] += 1
                        break
            denominator = computed_count if computed_count else 0
            for bin_id, bin_label, _ in FRAGMENT_BINS:
                count = bin_counts[bin_id]
                fragment_rows.append(
                    {
                        "sample_id": sample_id,
                        "group": group,
                        "length_bin": bin_id,
                        "length_bin_label": bin_label,
                        "fragment_count": count,
                        "total_valid_bed_rows": denominator,
                        "fragment_proportion": count / denominator if denominator else math.nan,
                    }
                )
        else:
            warnings.append(
                {
                    "sample_id": sample_id,
                    "group": group,
                    "issue": "missing_bed",
                    "input_count": input_count,
                    "computed_count": "",
                    "bed_path": str(bed_path) if bed_path is not None else "",
                }
            )
        final_count = computed_count if computed_count is not None else input_count
        final_epm = final_count / mapped_reads * 1_000_000 if final_count is not None and mapped_reads else input_epm
        if computed_count is not None and input_count is not None and int(input_count) != computed_count:
            warnings.append(
                {
                    "sample_id": sample_id,
                    "group": group,
                    "issue": "bed_count_differs_from_input_table",
                    "input_count": int(input_count),
                    "computed_count": computed_count,
                    "bed_path": str(bed_path),
                }
            )
        sample_rows.append(
            {
                "sample_id": sample_id,
                "group": group,
                "mapped_reads": mapped_reads,
                "eccdna_count_input_table": input_count,
                "epm_input_table": input_epm,
                "eccdna_count_qduh_bed": final_count,
                "epm_qduh_bed": final_epm,
                "log1p_epm_qduh_bed": math.log1p(final_epm) if final_epm is not None else math.nan,
                "bed_header_or_invalid_lines_skipped": skipped,
                "bed_empty_lines_skipped": empty,
                "bed_path": str(bed_path) if bed_path is not None else "",
            }
        )
    return sample_rows, fragment_rows, warnings


def add_metric_comparison(
    comparisons: list[dict[str, object]],
    *,
    rows: Iterable[dict[str, object]],
    value_key: str,
    family: str,
    fdr_family: str,
    source: str,
    metric: str,
    feature: str,
    unit: str,
    transform: str = "none",
    transform_function: Optional[Callable[[float], float]] = None,
    notes: str = "",
) -> None:
    grouped = {GROUP_COVID: [], GROUP_HC: []}
    for row in rows:
        group = normalize_group(row.get("group"))
        value = to_float(row.get(value_key))
        if group in grouped and value is not None:
            if transform_function is not None:
                value = transform_function(value)
            grouped[group].append(value)
    if not grouped[GROUP_COVID] and not grouped[GROUP_HC]:
        return
    comparisons.append(
        build_comparison(
            family=family,
            fdr_family=fdr_family,
            source=source,
            metric=metric,
            feature=feature,
            unit=unit,
            transform=transform,
            covid_values=grouped[GROUP_COVID],
            hc_values=grouped[GROUP_HC],
            notes=notes,
        )
    )


def add_feature_comparisons(
    comparisons: list[dict[str, object]],
    *,
    rows: list[dict[str, str]],
    feature_key: str,
    value_key: str,
    family: str,
    fdr_family: str,
    source: str,
    metric: str,
    unit: str,
    transform: str = "none",
    transform_function: Optional[Callable[[float], float]] = None,
    notes: str = "",
) -> None:
    by_feature: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        feature = canonical_feature_label(row.get(feature_key, ""))
        if feature:
            by_feature[feature].append(row)
    for feature in sorted(by_feature):
        add_metric_comparison(
            comparisons,
            rows=by_feature[feature],
            value_key=value_key,
            family=family,
            fdr_family=fdr_family,
            source=source,
            metric=metric,
            feature=feature,
            unit=unit,
            transform=transform,
            transform_function=transform_function,
            notes=notes,
        )


def write_analysis_plan(path: Path) -> None:
    rows = [
        {
            "manuscript_location": "Figure 1B; Table S2",
            "reported_metric": "Per-sample total eccDNA count",
            "redo_action": "Two-sided Wilcoxon rank-sum; median/IQR; Cliff's delta; BH within sample_burden family",
            "reviewer_link": "R1 Major 9; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 1C; Table S2",
            "reported_metric": "Total EPM",
            "redo_action": "Two-sided Wilcoxon rank-sum on raw EPM plus log1p(EPM) sensitivity; median/IQR; Cliff's delta",
            "reviewer_link": "R1 Major 9-10; R2 Minor 1-2; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 1E",
            "reported_metric": "Fragment-length bin proportions (<1 kb, 1-2 kb, >2 kb)",
            "redo_action": "Recomputed from QDUH BED intervals; two-sided Wilcoxon rank-sum by bin; BH within fragment bin family",
            "reviewer_link": "R2 Minor 1-2; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 2B; Table S3",
            "reported_metric": "Chromosome-level normalized distribution",
            "redo_action": "Two-sided Wilcoxon rank-sum per chromosome; BH-FDR within chromosome family; median/IQR; Cliff's delta",
            "reviewer_link": "R1 Major 9; R2 Minor 1-2; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 2C; Table S4",
            "reported_metric": "Repeat-class normalized mapping ratio",
            "redo_action": "Two-sided Wilcoxon rank-sum per repeat class; BH-FDR within repeat-class family; median/IQR; Cliff's delta",
            "reviewer_link": "R1 Major 9; R2 Minor 1-2; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 2D; Table S5/Table S6",
            "reported_metric": "Gene-feature ratios and enrichment scores",
            "redo_action": "Two-sided Wilcoxon rank-sum per gene feature; BH-FDR within each gene-feature metric family; median/IQR; Cliff's delta",
            "reviewer_link": "R1 Major 9; R2 Minor 1-2; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 4B; Table S8",
            "reported_metric": "ENCODE/Broad CD14+ histone-mark overlap EPM and enrichment metrics",
            "redo_action": "Two-sided Wilcoxon rank-sum per mark; BH-FDR within histone-mark metric families; raw EPM plus log1p(EPM) sensitivity; median/IQR; Cliff's delta",
            "reviewer_link": "R1 Major 9-10; R2 Major 2; R2 Minor 1-2; R3 Minor 2",
        },
        {
            "manuscript_location": "Figure 5B; Table S9",
            "reported_metric": "SARS-CoV-2 infected-cell H3K27ac/H3K27me3 overlap EPM and enrichment metrics",
            "redo_action": "Two-sided Wilcoxon rank-sum; raw EPM plus log1p(EPM) sensitivity; median/IQR; Cliff's delta",
            "reviewer_link": "R1 Major 9-10; R2 Major 2; R2 Minor 1-2; R3 Minor 2",
        },
    ]
    write_tsv(path, ["manuscript_location", "reported_metric", "redo_action", "reviewer_link"], rows)


def format_p(value: object) -> str:
    number = to_float(value)
    if number is None or math.isnan(number):
        return "NA"
    if number < 1e-4:
        return f"{number:.3e}"
    return f"{number:.4g}"


def write_summary(path: Path, comparisons: list[dict[str, object]]) -> None:
    by_key = {(row["family"], row["metric"], row["feature"], row["transform"]): row for row in comparisons}
    lines = [
        "COVID-19 versus HC EPM/basic-statistics reanalysis",
        f"Generated: {now_utc()}",
        "",
        "Primary burden comparisons:",
    ]
    for key in [
        ("sample_burden", "total_eccDNA_count", "all_samples", "none"),
        ("sample_burden", "total_EPM", "all_samples", "none"),
        ("sample_burden", "total_EPM", "all_samples", "log1p"),
    ]:
        row = by_key.get(key)
        if row:
            lines.append(
                "- {metric} ({transform}): COVID median {covid_median} [IQR {covid_iqr}] vs HC median {hc_median} [IQR {hc_iqr}], "
                "Cliff's delta {delta}, FDR {fdr}".format(
                    metric=row["metric"],
                    transform=row["transform"],
                    covid_median=format_p(row["covid_median"]),
                    covid_iqr=format_p(row["covid_iqr"]),
                    hc_median=format_p(row["hc_median"]),
                    hc_iqr=format_p(row["hc_iqr"]),
                    delta=format_p(row["cliffs_delta"]),
                    fdr=format_p(row["p_adj_bh"]),
                )
            )
    lines.extend(["", "FDR-significant feature counts by family (q < 0.05):"])
    family_counts: dict[str, tuple[int, int]] = {}
    for row in comparisons:
        family = str(row["fdr_family"])
        total, sig = family_counts.get(family, (0, 0))
        p_adj = to_float(row.get("p_adj_bh"))
        total += 1
        if p_adj is not None and p_adj < 0.05:
            sig += 1
        family_counts[family] = (total, sig)
    for family in sorted(family_counts):
        total, sig = family_counts[family]
        lines.append(f"- {family}: {sig}/{total}")
    path.write_text("\n".join(lines) + "\n")


def sample_id_from_eccbed(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if "/" in text:
        text = text.rsplit("/", 1)[-1]
    for suffix in ["_circle_site.bed", "_circle_site", ".bed", ".hg19.flt.bed", ".hg19.flt"]:
        if text.endswith(suffix):
            text = text[: -len(suffix)]
    return text


def read_epm_rows(path: Path, group: str, mark: str) -> dict[tuple[str, str], dict[str, str]]:
    rows = {}
    if not path.exists():
        return rows
    for row in read_tsv(path):
        sample_id = str(row.get("id", row.get("ID", ""))).strip()
        if not sample_id:
            sample_id = sample_id_from_eccbed(row.get("eccbed"))
        rows[(sample_id, mark)] = {
            "sample": sample_id,
            "group": group,
            "mark": mark,
            "mapped_reads_in_peaks": row.get("mapped_reads_in_peaks", ""),
            "eccdna_counts_in_peaks": row.get("eccdna_counts_in_peaks", ""),
            "epm": row.get("epm", row.get("EPM", "")),
        }
    return rows


def load_sars_cov2_histone_rows(input_dir: Path) -> list[dict[str, str]]:
    """Prefer QDUH 20251022 H3K27ac/H3K27me3 files when present.

    The manuscript workbook copy available locally contains only the H3K27ac
    rows in Table S9, whereas the QDUH 20251022 enrichment output used for
    Figure 5 contains both H3K27ac and H3K27me3.  If the QDUH-derived files are
    present in the input directory, this function merges them into the same
    column convention as Table S9.
    """

    enrichment_path = input_dir / "sars_cov2_20251022_eccDNA_enrichment.tsv"
    if not enrichment_path.exists():
        return read_tsv(input_dir / "histone_sars_cov2_overlap.tsv")

    epm_by_sample_mark: dict[tuple[str, str], dict[str, str]] = {}
    epm_inputs = [
        ("sars_cov2_H3K27ac_EPM_covid.tsv", GROUP_COVID, "H3K27ac"),
        ("sars_cov2_H3K27ac_EPM_hc.tsv", GROUP_HC, "H3K27ac"),
        ("sars_cov2_H3K27me3_EPM_covid.tsv", GROUP_COVID, "H3K27me3"),
        ("sars_cov2_H3K27me3_EPM_hc.tsv", GROUP_HC, "H3K27me3"),
    ]
    for filename, group, mark in epm_inputs:
        epm_by_sample_mark.update(read_epm_rows(input_dir / filename, group, mark))

    merged_rows: list[dict[str, str]] = []
    for row in read_tsv(enrichment_path):
        sample_id = sample_id_from_eccbed(row.get("sample"))
        group = normalize_group(row.get("group"))
        mark = canonical_feature_label(row.get("mark"))
        epm_row = epm_by_sample_mark.get((sample_id, mark), {})
        merged_rows.append(
            {
                "sample": sample_id,
                "group": group,
                "mark": mark,
                "total_eccdna_n": row.get("total_eccdna", ""),
                "observed_overlap_o": row.get("obs_overlap", ""),
                "background_total_b": row.get("bg_total", ""),
                "background_overlap_e": row.get("bg_overlap", ""),
                "observed_fraction": row.get("p_obs", ""),
                "expected_fraction": row.get("p_exp", ""),
                "enrichment_score": row.get("es_log2", row.get("ES_log2", "")),
                "mapped_reads_in_peaks": epm_row.get("mapped_reads_in_peaks", ""),
                "eccdna_counts_in_peaks": epm_row.get("eccdna_counts_in_peaks", ""),
                "epm": epm_row.get("epm", ""),
            }
        )
    return merged_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--covid-bed-dir", type=Path)
    parser.add_argument("--hc-bed-dir", type=Path)
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    results_dir = output_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    burden_rows = read_tsv(input_dir / "sample_burden.tsv")
    sample_rows, fragment_rows, warnings = recalculate_bed_counts_and_fragments(
        burden_rows, args.covid_bed_dir, args.hc_bed_dir
    )

    comparisons: list[dict[str, object]] = []
    add_metric_comparison(
        comparisons,
        rows=sample_rows,
        value_key="eccdna_count_qduh_bed",
        family="sample_burden",
        fdr_family="sample_burden",
        source="Table_S2 plus QDUH BED recount",
        metric="total_eccDNA_count",
        feature="all_samples",
        unit="count",
        notes="BED biological data rows were recounted on QDUH when BED files were available.",
    )
    add_metric_comparison(
        comparisons,
        rows=sample_rows,
        value_key="epm_qduh_bed",
        family="sample_burden",
        fdr_family="sample_burden",
        source="Table_S2 mapped reads plus QDUH BED recount",
        metric="total_EPM",
        feature="all_samples",
        unit="eccDNA per million mapped reads",
        notes="EPM = valid BED row count / mapped reads x 1e6.",
    )
    add_metric_comparison(
        comparisons,
        rows=sample_rows,
        value_key="epm_qduh_bed",
        family="sample_burden",
        fdr_family="sample_burden",
        source="Table_S2 mapped reads plus QDUH BED recount",
        metric="total_EPM",
        feature="all_samples",
        unit="log1p(EPM)",
        transform="log1p",
        transform_function=math.log1p,
        notes="Sensitivity analysis requested by Reviewer 1 Major Comment 10.",
    )

    chromosome_rows = read_tsv(input_dir / "chromosome_distribution.tsv")
    add_feature_comparisons(
        comparisons,
        rows=chromosome_rows,
        feature_key="chromosome",
        value_key="normalized_fraction_per_mb",
        family="chromosome_distribution",
        fdr_family="chromosome_distribution",
        source="Table_S3",
        metric="chromosome_normalized_fraction_per_mb",
        unit="fraction per Mb",
    )

    repeat_rows = read_tsv(input_dir / "repeat_class_enrichment.tsv")
    add_feature_comparisons(
        comparisons,
        rows=repeat_rows,
        feature_key="repeat_class",
        value_key="normalized_mapping_ratio",
        family="repeat_class_enrichment",
        fdr_family="repeat_class_enrichment",
        source="Table_S4",
        metric="repeat_normalized_mapping_ratio",
        unit="observed/genomic normalized mapping ratio",
    )

    gene_feature_rows = read_tsv(input_dir / "gene_feature_characteristics.tsv")
    # Table S5 contains a malformed Group column: 39 COVID-19 labels followed
    # by 39 HC labels rather than one group label per sample-feature row.  Use
    # the validated Table S2 sample-to-group mapping for all Table S5 rows.
    group_by_sample = {
        str(row.get("sample_id", "")).strip(): normalize_group(row.get("group"))
        for row in burden_rows
    }
    for row in gene_feature_rows:
        sample_id = str(row.get("sample_id", "")).strip()
        if sample_id in group_by_sample:
            row["group"] = group_by_sample[sample_id]
    add_feature_comparisons(
        comparisons,
        rows=gene_feature_rows,
        feature_key="element",
        value_key="eccdna_element_ratio",
        family="gene_feature",
        fdr_family="gene_feature_element_ratio",
        source="Table_S5",
        metric="gene_feature_eccDNA_element_ratio",
        unit="proportion",
    )
    add_feature_comparisons(
        comparisons,
        rows=gene_feature_rows,
        feature_key="element",
        value_key="nomalized_eccdna_ratio",
        family="gene_feature",
        fdr_family="gene_feature_normalized_ratio",
        source="Table_S5",
        metric="gene_feature_normalized_eccDNA_ratio",
        unit="observed/expected ratio",
        notes="Column name in the source workbook is misspelled as Nomalized EccDNA ratio.",
    )

    gene_enrichment_rows = read_tsv(input_dir / "gene_feature_enrichment_scores.tsv")
    add_feature_comparisons(
        comparisons,
        rows=gene_enrichment_rows,
        feature_key="element",
        value_key="enrichment_score",
        family="gene_feature",
        fdr_family="gene_feature_enrichment_score",
        source="Table_S6",
        metric="gene_feature_enrichment_score",
        unit="enrichment score",
    )

    histone_encode_rows = read_tsv(input_dir / "histone_encode_overlap.tsv")
    for value_key, metric, unit, transform, transform_function, fdr_family in [
        ("epm", "histone_mark_EPM", "eccDNA per million mapped reads in peaks", "none", None, "histone_encode_EPM"),
        ("epm", "histone_mark_EPM", "log1p(EPM)", "log1p", math.log1p, "histone_encode_log1p_EPM"),
        ("enrichment_score", "histone_mark_enrichment_score", "enrichment score", "none", None, "histone_encode_enrichment_score"),
        ("observed_fraction", "histone_mark_observed_fraction", "observed fraction", "none", None, "histone_encode_observed_fraction"),
    ]:
        add_feature_comparisons(
            comparisons,
            rows=histone_encode_rows,
            feature_key="mark",
            value_key=value_key,
            family="histone_mark_encode_cd14",
            fdr_family=fdr_family,
            source="Table_S8",
            metric=metric,
            unit=unit,
            transform=transform,
            transform_function=transform_function,
        )

    histone_sars_rows = load_sars_cov2_histone_rows(input_dir)
    for value_key, metric, unit, transform, transform_function, fdr_family in [
        ("epm", "sars_cov2_histone_EPM", "eccDNA per million mapped reads in peaks", "none", None, "histone_sars_cov2_EPM"),
        ("epm", "sars_cov2_histone_EPM", "log1p(EPM)", "log1p", math.log1p, "histone_sars_cov2_log1p_EPM"),
        ("enrichment_score", "sars_cov2_histone_enrichment_score", "enrichment score", "none", None, "histone_sars_cov2_enrichment_score"),
        ("observed_fraction", "sars_cov2_histone_observed_fraction", "observed fraction", "none", None, "histone_sars_cov2_observed_fraction"),
    ]:
        add_feature_comparisons(
            comparisons,
            rows=histone_sars_rows,
            feature_key="mark",
            value_key=value_key,
            family="histone_mark_sars_cov2",
            fdr_family=fdr_family,
            source="Table_S9",
            metric=metric,
            unit=unit,
            transform=transform,
            transform_function=transform_function,
        )

    if fragment_rows:
        add_feature_comparisons(
            comparisons,
            rows=[{key: str(value) for key, value in row.items()} for row in fragment_rows],
            feature_key="length_bin_label",
            value_key="fragment_proportion",
            family="fragment_length_bin",
            fdr_family="fragment_length_bin_proportion",
            source="QDUH BED intervals",
            metric="fragment_length_bin_proportion",
            unit="proportion of valid BED rows",
            notes="BED interval length = end - start; bins are <1 kb, 1-2 kb inclusive, and >2 kb.",
        )

    apply_fdr(comparisons)

    comparison_fields = [
        "family",
        "fdr_family",
        "source",
        "metric",
        "feature",
        "unit",
        "transform",
        "test",
        "n_covid",
        "n_hc",
        "covid_median",
        "covid_q1",
        "covid_q3",
        "covid_iqr",
        "hc_median",
        "hc_q1",
        "hc_q3",
        "hc_iqr",
        "median_difference_covid_minus_hc",
        "cliffs_delta",
        "rank_biserial_correlation",
        "mann_whitney_u_covid",
        "z_score",
        "p_value",
        "p_adj_bh",
        "direction",
        "notes",
    ]
    write_tsv(results_dir / "all_group_comparisons.tsv", comparison_fields, comparisons)

    significant_rows = [
        row for row in comparisons if (to_float(row.get("p_adj_bh")) is not None and to_float(row.get("p_adj_bh")) < 0.05)
    ]
    write_tsv(results_dir / "significant_group_comparisons_q005.tsv", comparison_fields, significant_rows)

    write_tsv(
        results_dir / "sample_burden_recomputed_from_qduh_bed.tsv",
        [
            "sample_id",
            "group",
            "mapped_reads",
            "eccdna_count_input_table",
            "epm_input_table",
            "eccdna_count_qduh_bed",
            "epm_qduh_bed",
            "log1p_epm_qduh_bed",
            "bed_header_or_invalid_lines_skipped",
            "bed_empty_lines_skipped",
            "bed_path",
        ],
        sample_rows,
    )
    if fragment_rows:
        write_tsv(
            results_dir / "fragment_length_bin_proportions.tsv",
            [
                "sample_id",
                "group",
                "length_bin",
                "length_bin_label",
                "fragment_count",
                "total_valid_bed_rows",
                "fragment_proportion",
            ],
            fragment_rows,
        )
    write_tsv(
        results_dir / "bed_count_warnings.tsv",
        ["sample_id", "group", "issue", "input_count", "computed_count", "bed_path"],
        warnings,
    )
    write_analysis_plan(results_dir / "statistics_to_redo_from_manuscript.tsv")
    write_summary(results_dir / "analysis_summary.txt", comparisons)

    input_manifest = []
    for path in sorted(input_dir.glob("*.tsv")):
        stat = path.stat()
        input_manifest.append(
            {
                "path": str(path),
                "size_bytes": stat.st_size,
                "mtime_epoch": int(stat.st_mtime),
            }
        )
    write_tsv(results_dir / "input_manifest.tsv", ["path", "size_bytes", "mtime_epoch"], input_manifest)

    write_tsv(
        results_dir / "software_versions.tsv",
        ["name", "version_or_value"],
        [
            {"name": "python", "version_or_value": sys.version.replace("\n", " ")},
            {"name": "platform", "version_or_value": platform.platform()},
            {"name": "script", "version_or_value": str(Path(__file__).resolve())},
            {"name": "run_utc", "version_or_value": now_utc()},
            {"name": "covid_bed_dir", "version_or_value": str(args.covid_bed_dir) if args.covid_bed_dir else ""},
            {"name": "hc_bed_dir", "version_or_value": str(args.hc_bed_dir) if args.hc_bed_dir else ""},
        ],
    )
    run_parameters = {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "covid_bed_dir": str(args.covid_bed_dir) if args.covid_bed_dir else None,
        "hc_bed_dir": str(args.hc_bed_dir) if args.hc_bed_dir else None,
        "fragment_bins": [{"id": bin_id, "label": label} for bin_id, label, _ in FRAGMENT_BINS],
        "statistical_test": "two-sided Wilcoxon rank-sum / Mann-Whitney U with tie-corrected normal approximation",
        "effect_size": "Cliff's delta; positive values indicate COVID-19 > HC",
        "multiple_testing": "Benjamini-Hochberg FDR within each fdr_family",
        "generated_utc": now_utc(),
    }
    (results_dir / "run_parameters.json").write_text(json.dumps(run_parameters, indent=2, ensure_ascii=False) + "\n")

    validation_lines = [
        f"Generated: {now_utc()}",
        f"PASS: read {len(burden_rows)} sample-burden rows",
        f"PASS: produced {len(comparisons)} group-comparison rows",
        "PASS: all group comparisons use two-sided Wilcoxon rank-sum tests",
        "PASS: BH-FDR was applied within every fdr_family",
    ]
    group_counts = defaultdict(int)
    for row in sample_rows:
        group_counts[normalize_group(row.get("group"))] += 1
    validation_lines.append(f"PASS: sample groups are COVID-19={group_counts[GROUP_COVID]}, HC={group_counts[GROUP_HC]}")
    if fragment_rows:
        validation_lines.append(f"PASS: computed fragment-length bins for {len(fragment_rows) // len(FRAGMENT_BINS)} samples")
    else:
        validation_lines.append("WARN: fragment-length bins were not computed because BED files were unavailable")
    if warnings:
        validation_lines.append(f"WARN: {len(warnings)} BED count/path warning(s); see bed_count_warnings.tsv")
    else:
        validation_lines.append("PASS: no BED count/path warnings")
    (results_dir / "validation_report.txt").write_text("\n".join(validation_lines) + "\n")


if __name__ == "__main__":
    main()
