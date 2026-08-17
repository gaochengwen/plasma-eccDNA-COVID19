#!/usr/bin/env python3
"""Objective, sample-weighted identification of eccDNA fragment-length peaks.

The analysis is intentionally discovery based: neither the number nor the
positions of peaks are supplied to the algorithm.  Unique Circle-Map intervals
are the observational units within each biological sample, while samples are
the bootstrap and group-averaging units.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import scipy
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import mannwhitneyu


GROUP_COVID = "COVID-19"
GROUP_HC = "HC"
GROUPS = (GROUP_COVID, GROUP_HC)
DATASETS = ("Overall", GROUP_COVID, GROUP_HC)

# Exact Nature scientific-illustration palette values.
COVID_COLOR = "#C93E3F"
COVID_LIGHT = "#FAD0CE"
HC_COLOR = "#0272B2"
HC_LIGHT = "#C8E7FB"
OVERALL_COLOR = "#49566D"
GREY = "#6F7B91"
LIGHT_GREY = "#C8CEDA"
BLACK = "#253247"

# Exact box/point colours extracted from the revised Figure 4 reference panels.
FIG4_COVID_BOX = "#FF8080"
FIG4_HC_BOX = "#8BABD3"
FIG4_COVID_POINT = "#D70000"
FIG4_HC_POINT = "#034E61"

DEFAULT_SEED = 20260729


@dataclass
class Sample:
    sample_id: str
    group: str
    bed_path: Path
    lengths: np.ndarray
    audit: dict[str, object]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--covid-bed-dir", required=True, type=Path)
    parser.add_argument("--hc-bed-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--subsample-repeats", type=int, default=100)
    parser.add_argument("--subsample-count", type=int, default=4000)
    parser.add_argument("--primary-low", type=int, default=50)
    parser.add_argument("--primary-high", type=int, default=1000)
    parser.add_argument("--sensitivity-high", type=int, default=2000)
    parser.add_argument("--primary-sigma", type=float, default=15.0)
    parser.add_argument("--sigmas", default="10,15,20,30")
    parser.add_argument("--search-low", type=int, default=100)
    parser.add_argument("--search-high", type=int, default=900)
    parser.add_argument("--peak-distance", type=int, default=100)
    parser.add_argument("--prominence-fraction", type=float, default=0.05)
    parser.add_argument("--minimum-width", type=float, default=10.0)
    parser.add_argument("--match-window", type=int, default=25)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def normalize_group(value: str) -> str:
    text = value.strip()
    lowered = text.lower()
    if lowered in {"covid", "covid-19", "covid_19", "case", "cases"}:
        return GROUP_COVID
    if lowered in {
        "hc",
        "healthy",
        "healthy control",
        "healthy controls",
        "normal",
        "control",
        "controls",
    }:
        return GROUP_HC
    raise ValueError(f"Unrecognized study group: {value!r}")


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if math.isnan(number):
            return "NA"
        if math.isinf(number):
            return "Inf" if number > 0 else "-Inf"
        return f"{number:.10g}"
    if isinstance(value, (np.bool_, bool)):
        return "TRUE" if bool(value) else "FALSE"
    return str(value)


def write_tsv(
    path: Path, rows: Iterable[dict[str, object]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        delimiter = "\t" if "\t" in first else ","
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header in metadata file: {path}")
        canonical = {
            field.strip().lower().replace(" ", "_"): field
            for field in reader.fieldnames
            if field is not None
        }
        sample_field = canonical.get("sample_id") or canonical.get("sample")
        group_field = canonical.get("group")
        if sample_field is None or group_field is None:
            raise ValueError("Metadata must contain sample_id and group columns")
        rows = []
        for row in reader:
            sample_id = (row.get(sample_field) or "").strip()
            group_text = (row.get(group_field) or "").strip()
            if sample_id:
                rows.append(
                    {"sample_id": sample_id, "group": normalize_group(group_text)}
                )

    sample_ids = [row["sample_id"] for row in rows]
    counts = Counter(row["group"] for row in rows)
    if len(rows) != 78 or len(set(sample_ids)) != 78:
        raise ValueError(
            f"Expected 78 unique samples; found {len(rows)} rows and "
            f"{len(set(sample_ids))} unique IDs"
        )
    if counts != Counter({GROUP_COVID: 39, GROUP_HC: 39}):
        raise ValueError(f"Expected 39 samples per group; found {dict(counts)}")
    return rows


def expected_bed_path(
    sample_id: str, group: str, covid_dir: Path, hc_dir: Path
) -> Path:
    directory = covid_dir if group == GROUP_COVID else hc_dir
    return directory / f"{sample_id}_circle_site.bed"


def validate_bed_inventory(
    metadata: Sequence[dict[str, str]], covid_dir: Path, hc_dir: Path
) -> None:
    expected = {
        expected_bed_path(row["sample_id"], row["group"], covid_dir, hc_dir)
        for row in metadata
    }
    missing = sorted(path for path in expected if not path.is_file())
    if missing:
        raise FileNotFoundError(
            "Missing expected BED files: " + ", ".join(str(path) for path in missing)
        )
    observed = set(covid_dir.glob("*_circle_site.bed")) | set(
        hc_dir.glob("*_circle_site.bed")
    )
    unexpected = sorted(observed - expected)
    if unexpected:
        raise ValueError(
            "Unexpected Circle-Map BED files in input directories: "
            + ", ".join(str(path) for path in unexpected)
        )


def parse_unique_lengths(path: Path) -> tuple[np.ndarray, dict[str, object]]:
    """Return one length per unique (chromosome, start, end) coordinate."""

    digest = hashlib.sha256()
    chromosome_codes: dict[bytes, int] = {}
    seen: set[tuple[int, int, int]] = set()
    lengths: list[int] = []
    physical_lines = 0
    blank_or_comment_lines = 0
    header_or_invalid_lines = 0
    negative_length_rows = 0
    zero_length_rows_retained = 0
    valid_coordinate_rows = 0
    duplicate_coordinate_rows = 0

    with path.open("rb") as handle:
        for raw_line in handle:
            digest.update(raw_line)
            physical_lines += 1
            stripped = raw_line.strip()
            if not stripped or stripped.startswith(b"#"):
                blank_or_comment_lines += 1
                continue
            fields = stripped.split(b"\t")
            if len(fields) < 3:
                header_or_invalid_lines += 1
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                header_or_invalid_lines += 1
                continue
            if end < start:
                negative_length_rows += 1
                continue
            if end == start:
                zero_length_rows_retained += 1
            chromosome = fields[0]
            chromosome_code = chromosome_codes.setdefault(
                chromosome, len(chromosome_codes)
            )
            valid_coordinate_rows += 1
            coordinate = (chromosome_code, start, end)
            if coordinate in seen:
                duplicate_coordinate_rows += 1
                continue
            seen.add(coordinate)
            lengths.append(end - start)

    length_array = np.asarray(lengths, dtype=np.int32)
    audit = {
        "bed_path": str(path),
        "size_bytes": path.stat().st_size,
        "mtime_epoch": path.stat().st_mtime,
        "sha256": digest.hexdigest(),
        "physical_lines": physical_lines,
        "blank_or_comment_lines": blank_or_comment_lines,
        "header_or_invalid_lines": header_or_invalid_lines,
        "negative_length_rows": negative_length_rows,
        "zero_length_rows_retained": zero_length_rows_retained,
        "valid_coordinate_rows_before_deduplication": valid_coordinate_rows,
        "duplicate_coordinate_rows_removed": duplicate_coordinate_rows,
        "unique_coordinate_rows": len(length_array),
        "minimum_length_bp": int(length_array.min()) if len(length_array) else None,
        "maximum_length_bp": int(length_array.max()) if len(length_array) else None,
        "n_50_1000_bp": int(
            np.count_nonzero((length_array >= 50) & (length_array <= 1000))
        ),
        "n_50_2000_bp": int(
            np.count_nonzero((length_array >= 50) & (length_array <= 2000))
        ),
    }
    return length_array, audit


def load_samples(
    metadata: Sequence[dict[str, str]], covid_dir: Path, hc_dir: Path
) -> list[Sample]:
    samples = []
    for index, row in enumerate(metadata, start=1):
        path = expected_bed_path(row["sample_id"], row["group"], covid_dir, hc_dir)
        print(
            f"[{index:02d}/{len(metadata)}] reading {row['sample_id']} ({row['group']})",
            flush=True,
        )
        lengths, audit = parse_unique_lengths(path)
        audit["sample_id"] = row["sample_id"]
        audit["group"] = row["group"]
        samples.append(
            Sample(
                sample_id=row["sample_id"],
                group=row["group"],
                bed_path=path,
                lengths=lengths,
                audit=audit,
            )
        )
    return samples


def smooth_normalized_counts(counts: np.ndarray, sigma: float) -> np.ndarray:
    total = float(np.sum(counts))
    if total <= 0:
        raise ValueError("Cannot estimate a distribution with zero eligible eccDNAs")
    probability = counts.astype(float) / total
    curve = gaussian_filter1d(
        probability, sigma=sigma, mode="reflect", truncate=4.0
    )
    area = float(np.trapezoid(curve, dx=1.0))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("Smoothed density has invalid area")
    return curve / area


def curve_from_lengths(
    lengths: np.ndarray, low: int, high: int, sigma: float
) -> tuple[np.ndarray, int]:
    eligible = lengths[(lengths >= low) & (lengths <= high)]
    counts = np.bincount(eligible - low, minlength=high - low + 1)
    return smooth_normalized_counts(counts, sigma), int(len(eligible))


def sample_curve_matrix(
    samples: Sequence[Sample], low: int, high: int, sigma: float
) -> tuple[np.ndarray, np.ndarray]:
    curves = []
    counts = []
    for sample in samples:
        curve, eligible_count = curve_from_lengths(sample.lengths, low, high, sigma)
        curves.append(curve)
        counts.append(eligible_count)
    return np.vstack(curves), np.asarray(counts, dtype=int)


def pooled_curve(
    samples: Sequence[Sample], low: int, high: int, sigma: float
) -> tuple[np.ndarray, int]:
    total_counts = np.zeros(high - low + 1, dtype=np.int64)
    n_eligible = 0
    for sample in samples:
        eligible = sample.lengths[
            (sample.lengths >= low) & (sample.lengths <= high)
        ]
        total_counts += np.bincount(
            eligible - low, minlength=high - low + 1
        ).astype(np.int64)
        n_eligible += len(eligible)
    return smooth_normalized_counts(total_counts, sigma), int(n_eligible)


def indices_for_dataset(samples: Sequence[Sample], dataset: str) -> np.ndarray:
    if dataset == "Overall":
        return np.arange(len(samples), dtype=int)
    return np.asarray(
        [index for index, sample in enumerate(samples) if sample.group == dataset],
        dtype=int,
    )


def detect_peaks(
    curve: np.ndarray,
    grid: np.ndarray,
    search_low: int,
    search_high: int,
    distance_bp: int,
    prominence_fraction: float,
    minimum_width_bp: float,
) -> list[dict[str, object]]:
    mask = (grid >= search_low) & (grid <= search_high)
    search_indices = np.flatnonzero(mask)
    search_curve = curve[mask]
    density_range = float(np.max(curve) - np.min(curve))
    prominence_threshold = prominence_fraction * density_range
    grid_step = float(np.median(np.diff(grid)))
    distance_points = max(1, int(math.ceil(distance_bp / grid_step)))
    width_points = minimum_width_bp / grid_step
    local_indices, properties = find_peaks(
        search_curve,
        distance=distance_points,
        prominence=prominence_threshold,
        width=width_points,
    )
    calls = []
    for rank, local_index in enumerate(local_indices):
        global_index = int(search_indices[local_index])
        calls.append(
            {
                "call_rank": rank + 1,
                "position_bp": int(grid[global_index]),
                "density": float(curve[global_index]),
                "prominence": float(properties["prominences"][rank]),
                "prominence_threshold": prominence_threshold,
                "prominence_fraction_of_density_range": (
                    float(properties["prominences"][rank]) / density_range
                    if density_range > 0
                    else math.nan
                ),
                "width_bp": float(properties["widths"][rank] * grid_step),
                "left_base_bp": int(
                    grid[int(search_indices[int(properties["left_bases"][rank])])]
                ),
                "right_base_bp": int(
                    grid[int(search_indices[int(properties["right_bases"][rank])])]
                ),
            }
        )
    return calls


def nearest_call(
    calls: Sequence[dict[str, object]], position: int, window: int
) -> dict[str, object] | None:
    matching = [
        call
        for call in calls
        if abs(int(call["position_bp"]) - int(position)) <= window
    ]
    if not matching:
        return None
    return min(
        matching,
        key=lambda call: (
            abs(int(call["position_bp"]) - int(position)),
            -float(call["prominence"]),
        ),
    )


def add_candidate_matches(
    rows: list[dict[str, object]],
    candidates: Sequence[dict[str, object]],
    match_window: int,
) -> None:
    for row in rows:
        row["matched_primary_peak_id"] = ""
        row["distance_to_primary_peak_bp"] = ""
        for candidate in candidates:
            distance = abs(
                int(row["position_bp"]) - int(candidate["primary_peak_position_bp"])
            )
            if distance <= match_window:
                row["matched_primary_peak_id"] = candidate["peak_id"]
                row["distance_to_primary_peak_bp"] = distance
                break


def bootstrap_dataset(
    curves: np.ndarray,
    grid: np.ndarray,
    dataset: str,
    candidates: Sequence[dict[str, object]],
    iterations: int,
    rng: np.random.Generator,
    args: argparse.Namespace,
) -> tuple[list[dict[str, object]], np.ndarray, np.ndarray]:
    records: list[dict[str, object]] = []
    bootstrap_curves = np.empty((iterations, curves.shape[1]), dtype=np.float64)
    for iteration in range(iterations):
        sampled_indices = rng.integers(0, curves.shape[0], size=curves.shape[0])
        mean_curve = np.mean(curves[sampled_indices], axis=0)
        bootstrap_curves[iteration] = mean_curve
        calls = detect_peaks(
            mean_curve,
            grid,
            args.search_low,
            args.search_high,
            args.peak_distance,
            args.prominence_fraction,
            args.minimum_width,
        )
        for candidate in candidates:
            matched = nearest_call(
                calls,
                int(candidate["primary_peak_position_bp"]),
                args.match_window,
            )
            records.append(
                {
                    "dataset": dataset,
                    "bootstrap_iteration": iteration + 1,
                    "peak_id": candidate["peak_id"],
                    "primary_peak_position_bp": candidate[
                        "primary_peak_position_bp"
                    ],
                    "detected_within_match_window": matched is not None,
                    "matched_position_bp": (
                        matched["position_bp"] if matched is not None else None
                    ),
                    "matched_prominence": (
                        matched["prominence"] if matched is not None else None
                    ),
                    "matched_width_bp": (
                        matched["width_bp"] if matched is not None else None
                    ),
                    "n_peaks_in_replicate": len(calls),
                }
            )
    lower, upper = np.percentile(bootstrap_curves, [2.5, 97.5], axis=0)
    return records, lower, upper


def percentile_or_na(values: Sequence[float], percentile: float) -> float:
    return float(np.percentile(np.asarray(values, dtype=float), percentile)) if values else math.nan


def summarize_bootstrap(
    records: Sequence[dict[str, object]], iterations: int
) -> dict[tuple[str, str], dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in records:
        grouped[(str(row["dataset"]), str(row["peak_id"]))].append(row)
    summaries = {}
    for key, rows in grouped.items():
        positions = [
            float(row["matched_position_bp"])
            for row in rows
            if row["matched_position_bp"] is not None
        ]
        summaries[key] = {
            "bootstrap_detection_count": len(positions),
            "bootstrap_detection_frequency": len(positions) / iterations,
            "bootstrap_position_median_bp": (
                float(np.median(positions)) if positions else math.nan
            ),
            "bootstrap_position_ci_low_bp": percentile_or_na(positions, 2.5),
            "bootstrap_position_ci_high_bp": percentile_or_na(positions, 97.5),
            "bootstrap_position_min_bp": min(positions) if positions else math.nan,
            "bootstrap_position_max_bp": max(positions) if positions else math.nan,
        }
    return summaries


def bh_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if not len(finite_indices):
        return adjusted
    finite_values = values[finite_indices]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    correction = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    correction = np.minimum.accumulate(correction[::-1])[::-1]
    restored = np.empty_like(correction)
    restored[order] = np.minimum(correction, 1.0)
    adjusted[finite_indices] = restored
    return adjusted


def group_summary(values: np.ndarray) -> dict[str, float]:
    return {
        "n": int(len(values)),
        "median": float(np.median(values)),
        "q1": float(np.percentile(values, 25)),
        "q3": float(np.percentile(values, 75)),
        "iqr": float(np.percentile(values, 75) - np.percentile(values, 25)),
    }


def peak_window_statistics(
    samples: Sequence[Sample],
    stable_candidates: Sequence[dict[str, object]],
    low: int,
    high: int,
    match_window: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sample_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    for candidate in stable_candidates:
        position = int(candidate["primary_peak_position_bp"])
        window_low = position - match_window
        window_high = position + match_window
        for sample in samples:
            eligible = sample.lengths[
                (sample.lengths >= low) & (sample.lengths <= high)
            ]
            window_count = int(
                np.count_nonzero(
                    (eligible >= window_low) & (eligible <= window_high)
                )
            )
            sample_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "group": sample.group,
                    "peak_id": candidate["peak_id"],
                    "consensus_peak_position_bp": position,
                    "window_low_bp_inclusive": window_low,
                    "window_high_bp_inclusive": window_high,
                    "window_eccdna_count": window_count,
                    "denominator_50_1000_bp": len(eligible),
                    "window_proportion": window_count / len(eligible),
                }
            )

    for candidate in stable_candidates:
        peak_id = str(candidate["peak_id"])
        covid = np.asarray(
            [
                float(row["window_proportion"])
                for row in sample_rows
                if row["peak_id"] == peak_id and row["group"] == GROUP_COVID
            ]
        )
        hc = np.asarray(
            [
                float(row["window_proportion"])
                for row in sample_rows
                if row["peak_id"] == peak_id and row["group"] == GROUP_HC
            ]
        )
        test = mannwhitneyu(covid, hc, alternative="two-sided", method="auto")
        delta = 2.0 * float(test.statistic) / (len(covid) * len(hc)) - 1.0
        covid_summary = group_summary(covid)
        hc_summary = group_summary(hc)
        comparison_rows.append(
            {
                "peak_id": peak_id,
                "consensus_peak_position_bp": candidate[
                    "primary_peak_position_bp"
                ],
                "window_definition": (
                    f"{int(candidate['primary_peak_position_bp']) - match_window}-"
                    f"{int(candidate['primary_peak_position_bp']) + match_window} bp, inclusive"
                ),
                "denominator_definition": "unique eccDNAs of 50-1,000 bp in the same sample",
                "n_covid": covid_summary["n"],
                "covid_median_proportion": covid_summary["median"],
                "covid_q1_proportion": covid_summary["q1"],
                "covid_q3_proportion": covid_summary["q3"],
                "covid_iqr_proportion": covid_summary["iqr"],
                "n_hc": hc_summary["n"],
                "hc_median_proportion": hc_summary["median"],
                "hc_q1_proportion": hc_summary["q1"],
                "hc_q3_proportion": hc_summary["q3"],
                "hc_iqr_proportion": hc_summary["iqr"],
                "median_difference_covid_minus_hc": (
                    covid_summary["median"] - hc_summary["median"]
                ),
                "mann_whitney_u_covid": float(test.statistic),
                "p_value_two_sided": float(test.pvalue),
                "p_adj_bh": math.nan,
                "cliffs_delta_covid_minus_hc": delta,
                "test": "two-sided Wilcoxon rank-sum (Mann-Whitney U)",
                "multiple_testing": "Benjamini-Hochberg across stable peak windows",
            }
        )
    if comparison_rows:
        adjusted = bh_adjust(
            [float(row["p_value_two_sided"]) for row in comparison_rows]
        )
        for row, q_value in zip(comparison_rows, adjusted):
            row["p_adj_bh"] = float(q_value)
    return sample_rows, comparison_rows


def run_equal_count_subsampling(
    samples: Sequence[Sample],
    candidates: Sequence[dict[str, object]],
    repeats: int,
    count: int,
    rng: np.random.Generator,
    args: argparse.Namespace,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if min(len(sample.lengths) for sample in samples) < count:
        smallest = min(samples, key=lambda sample: len(sample.lengths))
        raise ValueError(
            f"Cannot subsample {count} unique eccDNAs per sample: "
            f"{smallest.sample_id} has {len(smallest.lengths)}"
        )
    grid = np.arange(args.primary_low, args.primary_high + 1)
    detail_rows: list[dict[str, object]] = []
    sample_count_rows: list[dict[str, object]] = []
    dataset_indices = {
        dataset: indices_for_dataset(samples, dataset) for dataset in DATASETS
    }
    for repeat in range(repeats):
        curves = []
        for sample in samples:
            chosen_indices = rng.choice(
                len(sample.lengths), size=count, replace=False, shuffle=False
            )
            chosen_lengths = sample.lengths[chosen_indices]
            curve, eligible_count = curve_from_lengths(
                chosen_lengths,
                args.primary_low,
                args.primary_high,
                args.primary_sigma,
            )
            curves.append(curve)
            sample_count_rows.append(
                {
                    "subsample_repeat": repeat + 1,
                    "sample_id": sample.sample_id,
                    "group": sample.group,
                    "total_unique_eccdna_sampled": count,
                    "sampled_eccdna_50_1000_bp": eligible_count,
                }
            )
        matrix = np.vstack(curves)
        for dataset in DATASETS:
            mean_curve = np.mean(matrix[dataset_indices[dataset]], axis=0)
            calls = detect_peaks(
                mean_curve,
                grid,
                args.search_low,
                args.search_high,
                args.peak_distance,
                args.prominence_fraction,
                args.minimum_width,
            )
            for candidate in candidates:
                matched = nearest_call(
                    calls,
                    int(candidate["primary_peak_position_bp"]),
                    args.match_window,
                )
                detail_rows.append(
                    {
                        "subsample_repeat": repeat + 1,
                        "dataset": dataset,
                        "peak_id": candidate["peak_id"],
                        "primary_peak_position_bp": candidate[
                            "primary_peak_position_bp"
                        ],
                        "detected_within_match_window": matched is not None,
                        "matched_position_bp": (
                            matched["position_bp"] if matched is not None else None
                        ),
                        "matched_prominence": (
                            matched["prominence"] if matched is not None else None
                        ),
                        "matched_width_bp": (
                            matched["width_bp"] if matched is not None else None
                        ),
                        "n_peaks_in_subsample_curve": len(calls),
                    }
                )
        if (repeat + 1) % 10 == 0:
            print(
                f"[subsampling] completed {repeat + 1}/{repeats} repeats",
                flush=True,
            )
    return detail_rows, sample_count_rows


def summarize_subsampling(
    rows: Sequence[dict[str, object]], repeats: int
) -> dict[tuple[str, str], dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["dataset"]), str(row["peak_id"]))].append(row)
    result = {}
    for key, group_rows in grouped.items():
        positions = [
            float(row["matched_position_bp"])
            for row in group_rows
            if row["matched_position_bp"] is not None
        ]
        result[key] = {
            "subsample_detection_count": len(positions),
            "subsample_detection_frequency": len(positions) / repeats,
            "subsample_position_median_bp": (
                float(np.median(positions)) if positions else math.nan
            ),
            "subsample_position_ci_low_bp": percentile_or_na(positions, 2.5),
            "subsample_position_ci_high_bp": percentile_or_na(positions, 97.5),
            "subsample_position_min_bp": min(positions) if positions else math.nan,
            "subsample_position_max_bp": max(positions) if positions else math.nan,
        }
    return result


def compute_spacing_table(
    stable_candidates: Sequence[dict[str, object]],
    bootstrap_records: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    if len(stable_candidates) < 2:
        return []
    overall_records = [
        row for row in bootstrap_records if row["dataset"] == "Overall"
    ]
    positions_by_iteration: dict[int, dict[str, int]] = defaultdict(dict)
    for row in overall_records:
        if row["matched_position_bp"] is not None:
            positions_by_iteration[int(row["bootstrap_iteration"])][
                str(row["peak_id"])
            ] = int(row["matched_position_bp"])
    rows = []
    for left, right in zip(stable_candidates[:-1], stable_candidates[1:]):
        spacings = []
        for positions in positions_by_iteration.values():
            if left["peak_id"] in positions and right["peak_id"] in positions:
                spacings.append(positions[right["peak_id"]] - positions[left["peak_id"]])
        rows.append(
            {
                "left_peak_id": left["peak_id"],
                "right_peak_id": right["peak_id"],
                "left_primary_position_bp": left["primary_peak_position_bp"],
                "right_primary_position_bp": right["primary_peak_position_bp"],
                "primary_spacing_bp": (
                    int(right["primary_peak_position_bp"])
                    - int(left["primary_peak_position_bp"])
                ),
                "bootstrap_paired_detection_count": len(spacings),
                "bootstrap_spacing_median_bp": (
                    float(np.median(spacings)) if spacings else math.nan
                ),
                "bootstrap_spacing_ci_low_bp": percentile_or_na(spacings, 2.5),
                "bootstrap_spacing_ci_high_bp": percentile_or_na(spacings, 97.5),
            }
        )
    return rows


def compute_periodicity_summary(
    stable_candidates: Sequence[dict[str, object]],
    bootstrap_records: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    if len(stable_candidates) < 2:
        return []
    peak_ids = [str(candidate["peak_id"]) for candidate in stable_candidates]
    primary_positions = [
        int(candidate["primary_peak_position_bp"]) for candidate in stable_candidates
    ]
    primary_spacings = np.diff(primary_positions)

    positions_by_iteration: dict[int, dict[str, int]] = defaultdict(dict)
    for row in bootstrap_records:
        if row["dataset"] != "Overall" or row["matched_position_bp"] is None:
            continue
        positions_by_iteration[int(row["bootstrap_iteration"])][
            str(row["peak_id"])
        ] = int(row["matched_position_bp"])

    mean_spacings = []
    median_spacings = []
    for positions in positions_by_iteration.values():
        if not all(peak_id in positions for peak_id in peak_ids):
            continue
        ordered_positions = [positions[peak_id] for peak_id in peak_ids]
        replicate_spacings = np.diff(ordered_positions)
        mean_spacings.append(float(np.mean(replicate_spacings)))
        median_spacings.append(float(np.median(replicate_spacings)))

    return [
        {
            "stable_peak_ids": ",".join(peak_ids),
            "stable_peak_positions_bp": ",".join(
                str(position) for position in primary_positions
            ),
            "n_stable_peaks": len(stable_candidates),
            "n_adjacent_spacings": len(primary_spacings),
            "primary_adjacent_spacings_bp": ",".join(
                str(int(spacing)) for spacing in primary_spacings
            ),
            "primary_mean_spacing_bp": float(np.mean(primary_spacings)),
            "primary_median_spacing_bp": float(np.median(primary_spacings)),
            "bootstrap_complete_detection_count": len(mean_spacings),
            "bootstrap_mean_spacing_median_bp": (
                float(np.median(mean_spacings)) if mean_spacings else math.nan
            ),
            "bootstrap_mean_spacing_ci_low_bp": percentile_or_na(
                mean_spacings, 2.5
            ),
            "bootstrap_mean_spacing_ci_high_bp": percentile_or_na(
                mean_spacings, 97.5
            ),
            "bootstrap_median_spacing_median_bp": (
                float(np.median(median_spacings)) if median_spacings else math.nan
            ),
            "bootstrap_median_spacing_ci_low_bp": percentile_or_na(
                median_spacings, 2.5
            ),
            "bootstrap_median_spacing_ci_high_bp": percentile_or_na(
                median_spacings, 97.5
            ),
        }
    ]


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": 6,
            "axes.unicode_minus": False,
            "axes.labelsize": 6.5,
            "axes.linewidth": 0.6,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.fontsize": 5.5,
            "legend.frameon": False,
            "lines.linewidth": 1.0,
        }
    )


def save_figure(fig: plt.Figure, base_path: Path, dpi: int = 600) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base_path.with_suffix(".svg"))
    fig.savefig(base_path.with_suffix(".pdf"))
    fig.savefig(base_path.with_suffix(".png"), dpi=dpi)
    fig.savefig(
        base_path.with_suffix(".tiff"),
        dpi=dpi,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def make_main_figure(
    output_base: Path,
    grid: np.ndarray,
    covid_mean: np.ndarray,
    hc_mean: np.ndarray,
    covid_ci: tuple[np.ndarray, np.ndarray],
    hc_ci: tuple[np.ndarray, np.ndarray],
    stable_candidates: Sequence[dict[str, object]],
) -> None:
    configure_matplotlib()
    fig, ax = plt.subplots(figsize=(89 / 25.4, 62 / 25.4))
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.19, top=0.96)

    scale = 100.0
    ax.fill_between(
        grid,
        covid_ci[0] * scale,
        covid_ci[1] * scale,
        color=COVID_LIGHT,
        alpha=0.72,
        linewidth=0,
    )
    ax.fill_between(
        grid,
        hc_ci[0] * scale,
        hc_ci[1] * scale,
        color=HC_LIGHT,
        alpha=0.72,
        linewidth=0,
    )
    ax.plot(
        grid,
        covid_mean * scale,
        color=COVID_COLOR,
        label="COVID-19 (n = 39)",
        zorder=3,
    )
    ax.plot(
        grid, hc_mean * scale, color=HC_COLOR, label="HC (n = 39)", zorder=3
    )
    y_max = max(float(np.max(covid_ci[1])), float(np.max(hc_ci[1]))) * scale
    for index, candidate in enumerate(stable_candidates):
        position = int(candidate["primary_peak_position_bp"])
        ax.axvline(
            position,
            color=GREY,
            linestyle=(0, (2.5, 2.5)),
            linewidth=0.65,
            zorder=1,
        )
        text_y = y_max * (0.97 if index % 2 == 0 else 0.86)
        ax.text(
            position,
            text_y,
            f"{position} bp",
            ha="center",
            va="top",
            fontsize=5.2,
            color=BLACK,
        )
    ax.set_xlim(grid[0], grid[-1])
    ax.set_ylim(0, y_max * 1.03)
    ax.set_xlabel("eccDNA length (bp)")
    ax.set_ylabel("Probability density (% per bp)")
    ax.set_xticks([50, 250, 500, 750, 1000])
    ax.legend(loc="upper right", handlelength=2.4)
    save_figure(fig, output_base)


def deterministic_jitter(n: int, width: float = 0.065) -> np.ndarray:
    if n <= 1:
        return np.zeros(n)
    base = np.linspace(-width, width, n)
    order = np.argsort((np.arange(n) * 17) % n)
    return base[order]


def make_supplementary_figure(
    output_base: Path,
    grid: np.ndarray,
    sigma_curves: dict[float, np.ndarray],
    sigma_calls: dict[float, list[dict[str, object]]],
    candidates: Sequence[dict[str, object]],
    subsample_summary: dict[tuple[str, str], dict[str, object]],
    sample_window_rows: Sequence[dict[str, object]],
    comparison_rows: Sequence[dict[str, object]],
) -> None:
    configure_matplotlib()
    fig = plt.figure(figsize=(183 / 25.4, 132 / 25.4))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=(1.05, 1.0),
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.965,
        wspace=0.30,
        hspace=0.48,
    )
    top = outer[0, :].subgridspec(1, len(sigma_curves), wspace=0.20)
    sigma_axes = []
    maximum = max(float(np.max(curve)) for curve in sigma_curves.values()) * 100
    for index, sigma in enumerate(sorted(sigma_curves)):
        ax = fig.add_subplot(top[0, index], sharey=sigma_axes[0] if sigma_axes else None)
        sigma_axes.append(ax)
        plotted_curve = sigma_curves[sigma] * 100
        ax.plot(grid, plotted_curve, color=OVERALL_COLOR, linewidth=0.9)
        for call in sigma_calls[sigma]:
            position = int(call["position_bp"])
            peak_y = float(np.interp(position, grid, plotted_curve))
            ax.text(
                position,
                min(peak_y + maximum * 0.035, maximum * 1.08),
                str(position),
                ha="center",
                va="bottom",
                fontsize=4.8,
                color=BLACK,
            )
        ax.set_xlim(grid[0], grid[-1])
        ax.set_ylim(0, maximum * 1.14)
        ax.set_xticks([100, 500, 900])
        ax.set_title(f"sigma = {format_value(sigma)} bp", fontsize=6.2, pad=3)
        ax.set_xlabel("Length (bp)")
        if index == 0:
            ax.set_ylabel("Density (% per bp)")
        else:
            ax.tick_params(labelleft=False)

    ax_b = fig.add_subplot(outer[1, 0])
    peak_ids = [str(candidate["peak_id"]) for candidate in candidates]
    x = np.arange(len(peak_ids), dtype=float)
    width = 0.22
    display = [
        ("Overall", OVERALL_COLOR),
        (GROUP_COVID, COVID_COLOR),
        (GROUP_HC, HC_COLOR),
    ]
    for offset_index, (dataset, color) in enumerate(display):
        values = [
            100
            * float(
                subsample_summary.get((dataset, peak_id), {}).get(
                    "subsample_detection_frequency", math.nan
                )
            )
            for peak_id in peak_ids
        ]
        ax_b.bar(
            x + (offset_index - 1) * width,
            values,
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.35,
            label=dataset,
        )
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(
        [
            f"{candidate['peak_id']}\n({candidate['primary_peak_position_bp']} bp)"
            for candidate in candidates
        ]
    )
    ax_b.set_ylim(0, 105)
    ax_b.set_ylabel("Detection frequency (%)")
    ax_b.set_xlabel("Primary algorithmic candidate peak")
    ax_b.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=3,
        columnspacing=0.8,
        handlelength=1.4,
        borderaxespad=0,
    )

    ax_c = fig.add_subplot(outer[1, 1])
    stable_peak_ids = sorted(
        {str(row["peak_id"]) for row in sample_window_rows},
        key=lambda peak_id: int(peak_id.lstrip("P")),
    )
    comparison_by_peak = {str(row["peak_id"]): row for row in comparison_rows}
    if stable_peak_ids:
        centers = np.arange(len(stable_peak_ids), dtype=float)
        for peak_index, peak_id in enumerate(stable_peak_ids):
            for group_index, (group, box_color, point_color) in enumerate(
                (
                    (GROUP_COVID, FIG4_COVID_BOX, FIG4_COVID_POINT),
                    (GROUP_HC, FIG4_HC_BOX, FIG4_HC_POINT),
                )
            ):
                values = np.asarray(
                    [
                        100 * float(row["window_proportion"])
                        for row in sample_window_rows
                        if row["peak_id"] == peak_id and row["group"] == group
                    ]
                )
                position = centers[peak_index] + (-0.17 if group_index == 0 else 0.17)
                box = ax_c.boxplot(
                    [values],
                    positions=[position],
                    widths=0.28,
                    patch_artist=True,
                    showfliers=False,
                    whis=1.5,
                    medianprops={"color": "black", "linewidth": 0.75},
                    boxprops={"color": "black", "linewidth": 0.55},
                    whiskerprops={"color": "black", "linewidth": 0.55},
                    capprops={"color": "black", "linewidth": 0.55},
                )
                box["boxes"][0].set_facecolor(box_color)
                rng = np.random.default_rng(8400 + peak_index * 20 + group_index)
                ax_c.scatter(
                    position + rng.uniform(-0.075, 0.075, len(values)),
                    values,
                    s=5.2,
                    marker="o",
                    facecolor=point_color,
                    edgecolor="none",
                    alpha=0.88,
                    zorder=3,
                )
        y_top_by_peak = {}
        for peak_index, peak_id in enumerate(stable_peak_ids):
            peak_values = [
                100 * float(row["window_proportion"])
                for row in sample_window_rows
                if row["peak_id"] == peak_id
            ]
            y_top_by_peak[peak_id] = max(peak_values) if peak_values else 0.0
        y_axis_top = max(y_top_by_peak.values(), default=1.0) * 1.22
        ax_c.set_ylim(0, max(y_axis_top, 1.0))
        for peak_index, peak_id in enumerate(stable_peak_ids):
            q_value = float(comparison_by_peak[peak_id]["p_adj_bh"])
            q_label = f"q = {q_value:.2g}" if q_value >= 0.001 else "q < 0.001"
            left = centers[peak_index] - 0.17
            right = centers[peak_index] + 0.17
            bracket_y = y_top_by_peak[peak_id] + y_axis_top * 0.035
            bracket_h = y_axis_top * 0.025
            ax_c.plot(
                [left, left, right, right],
                [bracket_y, bracket_y + bracket_h, bracket_y + bracket_h, bracket_y],
                color="black",
                linewidth=0.55,
                clip_on=False,
            )
            ax_c.text(
                centers[peak_index],
                bracket_y + bracket_h + y_axis_top * 0.02,
                q_label,
                ha="center",
                va="bottom",
                fontsize=4.8,
            )
        ax_c.set_xticks(centers)
        ax_c.set_xticklabels(
            [
                (
                    f"{peak_id}\n"
                    f"({int(comparison_by_peak[peak_id]['consensus_peak_position_bp'])} bp)"
                )
                for peak_id in stable_peak_ids
            ]
        )
        handles = [
            Patch(
                facecolor=FIG4_COVID_BOX,
                edgecolor="black",
                linewidth=0.45,
                label="COVID-19",
            ),
            Patch(
                facecolor=FIG4_HC_BOX,
                edgecolor="black",
                linewidth=0.45,
                label="HC",
            ),
        ]
        ax_c.legend(
            handles=handles,
            loc="lower left",
            bbox_to_anchor=(0.0, 1.03),
            ncol=2,
            columnspacing=0.8,
            handletextpad=0.35,
            borderaxespad=0,
        )
    else:
        ax_c.text(
            0.5,
            0.5,
            "No peak met all prespecified\nstability criteria",
            transform=ax_c.transAxes,
            ha="center",
            va="center",
        )
        ax_c.set_xticks([])
    ax_c.set_ylabel("eccDNAs in peak window (%)")
    ax_c.set_xlabel("Stable consensus peak (+/-25 bp)")

    fig.text(0.012, 0.972, "a", fontsize=8, fontweight="bold", va="top")
    fig.text(0.012, 0.485, "b", fontsize=8, fontweight="bold", va="top")
    fig.text(0.515, 0.485, "c", fontsize=8, fontweight="bold", va="top")
    save_figure(fig, output_base)


def write_analysis_summary(
    path: Path,
    samples: Sequence[Sample],
    candidates: Sequence[dict[str, object]],
    stable_candidates: Sequence[dict[str, object]],
    spacing_rows: Sequence[dict[str, object]],
    periodicity_rows: Sequence[dict[str, object]],
    comparison_rows: Sequence[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    total_raw = sum(
        int(sample.audit["valid_coordinate_rows_before_deduplication"])
        for sample in samples
    )
    total_unique = sum(len(sample.lengths) for sample in samples)
    lines = [
        "# Objective eccDNA fragment-length peak analysis",
        "",
        f"Generated: {now_utc()}",
        "",
        "## Input audit",
        "",
        f"- Biological samples: {len(samples)} (39 COVID-19; 39 HC).",
        f"- Valid BED coordinate rows before within-sample deduplication: {total_raw:,}.",
        f"- Unique chromosome/start/end intervals analysed: {total_unique:,}.",
        f"- Exact coordinate duplicates removed: {total_raw - total_unique:,}.",
        (
            "- Smallest unique per-sample call set: "
            f"{min(len(sample.lengths) for sample in samples):,} eccDNAs."
        ),
        "",
        "## Fixed primary algorithm",
        "",
        (
            f"Each sample contributed an area-normalized 1-bp distribution over "
            f"{args.primary_low}-{args.primary_high:,} bp. Gaussian smoothing used "
            f"σ = {format_value(args.primary_sigma)} bp, reflect boundary mode and "
            "a four-sigma kernel truncation. Samples were averaged with equal weight."
        ),
        (
            f"`find_peaks` searched {args.search_low}-{args.search_high} bp with "
            f"distance ≥{args.peak_distance} bp, prominence ≥"
            f"{100 * args.prominence_fraction:g}% of the complete curve density range, "
            f"and width ≥{format_value(args.minimum_width)} bp."
        ),
        "",
        "## Algorithmic candidates and stability",
        "",
    ]
    if not candidates:
        lines.append("- The fixed primary algorithm detected no candidate peaks.")
    for candidate in candidates:
        lines.append(
            "- {peak_id}: {position} bp; overall bootstrap {overall:.1%}, "
            "COVID-19 {covid:.1%}, HC {hc:.1%}; bandwidth matches {bandwidth}/4; "
            "classification: {classification}.".format(
                peak_id=candidate["peak_id"],
                position=int(candidate["primary_peak_position_bp"]),
                overall=float(candidate["overall_bootstrap_detection_frequency"]),
                covid=float(candidate["covid_bootstrap_detection_frequency"]),
                hc=float(candidate["hc_bootstrap_detection_frequency"]),
                bandwidth=int(candidate["overall_bandwidth_detection_count"]),
                classification=(
                    "stable consensus peak"
                    if candidate["stable_consensus_peak"]
                    else "not stable under all prespecified criteria"
                ),
            )
        )
        if not candidate["stable_consensus_peak"]:
            lines.append(f"  Failure reason(s): {candidate['stability_failure_reasons']}.")

    lines.extend(["", "## Adjacent-peak spacing", ""])
    if spacing_rows:
        for row in spacing_rows:
            lines.append(
                f"- {row['left_peak_id']}–{row['right_peak_id']}: "
                f"{int(row['primary_spacing_bp'])} bp; bootstrap median "
                f"{format_value(row['bootstrap_spacing_median_bp'])} bp "
                f"(95% CI {format_value(row['bootstrap_spacing_ci_low_bp'])}–"
                f"{format_value(row['bootstrap_spacing_ci_high_bp'])} bp)."
            )
        if periodicity_rows:
            row = periodicity_rows[0]
            lines.append(
                "- Stable-peak spacing summary: primary mean "
                f"{format_value(row['primary_mean_spacing_bp'])} bp and median "
                f"{format_value(row['primary_median_spacing_bp'])} bp; bootstrap "
                "median of mean spacing "
                f"{format_value(row['bootstrap_mean_spacing_median_bp'])} bp "
                f"(95% CI {format_value(row['bootstrap_mean_spacing_ci_low_bp'])}–"
                f"{format_value(row['bootstrap_mean_spacing_ci_high_bp'])} bp)."
            )
    else:
        lines.append("- Fewer than two stable consensus peaks; spacing was not summarized.")

    lines.extend(["", "## Peak-window group comparisons", ""])
    if comparison_rows:
        for row in comparison_rows:
            lines.append(
                f"- {row['peak_id']} ({row['consensus_peak_position_bp']} ±25 bp): "
                f"COVID-19 median {100 * float(row['covid_median_proportion']):.2f}% "
                f"versus HC {100 * float(row['hc_median_proportion']):.2f}%; "
                f"BH q={float(row['p_adj_bh']):.3g}; "
                f"Cliff's δ={float(row['cliffs_delta_covid_minus_hc']):.3f}."
            )
    else:
        lines.append("- No stable consensus peak windows were available for comparison.")

    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            (
                f"The fixed algorithm detected {len(candidates)} primary candidate peak(s), "
                f"of which {len(stable_candidates)} met every prespecified consensus "
                "criterion. Peak counts and positions were not manually altered. "
                "Peak spacing may be described as nucleosome-associated or "
                "nucleosome-like periodic organization, but these data alone do not "
                "establish mono-, di-, tri- or tetranucleosomal origin."
            ),
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_text_drafts(
    output_dir: Path,
    candidates: Sequence[dict[str, object]],
    stable_candidates: Sequence[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    stable_positions = [
        str(int(candidate["primary_peak_position_bp"]))
        for candidate in stable_candidates
    ]
    position_phrase = (
        ", ".join(stable_positions[:-1]) + f", and {stable_positions[-1]} bp"
        if len(stable_positions) > 1
        else (f"{stable_positions[0]} bp" if stable_positions else "none")
    )
    methods = f"""**Objective identification of eccDNA fragment-length peaks**
EccDNA length was calculated as the BED end coordinate minus the start coordinate for each unique chromosome-start-end interval within each sample. The primary analysis was restricted to eccDNAs of {args.primary_low} to {args.primary_high:,} bp, with {args.primary_low} to {args.sensitivity_high:,} bp evaluated as a sensitivity range. A normalized 1-bp length histogram was generated separately for each sample and smoothed using a Gaussian kernel (σ = {format_value(args.primary_sigma)} bp; reflect boundary mode; kernel truncated at 4σ). Each smoothed sample-level curve was renormalized to unit area and the curves were averaged with equal sample weight. Candidate peaks were identified from the all-sample mean curve using `scipy.signal.find_peaks`, without prespecifying their number or positions. Peaks were required to occur between {args.search_low} and {args.search_high} bp, to have a minimum inter-peak distance of {args.peak_distance} bp, a prominence of at least {100 * args.prominence_fraction:g}% of the full smoothed-curve density range, and a width of at least {format_value(args.minimum_width)} bp at half prominence. The identical algorithm was applied independently to the COVID-19 and HC mean curves. Peak stability was evaluated in {args.bootstrap_iterations:,} sample-level bootstrap resamplings; positional confidence intervals were percentile intervals conditional on detection within ±{args.match_window} bp of the primary all-sample peak. A peak was classified as a stable consensus peak when it was detected in at least 90% of all-sample bootstrap replicates and at least 80% of bootstrap replicates in each group, had a corresponding peak within ±{args.match_window} bp in both group mean curves, and was recovered at the primary position window under at least three of four Gaussian bandwidths (σ = 10, 15, 20 and 30 bp). Sensitivity analyses also used the 50–2,000 bp range, molecule-pooled distributions, and {args.subsample_repeats} equal-count subsampling repeats of {args.subsample_count:,} unique eccDNAs per sample. For each stable peak, the sample-level proportion of 50–1,000-bp eccDNAs within ±{args.match_window} bp was compared between groups using a two-sided Wilcoxon rank-sum test, with Benjamini–Hochberg correction across peak windows and Cliff's δ as the effect size.
"""
    response = f"""**Response:** We thank the reviewer for identifying the insufficient description of the fragment-length peak analysis. We replaced the previous visual description with a prespecified, reproducible peak-detection workflow based on normalized sample-level distributions and equal sample weighting. The number and positions of peaks were not specified in advance. The fixed primary algorithm detected {len(candidates)} candidate peak(s); {len(stable_candidates)} met all prespecified overall, group-specific, bootstrap and bandwidth-stability criteria ({position_phrase}). We now report the Gaussian smoothing parameters, local-maximum definition, peak prominence and width, sample-level bootstrap detection frequencies and positional confidence intervals, together with 50–2,000-bp, pooled-molecule and equal-count subsampling sensitivity analyses. The revised figure shows group mean distributions with pointwise 95% sample-bootstrap intervals, and the complete objective peak metrics are provided in the supplementary table.
"""
    results = (
        f"The prespecified sample-weighted algorithm identified {len(candidates)} "
        f"candidate fragment-length peak(s) in the all-sample distribution. "
        f"{len(stable_candidates)} satisfied the complete stability definition"
    )
    if stable_candidates:
        results += f", at {position_phrase}"
    results += (
        ". Peak positions were estimated independently rather than adjusted to the "
        "previously reported values."
    )
    (output_dir / "methods_text_draft.md").write_text(methods, encoding="utf-8")
    (output_dir / "reviewer_response_draft.md").write_text(
        response, encoding="utf-8"
    )
    (output_dir / "results_text_draft.md").write_text(
        results + "\n", encoding="utf-8"
    )
    captions = f"""**Figure 1D. Objective eccDNA fragment-length peak detection.**
Sample-weighted fragment-length density curves for COVID-19 and healthy control (HC) samples were generated from unique high-confidence eccDNA intervals of {args.primary_low}-{args.primary_high:,} bp. Lines show equal-weight group means and shaded bands show pointwise 95% sample-level bootstrap intervals. Vertical dashed lines mark stable consensus peaks that satisfied all prespecified algorithmic, bootstrap and bandwidth-sensitivity criteria.

**Figure S6. Robustness of objective fragment-length peak detection.**
a, All-sample sample-weighted density curves smoothed with Gaussian bandwidths of 10, 15, 20 and 30 bp; labels denote automatically detected local maxima. b, Equal-count sensitivity analysis using {args.subsample_repeats} random subsampling repeats of {args.subsample_count:,} unique eccDNAs per sample. c, Sample-level proportions of 50-1,000-bp eccDNAs falling within each stable peak window (peak position +/- {args.match_window} bp); q values are from two-sided Wilcoxon rank-sum tests with Benjamini-Hochberg correction across stable peak windows.
"""
    (output_dir / "figure_caption_drafts.md").write_text(
        captions, encoding="utf-8"
    )


def main() -> int:
    args = parse_args()
    sigmas = tuple(float(value) for value in args.sigmas.split(","))
    if args.primary_sigma not in sigmas:
        raise ValueError("The primary sigma must be included in --sigmas")
    if not (0 < args.prominence_fraction < 1):
        raise ValueError("--prominence-fraction must be between 0 and 1")
    if args.search_low <= args.primary_low or args.search_high >= args.primary_high:
        raise ValueError("Peak search bounds must lie inside the primary range")
    if args.bootstrap_iterations < 1 or args.subsample_repeats < 1:
        raise ValueError("Bootstrap and subsampling repeat counts must be positive")

    output_dir = args.output_dir.resolve()
    tables_dir = output_dir / "tables"
    data_dir = output_dir / "source_data"
    figures_dir = output_dir / "figures"
    logs_dir = output_dir / "logs"
    for directory in (tables_dir, data_dir, figures_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    print(f"Analysis started: {now_utc()}", flush=True)
    metadata = read_metadata(args.metadata)
    validate_bed_inventory(metadata, args.covid_bed_dir, args.hc_bed_dir)
    samples = load_samples(metadata, args.covid_bed_dir, args.hc_bed_dir)

    audit_fields = [
        "sample_id",
        "group",
        "bed_path",
        "size_bytes",
        "mtime_epoch",
        "sha256",
        "physical_lines",
        "blank_or_comment_lines",
        "header_or_invalid_lines",
        "negative_length_rows",
        "zero_length_rows_retained",
        "valid_coordinate_rows_before_deduplication",
        "duplicate_coordinate_rows_removed",
        "unique_coordinate_rows",
        "minimum_length_bp",
        "maximum_length_bp",
        "n_50_1000_bp",
        "n_50_2000_bp",
    ]
    write_tsv(tables_dir / "input_bed_audit.tsv", [s.audit for s in samples], audit_fields)

    grid = np.arange(args.primary_low, args.primary_high + 1)
    curve_matrices: dict[float, np.ndarray] = {}
    eligible_counts = None
    for sigma in sigmas:
        print(f"[curves] primary range, sigma={sigma:g}", flush=True)
        matrix, counts = sample_curve_matrix(
            samples, args.primary_low, args.primary_high, sigma
        )
        curve_matrices[sigma] = matrix
        if eligible_counts is None:
            eligible_counts = counts
        elif not np.array_equal(eligible_counts, counts):
            raise AssertionError("Eligible counts unexpectedly changed with bandwidth")
    assert eligible_counts is not None
    primary_matrix = curve_matrices[args.primary_sigma]
    dataset_indices = {
        dataset: indices_for_dataset(samples, dataset) for dataset in DATASETS
    }
    dataset_means = {
        dataset: np.mean(primary_matrix[indices], axis=0)
        for dataset, indices in dataset_indices.items()
    }

    primary_calls = detect_peaks(
        dataset_means["Overall"],
        grid,
        args.search_low,
        args.search_high,
        args.peak_distance,
        args.prominence_fraction,
        args.minimum_width,
    )
    candidates: list[dict[str, object]] = []
    for index, call in enumerate(primary_calls, start=1):
        candidates.append(
            {
                "peak_id": f"P{index}",
                "primary_peak_position_bp": call["position_bp"],
                "primary_peak_density": call["density"],
                "primary_prominence": call["prominence"],
                "primary_prominence_threshold": call["prominence_threshold"],
                "primary_prominence_fraction_of_density_range": call[
                    "prominence_fraction_of_density_range"
                ],
                "primary_width_bp": call["width_bp"],
                "primary_left_base_bp": call["left_base_bp"],
                "primary_right_base_bp": call["right_base_bp"],
            }
        )
    print(
        "[peaks] primary candidates: "
        + (", ".join(str(c["primary_peak_position_bp"]) for c in candidates) or "none"),
        flush=True,
    )

    automatic_call_rows: list[dict[str, object]] = []
    sigma_call_lookup: dict[tuple[str, float], list[dict[str, object]]] = {}
    for sigma in sigmas:
        for dataset in DATASETS:
            curve = np.mean(
                curve_matrices[sigma][dataset_indices[dataset]], axis=0
            )
            calls = detect_peaks(
                curve,
                grid,
                args.search_low,
                args.search_high,
                args.peak_distance,
                args.prominence_fraction,
                args.minimum_width,
            )
            sigma_call_lookup[(dataset, sigma)] = calls
            for call in calls:
                automatic_call_rows.append(
                    {
                        "analysis": "sample_weighted_bandwidth",
                        "length_range_bp": "50-1000",
                        "search_range_bp": "100-900",
                        "weighting": "equal sample weight",
                        "dataset": dataset,
                        "gaussian_sigma_bp": sigma,
                        **call,
                    }
                )
    add_candidate_matches(automatic_call_rows, candidates, args.match_window)

    seed_sequence = np.random.SeedSequence(args.seed)
    bootstrap_sequences = seed_sequence.spawn(4)
    bootstrap_records: list[dict[str, object]] = []
    bootstrap_ci: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for dataset, child_seed in zip(DATASETS, bootstrap_sequences[:3]):
        print(
            f"[bootstrap] {dataset}: {args.bootstrap_iterations} sample-level replicates",
            flush=True,
        )
        records, lower, upper = bootstrap_dataset(
            primary_matrix[dataset_indices[dataset]],
            grid,
            dataset,
            candidates,
            args.bootstrap_iterations,
            np.random.default_rng(child_seed),
            args,
        )
        bootstrap_records.extend(records)
        bootstrap_ci[dataset] = (lower, upper)
    bootstrap_summary = summarize_bootstrap(
        bootstrap_records, args.bootstrap_iterations
    )

    # 50-2,000 bp sample-weighted sensitivity analysis.
    sensitivity_grid = np.arange(args.primary_low, args.sensitivity_high + 1)
    sensitivity_matrix, sensitivity_counts = sample_curve_matrix(
        samples,
        args.primary_low,
        args.sensitivity_high,
        args.primary_sigma,
    )
    sensitivity_means = {
        dataset: np.mean(sensitivity_matrix[indices], axis=0)
        for dataset, indices in dataset_indices.items()
    }
    sensitivity_restricted_calls: dict[str, list[dict[str, object]]] = {}
    sensitivity_full_calls: dict[str, list[dict[str, object]]] = {}
    sensitivity_call_rows: list[dict[str, object]] = []
    for dataset in DATASETS:
        restricted = detect_peaks(
            sensitivity_means[dataset],
            sensitivity_grid,
            args.search_low,
            args.search_high,
            args.peak_distance,
            args.prominence_fraction,
            args.minimum_width,
        )
        full = detect_peaks(
            sensitivity_means[dataset],
            sensitivity_grid,
            args.search_low,
            args.sensitivity_high - 100,
            args.peak_distance,
            args.prominence_fraction,
            args.minimum_width,
        )
        sensitivity_restricted_calls[dataset] = restricted
        sensitivity_full_calls[dataset] = full
        for search_label, calls in (
            ("100-900", restricted),
            (f"100-{args.sensitivity_high - 100}", full),
        ):
            for call in calls:
                sensitivity_call_rows.append(
                    {
                        "analysis": "sample_weighted_range_sensitivity",
                        "length_range_bp": f"50-{args.sensitivity_high}",
                        "search_range_bp": search_label,
                        "weighting": "equal sample weight",
                        "dataset": dataset,
                        "gaussian_sigma_bp": args.primary_sigma,
                        **call,
                    }
                )
    add_candidate_matches(sensitivity_call_rows, candidates, args.match_window)
    automatic_call_rows.extend(sensitivity_call_rows)

    # Molecule-pooled sensitivity analysis.
    pooled_curves: dict[str, np.ndarray] = {}
    pooled_counts: dict[str, int] = {}
    pooled_calls: dict[str, list[dict[str, object]]] = {}
    for dataset in DATASETS:
        subset = [samples[index] for index in dataset_indices[dataset]]
        curve, count = pooled_curve(
            subset, args.primary_low, args.primary_high, args.primary_sigma
        )
        calls = detect_peaks(
            curve,
            grid,
            args.search_low,
            args.search_high,
            args.peak_distance,
            args.prominence_fraction,
            args.minimum_width,
        )
        pooled_curves[dataset] = curve
        pooled_counts[dataset] = count
        pooled_calls[dataset] = calls
        for call in calls:
            automatic_call_rows.append(
                {
                    "analysis": "molecule_pooled_sensitivity",
                    "length_range_bp": "50-1000",
                    "search_range_bp": "100-900",
                    "weighting": "each unique eccDNA molecule",
                    "dataset": dataset,
                    "gaussian_sigma_bp": args.primary_sigma,
                    **call,
                }
            )
    add_candidate_matches(automatic_call_rows, candidates, args.match_window)

    # Build the prespecified consensus classification.
    for candidate in candidates:
        peak_id = str(candidate["peak_id"])
        covid_call = nearest_call(
            sigma_call_lookup[(GROUP_COVID, args.primary_sigma)],
            int(candidate["primary_peak_position_bp"]),
            args.match_window,
        )
        hc_call = nearest_call(
            sigma_call_lookup[(GROUP_HC, args.primary_sigma)],
            int(candidate["primary_peak_position_bp"]),
            args.match_window,
        )
        for prefix, call in (("covid", covid_call), ("hc", hc_call)):
            candidate[f"{prefix}_mean_peak_detected"] = call is not None
            candidate[f"{prefix}_mean_peak_position_bp"] = (
                call["position_bp"] if call else None
            )
            candidate[f"{prefix}_mean_peak_prominence"] = (
                call["prominence"] if call else None
            )
            candidate[f"{prefix}_mean_peak_width_bp"] = (
                call["width_bp"] if call else None
            )

        for dataset, prefix in (
            ("Overall", "overall"),
            (GROUP_COVID, "covid"),
            (GROUP_HC, "hc"),
        ):
            summary = bootstrap_summary[(dataset, peak_id)]
            for key, value in summary.items():
                candidate[f"{prefix}_{key}"] = value

        for dataset, prefix in (
            ("Overall", "overall"),
            (GROUP_COVID, "covid"),
            (GROUP_HC, "hc"),
        ):
            detections = 0
            positions = []
            for sigma in sigmas:
                call = nearest_call(
                    sigma_call_lookup[(dataset, sigma)],
                    int(candidate["primary_peak_position_bp"]),
                    args.match_window,
                )
                if call is not None:
                    detections += 1
                    positions.append(int(call["position_bp"]))
            candidate[f"{prefix}_bandwidth_detection_count"] = detections
            candidate[f"{prefix}_bandwidth_positions_bp"] = ",".join(
                str(position) for position in positions
            )

        for dataset, prefix in (
            ("Overall", "overall"),
            (GROUP_COVID, "covid"),
            (GROUP_HC, "hc"),
        ):
            sensitivity_call = nearest_call(
                sensitivity_restricted_calls[dataset],
                int(candidate["primary_peak_position_bp"]),
                args.match_window,
            )
            pooled_call = nearest_call(
                pooled_calls[dataset],
                int(candidate["primary_peak_position_bp"]),
                args.match_window,
            )
            candidate[f"{prefix}_range_50_2000_detected"] = (
                sensitivity_call is not None
            )
            candidate[f"{prefix}_range_50_2000_position_bp"] = (
                sensitivity_call["position_bp"] if sensitivity_call else None
            )
            candidate[f"{prefix}_pooled_detected"] = pooled_call is not None
            candidate[f"{prefix}_pooled_position_bp"] = (
                pooled_call["position_bp"] if pooled_call else None
            )

        failures = []
        if (
            float(candidate["overall_bootstrap_detection_frequency"]) < 0.90
        ):
            failures.append("overall bootstrap support <90%")
        if float(candidate["covid_bootstrap_detection_frequency"]) < 0.80:
            failures.append("COVID-19 bootstrap support <80%")
        if float(candidate["hc_bootstrap_detection_frequency"]) < 0.80:
            failures.append("HC bootstrap support <80%")
        if not candidate["covid_mean_peak_detected"]:
            failures.append("no COVID-19 mean-curve peak within ±25 bp")
        if not candidate["hc_mean_peak_detected"]:
            failures.append("no HC mean-curve peak within ±25 bp")
        if int(candidate["overall_bandwidth_detection_count"]) < 3:
            failures.append("detected under fewer than three bandwidths")
        candidate["stable_consensus_peak"] = not failures
        candidate["stability_failure_reasons"] = "; ".join(failures)

    stable_candidates = [
        candidate for candidate in candidates if candidate["stable_consensus_peak"]
    ]

    # Equal-count sensitivity analysis; it is reported but not used to redefine
    # the prespecified primary consensus classification.
    print(
        f"[subsampling] {args.subsample_repeats} repeats × "
        f"{args.subsample_count} unique eccDNAs/sample",
        flush=True,
    )
    subsample_records, subsample_sample_counts = run_equal_count_subsampling(
        samples,
        candidates,
        args.subsample_repeats,
        args.subsample_count,
        np.random.default_rng(bootstrap_sequences[3]),
        args,
    )
    subsample_summary = summarize_subsampling(
        subsample_records, args.subsample_repeats
    )
    for candidate in candidates:
        peak_id = str(candidate["peak_id"])
        for dataset, prefix in (
            ("Overall", "overall"),
            (GROUP_COVID, "covid"),
            (GROUP_HC, "hc"),
        ):
            summary = subsample_summary[(dataset, peak_id)]
            for key, value in summary.items():
                candidate[f"{prefix}_{key}"] = value

    peak_window_rows, comparison_rows = peak_window_statistics(
        samples,
        stable_candidates,
        args.primary_low,
        args.primary_high,
        args.match_window,
    )
    spacing_rows = compute_spacing_table(stable_candidates, bootstrap_records)
    periodicity_rows = compute_periodicity_summary(
        stable_candidates, bootstrap_records
    )
    for left_candidate, spacing in zip(stable_candidates[:-1], spacing_rows):
        left_candidate["adjacent_spacing_to_next_primary_peak_bp"] = spacing[
            "primary_spacing_bp"
        ]
        left_candidate["adjacent_spacing_bootstrap_median_bp"] = spacing[
            "bootstrap_spacing_median_bp"
        ]
        left_candidate["adjacent_spacing_bootstrap_ci_low_bp"] = spacing[
            "bootstrap_spacing_ci_low_bp"
        ]
        left_candidate["adjacent_spacing_bootstrap_ci_high_bp"] = spacing[
            "bootstrap_spacing_ci_high_bp"
        ]

    # Machine-readable tables.
    automatic_fields = [
        "analysis",
        "length_range_bp",
        "search_range_bp",
        "weighting",
        "dataset",
        "gaussian_sigma_bp",
        "call_rank",
        "position_bp",
        "density",
        "prominence",
        "prominence_threshold",
        "prominence_fraction_of_density_range",
        "width_bp",
        "left_base_bp",
        "right_base_bp",
        "matched_primary_peak_id",
        "distance_to_primary_peak_bp",
    ]
    write_tsv(
        tables_dir / "automatic_peak_calls.tsv",
        automatic_call_rows,
        automatic_fields,
    )
    bootstrap_fields = [
        "dataset",
        "bootstrap_iteration",
        "peak_id",
        "primary_peak_position_bp",
        "detected_within_match_window",
        "matched_position_bp",
        "matched_prominence",
        "matched_width_bp",
        "n_peaks_in_replicate",
    ]
    write_tsv(
        data_dir / "bootstrap_peak_detection_detail.tsv",
        bootstrap_records,
        bootstrap_fields,
    )
    subsample_fields = [
        "subsample_repeat",
        "dataset",
        "peak_id",
        "primary_peak_position_bp",
        "detected_within_match_window",
        "matched_position_bp",
        "matched_prominence",
        "matched_width_bp",
        "n_peaks_in_subsample_curve",
    ]
    write_tsv(
        data_dir / "equal_count_subsample_peak_detection.tsv",
        subsample_records,
        subsample_fields,
    )
    write_tsv(
        data_dir / "equal_count_subsample_eligible_counts.tsv",
        subsample_sample_counts,
        [
            "subsample_repeat",
            "sample_id",
            "group",
            "total_unique_eccdna_sampled",
            "sampled_eccdna_50_1000_bp",
        ],
    )

    curve_rows = []
    for index, position in enumerate(grid):
        curve_rows.append(
            {
                "length_bp": int(position),
                "overall_mean_density": dataset_means["Overall"][index],
                "overall_bootstrap_ci_low": bootstrap_ci["Overall"][0][index],
                "overall_bootstrap_ci_high": bootstrap_ci["Overall"][1][index],
                "covid_mean_density": dataset_means[GROUP_COVID][index],
                "covid_bootstrap_ci_low": bootstrap_ci[GROUP_COVID][0][index],
                "covid_bootstrap_ci_high": bootstrap_ci[GROUP_COVID][1][index],
                "hc_mean_density": dataset_means[GROUP_HC][index],
                "hc_bootstrap_ci_low": bootstrap_ci[GROUP_HC][0][index],
                "hc_bootstrap_ci_high": bootstrap_ci[GROUP_HC][1][index],
                "overall_pooled_density": pooled_curves["Overall"][index],
                "covid_pooled_density": pooled_curves[GROUP_COVID][index],
                "hc_pooled_density": pooled_curves[GROUP_HC][index],
            }
        )
    write_tsv(
        data_dir / "primary_group_curves_and_bootstrap_ci.tsv",
        curve_rows,
        list(curve_rows[0]),
    )

    bandwidth_curve_rows = []
    for sigma in sigmas:
        for index, position in enumerate(grid):
            bandwidth_curve_rows.append(
                {
                    "gaussian_sigma_bp": sigma,
                    "length_bp": int(position),
                    "overall_mean_density": float(
                        np.mean(curve_matrices[sigma], axis=0)[index]
                    ),
                    "covid_mean_density": float(
                        np.mean(
                            curve_matrices[sigma][dataset_indices[GROUP_COVID]],
                            axis=0,
                        )[index]
                    ),
                    "hc_mean_density": float(
                        np.mean(
                            curve_matrices[sigma][dataset_indices[GROUP_HC]],
                            axis=0,
                        )[index]
                    ),
                }
            )
    write_tsv(
        data_dir / "bandwidth_sensitivity_curves.tsv",
        bandwidth_curve_rows,
        list(bandwidth_curve_rows[0]),
    )

    range_curve_rows = []
    for index, position in enumerate(sensitivity_grid):
        range_curve_rows.append(
            {
                "length_bp": int(position),
                "overall_mean_density": sensitivity_means["Overall"][index],
                "covid_mean_density": sensitivity_means[GROUP_COVID][index],
                "hc_mean_density": sensitivity_means[GROUP_HC][index],
            }
        )
    write_tsv(
        data_dir / "range_50_2000_sensitivity_curves.tsv",
        range_curve_rows,
        list(range_curve_rows[0]),
    )

    sample_curve_rows = []
    for sample_index, sample in enumerate(samples):
        for grid_index, position in enumerate(grid):
            sample_curve_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "group": sample.group,
                    "length_bp": int(position),
                    "density_sigma15": primary_matrix[sample_index, grid_index],
                }
            )
    write_tsv(
        data_dir / "sample_level_primary_curves.tsv",
        sample_curve_rows,
        ["sample_id", "group", "length_bp", "density_sigma15"],
    )

    peak_window_fields = [
        "sample_id",
        "group",
        "peak_id",
        "consensus_peak_position_bp",
        "window_low_bp_inclusive",
        "window_high_bp_inclusive",
        "window_eccdna_count",
        "denominator_50_1000_bp",
        "window_proportion",
    ]
    write_tsv(
        data_dir / "peak_window_sample_proportions.tsv",
        peak_window_rows,
        peak_window_fields,
    )
    comparison_fields = [
        "peak_id",
        "consensus_peak_position_bp",
        "window_definition",
        "denominator_definition",
        "n_covid",
        "covid_median_proportion",
        "covid_q1_proportion",
        "covid_q3_proportion",
        "covid_iqr_proportion",
        "n_hc",
        "hc_median_proportion",
        "hc_q1_proportion",
        "hc_q3_proportion",
        "hc_iqr_proportion",
        "median_difference_covid_minus_hc",
        "mann_whitney_u_covid",
        "p_value_two_sided",
        "p_adj_bh",
        "cliffs_delta_covid_minus_hc",
        "test",
        "multiple_testing",
    ]
    write_tsv(
        tables_dir / "peak_window_group_comparisons.tsv",
        comparison_rows,
        comparison_fields,
    )
    spacing_fields = [
        "left_peak_id",
        "right_peak_id",
        "left_primary_position_bp",
        "right_primary_position_bp",
        "primary_spacing_bp",
        "bootstrap_paired_detection_count",
        "bootstrap_spacing_median_bp",
        "bootstrap_spacing_ci_low_bp",
        "bootstrap_spacing_ci_high_bp",
    ]
    write_tsv(
        tables_dir / "adjacent_peak_spacing.tsv", spacing_rows, spacing_fields
    )
    periodicity_fields = [
        "stable_peak_ids",
        "stable_peak_positions_bp",
        "n_stable_peaks",
        "n_adjacent_spacings",
        "primary_adjacent_spacings_bp",
        "primary_mean_spacing_bp",
        "primary_median_spacing_bp",
        "bootstrap_complete_detection_count",
        "bootstrap_mean_spacing_median_bp",
        "bootstrap_mean_spacing_ci_low_bp",
        "bootstrap_mean_spacing_ci_high_bp",
        "bootstrap_median_spacing_median_bp",
        "bootstrap_median_spacing_ci_low_bp",
        "bootstrap_median_spacing_ci_high_bp",
    ]
    write_tsv(
        tables_dir / "periodicity_summary.tsv",
        periodicity_rows,
        periodicity_fields,
    )

    candidate_fields = [
        "peak_id",
        "primary_peak_position_bp",
        "primary_peak_density",
        "primary_prominence",
        "primary_prominence_threshold",
        "primary_prominence_fraction_of_density_range",
        "primary_width_bp",
        "primary_left_base_bp",
        "primary_right_base_bp",
        "covid_mean_peak_detected",
        "covid_mean_peak_position_bp",
        "covid_mean_peak_prominence",
        "covid_mean_peak_width_bp",
        "hc_mean_peak_detected",
        "hc_mean_peak_position_bp",
        "hc_mean_peak_prominence",
        "hc_mean_peak_width_bp",
        "overall_bootstrap_detection_frequency",
        "overall_bootstrap_position_median_bp",
        "overall_bootstrap_position_ci_low_bp",
        "overall_bootstrap_position_ci_high_bp",
        "covid_bootstrap_detection_frequency",
        "covid_bootstrap_position_median_bp",
        "covid_bootstrap_position_ci_low_bp",
        "covid_bootstrap_position_ci_high_bp",
        "hc_bootstrap_detection_frequency",
        "hc_bootstrap_position_median_bp",
        "hc_bootstrap_position_ci_low_bp",
        "hc_bootstrap_position_ci_high_bp",
        "overall_bandwidth_detection_count",
        "overall_bandwidth_positions_bp",
        "covid_bandwidth_detection_count",
        "covid_bandwidth_positions_bp",
        "hc_bandwidth_detection_count",
        "hc_bandwidth_positions_bp",
        "overall_range_50_2000_detected",
        "overall_range_50_2000_position_bp",
        "covid_range_50_2000_detected",
        "covid_range_50_2000_position_bp",
        "hc_range_50_2000_detected",
        "hc_range_50_2000_position_bp",
        "overall_pooled_detected",
        "overall_pooled_position_bp",
        "covid_pooled_detected",
        "covid_pooled_position_bp",
        "hc_pooled_detected",
        "hc_pooled_position_bp",
        "overall_subsample_detection_frequency",
        "overall_subsample_position_median_bp",
        "overall_subsample_position_min_bp",
        "overall_subsample_position_max_bp",
        "covid_subsample_detection_frequency",
        "covid_subsample_position_median_bp",
        "covid_subsample_position_min_bp",
        "covid_subsample_position_max_bp",
        "hc_subsample_detection_frequency",
        "hc_subsample_position_median_bp",
        "hc_subsample_position_min_bp",
        "hc_subsample_position_max_bp",
        "adjacent_spacing_to_next_primary_peak_bp",
        "adjacent_spacing_bootstrap_median_bp",
        "adjacent_spacing_bootstrap_ci_low_bp",
        "adjacent_spacing_bootstrap_ci_high_bp",
        "stable_consensus_peak",
        "stability_failure_reasons",
    ]
    write_tsv(
        tables_dir / "objective_peak_identification_and_stability.tsv",
        candidates,
        candidate_fields,
    )

    # Publication figures.
    make_main_figure(
        figures_dir / "Figure_1D_fragment_length_peaks",
        grid,
        dataset_means[GROUP_COVID],
        dataset_means[GROUP_HC],
        bootstrap_ci[GROUP_COVID],
        bootstrap_ci[GROUP_HC],
        stable_candidates,
    )
    sigma_overall_curves = {
        sigma: np.mean(curve_matrices[sigma], axis=0) for sigma in sigmas
    }
    sigma_overall_calls = {
        sigma: sigma_call_lookup[("Overall", sigma)] for sigma in sigmas
    }
    make_supplementary_figure(
        figures_dir / "Figure_S6_fragment_length_peak_robustness",
        grid,
        sigma_overall_curves,
        sigma_overall_calls,
        candidates,
        subsample_summary,
        peak_window_rows,
        comparison_rows,
    )

    write_analysis_summary(
        output_dir / "analysis_summary.md",
        samples,
        candidates,
        stable_candidates,
        spacing_rows,
        periodicity_rows,
        comparison_rows,
        args,
    )
    write_text_drafts(output_dir, candidates, stable_candidates, args)

    manifest = {
        "analysis": "objective eccDNA fragment-length peak identification",
        "started_and_completed_utc": now_utc(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "metadata": str(args.metadata.resolve()),
        "covid_bed_dir": str(args.covid_bed_dir.resolve()),
        "hc_bed_dir": str(args.hc_bed_dir.resolve()),
        "output_dir": str(output_dir),
        "parameters": {
            "primary_length_range_bp_inclusive": [
                args.primary_low,
                args.primary_high,
            ],
            "sensitivity_length_range_bp_inclusive": [
                args.primary_low,
                args.sensitivity_high,
            ],
            "bin_width_bp": 1,
            "gaussian_sigma_primary_bp": args.primary_sigma,
            "gaussian_sigmas_sensitivity_bp": sigmas,
            "gaussian_boundary_mode": "reflect",
            "gaussian_truncate_sigma": 4.0,
            "evaluation_step_bp": 1,
            "peak_search_range_bp_inclusive": [
                args.search_low,
                args.search_high,
            ],
            "minimum_peak_distance_bp": args.peak_distance,
            "minimum_prominence_fraction_of_full_curve_density_range": (
                args.prominence_fraction
            ),
            "minimum_peak_width_at_half_prominence_bp": args.minimum_width,
            "peak_match_window_bp_inclusive": args.match_window,
            "bootstrap_iterations": args.bootstrap_iterations,
            "bootstrap_unit": "biological sample",
            "bootstrap_position_ci": (
                "2.5th and 97.5th percentiles conditional on peak detection "
                "within the matching window"
            ),
            "overall_bootstrap_stability_threshold": 0.90,
            "group_bootstrap_stability_threshold": 0.80,
            "subsample_repeats": args.subsample_repeats,
            "subsample_unique_eccdna_per_sample": args.subsample_count,
            "subsample_frame": (
                "all unique non-negative-length BED intervals before applying the "
                "50-1,000-bp primary range"
            ),
            "random_seed": args.seed,
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": mpl.__version__,
        },
        "sample_counts": {
            "overall": len(samples),
            "covid": int(len(dataset_indices[GROUP_COVID])),
            "hc": int(len(dataset_indices[GROUP_HC])),
        },
        "primary_candidate_positions_bp": [
            int(candidate["primary_peak_position_bp"]) for candidate in candidates
        ],
        "stable_consensus_positions_bp": [
            int(candidate["primary_peak_position_bp"])
            for candidate in stable_candidates
        ],
        "pooled_eligible_eccdna_counts": pooled_counts,
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    # Internal validation report.
    validation_lines = []

    def record_check(condition: bool, text: str) -> None:
        prefix = "PASS" if condition else "FAIL"
        validation_lines.append(f"{prefix}: {text}")
        if not condition:
            raise AssertionError(text)

    record_check(len(samples) == 78, "78 biological samples loaded")
    record_check(
        Counter(sample.group for sample in samples)
        == Counter({GROUP_COVID: 39, GROUP_HC: 39}),
        "39 COVID-19 and 39 HC samples loaded",
    )
    record_check(
        all(len(sample.lengths) >= args.subsample_count for sample in samples),
        f"every sample has at least {args.subsample_count} unique eccDNAs",
    )
    record_check(
        all(count > 0 for count in eligible_counts),
        "every sample has at least one unique 50-1,000-bp eccDNA",
    )
    record_check(
        all(
            abs(float(np.trapezoid(curve, dx=1.0)) - 1.0) < 1e-10
            for matrix in curve_matrices.values()
            for curve in matrix
        ),
        "every sample-level smoothed density integrates to one",
    )
    record_check(
        np.allclose(
            dataset_means["Overall"],
            (
                dataset_means[GROUP_COVID] + dataset_means[GROUP_HC]
            )
            / 2.0,
            rtol=0,
            atol=1e-14,
        ),
        "overall mean equals the mean of the two equally sized group means",
    )
    record_check(
        len(bootstrap_records)
        == len(candidates) * args.bootstrap_iterations * len(DATASETS),
        "bootstrap detail table has the expected number of candidate records",
    )
    record_check(
        len(subsample_records)
        == len(candidates) * args.subsample_repeats * len(DATASETS),
        "equal-count subsampling detail table has the expected number of candidate records",
    )
    record_check(
        all(
            (figures_dir / f"Figure_1D_fragment_length_peaks.{suffix}").is_file()
            for suffix in ("svg", "pdf", "png", "tiff")
        ),
        "main figure exported as SVG, PDF, PNG and TIFF",
    )
    record_check(
        all(
            (
                figures_dir
                / f"Figure_S6_fragment_length_peak_robustness.{suffix}"
            ).is_file()
            for suffix in ("svg", "pdf", "png", "tiff")
        ),
        "supplementary figure exported as SVG, PDF, PNG and TIFF",
    )
    validation_lines.append(
        f"INFO: fixed algorithm detected {len(candidates)} primary candidate peak(s)"
    )
    validation_lines.append(
        f"INFO: {len(stable_candidates)} candidate peak(s) met all stability criteria"
    )
    (output_dir / "validation_report.txt").write_text(
        "\n".join(validation_lines) + "\n", encoding="utf-8"
    )
    print(f"Analysis completed: {now_utc()}", flush=True)
    print(
        f"Stable consensus peaks: "
        + (
            ", ".join(
                str(candidate["primary_peak_position_bp"])
                for candidate in stable_candidates
            )
            if stable_candidates
            else "none"
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
