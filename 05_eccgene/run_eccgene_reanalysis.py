#!/usr/bin/env python3
"""Reanalyse eccGene abundance and detection frequency.

This script addresses reviewer concerns about treating continuous eccGene
abundance as a Fisher exact-test endpoint. It rebuilds per-sample gene-level
matrices from Circle-Map BED intervals, applies a prespecified gene assignment
rule, and runs separate analyses for continuous abundance and binary detection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import scipy
from scipy.stats import fisher_exact, mannwhitneyu
from scipy.stats.contingency import odds_ratio


CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
DEFINITIONS = ("junction", "interval", "midpoint")
GROUP_COVID = "COVID-19"
GROUP_HC = "HC"
PSEUDO_EA = 1e-6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--covid-bed-dir", required=True, type=Path)
    parser.add_argument("--hc-bed-dir", required=True, type=Path)
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--main-definition", choices=DEFINITIONS, default="junction")
    parser.add_argument("--bin-size", type=int, default=100_000)
    return parser.parse_args()


def normalize_group(raw: str) -> str:
    value = raw.strip().lower()
    if value in {"covid", "covid-19", "covid19", "patient", "patients"}:
        return GROUP_COVID
    if value in {"hc", "healthy", "healthy control", "healthy controls", "normal", "control"}:
        return GROUP_HC
    raise ValueError(f"Cannot normalize group label: {raw!r}")


def strip_bom(name: str) -> str:
    return name.lstrip("\ufeff")


def read_metadata(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in metadata file: {path}")
        fieldnames = {strip_bom(field): field for field in reader.fieldnames}
        sample_col = fieldnames.get("Sample ID") or fieldnames.get("sample") or fieldnames.get("sample_id")
        group_col = fieldnames.get("Group") or fieldnames.get("group")
        if sample_col is None or group_col is None:
            raise ValueError(f"Metadata must include sample and group columns: {path}")
        for row in reader:
            sample_id = row[sample_col].strip()
            if not sample_id:
                continue
            rows.append(
                {
                    "sample_id": sample_id,
                    "group": normalize_group(row[group_col]),
                    "sex": row.get(fieldnames.get("Gender", ""), row.get("sex", "")).strip()
                    if fieldnames.get("Gender", "") in row
                    else row.get("sex", "").strip(),
                    "age": row.get(fieldnames.get("Age", ""), row.get("age", "")).strip()
                    if fieldnames.get("Age", "") in row
                    else row.get("age", "").strip(),
                    "rca_ng_per_ul": row.get(
                        fieldnames.get("RCA concentration (ng/ul)", ""), ""
                    ).strip()
                    if fieldnames.get("RCA concentration (ng/ul)", "") in row
                    else "",
                }
            )
    return rows


def discover_samples(
    metadata_rows: list[dict[str, str]],
    covid_bed_dir: Path,
    hc_bed_dir: Path,
) -> tuple[list[dict[str, str]], list[str]]:
    samples: list[dict[str, str]] = []
    warnings: list[str] = []
    seen: set[str] = set()
    for row in metadata_rows:
        sample_id = row["sample_id"]
        group = row["group"]
        bed_dir = covid_bed_dir if group == GROUP_COVID else hc_bed_dir
        bed_path = bed_dir / f"{sample_id}_circle_site.bed"
        if sample_id in seen:
            warnings.append(f"Duplicate sample in metadata ignored after first occurrence: {sample_id}")
            continue
        seen.add(sample_id)
        if not bed_path.exists():
            warnings.append(f"Missing BED for metadata sample {sample_id}: {bed_path}")
            continue
        sample = dict(row)
        sample["bed_path"] = str(bed_path)
        samples.append(sample)

    for group, bed_dir in ((GROUP_COVID, covid_bed_dir), (GROUP_HC, hc_bed_dir)):
        for bed_path in sorted(bed_dir.glob("*.bed")):
            sample_id = bed_path.name.removesuffix("_circle_site.bed")
            if sample_id not in seen:
                warnings.append(f"BED not present in metadata ({group}): {bed_path}")

    covid_n = sum(1 for sample in samples if sample["group"] == GROUP_COVID)
    hc_n = sum(1 for sample in samples if sample["group"] == GROUP_HC)
    if covid_n == 0 or hc_n == 0:
        raise ValueError(f"Need both groups; found COVID-19={covid_n}, HC={hc_n}")
    return samples, warnings


def file_md5(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_gene_intervals(refgene: Path) -> tuple[list[tuple[str, int, int, str]], dict[str, int], dict[str, int]]:
    by_gene_chrom: dict[tuple[str, str], list[int]] = {}
    skipped = {"short": 0, "chrom": 0, "coordinate": 0, "empty_gene": 0}
    with refgene.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                skipped["short"] += 1
                continue
            chrom = fields[2]
            if chrom not in CANONICAL_CHROMS:
                skipped["chrom"] += 1
                continue
            gene = fields[12].strip()
            if not gene:
                skipped["empty_gene"] += 1
                continue
            try:
                start = int(fields[4])
                end = int(fields[5])
            except ValueError:
                skipped["coordinate"] += 1
                continue
            if end < start:
                start, end = end, start
            if end <= start:
                skipped["coordinate"] += 1
                continue
            key = (gene, chrom)
            if key not in by_gene_chrom:
                by_gene_chrom[key] = [start, end]
            else:
                by_gene_chrom[key][0] = min(by_gene_chrom[key][0], start)
                by_gene_chrom[key][1] = max(by_gene_chrom[key][1], end)

    intervals: list[tuple[str, int, int, str]] = []
    gene_lengths: dict[str, int] = defaultdict(int)
    gene_chrom_counts: dict[str, int] = defaultdict(int)
    for (gene, chrom), (start, end) in by_gene_chrom.items():
        intervals.append((chrom, start, end, gene))
        gene_lengths[gene] += end - start
        gene_chrom_counts[gene] += 1
    intervals.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    return intervals, dict(gene_lengths), dict(gene_chrom_counts)


def build_interval_bins(
    intervals: list[tuple[str, int, int, str]],
    bin_size: int,
) -> dict[str, dict[int, list[int]]]:
    bins: dict[str, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
    for idx, (chrom, start, end, _gene) in enumerate(intervals):
        first_bin = start // bin_size
        last_bin = (end - 1) // bin_size
        for bin_id in range(first_bin, last_bin + 1):
            bins[chrom][bin_id].append(idx)
    return {chrom: dict(chrom_bins) for chrom, chrom_bins in bins.items()}


def genes_for_point(
    chrom: str,
    point: int,
    intervals: list[tuple[str, int, int, str]],
    bins: dict[str, dict[int, list[int]]],
    bin_size: int,
) -> set[str]:
    genes: set[str] = set()
    for idx in bins.get(chrom, {}).get(point // bin_size, []):
        _chrom, start, end, gene = intervals[idx]
        if start <= point < end:
            genes.add(gene)
    return genes


def genes_for_interval(
    chrom: str,
    start: int,
    end: int,
    intervals: list[tuple[str, int, int, str]],
    bins: dict[str, dict[int, list[int]]],
    bin_size: int,
) -> set[str]:
    genes: set[str] = set()
    seen: set[int] = set()
    chrom_bins = bins.get(chrom, {})
    first_bin = start // bin_size
    last_bin = (end - 1) // bin_size
    for bin_id in range(first_bin, last_bin + 1):
        for idx in chrom_bins.get(bin_id, []):
            if idx in seen:
                continue
            seen.add(idx)
            _chrom, gene_start, gene_end, gene = intervals[idx]
            if gene_start < end and gene_end > start:
                genes.add(gene)
    return genes


def process_sample_bed(
    bed_path: Path,
    intervals: list[tuple[str, int, int, str]],
    bins: dict[str, dict[int, list[int]]],
    bin_size: int,
) -> tuple[dict[str, dict[str, int]], dict[str, int]]:
    counts: dict[str, dict[str, int]] = {definition: defaultdict(int) for definition in DEFINITIONS}
    seen_coords: set[tuple[str, int, int]] = set()
    metrics = {
        "bed_lines": 0,
        "valid_canonical_lines": 0,
        "unique_valid_canonical_eccdna": 0,
        "duplicate_unique_junction_lines": 0,
        "noncanonical_chrom_lines": 0,
        "invalid_coordinate_lines": 0,
    }
    with bed_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            metrics["bed_lines"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                metrics["invalid_coordinate_lines"] += 1
                continue
            chrom = fields[0]
            if chrom not in CANONICAL_CHROMS:
                metrics["noncanonical_chrom_lines"] += 1
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                metrics["invalid_coordinate_lines"] += 1
                continue
            if end <= start:
                metrics["invalid_coordinate_lines"] += 1
                continue
            metrics["valid_canonical_lines"] += 1
            coord = (chrom, start, end)
            if coord in seen_coords:
                metrics["duplicate_unique_junction_lines"] += 1
                continue
            seen_coords.add(coord)

            junction_genes = genes_for_point(chrom, start, intervals, bins, bin_size)
            for gene in junction_genes:
                counts["junction"][gene] += 1

            interval_genes = genes_for_interval(chrom, start, end, intervals, bins, bin_size)
            for gene in interval_genes:
                counts["interval"][gene] += 1

            midpoint = (start + end) // 2
            midpoint_genes = genes_for_point(chrom, midpoint, intervals, bins, bin_size)
            for gene in midpoint_genes:
                counts["midpoint"][gene] += 1

    metrics["unique_valid_canonical_eccdna"] = len(seen_coords)
    for definition in DEFINITIONS:
        metrics[f"{definition}_gene_assignments"] = sum(counts[definition].values())
        metrics[f"{definition}_genes_detected"] = len(counts[definition])
    return counts, metrics


def bh_adjust(p_values: Iterable[float]) -> list[float]:
    p_list = list(p_values)
    finite = [(idx, p) for idx, p in enumerate(p_list) if math.isfinite(p)]
    q_values = [math.nan] * len(p_list)
    if not finite:
        return q_values
    finite.sort(key=lambda item: item[1])
    m = len(finite)
    adjusted = [0.0] * m
    running = 1.0
    for rank_from_end, (idx, p) in enumerate(reversed(finite), start=1):
        rank = m - rank_from_end + 1
        running = min(running, p * m / rank)
        adjusted[rank - 1] = min(running, 1.0)
    for (idx, _p), q in zip(finite, adjusted):
        q_values[idx] = q
    return q_values


def mannwhitney_p(x: np.ndarray, y: np.ndarray) -> float:
    try:
        return float(mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").pvalue)
    except TypeError:
        return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except ValueError:
        return math.nan


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    if x.size == 0 or y.size == 0:
        return math.nan
    comparisons = x[:, None] - y[None, :]
    return float((np.sum(comparisons > 0) - np.sum(comparisons < 0)) / comparisons.size)


def conditional_odds_ratio_ci(table: list[list[int]]) -> tuple[float, float, float]:
    try:
        result = odds_ratio(table, kind="conditional")
        ci = result.confidence_interval(confidence_level=0.95)
        return float(result.statistic), float(ci.low), float(ci.high)
    except Exception:
        return math.nan, math.nan, math.nan


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "NA"
        if math.isinf(float(value)):
            return "Inf" if value > 0 else "-Inf"
        return f"{float(value):.10g}"
    return str(value)


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fieldnames})


def write_matrix(
    path: Path,
    genes: list[str],
    sample_ids: list[str],
    value_getter,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Gene", *sample_ids])
        for gene in genes:
            writer.writerow([gene, *[format_value(value_getter(gene, idx)) for idx in range(len(sample_ids))]])


def compute_definition_results(
    definition: str,
    counts_by_gene: dict[str, list[int]],
    gene_lengths: dict[str, int],
    samples: list[dict[str, str]],
    out_matrix_dir: Path,
    out_result_dir: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[float]]:
    sample_ids = [sample["sample_id"] for sample in samples]
    covid_idx = [idx for idx, sample in enumerate(samples) if sample["group"] == GROUP_COVID]
    hc_idx = [idx for idx, sample in enumerate(samples) if sample["group"] == GROUP_HC]
    genes = sorted(counts_by_gene)

    denominators = [0.0] * len(samples)
    for gene, counts in counts_by_gene.items():
        length = gene_lengths[gene]
        for idx, count in enumerate(counts):
            if count:
                denominators[idx] += count / length

    def count_value(gene: str, idx: int) -> int:
        return counts_by_gene[gene][idx]

    def detected_value(gene: str, idx: int) -> int:
        return 1 if counts_by_gene[gene][idx] > 0 else 0

    def ea_value(gene: str, idx: int) -> float:
        denominator = denominators[idx]
        if denominator <= 0:
            return 0.0
        return (counts_by_gene[gene][idx] / gene_lengths[gene]) / denominator * 1_000_000.0

    write_matrix(out_matrix_dir / f"{definition}_counts.tsv", genes, sample_ids, count_value)
    write_matrix(out_matrix_dir / f"{definition}_detected.tsv", genes, sample_ids, detected_value)
    write_matrix(out_matrix_dir / f"{definition}_EA.tsv", genes, sample_ids, ea_value)

    abundance_rows: list[dict[str, object]] = []
    log1p_rows: list[dict[str, object]] = []
    detection_rows: list[dict[str, object]] = []

    for gene in genes:
        length = gene_lengths[gene]
        counts = counts_by_gene[gene]
        ea = np.array([ea_value(gene, idx) for idx in range(len(samples))], dtype=float)
        covid = ea[covid_idx]
        hc = ea[hc_idx]
        p_wilcox = mannwhitney_p(covid, hc)
        delta = cliffs_delta(covid, hc)
        mean_covid = float(np.mean(covid))
        mean_hc = float(np.mean(hc))
        log2fc = math.log2((mean_covid + PSEUDO_EA) / (mean_hc + PSEUDO_EA))
        covid_detected = int(sum(1 for idx in covid_idx if counts[idx] > 0))
        hc_detected = int(sum(1 for idx in hc_idx if counts[idx] > 0))

        abundance_rows.append(
            {
                "Gene": gene,
                "gene_length_bp": length,
                "covid_n": len(covid_idx),
                "hc_n": len(hc_idx),
                "mean_EA_COVID": mean_covid,
                "mean_EA_HC": mean_hc,
                "median_EA_COVID": float(np.median(covid)),
                "median_EA_HC": float(np.median(hc)),
                "log2FC_mean_EA_COVID_vs_HC": log2fc,
                "cliffs_delta_COVID_vs_HC": delta,
                "covid_detected": covid_detected,
                "hc_detected": hc_detected,
                "wilcoxon_p": p_wilcox,
            }
        )

        covid_log1p = np.log1p(covid)
        hc_log1p = np.log1p(hc)
        log1p_rows.append(
            {
                "Gene": gene,
                "gene_length_bp": length,
                "mean_log1p_EA_COVID": float(np.mean(covid_log1p)),
                "mean_log1p_EA_HC": float(np.mean(hc_log1p)),
                "median_log1p_EA_COVID": float(np.median(covid_log1p)),
                "median_log1p_EA_HC": float(np.median(hc_log1p)),
                "mean_log1p_difference_COVID_minus_HC": float(np.mean(covid_log1p) - np.mean(hc_log1p)),
                "cliffs_delta_COVID_vs_HC": delta,
                "wilcoxon_p_log1p": mannwhitney_p(covid_log1p, hc_log1p),
            }
        )

        covid_not_detected = len(covid_idx) - covid_detected
        hc_not_detected = len(hc_idx) - hc_detected
        table = [[covid_detected, covid_not_detected], [hc_detected, hc_not_detected]]
        _sample_or, fisher_p = fisher_exact(table, alternative="two-sided")
        conditional_or, ci_low, ci_high = conditional_odds_ratio_ci(table)
        detection_rows.append(
            {
                "Gene": gene,
                "gene_length_bp": length,
                "covid_detected": covid_detected,
                "covid_not_detected": covid_not_detected,
                "hc_detected": hc_detected,
                "hc_not_detected": hc_not_detected,
                "covid_detection_frequency": covid_detected / len(covid_idx),
                "hc_detection_frequency": hc_detected / len(hc_idx),
                "frequency_difference_COVID_minus_HC": covid_detected / len(covid_idx)
                - hc_detected / len(hc_idx),
                "odds_ratio_conditional_COVID_vs_HC": conditional_or,
                "odds_ratio_95CI_low": ci_low,
                "odds_ratio_95CI_high": ci_high,
                "fisher_p": float(fisher_p),
            }
        )

    abundance_q = bh_adjust(row["wilcoxon_p"] for row in abundance_rows)
    for row, q_value in zip(abundance_rows, abundance_q):
        row["wilcoxon_q_BH"] = q_value
        row["direction"] = "COVID_up" if row["log2FC_mean_EA_COVID_vs_HC"] > 0 else "HC_up"

    log1p_q = bh_adjust(row["wilcoxon_p_log1p"] for row in log1p_rows)
    for row, q_value in zip(log1p_rows, log1p_q):
        row["wilcoxon_q_BH_log1p"] = q_value
        row["direction"] = "COVID_up" if row["mean_log1p_difference_COVID_minus_HC"] > 0 else "HC_up"

    detection_q = bh_adjust(row["fisher_p"] for row in detection_rows)
    for row, q_value in zip(detection_rows, detection_q):
        row["fisher_q_BH"] = q_value
        row["direction"] = "COVID_up" if row["frequency_difference_COVID_minus_HC"] > 0 else "HC_up"

    abundance_fields = [
        "Gene",
        "gene_length_bp",
        "covid_n",
        "hc_n",
        "mean_EA_COVID",
        "mean_EA_HC",
        "median_EA_COVID",
        "median_EA_HC",
        "log2FC_mean_EA_COVID_vs_HC",
        "cliffs_delta_COVID_vs_HC",
        "covid_detected",
        "hc_detected",
        "wilcoxon_p",
        "wilcoxon_q_BH",
        "direction",
    ]
    log1p_fields = [
        "Gene",
        "gene_length_bp",
        "mean_log1p_EA_COVID",
        "mean_log1p_EA_HC",
        "median_log1p_EA_COVID",
        "median_log1p_EA_HC",
        "mean_log1p_difference_COVID_minus_HC",
        "cliffs_delta_COVID_vs_HC",
        "wilcoxon_p_log1p",
        "wilcoxon_q_BH_log1p",
        "direction",
    ]
    detection_fields = [
        "Gene",
        "gene_length_bp",
        "covid_detected",
        "covid_not_detected",
        "hc_detected",
        "hc_not_detected",
        "covid_detection_frequency",
        "hc_detection_frequency",
        "frequency_difference_COVID_minus_HC",
        "odds_ratio_conditional_COVID_vs_HC",
        "odds_ratio_95CI_low",
        "odds_ratio_95CI_high",
        "fisher_p",
        "fisher_q_BH",
        "direction",
    ]

    abundance_rows_sorted = sorted(
        abundance_rows,
        key=lambda row: (
            row["wilcoxon_q_BH"] if math.isfinite(row["wilcoxon_q_BH"]) else math.inf,
            row["wilcoxon_p"] if math.isfinite(row["wilcoxon_p"]) else math.inf,
            row["Gene"],
        ),
    )
    log1p_rows_sorted = sorted(
        log1p_rows,
        key=lambda row: (
            row["wilcoxon_q_BH_log1p"] if math.isfinite(row["wilcoxon_q_BH_log1p"]) else math.inf,
            row["wilcoxon_p_log1p"] if math.isfinite(row["wilcoxon_p_log1p"]) else math.inf,
            row["Gene"],
        ),
    )
    detection_rows_sorted = sorted(
        detection_rows,
        key=lambda row: (
            row["fisher_q_BH"] if math.isfinite(row["fisher_q_BH"]) else math.inf,
            row["fisher_p"] if math.isfinite(row["fisher_p"]) else math.inf,
            row["Gene"],
        ),
    )
    write_tsv(out_result_dir / f"{definition}_abundance_wilcoxon.tsv", abundance_fields, abundance_rows_sorted)
    write_tsv(
        out_result_dir / f"{definition}_abundance_log1p_sensitivity.tsv",
        log1p_fields,
        log1p_rows_sorted,
    )
    write_tsv(out_result_dir / f"{definition}_detection_fisher.tsv", detection_fields, detection_rows_sorted)

    abundance_sig = [row for row in abundance_rows_sorted if row["wilcoxon_q_BH"] < 0.05]
    detection_sig = [row for row in detection_rows_sorted if row["fisher_q_BH"] < 0.05]
    write_tsv(out_result_dir / f"{definition}_significant_abundance_q005.tsv", abundance_fields, abundance_sig)
    write_tsv(out_result_dir / f"{definition}_significant_detection_q005.tsv", detection_fields, detection_sig)

    return abundance_rows, log1p_rows, detection_rows, abundance_sig, denominators


def write_sample_metadata(path: Path, samples: list[dict[str, str]]) -> None:
    fields = ["sample_id", "group", "sex", "age", "rca_ng_per_ul", "bed_path"]
    write_tsv(path, fields, samples)


def write_annotation_outputs(
    annotation_dir: Path,
    intervals: list[tuple[str, int, int, str]],
    gene_lengths: dict[str, int],
    gene_chrom_counts: dict[str, int],
) -> None:
    interval_rows = [
        {
            "chrom": chrom,
            "start_0based": start,
            "end_0based_exclusive": end,
            "gene": gene,
            "interval_length_bp": end - start,
            "gene_length_bp": gene_lengths[gene],
            "gene_chrom_interval_count": gene_chrom_counts[gene],
        }
        for chrom, start, end, gene in intervals
    ]
    write_tsv(
        annotation_dir / "gene_intervals_refGene_merged.tsv",
        [
            "chrom",
            "start_0based",
            "end_0based_exclusive",
            "gene",
            "interval_length_bp",
            "gene_length_bp",
            "gene_chrom_interval_count",
        ],
        interval_rows,
    )
    length_rows = [
        {
            "Gene": gene,
            "gene_length_bp": length,
            "gene_chrom_interval_count": gene_chrom_counts[gene],
        }
        for gene, length in sorted(gene_lengths.items())
    ]
    write_tsv(
        annotation_dir / "gene_lengths_refGene_merged.tsv",
        ["Gene", "gene_length_bp", "gene_chrom_interval_count"],
        length_rows,
    )


def summarize_definitions(
    result_dir: Path,
    counts_all: dict[str, dict[str, list[int]]],
    abundance_all: dict[str, list[dict[str, object]]],
    detection_all: dict[str, list[dict[str, object]]],
    denominators_all: dict[str, list[float]],
    samples: list[dict[str, str]],
) -> None:
    rows: list[dict[str, object]] = []
    for definition in DEFINITIONS:
        abundance_rows = abundance_all[definition]
        detection_rows = detection_all[definition]
        denom = np.array(denominators_all[definition], dtype=float)
        rows.append(
            {
                "definition": definition,
                "genes_detected_any": len(counts_all[definition]),
                "total_gene_assignments": int(sum(sum(values) for values in counts_all[definition].values())),
                "samples": len(samples),
                "covid_samples": sum(1 for sample in samples if sample["group"] == GROUP_COVID),
                "hc_samples": sum(1 for sample in samples if sample["group"] == GROUP_HC),
                "ea_denominator_min": float(np.min(denom)),
                "ea_denominator_median": float(np.median(denom)),
                "ea_denominator_max": float(np.max(denom)),
                "significant_abundance_q005": sum(1 for row in abundance_rows if row["wilcoxon_q_BH"] < 0.05),
                "significant_abundance_q005_COVID_up": sum(
                    1
                    for row in abundance_rows
                    if row["wilcoxon_q_BH"] < 0.05 and row["direction"] == "COVID_up"
                ),
                "significant_detection_q005": sum(1 for row in detection_rows if row["fisher_q_BH"] < 0.05),
                "significant_detection_q005_COVID_up": sum(
                    1
                    for row in detection_rows
                    if row["fisher_q_BH"] < 0.05 and row["direction"] == "COVID_up"
                ),
            }
        )
    write_tsv(
        result_dir / "definition_summary.tsv",
        [
            "definition",
            "genes_detected_any",
            "total_gene_assignments",
            "samples",
            "covid_samples",
            "hc_samples",
            "ea_denominator_min",
            "ea_denominator_median",
            "ea_denominator_max",
            "significant_abundance_q005",
            "significant_abundance_q005_COVID_up",
            "significant_detection_q005",
            "significant_detection_q005_COVID_up",
        ],
        rows,
    )


def significant_set(rows: list[dict[str, object]], q_key: str, direction: str | None = None) -> set[str]:
    result = set()
    for row in rows:
        if row[q_key] < 0.05 and (direction is None or row["direction"] == direction):
            result.add(str(row["Gene"]))
    return result


def write_overlap_summary(
    result_dir: Path,
    abundance_all: dict[str, list[dict[str, object]]],
    detection_all: dict[str, list[dict[str, object]]],
) -> None:
    rows: list[dict[str, object]] = []
    for endpoint, data, q_key in (
        ("abundance_wilcoxon", abundance_all, "wilcoxon_q_BH"),
        ("detection_fisher", detection_all, "fisher_q_BH"),
    ):
        for set_filter, direction in (("all_significant", None), ("covid_up_significant", "COVID_up")):
            sets = {definition: significant_set(data[definition], q_key, direction) for definition in DEFINITIONS}
            for i, definition_a in enumerate(DEFINITIONS):
                for definition_b in DEFINITIONS[i + 1 :]:
                    set_a = sets[definition_a]
                    set_b = sets[definition_b]
                    union = set_a | set_b
                    intersection = set_a & set_b
                    rows.append(
                        {
                            "endpoint": endpoint,
                            "set_filter": set_filter,
                            "definition_a": definition_a,
                            "definition_b": definition_b,
                            "n_a": len(set_a),
                            "n_b": len(set_b),
                            "n_intersection": len(intersection),
                            "n_union": len(union),
                            "jaccard": len(intersection) / len(union) if union else math.nan,
                        }
                    )
    write_tsv(
        result_dir / "definition_overlap_summary.tsv",
        [
            "endpoint",
            "set_filter",
            "definition_a",
            "definition_b",
            "n_a",
            "n_b",
            "n_intersection",
            "n_union",
            "jaccard",
        ],
        rows,
    )


def write_candidate_stability(
    result_dir: Path,
    abundance_all: dict[str, list[dict[str, object]]],
    detection_all: dict[str, list[dict[str, object]]],
    main_definition: str,
) -> None:
    abundance_by_def = {
        definition: {str(row["Gene"]): row for row in rows} for definition, rows in abundance_all.items()
    }
    detection_by_def = {
        definition: {str(row["Gene"]): row for row in rows} for definition, rows in detection_all.items()
    }
    main_candidates: dict[str, set[str]] = defaultdict(set)
    for row in abundance_all[main_definition]:
        if row["wilcoxon_q_BH"] < 0.05:
            main_candidates[str(row["Gene"])].add("abundance")
    for row in detection_all[main_definition]:
        if row["fisher_q_BH"] < 0.05:
            main_candidates[str(row["Gene"])].add("detection")

    rows: list[dict[str, object]] = []
    for gene in sorted(main_candidates):
        row: dict[str, object] = {
            "Gene": gene,
            "main_definition": main_definition,
            "main_significant_endpoints": ",".join(sorted(main_candidates[gene])),
        }
        for definition in DEFINITIONS:
            abundance = abundance_by_def[definition].get(gene)
            detection = detection_by_def[definition].get(gene)
            row[f"{definition}_abundance_q"] = abundance.get("wilcoxon_q_BH") if abundance else math.nan
            row[f"{definition}_abundance_log2FC"] = (
                abundance.get("log2FC_mean_EA_COVID_vs_HC") if abundance else math.nan
            )
            row[f"{definition}_abundance_direction"] = abundance.get("direction") if abundance else "not_detected"
            row[f"{definition}_detection_q"] = detection.get("fisher_q_BH") if detection else math.nan
            row[f"{definition}_detection_OR"] = (
                detection.get("odds_ratio_conditional_COVID_vs_HC") if detection else math.nan
            )
            row[f"{definition}_detection_direction"] = detection.get("direction") if detection else "not_detected"

        main_abundance_direction = row.get(f"{main_definition}_abundance_direction")
        main_detection_direction = row.get(f"{main_definition}_detection_direction")
        row["abundance_significant_same_direction_all_definitions"] = all(
            row.get(f"{definition}_abundance_q", math.nan) < 0.05
            and row.get(f"{definition}_abundance_direction") == main_abundance_direction
            for definition in DEFINITIONS
        )
        row["detection_significant_same_direction_all_definitions"] = all(
            row.get(f"{definition}_detection_q", math.nan) < 0.05
            and row.get(f"{definition}_detection_direction") == main_detection_direction
            for definition in DEFINITIONS
        )
        rows.append(row)

    fields = [
        "Gene",
        "main_definition",
        "main_significant_endpoints",
        "junction_abundance_q",
        "junction_abundance_log2FC",
        "junction_abundance_direction",
        "junction_detection_q",
        "junction_detection_OR",
        "junction_detection_direction",
        "interval_abundance_q",
        "interval_abundance_log2FC",
        "interval_abundance_direction",
        "interval_detection_q",
        "interval_detection_OR",
        "interval_detection_direction",
        "midpoint_abundance_q",
        "midpoint_abundance_log2FC",
        "midpoint_abundance_direction",
        "midpoint_detection_q",
        "midpoint_detection_OR",
        "midpoint_detection_direction",
        "abundance_significant_same_direction_all_definitions",
        "detection_significant_same_direction_all_definitions",
    ]
    write_tsv(result_dir / "candidate_stability.tsv", fields, rows)


def write_methods_and_validation(
    output_dir: Path,
    args: argparse.Namespace,
    samples: list[dict[str, str]],
    metadata_warnings: list[str],
    intervals: list[tuple[str, int, int, str]],
    gene_lengths: dict[str, int],
    gene_chrom_counts: dict[str, int],
    sample_metric_rows: list[dict[str, object]],
    elapsed_seconds: float,
) -> None:
    result_dir = output_dir / "results"
    method_lines = [
        "Revised eccGene definition and differential analysis",
        "",
        "Locked annotation and assignment rules:",
        f"- Reference annotation: UCSC refGene table on hg38 from {args.refgene}.",
        f"- refGene MD5: {file_md5(args.refgene)}.",
        "- Canonical chromosomes retained: chr1-chr22, chrX, chrY, chrM.",
        "- Transcript handling: all refGene transcript txStart-txEnd spans were collapsed to one",
        "  bounding gene-body interval per gene symbol per chromosome.",
        "- Gene length: sum of collapsed gene-body interval lengths across retained chromosomes",
        "  for a gene symbol.",
        "- Per-sample eccDNA de-duplication: unique (chrom, start, end) BED coordinates.",
        "- Main eccGene definition: junction. A Circle-Map BED start coordinate is represented",
        "  as a 1 bp half-open interval [start, start+1) and assigned to genes whose merged",
        "  gene body contains that coordinate.",
        "- Sensitivity definitions: interval, any full eccDNA interval/gene-body intersection;",
        "  midpoint, floor((start+end)/2) contained in gene body.",
        "- Multi-gene assignment: if one eccDNA maps to multiple genes under a definition,",
        "  it contributes one isoform count to each assigned gene. Counts are not fractional.",
        "- Multiple eccDNA isoforms from the same gene: each distinct per-sample (chrom,start,end)",
        "  assigned to that gene contributes one count.",
        "- EA formula for sample s, gene g and definition d:",
        "  EA_sg = ((X_sg / L_g) / sum_h(X_sh / L_h)) * 1e6, where X is the unique",
        "  assigned eccDNA isoform count and L is the merged gene length. The denominator",
        "  is the per-sample sum across all genes detected under the same definition.",
        "- Abundance endpoint: per-sample EA values, compared between COVID-19 and HC by",
        "  two-sided Wilcoxon rank-sum test; BH-FDR reported across genes. Effect sizes:",
        "  log2 fold change of mean EA with pseudocount 1e-6 EA units, and Cliff's delta",
        "  for COVID-19 versus HC.",
        "- log1p(EA) sensitivity: Wilcoxon rank-sum test on log1p(EA), with BH-FDR. Because",
        "  log1p is monotonic, the rank-test p-values are expected to be rank-equivalent",
        "  to the raw EA test; group log1p means/medians are reported separately.",
        "- Detection endpoint: detected=1 if X_sg>0 in a sample, otherwise 0. Fisher's exact",
        "  test is applied only to this binary endpoint. Conditional odds ratio and exact",
        "  95% CI are reported with BH-FDR across genes.",
        "",
        f"Samples analysed: {len(samples)} total; "
        f"{sum(1 for sample in samples if sample['group'] == GROUP_COVID)} COVID-19 and "
        f"{sum(1 for sample in samples if sample['group'] == GROUP_HC)} HC.",
        f"Collapsed gene intervals: {len(intervals)} intervals for {len(gene_lengths)} gene symbols.",
        f"Genes on more than one retained chromosome: {sum(1 for count in gene_chrom_counts.values() if count > 1)}.",
        f"Elapsed runtime: {elapsed_seconds:.1f} seconds.",
    ]
    (result_dir / "eccgene_analysis_summary.txt").write_text("\n".join(method_lines) + "\n", encoding="utf-8")

    validation_lines = []
    expected_samples = 78
    covid_n = sum(1 for sample in samples if sample["group"] == GROUP_COVID)
    hc_n = sum(1 for sample in samples if sample["group"] == GROUP_HC)
    validation_lines.append(
        "PASS: both groups present"
        if covid_n > 0 and hc_n > 0
        else "FAIL: one or both groups are absent"
    )
    validation_lines.append(
        f"PASS: analysed {expected_samples} samples"
        if len(samples) == expected_samples
        else f"WARN: analysed {len(samples)} samples, expected {expected_samples}"
    )
    validation_lines.append(
        "PASS: group sizes are 39 COVID-19 and 39 HC"
        if covid_n == 39 and hc_n == 39
        else f"WARN: group sizes are COVID-19={covid_n}, HC={hc_n}"
    )
    validation_lines.append(
        "PASS: no metadata/BED discovery warnings"
        if not metadata_warnings
        else f"WARN: {len(metadata_warnings)} metadata/BED discovery warnings"
    )
    total_unique = sum(int(row["unique_valid_canonical_eccdna"]) for row in sample_metric_rows)
    validation_lines.append(
        "PASS: unique canonical eccDNA records were parsed"
        if total_unique > 0
        else "FAIL: no unique canonical eccDNA records were parsed"
    )
    validation_lines.extend(f"WARN: {warning}" for warning in metadata_warnings)
    (result_dir / "validation_report.txt").write_text("\n".join(validation_lines) + "\n", encoding="utf-8")


def write_run_metadata(
    output_dir: Path,
    args: argparse.Namespace,
    samples: list[dict[str, str]],
    intervals: list[tuple[str, int, int, str]],
    gene_lengths: dict[str, int],
) -> None:
    result_dir = output_dir / "results"
    software_rows = [
        {"name": "python", "version": sys.version.split()[0]},
        {"name": "numpy", "version": np.__version__},
        {"name": "scipy", "version": scipy.__version__},
        {"name": "platform", "version": platform.platform()},
    ]
    write_tsv(result_dir / "software_versions.tsv", ["name", "version"], software_rows)

    manifest_rows = [
        {
            "file_type": "metadata",
            "path": str(args.metadata),
            "size_bytes": args.metadata.stat().st_size,
            "mtime_epoch": args.metadata.stat().st_mtime,
            "md5": file_md5(args.metadata),
        },
        {
            "file_type": "refgene",
            "path": str(args.refgene),
            "size_bytes": args.refgene.stat().st_size,
            "mtime_epoch": args.refgene.stat().st_mtime,
            "md5": file_md5(args.refgene),
        },
    ]
    for sample in samples:
        bed_path = Path(sample["bed_path"])
        stat = bed_path.stat()
        manifest_rows.append(
            {
                "file_type": "circlemap_bed",
                "path": str(bed_path),
                "sample_id": sample["sample_id"],
                "group": sample["group"],
                "size_bytes": stat.st_size,
                "mtime_epoch": stat.st_mtime,
                "md5": "not_computed_large_input",
            }
        )
    write_tsv(
        result_dir / "input_manifest.tsv",
        ["file_type", "sample_id", "group", "path", "size_bytes", "mtime_epoch", "md5"],
        manifest_rows,
    )

    parameters = {
        "metadata": str(args.metadata),
        "covid_bed_dir": str(args.covid_bed_dir),
        "hc_bed_dir": str(args.hc_bed_dir),
        "refgene": str(args.refgene),
        "refgene_md5": file_md5(args.refgene),
        "output_dir": str(args.output_dir),
        "main_definition": args.main_definition,
        "sensitivity_definitions": list(DEFINITIONS),
        "bin_size": args.bin_size,
        "pseudo_ea_for_log2fc": PSEUDO_EA,
        "sample_count": len(samples),
        "gene_interval_count": len(intervals),
        "gene_symbol_count": len(gene_lengths),
    }
    (result_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    start_time = time.time()
    args = parse_args()
    output_dir = args.output_dir
    annotation_dir = output_dir / "annotation"
    matrix_dir = output_dir / "matrices"
    result_dir = output_dir / "results"
    for path in (annotation_dir, matrix_dir, result_dir, output_dir / "logs"):
        path.mkdir(parents=True, exist_ok=True)

    metadata_rows = read_metadata(args.metadata)
    samples, metadata_warnings = discover_samples(metadata_rows, args.covid_bed_dir, args.hc_bed_dir)
    write_sample_metadata(output_dir / "sample_metadata.tsv", samples)

    intervals, gene_lengths, gene_chrom_counts = build_gene_intervals(args.refgene)
    if not intervals:
        raise ValueError("No gene intervals were built from refGene")
    write_annotation_outputs(annotation_dir, intervals, gene_lengths, gene_chrom_counts)
    bins = build_interval_bins(intervals, args.bin_size)

    sample_ids = [sample["sample_id"] for sample in samples]
    counts_all: dict[str, dict[str, list[int]]] = {definition: {} for definition in DEFINITIONS}
    sample_metric_rows: list[dict[str, object]] = []

    for sample_idx, sample in enumerate(samples):
        bed_path = Path(sample["bed_path"])
        print(
            f"[{sample_idx + 1}/{len(samples)}] processing {sample['sample_id']} "
            f"({sample['group']}) from {bed_path}",
            flush=True,
        )
        sample_counts, metrics = process_sample_bed(bed_path, intervals, bins, args.bin_size)
        metric_row: dict[str, object] = {
            "sample_id": sample["sample_id"],
            "group": sample["group"],
            "bed_path": str(bed_path),
        }
        metric_row.update(metrics)
        sample_metric_rows.append(metric_row)
        for definition in DEFINITIONS:
            for gene, count in sample_counts[definition].items():
                if gene not in counts_all[definition]:
                    counts_all[definition][gene] = [0] * len(samples)
                counts_all[definition][gene][sample_idx] = count

    metric_fields = [
        "sample_id",
        "group",
        "bed_path",
        "bed_lines",
        "valid_canonical_lines",
        "unique_valid_canonical_eccdna",
        "duplicate_unique_junction_lines",
        "noncanonical_chrom_lines",
        "invalid_coordinate_lines",
        "junction_gene_assignments",
        "junction_genes_detected",
        "interval_gene_assignments",
        "interval_genes_detected",
        "midpoint_gene_assignments",
        "midpoint_genes_detected",
    ]
    write_tsv(result_dir / "sample_eccgene_metrics.tsv", metric_fields, sample_metric_rows)

    abundance_all: dict[str, list[dict[str, object]]] = {}
    log1p_all: dict[str, list[dict[str, object]]] = {}
    detection_all: dict[str, list[dict[str, object]]] = {}
    denominators_all: dict[str, list[float]] = {}

    for definition in DEFINITIONS:
        print(f"[stats] computing {definition}", flush=True)
        abundance_rows, log1p_rows, detection_rows, _abundance_sig, denominators = compute_definition_results(
            definition,
            counts_all[definition],
            gene_lengths,
            samples,
            matrix_dir,
            result_dir,
        )
        abundance_all[definition] = abundance_rows
        log1p_all[definition] = log1p_rows
        detection_all[definition] = detection_rows
        denominators_all[definition] = denominators

    summarize_definitions(result_dir, counts_all, abundance_all, detection_all, denominators_all, samples)
    write_overlap_summary(result_dir, abundance_all, detection_all)
    write_candidate_stability(result_dir, abundance_all, detection_all, args.main_definition)
    write_run_metadata(output_dir, args, samples, intervals, gene_lengths)
    write_methods_and_validation(
        output_dir,
        args,
        samples,
        metadata_warnings,
        intervals,
        gene_lengths,
        gene_chrom_counts,
        sample_metric_rows,
        time.time() - start_time,
    )
    print(f"[done] revised eccGene analysis written to {output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
