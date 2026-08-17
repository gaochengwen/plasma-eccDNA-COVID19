#!/usr/bin/env python3
"""Build methods-consistent and strict Circle-Map sensitivity callsets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


SPLIT_THRESHOLDS = (3, 5, 10)
SCORE_THRESHOLDS = (200.0, 500.0, 1000.0)
OUTPUT_CALLSETS = {
    "methods": (3, 200.0),
    "strict": (10, 1000.0),
    "very_strict": (20, 2000.0),
}
LENGTH_BINS = (
    ("lt200", 0, 200),
    ("200_399", 200, 400),
    ("400_599", 400, 600),
    ("600_999", 600, 1000),
    ("1000_1999", 1000, 2000),
    ("ge2000", 2000, math.inf),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def length_bin(length: int) -> str:
    for name, left, right in LENGTH_BINS:
        if left <= length < right:
            return name
    raise AssertionError(length)


def coverage_pass(values: list[str]) -> bool:
    mean_coverage = float(values[6])
    coverage_sd = float(values[7])
    start_increase = float(values[8])
    end_increase = float(values[9])
    uncovered_fraction = float(values[10])
    return (
        mean_coverage > coverage_sd
        and start_increase >= 0.3
        and end_increase >= 0.3
        and uncovered_fraction < 0.1
    )


def parse_call(line: str) -> tuple[list[str], int, float, int] | None:
    values = line.rstrip("\n").split("\t")
    if len(values) < 11 or values[0].lower() in {"chrom", "chromosome"}:
        return None
    try:
        start = int(values[1])
        end = int(values[2])
        split_reads = int(float(values[4]))
        score = float(values[5])
    except ValueError:
        return None
    if start < 0 or end <= start:
        return None
    return values, split_reads, score, end - start


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_id", "group", "eccdna_bed", "mapped_alignments_idxstats"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Manifest lacks required fields: {path}")
    if len(rows) != 78 or len({row["sample_id"] for row in rows}) != 78:
        raise ValueError("Expected exactly 78 unique samples")
    return rows


def process_sample(row: dict[str, str], outdir: Path) -> tuple[list[dict], list[dict]]:
    sample = row["sample_id"]
    group = row["group"]
    bed = Path(row["eccdna_bed"])
    mapped = int(row["mapped_alignments_idxstats"])
    if not bed.is_file():
        raise FileNotFoundError(bed)

    handles = {}
    output_paths = {}
    for name in OUTPUT_CALLSETS:
        directory = outdir / "callsets" / f"circlemap_{name}"
        directory.mkdir(parents=True, exist_ok=True)
        output_path = directory / f"{sample}_circle_site.bed"
        handles[name] = output_path.open("w")
        output_paths[name] = output_path

    current_count = 0
    invalid_or_header = 0
    coverage_eligible_count = 0
    threshold_counts: Counter[tuple[int, float]] = Counter()
    threshold_lengths: Counter[tuple[int, float, str]] = Counter()
    output_counts: Counter[str] = Counter()
    output_lengths: Counter[tuple[str, str]] = Counter()

    try:
        with bed.open() as source:
            for line in source:
                parsed = parse_call(line)
                if parsed is None:
                    invalid_or_header += 1
                    continue
                values, split_reads, score, length = parsed
                current_count += 1
                try:
                    passes_coverage = coverage_pass(values)
                except ValueError:
                    invalid_or_header += 1
                    current_count -= 1
                    continue
                if not passes_coverage:
                    continue
                coverage_eligible_count += 1
                bin_name = length_bin(length)

                for split_min in SPLIT_THRESHOLDS:
                    for score_gt in SCORE_THRESHOLDS:
                        if split_reads >= split_min and score > score_gt:
                            threshold_counts[(split_min, score_gt)] += 1
                            threshold_lengths[(split_min, score_gt, bin_name)] += 1

                for name, (split_min, score_gt) in OUTPUT_CALLSETS.items():
                    if split_reads >= split_min and score > score_gt:
                        handles[name].write(line if line.endswith("\n") else line + "\n")
                        output_counts[name] += 1
                        output_lengths[(name, bin_name)] += 1
    finally:
        for handle in handles.values():
            handle.close()

    threshold_rows = []
    for split_min in SPLIT_THRESHOLDS:
        for score_gt in SCORE_THRESHOLDS:
            count = threshold_counts[(split_min, score_gt)]
            threshold_rows.append(
                {
                    "sample_id": sample,
                    "group": group,
                    "split_reads_min": split_min,
                    "score_strictly_greater_than": f"{score_gt:g}",
                    "coverage_filters_applied": 1,
                    "call_count": count,
                    "mapped_alignments": mapped,
                    "epm": count / mapped * 1_000_000,
                    **{
                        f"length_{name}_count": threshold_lengths[
                            (split_min, score_gt, name)
                        ]
                        for name, _, _ in LENGTH_BINS
                    },
                }
            )

    sample_rows = []
    for name in OUTPUT_CALLSETS:
        count = output_counts[name]
        sample_rows.append(
            {
                "sample_id": sample,
                "group": group,
                "callset": f"circlemap_{name}",
                "source_current_valid_rows": current_count,
                "coverage_eligible_rows": coverage_eligible_count,
                "call_count": count,
                "retention_fraction_of_current": (
                    count / current_count if current_count else float("nan")
                ),
                "mapped_alignments": mapped,
                "epm": count / mapped * 1_000_000,
                "invalid_or_header_rows": invalid_or_header,
                **{
                    f"length_{bin_name}_count": output_lengths[(name, bin_name)]
                    for bin_name, _, _ in LENGTH_BINS
                },
                "output_bed": str(output_paths[name]),
                "output_sha256": sha256(output_paths[name]),
            }
        )
    return threshold_rows, sample_rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = load_manifest(args.manifest)
    all_threshold_rows = []
    all_sample_rows = []
    for index, row in enumerate(manifest_rows, start=1):
        threshold_rows, sample_rows = process_sample(row, args.output_dir)
        all_threshold_rows.extend(threshold_rows)
        all_sample_rows.extend(sample_rows)
        print(
            f"[{index:02d}/78] {row['sample_id']} complete",
            file=sys.stderr,
            flush=True,
        )

    tables = args.output_dir / "tables"
    write_tsv(tables / "circlemap_threshold_grid_sample_metrics.tsv", all_threshold_rows)
    write_tsv(tables / "circlemap_output_callset_sample_metrics.tsv", all_sample_rows)

    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(args.manifest),
        "manifest_sha256": sha256(args.manifest),
        "sample_count": len(manifest_rows),
        "split_thresholds": SPLIT_THRESHOLDS,
        "score_thresholds_strictly_greater_than": SCORE_THRESHOLDS,
        "output_callsets": OUTPUT_CALLSETS,
        "coverage_filters": {
            "mean_coverage_greater_than_sd": True,
            "start_coverage_increase_min": 0.3,
            "end_coverage_increase_min": 0.3,
            "uncovered_fraction_strictly_less_than": 0.1,
        },
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.output_dir / "circlemap_sensitivity_run_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")

    expected_sample_rows = len(manifest_rows) * len(OUTPUT_CALLSETS)
    expected_threshold_rows = (
        len(manifest_rows) * len(SPLIT_THRESHOLDS) * len(SCORE_THRESHOLDS)
    )
    if len(all_sample_rows) != expected_sample_rows:
        raise RuntimeError("Unexpected output callset sample-row count")
    if len(all_threshold_rows) != expected_threshold_rows:
        raise RuntimeError("Unexpected threshold-grid sample-row count")
    if any(not Path(row["output_bed"]).is_file() for row in all_sample_rows):
        raise RuntimeError("One or more output BED files are missing")

    with (args.output_dir / "circlemap_sensitivity_validation.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write(f"samples={len(manifest_rows)}\n")
        handle.write(f"callset_sample_rows={len(all_sample_rows)}\n")
        handle.write(f"threshold_grid_sample_rows={len(all_threshold_rows)}\n")


if __name__ == "__main__":
    main()

