#!/usr/bin/env python3
"""Apply prespecified breakpoint-level artifact masks to Circle-Map calls."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import platform
import sys
from array import array
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
WINDOW_BP = 100
MIN_UMAP_FRACTION = 0.90
HIGH_SUPPORT_SPLIT_MIN = 10
HIGH_SUPPORT_SCORE_GT = 1000.0
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


class IntervalIndex:
    def __init__(self, path: Path) -> None:
        starts: dict[str, array] = {}
        ends: dict[str, array] = {}
        previous: dict[str, tuple[int, int]] = {}
        with path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                start, end = int(fields[1]), int(fields[2])
                if chrom not in CANONICAL or end <= start:
                    continue
                if chrom in previous:
                    prev_start, prev_end = previous[chrom]
                    if start < prev_start or start < prev_end:
                        raise ValueError(f"Intervals not sorted and disjoint: {path}")
                starts.setdefault(chrom, array("q")).append(start)
                ends.setdefault(chrom, array("q")).append(end)
                previous[chrom] = (start, end)
        self.starts = starts
        self.ends = ends
        self.path = path

    def overlap_bp(self, chrom: str, start: int, end: int) -> int:
        if end <= start or chrom not in self.starts:
            return 0
        starts = self.starts[chrom]
        ends = self.ends[chrom]
        index = max(0, bisect.bisect_right(ends, start))
        overlap = 0
        while index < len(starts) and starts[index] < end:
            overlap += max(0, min(end, ends[index]) - max(start, starts[index]))
            index += 1
        return overlap

    def any_overlap(self, chrom: str, start: int, end: int) -> bool:
        if end <= start or chrom not in self.starts:
            return False
        ends = self.ends[chrom]
        index = bisect.bisect_right(ends, start)
        return index < len(ends) and self.starts[chrom][index] < end


def length_bin(length: int) -> str:
    for name, left, right in LENGTH_BINS:
        if left <= length < right:
            return name
    raise AssertionError(length)


def window(center: int) -> tuple[int, int]:
    return max(0, center - WINDOW_BP), center + WINDOW_BP


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_id", "group", "mapped_alignments_idxstats"}
    if len(rows) != 78 or not rows or not required.issubset(rows[0]):
        raise ValueError(f"Expected 78-sample manifest with required fields: {path}")
    return rows


def process_sample(
    row: dict[str, str],
    root: Path,
    forbidden: IntervalIndex,
    umap: IntervalIndex,
    segdup: IntervalIndex,
    repeats: IntervalIndex,
) -> tuple[list[dict], list[dict]]:
    sample = row["sample_id"]
    source = root / "callsets" / "circlemap_methods" / f"{sample}_circle_site.bed"
    if not source.is_file():
        raise FileNotFoundError(source)
    output_paths = {
        "circlemap_artifact_masked": root / "callsets" / "circlemap_artifact_masked" / source.name,
        "circlemap_high_support_masked": root
        / "callsets"
        / "circlemap_high_support_masked"
        / source.name,
    }
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    handles = {name: path.open("w") for name, path in output_paths.items()}

    stages = Counter()
    output_counts = Counter()
    output_lengths = Counter()
    try:
        with source.open() as source_handle:
            for line in source_handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    stages["invalid"] += 1
                    continue
                try:
                    chrom = fields[0]
                    start, end = int(fields[1]), int(fields[2])
                    split_reads = int(float(fields[4]))
                    score = float(fields[5])
                except ValueError:
                    stages["invalid"] += 1
                    continue
                stages["source_methods"] += 1
                if chrom not in CANONICAL or start < 0 or end <= start:
                    stages["failed_canonical_or_coordinate"] += 1
                    continue
                stages["canonical_coordinate"] += 1
                start_window = window(start)
                end_window = window(end)
                if (
                    forbidden.any_overlap(chrom, start, end)
                    or forbidden.any_overlap(chrom, *start_window)
                    or forbidden.any_overlap(chrom, *end_window)
                ):
                    stages["failed_blacklist_gap"] += 1
                    continue
                stages["blacklist_gap_pass"] += 1
                start_umap = umap.overlap_bp(chrom, *start_window) / (
                    start_window[1] - start_window[0]
                )
                end_umap = umap.overlap_bp(chrom, *end_window) / (
                    end_window[1] - end_window[0]
                )
                if min(start_umap, end_umap) < MIN_UMAP_FRACTION:
                    stages["failed_umap"] += 1
                    continue
                stages["umap_pass"] += 1
                if segdup.any_overlap(chrom, *start_window) or segdup.any_overlap(
                    chrom, *end_window
                ):
                    stages["failed_segmental_duplication"] += 1
                    continue
                stages["segmental_duplication_pass"] += 1
                if repeats.any_overlap(chrom, *start_window) or repeats.any_overlap(
                    chrom, *end_window
                ):
                    stages["failed_repeatmasker"] += 1
                    continue
                stages["repeatmasker_pass"] += 1

                handles["circlemap_artifact_masked"].write(line)
                output_counts["circlemap_artifact_masked"] += 1
                output_lengths[
                    ("circlemap_artifact_masked", length_bin(end - start))
                ] += 1
                if split_reads >= HIGH_SUPPORT_SPLIT_MIN and score > HIGH_SUPPORT_SCORE_GT:
                    handles["circlemap_high_support_masked"].write(line)
                    output_counts["circlemap_high_support_masked"] += 1
                    output_lengths[
                        ("circlemap_high_support_masked", length_bin(end - start))
                    ] += 1
    finally:
        for handle in handles.values():
            handle.close()

    stage_rows = [
        {
            "sample_id": sample,
            "group": row["group"],
            "stage_or_failure": stage,
            "call_count": count,
        }
        for stage, count in sorted(stages.items())
    ]
    sample_rows = []
    mapped = int(row["mapped_alignments_idxstats"])
    for callset, output_path in output_paths.items():
        count = output_counts[callset]
        sample_rows.append(
            {
                "sample_id": sample,
                "group": row["group"],
                "callset": callset,
                "call_count": count,
                "retention_fraction_of_methods": (
                    count / stages["source_methods"] if stages["source_methods"] else math.nan
                ),
                "mapped_alignments": mapped,
                "epm": count / mapped * 1_000_000,
                **{
                    f"length_{name}_count": output_lengths[(callset, name)]
                    for name, _, _ in LENGTH_BINS
                },
                "output_bed": str(output_path),
                "output_sha256": sha256(output_path),
            }
        )
    return stage_rows, sample_rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--forbidden-bed", required=True, type=Path)
    parser.add_argument("--umap-bed", required=True, type=Path)
    parser.add_argument("--segdup-bed", required=True, type=Path)
    parser.add_argument("--repeat-bed", required=True, type=Path)
    args = parser.parse_args()

    indexes = {
        "forbidden": IntervalIndex(args.forbidden_bed),
        "umap": IntervalIndex(args.umap_bed),
        "segdup": IntervalIndex(args.segdup_bed),
        "repeats": IntervalIndex(args.repeat_bed),
    }
    stage_rows = []
    sample_rows = []
    manifest = load_manifest(args.manifest)
    for index, row in enumerate(manifest, start=1):
        sample_stages, sample_metrics = process_sample(
            row,
            args.root,
            indexes["forbidden"],
            indexes["umap"],
            indexes["segdup"],
            indexes["repeats"],
        )
        stage_rows.extend(sample_stages)
        sample_rows.extend(sample_metrics)
        print(f"[{index:02d}/78] {row['sample_id']} masked", file=sys.stderr, flush=True)

    write_tsv(args.root / "tables" / "artifact_filter_stage_counts.tsv", stage_rows)
    write_tsv(args.root / "tables" / "artifact_masked_sample_metrics.tsv", sample_rows)
    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_system": "GRCh38/hg38, BED 0-based half-open",
        "breakpoint_window_bp_each_side": WINDOW_BP,
        "minimum_umap_k100_unique_fraction_each_breakpoint_window": MIN_UMAP_FRACTION,
        "blacklist_gap_rule": "no full-circle or breakpoint-window overlap",
        "segmental_duplication_rule": "no breakpoint-window overlap",
        "repeatmasker_rule": "no breakpoint-window overlap; internal repeats retained",
        "high_support_rule": {
            "split_reads_min": HIGH_SUPPORT_SPLIT_MIN,
            "score_strictly_greater_than": HIGH_SUPPORT_SCORE_GT,
        },
        "resources": {
            key: {"path": str(value.path), "sha256": sha256(value.path)}
            for key, value in indexes.items()
        },
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.root / "artifact_filter_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")
    if len(sample_rows) != 156:
        raise RuntimeError(f"Expected 156 sample metric rows, observed {len(sample_rows)}")
    with (args.root / "artifact_filter_validation.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write("samples=78\n")
        handle.write(f"sample_metric_rows={len(sample_rows)}\n")


if __name__ == "__main__":
    main()
