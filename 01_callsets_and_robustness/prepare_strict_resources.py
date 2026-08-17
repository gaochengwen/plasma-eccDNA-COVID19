#!/usr/bin/env python3
"""Download, normalize, merge, and checksum GRCh38 robustness masks."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import sys
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
SOURCES = {
    "umap_k100": {
        "url": "https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz",
        "filename": "k100.umap.bed.gz",
        "coordinate_columns": (0, 1, 2),
        "description": "Umap single-read unique mappability intervals for 100-mers",
    },
    "segmental_duplication": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
        "filename": "genomicSuperDups.txt.gz",
        "coordinate_columns": (1, 2, 3),
        "description": "UCSC GRCh38 genomicSuperDups",
    },
    "repeatmasker": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
        "filename": "rmsk.txt.gz",
        "coordinate_columns": (5, 6, 7),
        "description": "UCSC GRCh38 RepeatMasker table",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    if path.is_file() and path.stat().st_size > 0:
        return
    temporary = path.with_suffix(path.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "eccDNA-robustness/1.0"})
    with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
    if temporary.stat().st_size == 0:
        raise RuntimeError(f"Empty download: {url}")
    temporary.replace(path)


def normalize_and_merge(
    source: Path, output: Path, coordinate_columns: tuple[int, int, int]
) -> tuple[int, int, int]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    raw_valid = 0
    skipped = 0
    with gzip.open(source, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                chrom = fields[coordinate_columns[0]]
                start = int(fields[coordinate_columns[1]])
                end = int(fields[coordinate_columns[2]])
            except (IndexError, ValueError):
                skipped += 1
                continue
            if chrom not in CANONICAL or start < 0 or end <= start:
                skipped += 1
                continue
            intervals[chrom].append((start, end))
            raw_valid += 1

    merged_count = 0
    with output.open("w") as out:
        for chrom in sorted(
            intervals,
            key=lambda value: (
                23 if value == "chrX" else 24 if value == "chrY" else 25 if value == "chrM" else int(value[3:])
            ),
        ):
            values = sorted(intervals[chrom])
            merged_start, merged_end = values[0]
            for start, end in values[1:]:
                if start <= merged_end:
                    merged_end = max(merged_end, end)
                else:
                    out.write(f"{chrom}\t{merged_start}\t{merged_end}\n")
                    merged_count += 1
                    merged_start, merged_end = start, end
            out.write(f"{chrom}\t{merged_start}\t{merged_end}\n")
            merged_count += 1
    return raw_valid, merged_count, skipped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    raw_dir = args.output_dir / "raw"
    normalized_dir = args.output_dir / "normalized"
    raw_dir.mkdir(parents=True, exist_ok=True)
    normalized_dir.mkdir(parents=True, exist_ok=True)

    for key in (
        "http_proxy",
        "https_proxy",
        "all_proxy",
        "HTTP_PROXY",
        "HTTPS_PROXY",
        "ALL_PROXY",
    ):
        os.environ.pop(key, None)

    rows = []
    for resource_id, metadata in SOURCES.items():
        raw_path = raw_dir / str(metadata["filename"])
        output_path = normalized_dir / f"{resource_id}.canonical.merged.bed"
        download(str(metadata["url"]), raw_path)
        raw_count, merged_count, skipped = normalize_and_merge(
            raw_path, output_path, metadata["coordinate_columns"]  # type: ignore[arg-type]
        )
        rows.append(
            {
                "resource_id": resource_id,
                "description": metadata["description"],
                "source_url": metadata["url"],
                "downloaded_file": str(raw_path),
                "downloaded_bytes": raw_path.stat().st_size,
                "downloaded_sha256": sha256(raw_path),
                "normalized_file": str(output_path),
                "normalized_sha256": sha256(output_path),
                "raw_valid_canonical_intervals": raw_count,
                "merged_canonical_intervals": merged_count,
                "skipped_rows": skipped,
            }
        )

    manifest = args.output_dir / "strict_resource_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_system": "GRCh38/hg38, BED 0-based half-open",
        "canonical_chromosomes": sorted(CANONICAL),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.output_dir / "strict_resource_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")

    if len(rows) != 3 or any(int(row["merged_canonical_intervals"]) == 0 for row in rows):
        raise RuntimeError("Strict resource validation failed")
    with (args.output_dir / "strict_resource_validation.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write(f"resources={len(rows)}\n")


if __name__ == "__main__":
    main()
