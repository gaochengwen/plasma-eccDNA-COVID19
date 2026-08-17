#!/usr/bin/env python3
"""Create group-specific symlinks for an existing 78-sample callset."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--callset-dir", required=True, type=Path)
    parser.add_argument("--covid-dir", required=True, type=Path)
    parser.add_argument("--hc-dir", required=True, type=Path)
    args = parser.parse_args()

    args.covid_dir.mkdir(parents=True, exist_ok=True)
    args.hc_dir.mkdir(parents=True, exist_ok=True)
    counts = {"COVID-19": 0, "HC": 0}
    seen: set[str] = set()
    with args.metadata.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = row["sample_id"].strip()
            group = row["group"].strip()
            if group not in counts:
                raise ValueError(f"Unexpected group {group!r} for {sample}")
            if sample in seen:
                raise ValueError(f"Duplicate sample {sample}")
            seen.add(sample)
            source = args.callset_dir / f"{sample}_circle_site.bed"
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(source)
            destination_dir = args.covid_dir if group == "COVID-19" else args.hc_dir
            destination = destination_dir / source.name
            if destination.is_symlink():
                if destination.resolve() != source.resolve():
                    destination.unlink()
            elif destination.exists():
                raise FileExistsError(destination)
            if not destination.exists():
                os.symlink(source, destination)
            counts[group] += 1
    if counts != {"COVID-19": 39, "HC": 39}:
        raise ValueError(f"Expected 39+39 samples; observed {counts}")
    print(f"Prepared artifact callset links: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
