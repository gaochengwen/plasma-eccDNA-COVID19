#!/usr/bin/env python3
"""Audit four validation targets against all prespecified artifact masks."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

from build_masked_callsets import IntervalIndex, WINDOW_BP


TARGETS = {
    "CTNNA2": ("chr2", 79_885_063, 79_885_444),
    "CAB39": ("chr2", 230_727_472, 230_727_831),
    "chr22_target": ("chr22", 42_359_179, 42_359_558),
    "SDK1": ("chr7", 3_619_742, 3_620_122),
    "TCF7L1": ("chr2", 85_213_855, 85_214_235),
}


def window(position: int) -> tuple[int, int]:
    return max(0, position - WINDOW_BP), position + WINDOW_BP


def repeatmasker_annotations(path: Path) -> dict[tuple[str, str], set[str]]:
    output = {(target, side): set() for target in TARGETS for side in ("start", "end")}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            chrom = fields[5]
            try:
                start, end = int(fields[6]), int(fields[7])
            except ValueError:
                continue
            for target, (target_chrom, target_start, target_end) in TARGETS.items():
                if chrom != target_chrom:
                    continue
                for side, position in (("start", target_start), ("end", target_end)):
                    left, right = window(position)
                    if start < right and end > left:
                        output[(target, side)].add(
                            f"{fields[10]}|{fields[11]}|{fields[12]}:{start}-{end}"
                        )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--forbidden-bed", required=True, type=Path)
    parser.add_argument("--umap-bed", required=True, type=Path)
    parser.add_argument("--segdup-bed", required=True, type=Path)
    parser.add_argument("--repeat-bed", required=True, type=Path)
    parser.add_argument("--repeatmasker-raw", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    indexes = {
        "blacklist_gap": IntervalIndex(args.forbidden_bed),
        "umap_k100_unique": IntervalIndex(args.umap_bed),
        "segmental_duplication": IntervalIndex(args.segdup_bed),
        "repeatmasker": IntervalIndex(args.repeat_bed),
    }
    annotations = repeatmasker_annotations(args.repeatmasker_raw)
    rows = []
    for target, (chrom, start, end) in TARGETS.items():
        for side, position in (("start", start), ("end", end)):
            left, right = window(position)
            row = {
                "target": target,
                "chrom": chrom,
                "target_start": start,
                "target_end": end,
                "breakpoint_side": side,
                "breakpoint_position": position,
                "window_start": left,
                "window_end": right,
                "window_length_bp": right - left,
            }
            for resource, index in indexes.items():
                overlap = index.overlap_bp(chrom, left, right)
                row[f"{resource}_overlap_bp"] = overlap
                row[f"{resource}_overlap_fraction"] = overlap / (right - left)
            row["repeatmasker_annotations"] = ";".join(
                sorted(annotations[(target, side)])
            )
            rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
