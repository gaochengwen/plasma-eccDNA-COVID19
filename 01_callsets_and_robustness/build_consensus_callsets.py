#!/usr/bin/env python3
"""Match methods-consistent Circle-Map calls to Circle_finder breakpoints."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import platform
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


TOLERANCES = (2, 5, 10, 20)
MIN_CIRCLE_FINDER_SUPPORT = 2
MATCH_FIELDS = (
    "sample_id",
    "group",
    "tolerance_bp",
    "chrom",
    "circlemap_start",
    "circlemap_end",
    "circlefinder_start",
    "circlefinder_end",
    "start_difference_cf_minus_cm",
    "end_difference_cf_minus_cm",
    "total_absolute_breakpoint_distance",
    "circlemap_discordant_reads",
    "circlemap_split_reads",
    "circlemap_score",
    "circlefinder_junction_support",
)


@dataclass(frozen=True)
class Call:
    chrom: str
    start: int
    end: int
    fields: tuple[str, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "sample_id" not in rows[0] or "group" not in rows[0]:
        raise ValueError(f"Invalid manifest: {path}")
    return rows


def load_sample_filter(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError(f"Sample filter has no sample_id: {path}")
        return {row["sample_id"] for row in reader}


def load_circlemap(path: Path) -> list[Call]:
    calls = []
    with path.open() as handle:
        for line in handle:
            fields = tuple(line.rstrip("\n").split("\t"))
            if len(fields) < 11:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            calls.append(Call(fields[0], start, end, fields))
    calls.sort(key=lambda call: (call.chrom, call.start, call.end))
    return calls


def load_circlefinder(path: Path) -> tuple[list[Call], int]:
    calls = []
    below_support = 0
    with path.open() as handle:
        for line in handle:
            fields = tuple(line.rstrip("\n").split("\t"))
            if len(fields) < 4:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
                support = int(fields[3])
            except ValueError:
                continue
            if support < MIN_CIRCLE_FINDER_SUPPORT:
                below_support += 1
                continue
            calls.append(Call(fields[0], start, end, fields))
    calls.sort(key=lambda call: (call.chrom, call.start, call.end))
    return calls, below_support


def candidate_edges(
    circlemap: list[Call], circlefinder: list[Call], tolerance: int
) -> tuple[list[tuple], Counter[int], Counter[int]]:
    cm_by_chrom: dict[str, list[tuple[int, Call]]] = defaultdict(list)
    for index, call in enumerate(circlemap):
        cm_by_chrom[call.chrom].append((index, call))
    cm_starts = {
        chrom: [call.start for _index, call in calls]
        for chrom, calls in cm_by_chrom.items()
    }

    edges = []
    cm_degrees: Counter[int] = Counter()
    cf_degrees: Counter[int] = Counter()
    for cf_index, cf_call in enumerate(circlefinder):
        chromosome_calls = cm_by_chrom.get(cf_call.chrom, [])
        starts = cm_starts.get(cf_call.chrom, [])
        left = bisect.bisect_left(starts, cf_call.start - tolerance)
        right = bisect.bisect_right(starts, cf_call.start + tolerance)
        for cm_index, cm_call in chromosome_calls[left:right]:
            start_diff = cf_call.start - cm_call.start
            end_diff = cf_call.end - cm_call.end
            if abs(end_diff) > tolerance:
                continue
            distance = abs(start_diff) + abs(end_diff)
            edges.append(
                (
                    distance,
                    abs(start_diff),
                    abs(end_diff),
                    cm_call.start,
                    cm_call.end,
                    cf_call.start,
                    cf_call.end,
                    cm_index,
                    cf_index,
                    start_diff,
                    end_diff,
                )
            )
            cm_degrees[cm_index] += 1
            cf_degrees[cf_index] += 1
    return edges, cm_degrees, cf_degrees


def greedy_match(
    circlemap: list[Call], circlefinder: list[Call], tolerance: int
) -> tuple[list[tuple[int, int, int, int, int]], dict[str, int]]:
    edges, cm_degrees, cf_degrees = candidate_edges(
        circlemap, circlefinder, tolerance
    )
    edges.sort()
    used_cm = set()
    used_cf = set()
    matches = []
    for edge in edges:
        (
            distance,
            _abs_start,
            _abs_end,
            _cm_start,
            _cm_end,
            _cf_start,
            _cf_end,
            cm_index,
            cf_index,
            start_diff,
            end_diff,
        ) = edge
        if cm_index in used_cm or cf_index in used_cf:
            continue
        used_cm.add(cm_index)
        used_cf.add(cf_index)
        matches.append((cm_index, cf_index, distance, start_diff, end_diff))
    matches.sort(key=lambda item: item[0])
    diagnostics = {
        "candidate_edge_count": len(edges),
        "ambiguous_circlemap_call_count": sum(
            degree > 1 for degree in cm_degrees.values()
        ),
        "ambiguous_circlefinder_call_count": sum(
            degree > 1 for degree in cf_degrees.values()
        ),
    }
    return matches, diagnostics


def write_consensus_bed(path: Path, circlemap: list[Call], matches: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for cm_index, _cf_index, _distance, _start_diff, _end_diff in matches:
            handle.write("\t".join(circlemap[cm_index].fields) + "\n")


def process_sample(
    sample: str,
    group: str,
    root: Path,
    match_writer: csv.DictWriter,
) -> tuple[list[dict], int]:
    cm_path = root / "callsets" / "circlemap_methods" / f"{sample}_circle_site.bed"
    cf_path = (
        root
        / "results"
        / "circlefinder_raw"
        / sample
        / f"{sample}.circle_finder.deduplicated.bed"
    )
    if not cm_path.is_file():
        raise FileNotFoundError(cm_path)
    if not cf_path.is_file():
        raise FileNotFoundError(cf_path)

    circlemap = load_circlemap(cm_path)
    circlefinder, cf_below_support = load_circlefinder(cf_path)
    summary_rows = []
    match_row_count = 0
    for tolerance in TOLERANCES:
        matches, diagnostics = greedy_match(circlemap, circlefinder, tolerance)
        out_path = (
            root
            / "callsets"
            / f"consensus_t{tolerance}"
            / f"{sample}_circle_site.bed"
        )
        write_consensus_bed(out_path, circlemap, matches)
        match_count = len(matches)
        union_count = len(circlemap) + len(circlefinder) - match_count
        summary_rows.append(
            {
                "sample_id": sample,
                "group": group,
                "tolerance_bp": tolerance,
                "circlemap_methods_count": len(circlemap),
                "circlefinder_support_ge2_count": len(circlefinder),
                "circlefinder_support_lt2_excluded": cf_below_support,
                "consensus_count": match_count,
                "circlemap_retention_fraction": (
                    match_count / len(circlemap) if circlemap else float("nan")
                ),
                "circlefinder_retention_fraction": (
                    match_count / len(circlefinder)
                    if circlefinder
                    else float("nan")
                ),
                "jaccard_index": (
                    match_count / union_count if union_count else float("nan")
                ),
                **diagnostics,
                "consensus_bed": str(out_path),
                "consensus_bed_sha256": sha256(out_path),
            }
        )
        for cm_index, cf_index, distance, start_diff, end_diff in matches:
            cm = circlemap[cm_index]
            cf = circlefinder[cf_index]
            match_writer.writerow(
                {
                    "sample_id": sample,
                    "group": group,
                    "tolerance_bp": tolerance,
                    "chrom": cm.chrom,
                    "circlemap_start": cm.start,
                    "circlemap_end": cm.end,
                    "circlefinder_start": cf.start,
                    "circlefinder_end": cf.end,
                    "start_difference_cf_minus_cm": start_diff,
                    "end_difference_cf_minus_cm": end_diff,
                    "total_absolute_breakpoint_distance": distance,
                    "circlemap_discordant_reads": cm.fields[3],
                    "circlemap_split_reads": cm.fields[4],
                    "circlemap_score": cm.fields[5],
                    "circlefinder_junction_support": cf.fields[3],
                }
            )
            match_row_count += 1
    return summary_rows, match_row_count


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        handle_context = gzip.open(path, mode="wt", newline="")
    else:
        handle_context = path.open(mode="w", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-filter", type=Path)
    parser.add_argument("--label", default="all")
    args = parser.parse_args()

    selected = load_sample_filter(args.sample_filter)
    manifest = load_manifest(args.manifest)
    if selected is not None:
        manifest = [row for row in manifest if row["sample_id"] in selected]
        missing = selected - {row["sample_id"] for row in manifest}
        if missing:
            raise ValueError(f"Samples absent from manifest: {sorted(missing)}")

    tables = args.root / "tables"
    summary_path = tables / f"caller_concordance_{args.label}.tsv"
    matches_path = tables / f"caller_matches_{args.label}.tsv.gz"
    tables.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    match_row_count = 0
    with gzip.open(matches_path, mode="wt", newline="") as match_handle:
        match_writer = csv.DictWriter(
            match_handle,
            fieldnames=MATCH_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        match_writer.writeheader()
        for index, row in enumerate(manifest, start=1):
            sample_summary, sample_match_count = process_sample(
                row["sample_id"], row["group"], args.root, match_writer
            )
            summary_rows.extend(sample_summary)
            match_row_count += sample_match_count
            print(
                f"[{index}/{len(manifest)}] {row['sample_id']} matched",
                file=sys.stderr,
                flush=True,
            )

    write_tsv(summary_path, summary_rows)

    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(manifest),
        "sample_filter": str(args.sample_filter) if args.sample_filter else None,
        "tolerances_bp": TOLERANCES,
        "circlefinder_min_junction_support": MIN_CIRCLE_FINDER_SUPPORT,
        "matching_algorithm": (
            "deterministic greedy one-to-one matching ordered by total absolute "
            "breakpoint distance and then coordinate tie-breakers"
        ),
        "circlemap_input": "methods-consistent callset",
        "consensus_coordinates": "Circle-Map coordinates",
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.root / f"consensus_run_parameters_{args.label}.json").open("w") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")

    expected = len(manifest) * len(TOLERANCES)
    if len(summary_rows) != expected:
        raise RuntimeError(f"Expected {expected} summary rows, observed {len(summary_rows)}")
    with (args.root / f"consensus_validation_{args.label}.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write(f"samples={len(manifest)}\n")
        handle.write(f"summary_rows={len(summary_rows)}\n")
        handle.write(f"match_rows={match_row_count}\n")


if __name__ == "__main__":
    main()
