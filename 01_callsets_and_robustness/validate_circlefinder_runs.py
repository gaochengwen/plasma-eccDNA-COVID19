#!/usr/bin/env python3
"""Validate completed Circle_finder outputs and record auditable QC metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_fastq_recovery(path: Path) -> dict[str, int | str]:
    text = path.read_text(errors="replace") if path.is_file() else ""
    result: dict[str, int | str] = {"fastq_recovery_log_present": int(bool(text))}
    patterns = {
        "fastq_reads_processed": r"processed\s+(\d+)\s+reads",
        "fastq_singletons_discarded": r"discarded\s+(\d+)\s+singletons",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.IGNORECASE)
        result[key] = int(match.group(1)) if match else ""
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-filter", type=Path)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()

    manifest = read_tsv(args.manifest)
    if args.sample_filter:
        selected = {row["sample_id"] for row in read_tsv(args.sample_filter)}
        manifest = [row for row in manifest if row["sample_id"] in selected]
        missing = selected - {row["sample_id"] for row in manifest}
        if missing:
            raise ValueError(f"Samples absent from manifest: {sorted(missing)}")
    rows = []
    failures = []
    for metadata in manifest:
        sample = metadata["sample_id"]
        result_dir = args.root / "results" / "circlefinder_raw" / sample
        bed = result_dir / f"{sample}.circle_finder.deduplicated.bed"
        raw = result_dir / f"{sample}.circle_finder.raw.tsv.gz"
        summary = result_dir / f"{sample}.run_summary.tsv"
        complete = result_dir / ".complete"
        checksum_file = result_dir / f"{sample}.sha256"
        required = (bed, raw, summary, complete, checksum_file)
        missing = [str(path) for path in required if not path.exists()]
        status = "PASS"
        reason = ""
        raw_rows = -1
        bed_rows = -1
        checksum_pass = 0
        if missing:
            status = "FAIL"
            reason = "missing:" + ",".join(missing)
        else:
            try:
                with gzip.open(raw, "rt", errors="replace") as handle:
                    raw_rows = sum(1 for _ in handle)
                with bed.open() as handle:
                    bed_rows = sum(1 for line in handle if line.strip())
                expected = {}
                with checksum_file.open() as handle:
                    for line in handle:
                        fields = line.rstrip("\n").split()
                        if len(fields) >= 2:
                            expected[Path(fields[-1]).name] = fields[0]
                checksum_pass = int(
                    expected.get(bed.name) == sha256(bed)
                    and expected.get(raw.name) == sha256(raw)
                    and expected.get(summary.name) == sha256(summary)
                )
                if not checksum_pass:
                    status = "FAIL"
                    reason = "checksum_mismatch"
                elif raw_rows <= 0 or bed_rows <= 0:
                    status = "FAIL"
                    reason = f"nonpositive_calls:raw={raw_rows},deduplicated={bed_rows}"
            except (OSError, EOFError, ValueError) as error:
                status = "FAIL"
                reason = f"read_error:{error}"
        recovery = parse_fastq_recovery(result_dir / "fastq_recovery.stderr.log")
        rows.append(
            {
                "sample_id": sample,
                "group": metadata["group"],
                "status": status,
                "failure_reason": reason,
                "raw_rows": raw_rows,
                "deduplicated_calls": bed_rows,
                "checksum_pass": checksum_pass,
                "raw_gzip_bytes": raw.stat().st_size if raw.is_file() else -1,
                "deduplicated_bed_bytes": bed.stat().st_size if bed.is_file() else -1,
                **recovery,
            }
        )
        if status != "PASS":
            failures.append(sample)

    table = args.root / "tables" / f"circlefinder_validation_{args.label}.tsv"
    table.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "label": args.label,
        "sample_count": len(rows),
        "positive_raw_and_deduplicated_calls_required": True,
        "sha256_validation_required": True,
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.root / f"circlefinder_validation_parameters_{args.label}.json").open(
        "w"
    ) as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")
    validation = args.root / f"circlefinder_validation_{args.label}.txt"
    with validation.open("w") as handle:
        handle.write("PASS\n" if not failures else "FAIL\n")
        handle.write(f"samples={len(rows)}\n")
        handle.write(f"passed={len(rows) - len(failures)}\n")
        handle.write(f"failed={len(failures)}\n")
        if failures:
            handle.write("failed_samples=" + ",".join(failures) + "\n")
    if failures:
        raise RuntimeError(f"Circle_finder validation failed: {failures}")


if __name__ == "__main__":
    main()
