#!/usr/bin/env python3
"""Calculate eccDNA EPM inside a peak BED using BAM indexes.

This is a drop-in reimplementation of the existing shell definition:

    EPM = unique eccDNA intervals overlapping peaks / mapped alignments in peaks * 1e6

The numerator is computed with `bedtools intersect -u`.  The denominator uses
`samtools view -L peaks -c -F 3332`, which counts alignments overlapping the
same peak BED through the BAM index and avoids streaming the entire BAM through
bedtools.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path


EXCLUDE_FLAGS = "3332"


def run_stdout(args: list[str]) -> str:
    completed = subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return completed.stdout


def run(args: list[str]) -> None:
    subprocess.run(args, check=True)


def sort_bed(input_path: Path, output_path: Path) -> None:
    with output_path.open("w") as out_handle:
        subprocess.run(["sort", "-k1,1", "-k2,2n", str(input_path)], check=True, stdout=out_handle)


def filter_and_sort_eccbed(input_path: Path, output_path: Path) -> None:
    awk_program = 'BEGIN{OFS="\\t"} $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print}'
    with output_path.open("w") as out_handle:
        awk = subprocess.Popen(["awk", awk_program, str(input_path)], stdout=subprocess.PIPE)
        sort = subprocess.run(["sort", "-k1,1", "-k2,2n"], stdin=awk.stdout, stdout=out_handle, check=True)
        if awk.stdout is not None:
            awk.stdout.close()
        awk_return = awk.wait()
        if awk_return != 0:
            raise subprocess.CalledProcessError(awk_return, ["awk", awk_program, str(input_path)])


def count_ecc_overlaps(bedtools: str, ecc_sorted: Path, peaks_sorted: Path) -> int:
    process = subprocess.Popen(
        [bedtools, "intersect", "-a", str(ecc_sorted), "-b", str(peaks_sorted), "-u"],
        stdout=subprocess.PIPE,
        text=False,
    )
    count = 0
    assert process.stdout is not None
    for _ in process.stdout:
        count += 1
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, [bedtools, "intersect"])
    return count


def count_mapped_reads_in_peaks(samtools: str, bam: Path, peaks_sorted: Path) -> int:
    output = run_stdout([samtools, "view", "-c", "-F", EXCLUDE_FLAGS, "-L", str(peaks_sorted), str(bam)])
    return int(output.strip())


def sample_id_from_bam(path: Path) -> str:
    name = path.name
    if name.startswith("sorted_"):
        name = name[len("sorted_") :]
    if name.endswith("_circle.bam"):
        name = name[: -len("_circle.bam")]
    return name


def process_group(
    *,
    peaks: Path,
    ecc_dir: Path,
    bam_dir: Path,
    output: Path,
    samtools: str,
    bedtools: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="epm_peaks_") as tmp_name:
        tmp_dir = Path(tmp_name)
        peaks_sorted = tmp_dir / "peaks.sorted.bed"
        sort_bed(peaks, peaks_sorted)
        with output.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["ID", "bam", "eccbed", "mapped_reads_in_peaks", "eccDNA_counts_in_peaks", "EPM"])
            for bam in sorted(bam_dir.glob("sorted_*_circle.bam")):
                sample_id = sample_id_from_bam(bam)
                eccbed = ecc_dir / f"{sample_id}_circle_site.bed"
                if not eccbed.exists():
                    continue
                ecc_sorted = tmp_dir / f"{sample_id}.ecc.sorted.bed"
                filter_and_sort_eccbed(eccbed, ecc_sorted)
                ecc_count = count_ecc_overlaps(bedtools, ecc_sorted, peaks_sorted)
                mapped_reads = count_mapped_reads_in_peaks(samtools, bam, peaks_sorted)
                epm = ecc_count * 1_000_000 / mapped_reads if mapped_reads else float("nan")
                writer.writerow([sample_id, str(bam), str(eccbed), mapped_reads, ecc_count, f"{epm:.6f}"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True, type=Path)
    parser.add_argument("--ecc-dir", required=True, type=Path)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samtools", default=shutil.which("samtools") or "/opt/samtools-1.9/samtools")
    parser.add_argument("--bedtools", default=shutil.which("bedtools") or "/opt/bedtools2-2.28.0/bin/bedtools")
    args = parser.parse_args()

    process_group(
        peaks=args.peaks,
        ecc_dir=args.ecc_dir,
        bam_dir=args.bam_dir,
        output=args.output,
        samtools=args.samtools,
        bedtools=args.bedtools,
    )


if __name__ == "__main__":
    main()
