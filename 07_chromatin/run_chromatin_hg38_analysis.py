#!/usr/bin/env python3
"""Redo chromatin relative-enrichment analysis in a single hg38 coordinate system.

The workflow keeps eccDNA BED files in hg38, lifts legacy hg19 ChIP-seq peak BEDs
to hg38, excludes hg38 blacklist/gap regions, and computes both absolute burden
and relative enrichment after controlling for each sample's total eligible
eccDNA count.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


CANON_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CANON_SET = set(CANON_CHROMS)
MODES = ("full_interval", "junction_start_end", "midpoint")
SEED = 20260723

REFERENCE_URLS = {
    "hg19ToHg38.over.chain.gz": [
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    ],
    "hg38.chrom.sizes": [
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
        "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    ],
    "hg38-blacklist.v2.bed.gz": [
        "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
    ],
    "hg38.gap.txt.gz": [
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz",
        "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz",
    ],
}


PEAK_SOURCES = [
    {
        "peak_set_id": "CD14_H3K27ac",
        "dataset": "CD14_reference",
        "mark": "H3K27ac",
        "source_build": "hg19",
        "source_path": "enhancer/20260106/H3k27ac.th3.bed",
        "source_note": "Existing manuscript ChIP-seq peak BED.",
        "analysis_tier": "primary",
    },
    {
        "peak_set_id": "CD14_H3K27me3",
        "dataset": "CD14_reference",
        "mark": "H3K27me3",
        "source_build": "hg19",
        "source_path": "enhancer/20260106/H3k27me3.th3.bed",
        "source_note": "Existing manuscript ChIP-seq peak BED.",
        "analysis_tier": "primary",
    },
    {
        "peak_set_id": "CD14_H3K4me1",
        "dataset": "CD14_reference",
        "mark": "H3K4me1",
        "source_build": "hg19",
        "source_path": "enhancer/20260106/H3k4me1.th3.bed",
        "source_note": "Existing manuscript ChIP-seq peak BED.",
        "analysis_tier": "primary",
    },
    {
        "peak_set_id": "CD14_H3K4me3",
        "dataset": "CD14_reference",
        "mark": "H3K4me3",
        "source_build": "hg19",
        "source_path": "enhancer/20260106/H3k4me3.th3.bed",
        "source_note": "Existing manuscript ChIP-seq peak BED.",
        "analysis_tier": "primary",
    },
    {
        "peak_set_id": "CD14_H3K9ac",
        "dataset": "CD14_reference",
        "mark": "H3K9ac",
        "source_build": "hg19",
        "source_path": "enhancer/20260106/H3k9ac.th3.bed",
        "source_note": "Existing manuscript ChIP-seq peak BED.",
        "analysis_tier": "primary",
    },
    {
        "peak_set_id": "CD14_H3K9me3",
        "dataset": "CD14_reference",
        "mark": "H3K9me3",
        "source_build": "hg19",
        "source_path": "enhancer/20260106/H3k9me3.th3.bed",
        "source_note": "Existing manuscript ChIP-seq peak BED.",
        "analysis_tier": "primary",
    },
    {
        "peak_set_id": "A549_GSE179184_H3K27ac_official",
        "dataset": "A549_GSE179184_official",
        "mark": "H3K27ac",
        "source_build": "hg19",
        "source_path": "enhancer/GSE179184_A549-ACE2-H3K27ac-ChIP-Seq.merged.peak.bed",
        "source_note": "GSE179184 processed merged peak BED present in project.",
        "analysis_tier": "secondary",
    },
    {
        "peak_set_id": "A549_GSE179184_H3K4me3_official",
        "dataset": "A549_GSE179184_official",
        "mark": "H3K4me3",
        "source_build": "hg19",
        "source_path": "enhancer/GSE179184_A549-ACE2-H3K4me3-ChIP-Seq.merged.peak.bed",
        "source_note": "GSE179184 processed merged peak BED present in project.",
        "analysis_tier": "secondary",
    },
    {
        "peak_set_id": "A549_legacy_H3K27ac_thresholded",
        "dataset": "A549_legacy_thresholded",
        "mark": "H3K27ac",
        "source_build": "hg19",
        "source_path": "enhancer/20251022/H3K27ac.th3.bed",
        "source_note": "Legacy thresholded BED used for compatibility with previous figure analysis; not treated as the primary public processed peak set.",
        "analysis_tier": "legacy_sensitivity",
    },
    {
        "peak_set_id": "A549_legacy_H3K27me3_thresholded",
        "dataset": "A549_legacy_thresholded",
        "mark": "H3K27me3",
        "source_build": "hg19",
        "source_path": "enhancer/20251022/H3K27me3.th3.bed",
        "source_note": "Legacy thresholded BED used for compatibility with previous figure analysis; not treated as the primary public processed peak set.",
        "analysis_tier": "legacy_sensitivity",
    },
]


@dataclass
class Dirs:
    outdir: Path
    reference: Path
    inputs: Path
    liftover: Path
    peaks: Path
    tables: Path
    figures: Path
    logs: Path


@dataclass
class Track:
    peak_set_id: str
    dataset: str
    mark: str
    analysis_tier: str
    intervals: Dict[str, Tuple[np.ndarray, np.ndarray]]
    free_components: Dict[str, List[Tuple[int, int]]]
    prefix: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]
    merged_interval_count: int
    coverage_bp: int


def log(message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {message}", flush=True)


def ensure_dirs(outdir: Path) -> Dirs:
    dirs = Dirs(
        outdir=outdir,
        reference=outdir / "reference",
        inputs=outdir / "inputs",
        liftover=outdir / "liftover",
        peaks=outdir / "peaks_hg38",
        tables=outdir / "tables",
        figures=outdir / "figures",
        logs=outdir / "logs",
    )
    for d in dirs.__dict__.values():
        Path(d).mkdir(parents=True, exist_ok=True)
    return dirs


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def clear_proxy_env() -> None:
    for key in ("http_proxy", "https_proxy", "all_proxy", "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY"):
        os.environ.pop(key, None)


def download_one(name: str, urls: Sequence[str], outpath: Path, force: bool = False) -> str:
    if outpath.exists() and outpath.stat().st_size > 0 and not force:
        return "existing"
    clear_proxy_env()
    last_error = None
    tmp = outpath.with_suffix(outpath.suffix + ".tmp")
    for url in urls:
        try:
            log(f"Downloading {name} from {url}")
            with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as out:
                shutil.copyfileobj(response, out)
            if tmp.stat().st_size == 0:
                raise RuntimeError("downloaded file is empty")
            tmp.replace(outpath)
            (outpath.with_suffix(outpath.suffix + ".url")).write_text(url + "\n", encoding="utf-8")
            return url
        except Exception as exc:  # noqa: BLE001 - retry all URL failures.
            last_error = exc
            if tmp.exists():
                tmp.unlink()
    raise RuntimeError(f"Unable to download {name}: {last_error}")


def prepare_reference_files(dirs: Dirs, force_download: bool = False) -> pd.DataFrame:
    rows = []
    for name, urls in REFERENCE_URLS.items():
        outpath = dirs.reference / name
        source = download_one(name, urls, outpath, force=force_download)
        rows.append(
            {
                "filename": name,
                "source_url": source if source != "existing" else read_url_sidecar(outpath),
                "bytes": outpath.stat().st_size,
                "sha256": sha256_file(outpath),
                "server_file_path": str(outpath),
                "status": "present",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(dirs.reference / "reference_manifest.tsv", sep="\t", index=False)
    return manifest


def read_url_sidecar(path: Path) -> str:
    sidecar = path.with_suffix(path.suffix + ".url")
    if sidecar.exists():
        return sidecar.read_text(encoding="utf-8").strip()
    return "existing_without_url_sidecar"


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    sizes = {}
    with open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.rstrip("\n").split("\t")[:2]
            if chrom in CANON_SET:
                sizes[chrom] = int(size)
    missing = [chrom for chrom in CANON_CHROMS if chrom not in sizes]
    if missing:
        raise RuntimeError(f"Missing canonical chromosome sizes: {missing}")
    return sizes


def read_bed3(path: Path, chrom_sizes: Optional[Mapping[str, int]] = None) -> pd.DataFrame:
    rows = []
    with open_text_maybe_gzip(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(float(parts[1]))
                end = int(float(parts[2]))
            except ValueError:
                continue
            rows.append((chrom, start, end))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    if df.empty:
        return df
    df["length"] = df["end"] - df["start"]
    df["valid_interval"] = df["length"] > 0
    df["canonical"] = df["chrom"].isin(CANON_SET)
    if chrom_sizes is not None:
        df["within_chrom_bounds"] = [
            bool(row.canonical and row.start >= 0 and row.end <= chrom_sizes.get(row.chrom, -1))
            for row in df.itertuples(index=False)
        ]
    else:
        df["within_chrom_bounds"] = True
    return df


def write_bed(path: Path, intervals: Mapping[str, Sequence[Tuple[int, int]]], extra_name: Optional[str] = None) -> None:
    with open(path, "wt", encoding="utf-8") as out:
        for chrom in CANON_CHROMS:
            for start, end in intervals.get(chrom, []):
                if extra_name is None:
                    out.write(f"{chrom}\t{start}\t{end}\n")
                else:
                    out.write(f"{chrom}\t{start}\t{end}\t{extra_name}\n")


def merge_intervals_from_df(df: pd.DataFrame) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {chrom: [] for chrom in CANON_CHROMS}
    if df.empty:
        return intervals
    work = df[df["chrom"].isin(CANON_SET) & (df["end"] > df["start"])][["chrom", "start", "end"]].copy()
    if work.empty:
        return intervals
    work = work.sort_values(["chrom", "start", "end"])
    for chrom, sub in work.groupby("chrom", sort=False):
        merged: List[Tuple[int, int]] = []
        cur_s = None
        cur_e = None
        for start, end in sub[["start", "end"]].itertuples(index=False, name=None):
            start = int(start)
            end = int(end)
            if cur_s is None:
                cur_s, cur_e = start, end
            elif start <= cur_e:
                cur_e = max(cur_e, end)
            else:
                merged.append((cur_s, cur_e))
                cur_s, cur_e = start, end
        if cur_s is not None:
            merged.append((cur_s, cur_e))
        intervals[str(chrom)] = merged
    return intervals


def interval_dict_to_arrays(intervals: Mapping[str, Sequence[Tuple[int, int]]]) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    arrays = {}
    for chrom in CANON_CHROMS:
        vals = intervals.get(chrom, [])
        if vals:
            arr = np.asarray(vals, dtype=np.int64)
            arrays[chrom] = (arr[:, 0], arr[:, 1])
        else:
            arrays[chrom] = (np.asarray([], dtype=np.int64), np.asarray([], dtype=np.int64))
    return arrays


def prefix_tracks(intervals: Mapping[str, Tuple[np.ndarray, np.ndarray]]) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    prefix = {}
    for chrom, (starts, ends) in intervals.items():
        lengths = ends - starts
        prefix[chrom] = (starts, ends, np.concatenate([[0], np.cumsum(lengths, dtype=np.int64)]))
    return prefix


def coverage_to(prefix: Tuple[np.ndarray, np.ndarray, np.ndarray], x: int) -> int:
    starts, ends, csum = prefix
    if len(starts) == 0:
        return 0
    idx = int(np.searchsorted(ends, x, side="right"))
    total = int(csum[idx])
    if idx < len(starts) and starts[idx] < x:
        total += max(0, int(x - starts[idx]))
    return total


def coverage_in_range(prefix: Tuple[np.ndarray, np.ndarray, np.ndarray], start: int, end: int) -> int:
    if end <= start:
        return 0
    return coverage_to(prefix, end) - coverage_to(prefix, start)


def subtract_intervals(
    bases: Sequence[Tuple[int, int]],
    blockers: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    if not bases:
        return []
    if not blockers:
        return list(bases)
    out: List[Tuple[int, int]] = []
    j = 0
    nb = len(blockers)
    for base_s, base_e in bases:
        cursor = base_s
        while j < nb and blockers[j][1] <= base_s:
            j += 1
        k = j
        while k < nb and blockers[k][0] < base_e:
            blk_s, blk_e = blockers[k]
            if blk_s > cursor:
                out.append((cursor, min(blk_s, base_e)))
            cursor = max(cursor, blk_e)
            if cursor >= base_e:
                break
            k += 1
        if cursor < base_e:
            out.append((cursor, base_e))
    return [(s, e) for s, e in out if e > s]


def build_forbidden_and_allowed(dirs: Dirs, chrom_sizes: Mapping[str, int]) -> Dict[str, List[Tuple[int, int]]]:
    blacklist = read_bed3(dirs.reference / "hg38-blacklist.v2.bed.gz", chrom_sizes)
    gap_rows = []
    gap_path = dirs.reference / "hg38.gap.txt.gz"
    with gzip.open(gap_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[1]
            if chrom not in CANON_SET:
                continue
            gap_rows.append((chrom, int(parts[2]), int(parts[3])))
    gaps = pd.DataFrame(gap_rows, columns=["chrom", "start", "end"])
    for df in (blacklist, gaps):
        if "length" not in df.columns and not df.empty:
            df["length"] = df["end"] - df["start"]
    merged_forbidden = merge_intervals_from_df(pd.concat([blacklist[["chrom", "start", "end"]], gaps], ignore_index=True))
    allowed = {}
    whole = {chrom: [(0, int(chrom_sizes[chrom]))] for chrom in CANON_CHROMS}
    for chrom in CANON_CHROMS:
        allowed[chrom] = subtract_intervals(whole[chrom], merged_forbidden.get(chrom, []))
    write_bed(dirs.reference / "hg38_forbidden.blacklist_gaps.canonical.merged.bed", merged_forbidden)
    write_bed(dirs.reference / "hg38_allowed_components.canonical.bed", allowed)
    summary = {
        "blacklist_intervals_raw": int(len(blacklist)),
        "gap_intervals_raw": int(len(gaps)),
        "forbidden_merged_intervals_canonical": int(sum(len(v) for v in merged_forbidden.values())),
        "forbidden_bp_canonical": int(sum(e - s for vals in merged_forbidden.values() for s, e in vals)),
        "allowed_components_canonical": int(sum(len(v) for v in allowed.values())),
        "allowed_bp_canonical": int(sum(e - s for vals in allowed.values() for s, e in vals)),
    }
    pd.DataFrame([summary]).to_csv(dirs.reference / "hg38_blacklist_gap_summary.tsv", sep="\t", index=False)
    return allowed


def sanitize_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def make_peak_manifest(project_dir: Path, dirs: Dirs) -> pd.DataFrame:
    rows = []
    for source in PEAK_SOURCES:
        source_path = project_dir / source["source_path"]
        rows.append({**source, "source_abs_path": str(source_path), "source_exists": source_path.exists()})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(dirs.inputs / "chipseq_peak_sources.manifest.tsv", sep="\t", index=False)
    missing = manifest[~manifest["source_exists"]]
    if not missing.empty:
        missing_ids = ", ".join(missing["peak_set_id"].tolist())
        raise FileNotFoundError(f"Missing peak source files: {missing_ids}")
    return manifest


def make_sample_manifest(project_dir: Path, dirs: Dirs, samtools: str) -> pd.DataFrame:
    rows = []
    for group, bed_dir, bam_dir in (
        ("COVID", project_dir / "covid", project_dir / "covid_bam"),
        ("HC", project_dir / "normal", project_dir / "normal_bam"),
    ):
        for bed_path in sorted(bed_dir.glob("*_circle_site.bed")):
            sample = bed_path.name.replace("_circle_site.bed", "")
            bam_path = bam_dir / f"sorted_{sample}_circle.bam"
            mapped = mapped_read_count(bam_path, samtools)
            rows.append(
                {
                    "sample_id": sample,
                    "group": group,
                    "eccdna_bed": str(bed_path),
                    "bam": str(bam_path),
                    "bam_exists": bam_path.exists(),
                    "mapped_alignments_idxstats": mapped,
                }
            )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(dirs.inputs / "eccdna_sample_manifest.tsv", sep="\t", index=False)
    if manifest.empty:
        raise RuntimeError("No eccDNA sample BED files were found.")
    if not manifest["bam_exists"].all():
        raise FileNotFoundError("One or more sample BAM files are missing. See eccdna_sample_manifest.tsv.")
    return manifest


def mapped_read_count(bam_path: Path, samtools: str) -> int:
    if not bam_path.exists():
        return -1
    try:
        proc = subprocess.run([samtools, "idxstats", str(bam_path)], check=True, text=True, capture_output=True)
        total = 0
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                total += int(parts[2])
        return int(total)
    except Exception:
        proc = subprocess.run([samtools, "view", "-c", "-F", "0x4", str(bam_path)], check=True, text=True, capture_output=True)
        return int(proc.stdout.strip())


def standardize_peak_bed(source_path: Path, outpath: Path) -> Tuple[int, int]:
    original = 0
    valid = 0
    with open_text_maybe_gzip(source_path) as inp, open(outpath, "wt", encoding="utf-8") as out:
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            original += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            try:
                start = int(float(parts[1]))
                end = int(float(parts[2]))
            except ValueError:
                continue
            if end <= start:
                continue
            peak_id = f"peak_{original:09d}"
            out.write(f"{parts[0]}\t{start}\t{end}\t{peak_id}\t{end - start}\n")
            valid += 1
    return original, valid


def run_liftover_for_peaks(
    peak_manifest: pd.DataFrame,
    dirs: Dirs,
    chain: Path,
    liftover: str,
    chrom_sizes: Mapping[str, int],
) -> Tuple[pd.DataFrame, Dict[str, Track]]:
    rows = []
    tracks: Dict[str, Track] = {}
    for source in peak_manifest.to_dict(orient="records"):
        peak_id = source["peak_set_id"]
        stem = sanitize_id(peak_id)
        source_path = Path(source["source_abs_path"])
        std_bed = dirs.liftover / f"{stem}.source.hg19.with_id.bed"
        lifted = dirs.liftover / f"{stem}.lifted.hg38.raw.bed"
        unmapped = dirs.liftover / f"{stem}.lifted.hg38.unmapped.txt"
        lifted_multiple = dirs.liftover / f"{stem}.lifted.hg38.multiple.raw.bed"
        unmapped_multiple = dirs.liftover / f"{stem}.lifted.hg38.multiple.unmapped.txt"
        original_count, valid_input_count = standardize_peak_bed(source_path, std_bed)

        log(f"liftOver {peak_id}: {original_count} source rows")
        subprocess.run([liftover, str(std_bed), str(chain), str(lifted), str(unmapped)], check=True)
        subprocess.run(
            [liftover, "-multiple", str(std_bed), str(chain), str(lifted_multiple), str(unmapped_multiple)],
            check=True,
        )

        lifted_df = read_lifted_peak_bed(lifted, chrom_sizes)
        multi_df = read_lifted_peak_bed(lifted_multiple, chrom_sizes)
        success_ids = lifted_df["peak_original_id"].nunique() if not lifted_df.empty else 0
        multi_counts = multi_df.groupby("peak_original_id").size() if not multi_df.empty else pd.Series(dtype=int)
        multi_mapped_ids = int((multi_counts > 1).sum())

        if lifted_df.empty:
            filtered = lifted_df
        else:
            filtered = lifted_df[
                lifted_df["canonical"]
                & lifted_df["within_chrom_bounds"]
                & (lifted_df["length"] > 0)
                & (lifted_df["length_ratio"] >= 0.5)
                & (lifted_df["length_ratio"] <= 2.0)
            ].copy()

        merged = merge_intervals_from_df(filtered)
        merged_arrays = interval_dict_to_arrays(merged)
        merged_prefix = prefix_tracks(merged_arrays)
        peak_free = {}
        # Filled by the caller after allowed components are known.
        unmerged_bed = dirs.peaks / f"{stem}.liftover.hg38.canonical.filtered.bed"
        merged_bed = dirs.peaks / f"{stem}.liftover.hg38.canonical.filtered.merged.bed"
        write_unmerged_peak_bed(unmerged_bed, filtered, peak_id)
        write_bed(merged_bed, merged, extra_name=peak_id)

        row = {
            **{k: source[k] for k in source if k not in {"source_exists"}},
            "source_abs_path": str(source_path),
            "original_peak_count": int(original_count),
            "valid_source_interval_count": int(valid_input_count),
            "liftover_success_count": int(success_ids),
            "liftover_failure_count": int(original_count - success_ids),
            "liftover_success_rate": float(success_ids / original_count) if original_count else math.nan,
            "multiple_mapping_input_count": int(multi_mapped_ids),
            "successful_rows_raw": int(len(lifted_df)),
            "noncanonical_success_count": int((~lifted_df["canonical"]).sum()) if not lifted_df.empty else 0,
            "out_of_bounds_success_count": int((~lifted_df["within_chrom_bounds"]).sum()) if not lifted_df.empty else 0,
            "zero_or_negative_success_count": int((lifted_df["length"] <= 0).sum()) if not lifted_df.empty else 0,
            "abnormal_length_ratio_count_0.5_2": int(
                ((lifted_df["length_ratio"] < 0.5) | (lifted_df["length_ratio"] > 2.0)).sum()
            )
            if not lifted_df.empty
            else 0,
            "length_changed_gt10pct_count": int((np.abs(lifted_df["length_ratio"] - 1.0) > 0.10).sum())
            if not lifted_df.empty
            else 0,
            "analysis_unmerged_peak_count": int(len(filtered)),
            "analysis_merged_peak_count": int(sum(len(v) for v in merged.values())),
            "analysis_merged_coverage_bp": int(sum(e - s for vals in merged.values() for s, e in vals)),
            "filtered_peak_bed": str(unmerged_bed),
            "merged_peak_bed": str(merged_bed),
        }
        rows.append(row)
        tracks[peak_id] = Track(
            peak_set_id=peak_id,
            dataset=source["dataset"],
            mark=source["mark"],
            analysis_tier=source["analysis_tier"],
            intervals=merged_arrays,
            free_components=peak_free,
            prefix=merged_prefix,
            merged_interval_count=row["analysis_merged_peak_count"],
            coverage_bp=row["analysis_merged_coverage_bp"],
        )
    qc = pd.DataFrame(rows)
    qc.to_csv(dirs.tables / "chipseq_liftover_qc.tsv", sep="\t", index=False)
    return qc, tracks


def read_lifted_peak_bed(path: Path, chrom_sizes: Mapping[str, int]) -> pd.DataFrame:
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(
            columns=[
                "chrom",
                "start",
                "end",
                "peak_original_id",
                "original_length",
                "length",
                "length_ratio",
                "canonical",
                "within_chrom_bounds",
            ]
        )
    with open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, start_s, end_s, peak_original_id, orig_len_s = parts[:5]
            try:
                start = int(start_s)
                end = int(end_s)
                orig_len = int(float(orig_len_s))
            except ValueError:
                continue
            length = end - start
            canonical = chrom in CANON_SET
            within = bool(canonical and start >= 0 and end <= chrom_sizes.get(chrom, -1))
            ratio = float(length / orig_len) if orig_len else math.nan
            rows.append((chrom, start, end, peak_original_id, orig_len, length, ratio, canonical, within))
    return pd.DataFrame(
        rows,
        columns=[
            "chrom",
            "start",
            "end",
            "peak_original_id",
            "original_length",
            "length",
            "length_ratio",
            "canonical",
            "within_chrom_bounds",
        ],
    )


def write_unmerged_peak_bed(path: Path, df: pd.DataFrame, peak_set_id: str) -> None:
    with open(path, "wt", encoding="utf-8") as out:
        if df.empty:
            return
        work = df.sort_values(["chrom", "start", "end"])
        for i, row in enumerate(work.itertuples(index=False), 1):
            out.write(f"{row.chrom}\t{int(row.start)}\t{int(row.end)}\t{peak_set_id}_{i}\n")


def finalize_track_free_components(tracks: Dict[str, Track], allowed: Mapping[str, Sequence[Tuple[int, int]]]) -> None:
    for track in tracks.values():
        free = {}
        for chrom in CANON_CHROMS:
            starts, ends = track.intervals[chrom]
            blockers = list(zip(starts.tolist(), ends.tolist()))
            free[chrom] = subtract_intervals(allowed.get(chrom, []), blockers)
        track.free_components = free


def interval_hits(starts: np.ndarray, ends: np.ndarray, q_starts: np.ndarray, q_ends: np.ndarray) -> np.ndarray:
    if len(starts) == 0 or len(q_starts) == 0:
        return np.zeros(len(q_starts), dtype=bool)
    idx = np.searchsorted(starts, q_ends, side="left") - 1
    out = np.zeros(len(q_starts), dtype=bool)
    valid = idx >= 0
    if valid.any():
        out[valid] = ends[idx[valid]] > q_starts[valid]
    return out


def point_hits(starts: np.ndarray, ends: np.ndarray, points: np.ndarray) -> np.ndarray:
    if len(starts) == 0 or len(points) == 0:
        return np.zeros(len(points), dtype=bool)
    idx = np.searchsorted(starts, points, side="right") - 1
    out = np.zeros(len(points), dtype=bool)
    valid = idx >= 0
    if valid.any():
        out[valid] = ends[idx[valid]] > points[valid]
    return out


def load_sample_eccdna(path: Path, chrom_sizes: Mapping[str, int], forbidden: Mapping[str, Tuple[np.ndarray, np.ndarray]]) -> Tuple[pd.DataFrame, dict]:
    df = read_bed3(path, chrom_sizes)
    if df.empty:
        return df, {
            "original_eccdna_count": 0,
            "valid_coordinate_count": 0,
            "canonical_count": 0,
            "within_chrom_bounds_count": 0,
            "blacklist_gap_overlap_count": 0,
            "eligible_eccdna_count": 0,
        }
    valid_coordinate = df["valid_interval"]
    canonical = valid_coordinate & df["canonical"]
    within = canonical & df["within_chrom_bounds"]
    eligible_mask = within.to_numpy().copy()
    forbidden_hits = np.zeros(len(df), dtype=bool)
    for chrom in CANON_CHROMS:
        mask = (df["chrom"].to_numpy() == chrom) & within.to_numpy()
        if not mask.any():
            continue
        starts, ends = forbidden[chrom]
        idx = np.where(mask)[0]
        forbidden_hits[idx] = interval_hits(
            starts,
            ends,
            df["start"].to_numpy(dtype=np.int64)[idx],
            df["end"].to_numpy(dtype=np.int64)[idx],
        )
    eligible_mask &= ~forbidden_hits
    eligible = df.loc[eligible_mask, ["chrom", "start", "end", "length"]].copy()
    eligible["start"] = eligible["start"].astype(np.int64)
    eligible["end"] = eligible["end"].astype(np.int64)
    eligible["length"] = eligible["length"].astype(np.int64)
    qc = {
        "original_eccdna_count": int(len(df)),
        "valid_coordinate_count": int(valid_coordinate.sum()),
        "canonical_count": int(canonical.sum()),
        "within_chrom_bounds_count": int(within.sum()),
        "blacklist_gap_overlap_count": int((forbidden_hits & within.to_numpy()).sum()),
        "eligible_eccdna_count": int(len(eligible)),
    }
    return eligible, qc


class ProbabilityCalculator:
    def __init__(self, allowed: Mapping[str, Sequence[Tuple[int, int]]]):
        self.allowed = allowed
        self.cache: Dict[Tuple[str, str, int], dict] = {}

    def probabilities(self, track: Track, chrom: str, length: int) -> dict:
        key = (track.peak_set_id, chrom, int(length))
        if key in self.cache:
            return self.cache[key]
        comps = self.allowed.get(chrom, [])
        denom = 0
        midpoint_hits = 0
        start_hits = 0
        end_hits = 0
        for comp_s, comp_e in comps:
            span = comp_e - comp_s
            if span < length:
                continue
            placements = span - length + 1
            denom += placements
            k = length // 2
            midpoint_hits += coverage_in_range(track.prefix[chrom], comp_s + k, comp_e - length + k + 1)
            start_hits += coverage_in_range(track.prefix[chrom], comp_s, comp_e - length + 1)
            end_hits += coverage_in_range(track.prefix[chrom], comp_s + length - 1, comp_e)

        no_full_hits = 0
        for free_s, free_e in track.free_components.get(chrom, []):
            span = free_e - free_s
            if span >= length:
                no_full_hits += span - length + 1

        if denom <= 0:
            out = {
                "placement_count": 0,
                "p_full_interval": math.nan,
                "p_midpoint": math.nan,
                "p_junction_start": math.nan,
                "p_junction_end": math.nan,
                "p_junction_unit": math.nan,
            }
        else:
            full_hits = max(0, denom - no_full_hits)
            out = {
                "placement_count": int(denom),
                "p_full_interval": float(full_hits / denom),
                "p_midpoint": float(midpoint_hits / denom),
                "p_junction_start": float(start_hits / denom),
                "p_junction_end": float(end_hits / denom),
                "p_junction_unit": float((start_hits + end_hits) / (2 * denom)),
            }
        self.cache[key] = out
        return out


def aggregate_expected(
    length_bins: pd.DataFrame,
    eligible_count: int,
    track: Track,
    calculator: ProbabilityCalculator,
    mode: str,
) -> Tuple[float, float, int, str]:
    if eligible_count == 0:
        return math.nan, math.nan, 0, "no_eligible_eccdna"
    n = int(eligible_count)
    denom = 2 * n if mode == "junction_start_end" else n
    expected_hits = 0.0
    variance_hits = 0.0
    placement_bins = 0
    for row in length_bins.itertuples(index=False):
        chrom = str(row.chrom)
        length = int(row.length)
        count = int(row.count)
        probs = calculator.probabilities(track, chrom, length)
        placement_bins += 1
        if mode == "full_interval":
            p = probs["p_full_interval"]
            unit_count = count
        elif mode == "midpoint":
            p = probs["p_midpoint"]
            unit_count = count
        else:
            p = probs["p_junction_unit"]
            unit_count = 2 * count
        if math.isnan(p):
            return math.nan, math.nan, placement_bins, "no_valid_matched_placement"
        expected_hits += unit_count * p
        variance_hits += unit_count * p * (1.0 - p)
    return float(expected_hits / denom), float(math.sqrt(max(0.0, variance_hits)) / denom), placement_bins, "ok"


def aggregate_expected_all_modes(
    length_bins: pd.DataFrame,
    eligible_count: int,
    track: Track,
    calculator: ProbabilityCalculator,
) -> Dict[str, Tuple[float, float, int, str]]:
    if eligible_count == 0:
        return {mode: (math.nan, math.nan, 0, "no_eligible_eccdna") for mode in MODES}

    n = int(eligible_count)
    accum = {
        "full_interval": {"expected": 0.0, "variance": 0.0, "denom": n, "bins": 0, "status": "ok"},
        "midpoint": {"expected": 0.0, "variance": 0.0, "denom": n, "bins": 0, "status": "ok"},
        "junction_start_end": {"expected": 0.0, "variance": 0.0, "denom": 2 * n, "bins": 0, "status": "ok"},
    }

    for row in length_bins.itertuples(index=False):
        chrom = str(row.chrom)
        length = int(row.length)
        count = int(row.count)
        probs = calculator.probabilities(track, chrom, length)
        per_mode = {
            "full_interval": (count, probs["p_full_interval"]),
            "midpoint": (count, probs["p_midpoint"]),
            "junction_start_end": (2 * count, probs["p_junction_unit"]),
        }
        for mode, (unit_count, p) in per_mode.items():
            accum[mode]["bins"] += 1
            if math.isnan(p):
                accum[mode]["status"] = "no_valid_matched_placement"
                continue
            accum[mode]["expected"] += unit_count * p
            accum[mode]["variance"] += unit_count * p * (1.0 - p)

    out = {}
    for mode, vals in accum.items():
        if vals["status"] != "ok":
            out[mode] = (math.nan, math.nan, int(vals["bins"]), vals["status"])
        else:
            denom = float(vals["denom"])
            out[mode] = (
                float(vals["expected"] / denom),
                float(math.sqrt(max(0.0, vals["variance"])) / denom),
                int(vals["bins"]),
                "ok",
            )
    return out


def observed_for_track(ecc: pd.DataFrame, track: Track) -> Dict[str, dict]:
    out = {}
    n = int(len(ecc))
    for mode in MODES:
        out[mode] = {"denominator_units": 2 * n if mode == "junction_start_end" else n}
    if n == 0:
        for mode in MODES:
            out[mode].update(
                {
                    "observed_overlap_units": 0,
                    "observed_fraction": math.nan,
                    "junction_any_eccdna_count": math.nan,
                }
            )
        return out

    chrom_values = ecc["chrom"].to_numpy()
    starts_all = ecc["start"].to_numpy(dtype=np.int64)
    ends_all = ecc["end"].to_numpy(dtype=np.int64)
    lengths_all = ecc["length"].to_numpy(dtype=np.int64)
    full_hits = np.zeros(n, dtype=bool)
    midpoint_hits_arr = np.zeros(n, dtype=bool)
    start_hits_arr = np.zeros(n, dtype=bool)
    end_hits_arr = np.zeros(n, dtype=bool)
    midpoint_all = starts_all + (lengths_all // 2)
    for chrom in CANON_CHROMS:
        idx = np.where(chrom_values == chrom)[0]
        if len(idx) == 0:
            continue
        peak_starts, peak_ends = track.intervals[chrom]
        full_hits[idx] = interval_hits(peak_starts, peak_ends, starts_all[idx], ends_all[idx])
        midpoint_hits_arr[idx] = point_hits(peak_starts, peak_ends, midpoint_all[idx])
        start_hits_arr[idx] = point_hits(peak_starts, peak_ends, starts_all[idx])
        end_hits_arr[idx] = point_hits(peak_starts, peak_ends, ends_all[idx] - 1)

    junction_units = int(start_hits_arr.sum() + end_hits_arr.sum())
    junction_any = int((start_hits_arr | end_hits_arr).sum())
    out["full_interval"].update(
        {
            "observed_overlap_units": int(full_hits.sum()),
            "observed_fraction": float(full_hits.mean()),
            "junction_any_eccdna_count": math.nan,
        }
    )
    out["midpoint"].update(
        {
            "observed_overlap_units": int(midpoint_hits_arr.sum()),
            "observed_fraction": float(midpoint_hits_arr.mean()),
            "junction_any_eccdna_count": math.nan,
        }
    )
    out["junction_start_end"].update(
        {
            "observed_overlap_units": junction_units,
            "observed_fraction": float(junction_units / (2 * n)),
            "junction_any_eccdna_count": junction_any,
        }
    )
    return out


def analyze_samples(
    sample_manifest: pd.DataFrame,
    tracks: Dict[str, Track],
    dirs: Dirs,
    chrom_sizes: Mapping[str, int],
    forbidden_intervals: Mapping[str, Sequence[Tuple[int, int]]],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    forbidden_arrays = interval_dict_to_arrays(forbidden_intervals)
    calculator = ProbabilityCalculator(allowed)
    sample_qc_rows = []
    relative_rows = []
    absolute_rows = []
    progress_path = dirs.logs / "sample_analysis_progress.tsv"
    with open(progress_path, "wt", encoding="utf-8") as progress:
        progress.write("timestamp\tsample_id\tgroup\teligible_eccdna_count\tstatus\n")
    for sample in sample_manifest.sort_values(["group", "sample_id"]).to_dict(orient="records"):
        log(f"Analyzing sample {sample['sample_id']} ({sample['group']})")
        ecc, qc = load_sample_eccdna(Path(sample["eccdna_bed"]), chrom_sizes, forbidden_arrays)
        if ecc.empty:
            length_bins = pd.DataFrame(columns=["chrom", "length", "count"])
        else:
            length_bins = ecc.groupby(["chrom", "length"], sort=False).size().reset_index(name="count")
        qc_row = {**sample, **qc}
        sample_qc_rows.append(qc_row)
        mapped = int(sample["mapped_alignments_idxstats"])
        for track in tracks.values():
            obs_by_mode = observed_for_track(ecc, track)
            expected_by_mode = aggregate_expected_all_modes(length_bins, len(ecc), track, calculator)
            for mode in MODES:
                expected_fraction, null_sd, placement_bins, status = expected_by_mode[mode]
                obs = obs_by_mode[mode]
                observed_fraction = obs["observed_fraction"]
                if expected_fraction and not math.isnan(expected_fraction) and expected_fraction > 0:
                    ratio = observed_fraction / expected_fraction if not math.isnan(observed_fraction) else math.nan
                    log2_ratio = math.log2(ratio) if ratio > 0 else -math.inf
                else:
                    ratio = math.nan
                    log2_ratio = math.nan
                common = {
                    "sample_id": sample["sample_id"],
                    "group": sample["group"],
                    "peak_set_id": track.peak_set_id,
                    "dataset": track.dataset,
                    "mark": track.mark,
                    "analysis_tier": track.analysis_tier,
                    "mode": mode,
                    "eligible_eccdna_count": int(len(ecc)),
                    "denominator_units": int(obs["denominator_units"]),
                    "mapped_alignments_idxstats": mapped,
                    "observed_overlap_units": int(obs["observed_overlap_units"]),
                    "observed_fraction": observed_fraction,
                    "mean_shuffled_overlap_fraction": expected_fraction,
                    "shuffled_fraction_sd_exact_or_approx": null_sd,
                    "enrichment_ratio": ratio,
                    "log2_enrichment_ratio": log2_ratio,
                    "matched_control_method": "exact_all_valid_chromosome_and_length_matched_placements_excluding_hg38_blacklist_gaps",
                    "placement_bins": int(placement_bins),
                    "matched_control_status": status,
                }
                relative_rows.append(common)
                epm = (obs["observed_overlap_units"] / mapped * 1_000_000.0) if mapped > 0 else math.nan
                abs_row = {
                    **common,
                    "absolute_overlap_count_or_units": int(obs["observed_overlap_units"]),
                    "absolute_epm_per_total_mapped_alignments": epm,
                    "absolute_metric_note": "junction_start_end is breakpoint-unit count; full_interval and midpoint are eccDNA counts",
                    "junction_any_eccdna_count": obs.get("junction_any_eccdna_count", math.nan),
                    "junction_any_eccdna_epm": (
                        obs.get("junction_any_eccdna_count", math.nan) / mapped * 1_000_000.0
                        if mode == "junction_start_end" and mapped > 0
                        else math.nan
                    ),
                }
                absolute_rows.append(abs_row)
        with open(progress_path, "a", encoding="utf-8") as progress:
            progress.write(
                f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{sample['sample_id']}\t{sample['group']}\t{len(ecc)}\tok\n"
            )
    sample_qc = pd.DataFrame(sample_qc_rows)
    relative = pd.DataFrame(relative_rows)
    absolute = pd.DataFrame(absolute_rows)
    sample_qc.to_csv(dirs.tables / "eccdna_sample_qc.tsv", sep="\t", index=False)
    relative.to_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t", index=False)
    absolute.to_csv(dirs.tables / "sample_absolute_burden.tsv", sep="\t", index=False)
    return sample_qc, relative, absolute


def bh_fdr(pvalues: Sequence[float]) -> List[float]:
    p = np.asarray([np.nan if v is None else v for v in pvalues], dtype=float)
    out = np.full(len(p), np.nan, dtype=float)
    valid = np.where(~np.isnan(p))[0]
    if len(valid) == 0:
        return out.tolist()
    order = valid[np.argsort(p[valid])]
    ranked = p[order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.minimum(q, 1.0)
    return out.tolist()


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) == 0 or len(y) == 0:
        return math.nan
    # For n=39 per group this direct computation is clearer and inexpensive.
    comp = x[:, None] - y[None, :]
    return float((np.sum(comp > 0) - np.sum(comp < 0)) / comp.size)


def group_stats(df: pd.DataFrame, value_col: str, outpath: Path) -> pd.DataFrame:
    try:
        from scipy import stats
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("scipy is required for Mann-Whitney group statistics in the QDUH Miniforge environment.") from exc

    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode"]
    for key_vals, sub in df.groupby(keys, dropna=False):
        covid = sub.loc[sub["group"] == "COVID", value_col].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
        hc = sub.loc[sub["group"] == "HC", value_col].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
        if len(covid) > 0 and len(hc) > 0:
            test = stats.mannwhitneyu(covid, hc, alternative="two-sided")
            p = float(test.pvalue)
            delta = cliffs_delta(covid, hc)
        else:
            p = math.nan
            delta = math.nan
        rows.append(
            {
                **dict(zip(keys, key_vals)),
                "metric": value_col,
                "n_covid": int(len(covid)),
                "n_hc": int(len(hc)),
                "median_covid": float(np.median(covid)) if len(covid) else math.nan,
                "median_hc": float(np.median(hc)) if len(hc) else math.nan,
                "mean_covid": float(np.mean(covid)) if len(covid) else math.nan,
                "mean_hc": float(np.mean(hc)) if len(hc) else math.nan,
                "median_difference_covid_minus_hc": (
                    float(np.median(covid) - np.median(hc)) if len(covid) and len(hc) else math.nan
                ),
                "mannwhitneyu_p": p,
                "cliffs_delta_covid_vs_hc": delta,
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df["bh_fdr_all_tests"] = bh_fdr(stats_df["mannwhitneyu_p"].tolist())
    stats_df.to_csv(outpath, sep="\t", index=False)
    return stats_df


def run_stats(dirs: Dirs, relative: pd.DataFrame, absolute: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    relative_stats = group_stats(relative, "log2_enrichment_ratio", dirs.tables / "relative_enrichment_group_stats.tsv")
    absolute_stats = group_stats(
        absolute,
        "absolute_epm_per_total_mapped_alignments",
        dirs.tables / "absolute_burden_group_stats.tsv",
    )
    observed_stats = group_stats(relative, "observed_fraction", dirs.tables / "observed_fraction_group_stats.tsv")
    return relative_stats, absolute_stats, observed_stats


def setup_matplotlib():
    import matplotlib as mpl
    from matplotlib import font_manager

    available = {font.name for font in font_manager.fontManager.ttflist}
    font_family = next(
        (name for name in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"] if name in available),
        "DejaVu Sans",
    )

    mpl.rcParams.update(
        {
            "font.family": font_family,
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.55,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "font.size": 6,
            "axes.labelsize": 6,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
        }
    )


def make_box_panel(ax, df: pd.DataFrame, value_col: str, marks: Sequence[str], ylabel: str, ylim: Optional[Tuple[float, float]] = None) -> None:
    colors = {"HC": "#0272B2", "COVID": "#EC6F00"}
    rng = np.random.default_rng(SEED)
    width = 0.28
    offsets = {"HC": -width / 1.5, "COVID": width / 1.5}
    for i, mark in enumerate(marks):
        for group in ("HC", "COVID"):
            values = (
                df[(df["mark"] == mark) & (df["group"] == group)][value_col]
                .replace([np.inf, -np.inf], np.nan)
                .dropna()
                .to_numpy(float)
            )
            pos = i + offsets[group]
            if len(values) == 0:
                continue
            ax.boxplot(
                values,
                positions=[pos],
                widths=width,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "#000000", "linewidth": 0.6},
                boxprops={"facecolor": colors[group], "edgecolor": "#253247", "linewidth": 0.45, "alpha": 0.55},
                whiskerprops={"color": "#253247", "linewidth": 0.45},
                capprops={"color": "#253247", "linewidth": 0.45},
            )
            jitter = rng.normal(0, 0.025, len(values))
            ax.scatter(
                np.full(len(values), pos) + jitter,
                values,
                s=5,
                facecolors=colors[group],
                edgecolors="#253247",
                linewidths=0.2,
                alpha=0.75,
                zorder=3,
            )
    ax.set_xticks(range(len(marks)))
    ax.set_xticklabels(marks, rotation=35, ha="right")
    ax.set_ylabel(ylabel)
    ax.axhline(0, color="#99A3B4", linewidth=0.55, zorder=0)
    ax.grid(axis="y", color="#E6E6ED", linewidth=0.45)
    if ylim is not None:
        ax.set_ylim(*ylim)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def make_heatmap_panel(ax, stats_df: pd.DataFrame, marks: Sequence[str], modes: Sequence[str]) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    mode_labels = {"full_interval": "Full", "junction_start_end": "Junction", "midpoint": "Midpoint"}
    mat = np.full((len(modes), len(marks)), np.nan)
    for i, mode in enumerate(modes):
        for j, mark in enumerate(marks):
            row = stats_df[(stats_df["mode"] == mode) & (stats_df["mark"] == mark)]
            if not row.empty:
                mat[i, j] = row.iloc[0]["median_difference_covid_minus_hc"]
    vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
    vmax = max(vmax, 0.25)
    cmap = LinearSegmentedColormap.from_list("nature_blue_red", ["#0272B2", "#FFFFFF", "#C93E3F"])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(marks)))
    ax.set_xticklabels(marks, rotation=35, ha="right")
    ax.set_yticks(range(len(modes)))
    ax.set_yticklabels([mode_labels.get(m, m) for m in modes])
    for i in range(len(modes)):
        for j in range(len(marks)):
            val = mat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5, color="#000000")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Median delta log2 enrichment", fontsize=5.5)
    cbar.ax.tick_params(labelsize=5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.45)


def make_figures(dirs: Dirs) -> None:
    import matplotlib.pyplot as plt

    setup_matplotlib()
    relative = pd.read_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t")
    absolute = pd.read_csv(dirs.tables / "sample_absolute_burden.tsv", sep="\t")
    rel_stats = pd.read_csv(dirs.tables / "relative_enrichment_group_stats.tsv", sep="\t")

    cd14_rel = relative[(relative["dataset"] == "CD14_reference") & (relative["mode"] == "full_interval")].copy()
    cd14_abs = absolute[(absolute["dataset"] == "CD14_reference") & (absolute["mode"] == "full_interval")].copy()
    marks = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]
    cd14_stats = rel_stats[rel_stats["dataset"] == "CD14_reference"].copy()

    width = 183 / 25.4
    height = 128 / 25.4
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1.0], height_ratios=[1.0, 1.0], wspace=0.38, hspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    make_box_panel(ax_a, cd14_rel, "log2_enrichment_ratio", marks, "Relative enrichment\nlog2(observed/expected)")
    make_heatmap_panel(ax_b, cd14_stats, marks, MODES)
    abs_plot = cd14_abs.copy()
    abs_plot["log10_epm_plus1"] = np.log10(abs_plot["absolute_epm_per_total_mapped_alignments"] + 1.0)
    make_box_panel(ax_c, abs_plot, "log10_epm_plus1", marks, "Absolute burden\nlog10(EPM + 1)")

    secondary = relative[
        (relative["dataset"] == "A549_GSE179184_official") & (relative["mode"] == "full_interval")
    ].copy()
    secondary_marks = [m for m in ["H3K27ac", "H3K4me3"] if m in set(secondary["mark"])]
    make_box_panel(
        ax_d,
        secondary,
        "log2_enrichment_ratio",
        secondary_marks,
        "A549 official peaks\nlog2 enrichment",
    )

    for label, ax in zip(("a", "b", "c", "d"), (ax_a, ax_b, ax_c, ax_d)):
        ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top", ha="left")

    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#0272B2", markeredgecolor="#253247", markersize=4, label="HC"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="#EC6F00", markeredgecolor="#253247", markersize=4, label="COVID-19"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.52, 0.995), ncol=2, frameon=False)

    for ext in ("pdf", "svg", "png"):
        out = dirs.figures / f"chromatin_hg38_relative_and_absolute.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_run_summary(
    dirs: Dirs,
    args: argparse.Namespace,
    sample_qc: Optional[pd.DataFrame] = None,
    liftover_qc: Optional[pd.DataFrame] = None,
    rel_stats: Optional[pd.DataFrame] = None,
    abs_stats: Optional[pd.DataFrame] = None,
) -> None:
    if sample_qc is None and (dirs.tables / "eccdna_sample_qc.tsv").exists():
        sample_qc = pd.read_csv(dirs.tables / "eccdna_sample_qc.tsv", sep="\t")
    if liftover_qc is None and (dirs.tables / "chipseq_liftover_qc.tsv").exists():
        liftover_qc = pd.read_csv(dirs.tables / "chipseq_liftover_qc.tsv", sep="\t")
    if rel_stats is None and (dirs.tables / "relative_enrichment_group_stats.tsv").exists():
        rel_stats = pd.read_csv(dirs.tables / "relative_enrichment_group_stats.tsv", sep="\t")
    if abs_stats is None and (dirs.tables / "absolute_burden_group_stats.tsv").exists():
        abs_stats = pd.read_csv(dirs.tables / "absolute_burden_group_stats.tsv", sep="\t")

    lines = []
    lines.append("# Chromatin relative enrichment redo (hg38)")
    lines.append("")
    lines.append(f"Run time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Host: {socket.gethostname()}")
    lines.append(f"Project directory: `{args.project_dir}`")
    lines.append(f"Output directory: `{args.outdir}`")
    lines.append("")
    lines.append("## Coordinate system")
    lines.append("")
    lines.append("- eccDNA BED inputs were treated as hg38 and filtered to canonical chr1-22/X/Y.")
    lines.append("- All ChIP-seq peak BEDs in the manifest were lifted from hg19 to hg38 with UCSC liftOver.")
    lines.append("- Primary analysis excluded hg38 blacklist and assembly-gap intervals from both observed eccDNA denominators and the matched placement space.")
    lines.append("")
    lines.append("## Relative enrichment definition")
    lines.append("")
    lines.append("- Observed fraction = overlapping eligible eccDNA units / total eligible eccDNA units for each sample.")
    lines.append("- Expected fraction = exact mean over all chromosome- and length-matched placements in hg38 allowed regions.")
    lines.append("- Enrichment ratio = observed fraction / expected fraction; values above 1 indicate relative enrichment after controlling for total eccDNA burden.")
    lines.append("- Modes: full interval, junction/start-end breakpoint units, and midpoint.")
    lines.append("")
    lines.append("## Output tables")
    lines.append("")
    for rel in [
        "inputs/eccdna_sample_manifest.tsv",
        "inputs/chipseq_peak_sources.manifest.tsv",
        "reference/reference_manifest.tsv",
        "reference/hg38_blacklist_gap_summary.tsv",
        "tables/chipseq_liftover_qc.tsv",
        "tables/eccdna_sample_qc.tsv",
        "tables/sample_relative_enrichment.tsv",
        "tables/sample_absolute_burden.tsv",
        "tables/relative_enrichment_group_stats.tsv",
        "tables/absolute_burden_group_stats.tsv",
        "figures/chromatin_hg38_relative_and_absolute.pdf",
        "figures/chromatin_hg38_relative_and_absolute.svg",
        "figures/chromatin_hg38_relative_and_absolute.png",
    ]:
        lines.append(f"- `{rel}`")
    lines.append("")
    if liftover_qc is not None and not liftover_qc.empty:
        lines.append("## liftOver QC overview")
        lines.append("")
        cols = [
            "peak_set_id",
            "original_peak_count",
            "liftover_success_count",
            "liftover_failure_count",
            "liftover_success_rate",
            "multiple_mapping_input_count",
            "abnormal_length_ratio_count_0.5_2",
            "analysis_merged_peak_count",
        ]
        lines.append(markdown_table(liftover_qc[cols]))
        lines.append("")
    if sample_qc is not None and not sample_qc.empty:
        lines.append("## eccDNA QC overview")
        lines.append("")
        grp = sample_qc.groupby("group").agg(
            samples=("sample_id", "count"),
            original_eccdna=("original_eccdna_count", "sum"),
            eligible_eccdna=("eligible_eccdna_count", "sum"),
            blacklist_gap_overlap=("blacklist_gap_overlap_count", "sum"),
        )
        lines.append(markdown_table(grp.reset_index()))
        lines.append("")
    if rel_stats is not None and not rel_stats.empty:
        lines.append("## Main CD14 relative enrichment statistics")
        lines.append("")
        main = rel_stats[(rel_stats["dataset"] == "CD14_reference") & (rel_stats["mode"] == "full_interval")].copy()
        if not main.empty:
            keep = [
                "mark",
                "median_covid",
                "median_hc",
                "median_difference_covid_minus_hc",
                "mannwhitneyu_p",
                "bh_fdr_all_tests",
            ]
            lines.append(markdown_table(main[keep]))
            lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- Absolute burden is reported separately as peak-overlapping count/EPM per total mapped alignments and is not used as evidence for preferential enrichment.")
    lines.append("- GC, mappability, and repeat-class matching were not added because no validated matched hg38 resources were present in the project directory; chromosome and exact-length matching plus blacklist/gap exclusion are the primary reviewer-requested controls.")
    lines.append("- The A549 legacy thresholded peak sets are retained only as sensitivity outputs; official GSE179184 processed peaks are reported separately.")
    (dirs.outdir / "README_chromatin_hg38.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    show = df.copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "NA" if pd.isna(x) else f"{x:.4g}")
        else:
            show[col] = show[col].map(lambda x: "NA" if pd.isna(x) else str(x))
    headers = list(show.columns)
    rows = show.values.tolist()
    widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            widths[i] = max(widths[i], len(str(val)))
    sep = "| " + " | ".join("-" * w for w in widths) + " |"
    header = "| " + " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers)) + " |"
    body = ["| " + " | ".join(str(v).ljust(widths[i]) for i, v in enumerate(row)) + " |" for row in rows]
    return "\n".join([header, sep] + body)


def run_prepare(args: argparse.Namespace) -> None:
    dirs = ensure_dirs(Path(args.outdir))
    prepare_reference_files(dirs, force_download=args.force_download)
    chrom_sizes = read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    build_forbidden_and_allowed(dirs, chrom_sizes)
    make_peak_manifest(Path(args.project_dir), dirs)
    write_run_summary(dirs, args)


def run_analyze(args: argparse.Namespace) -> None:
    dirs = ensure_dirs(Path(args.outdir))
    prepare_reference_files(dirs, force_download=False)
    chrom_sizes = read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    allowed = build_forbidden_and_allowed(dirs, chrom_sizes)
    forbidden_df = read_bed3(dirs.reference / "hg38_forbidden.blacklist_gaps.canonical.merged.bed", chrom_sizes)
    forbidden = merge_intervals_from_df(forbidden_df)
    peak_manifest = make_peak_manifest(Path(args.project_dir), dirs)
    sample_manifest = make_sample_manifest(Path(args.project_dir), dirs, args.samtools)
    liftover_qc, tracks = run_liftover_for_peaks(
        peak_manifest,
        dirs,
        dirs.reference / "hg19ToHg38.over.chain.gz",
        args.liftover,
        chrom_sizes,
    )
    finalize_track_free_components(tracks, allowed)
    sample_qc, relative, absolute = analyze_samples(sample_manifest, tracks, dirs, chrom_sizes, forbidden, allowed)
    rel_stats, abs_stats, _ = run_stats(dirs, relative, absolute)
    make_figures(dirs)
    write_run_summary(dirs, args, sample_qc=sample_qc, liftover_qc=liftover_qc, rel_stats=rel_stats, abs_stats=abs_stats)
    manifest = {
        "project_dir": str(args.project_dir),
        "outdir": str(args.outdir),
        "stage": args.stage,
        "seed": SEED,
        "python": sys.version,
        "hostname": socket.gethostname(),
        "liftover": args.liftover,
        "samtools": args.samtools,
        "control_method": "exact chromosome- and exact-length-matched placement mean excluding hg38 blacklist/gaps",
    }
    (dirs.outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_figures_only(args: argparse.Namespace) -> None:
    dirs = ensure_dirs(Path(args.outdir))
    make_figures(dirs)
    write_run_summary(dirs, args)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default="/gpfs/data/gao/Covid-eccDNA", help="QDUH Covid-eccDNA project directory.")
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin", help="Output directory.")
    parser.add_argument(
        "--stage",
        choices=["prepare", "analyze", "figures"],
        default="analyze",
        help="prepare downloads/references only; analyze runs full workflow; figures rebuilds plots from tables.",
    )
    parser.add_argument("--liftover", default="/gpfs/softwares/ricopili202511/dependencies/liftover/liftOver")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--force-download", action="store_true", help="Redownload small reference files.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.stage == "prepare":
        run_prepare(args)
    elif args.stage == "analyze":
        run_analyze(args)
    elif args.stage == "figures":
        run_figures_only(args)
    else:
        raise AssertionError(args.stage)
    log("Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
