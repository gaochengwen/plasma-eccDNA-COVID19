#!/usr/bin/env python3
"""Discover and download ENCODE histone ChIP-seq peak files already released on GRCh38.

Reviewer 3 major point 4 asks whether the chromatin association holds in cell
types other than CD14+ monocytes. Using ENCODE files that are natively GRCh38
avoids a second liftOver step, so the cross-cell-type comparison does not inherit
any conversion artefact from the primary analysis.

For each requested (biosample, target) pair this selects the experiment's
preferred default peak BED, preferring replicated over pseudoreplicated peaks,
and writes a manifest that `run_multicell_chromatin.py` consumes.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Sequence

ENCODE_BASE = "https://www.encodeproject.org"

# Cell types chosen for relevance to COVID-19: circulating innate and adaptive
# immune cells, the respiratory epithelium, and vascular endothelium. CD14+
# monocytes are included as a GRCh38-native counterpart to the lifted hg19
# reference set used in the primary analysis.
DEFAULT_BIOSAMPLES = [
    "CD14-positive monocyte",
    "neutrophil",
    "CD4-positive, alpha-beta T cell",
    "B cell",
    "upper lobe of left lung",
    "endothelial cell of umbilical vein",
]

DEFAULT_TARGETS = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3", "H3K9ac"]


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def api_get(path: str, params: Dict[str, object]) -> dict:
    query = []
    for key, value in params.items():
        if isinstance(value, (list, tuple)):
            query.extend((key, str(v)) for v in value)
        else:
            query.append((key, str(value)))
    url = f"{ENCODE_BASE}{path}?{urllib.parse.urlencode(query)}"
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pick_peak_file(experiment_accession: str) -> Optional[dict]:
    """Return the best GRCh38 peak BED for an experiment, or None."""
    detail = api_get(f"/experiments/{experiment_accession}/", {"frame": "embedded"})
    candidates = []
    for file_record in detail.get("files", []):
        if file_record.get("status") != "released":
            continue
        if file_record.get("assembly") != "GRCh38":
            continue
        if file_record.get("file_format") != "bed":
            continue
        output_type = file_record.get("output_type", "")
        if "peaks" not in output_type:
            continue
        # Prefer the file ENCODE itself flags as the default representation,
        # then replicated over pseudoreplicated peaks.
        score = 0
        if file_record.get("preferred_default"):
            score += 100
        if output_type == "replicated peaks":
            score += 50
        elif output_type == "pseudoreplicated peaks":
            score += 30
        elif output_type == "stable peaks":
            score += 40
        candidates.append((score, file_record))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def discover(biosamples: Sequence[str], targets: Sequence[str]) -> List[dict]:
    rows: List[dict] = []
    for biosample in biosamples:
        for target in targets:
            params = {
                "type": "Experiment",
                "assay_title": "Histone ChIP-seq",
                "biosample_ontology.term_name": biosample,
                "target.label": target,
                "status": "released",
                "assembly": "GRCh38",
                "limit": 25,
                "format": "json",
            }
            try:
                result = api_get("/search/", params)
            except Exception as exc:  # noqa: BLE001
                log(f"query failed for {biosample} / {target}: {exc}")
                continue
            experiments = result.get("@graph", [])
            if not experiments:
                log(f"no GRCh38 experiment: {biosample} / {target}")
                continue
            chosen = None
            for experiment in experiments:
                accession = experiment["accession"]
                try:
                    peak_file = pick_peak_file(accession)
                except Exception as exc:  # noqa: BLE001
                    log(f"file lookup failed for {accession}: {exc}")
                    continue
                if peak_file is not None:
                    chosen = (experiment, peak_file)
                    break
            if chosen is None:
                log(f"no GRCh38 peak BED: {biosample} / {target}")
                continue
            experiment, peak_file = chosen
            rows.append(
                {
                    "biosample": biosample,
                    "target": target,
                    "experiment_accession": experiment["accession"],
                    "file_accession": peak_file["accession"],
                    "output_type": peak_file.get("output_type", ""),
                    "file_type": peak_file.get("file_type", ""),
                    "assembly": peak_file.get("assembly", ""),
                    "href": ENCODE_BASE + peak_file["href"],
                    "encode_md5": peak_file.get("md5sum", ""),
                    "biological_replicates": ";".join(
                        str(v) for v in peak_file.get("biological_replicates", [])
                    ),
                }
            )
            log(
                f"selected {biosample} / {target}: {experiment['accession']} "
                f"{peak_file['accession']} ({peak_file.get('output_type')})"
            )
    return rows


def download(rows: List[dict], outdir: Path) -> List[dict]:
    outdir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        suffix = ".bed.gz"
        stem = f"ENCODE_{row['biosample'].replace(' ', '_').replace(',', '')}_{row['target']}"
        stem = "".join(ch if (ch.isalnum() or ch in "_-") else "_" for ch in stem)
        gz_path = outdir / f"{stem}.{row['file_accession']}{suffix}"
        bed_path = outdir / f"{stem}.{row['file_accession']}.bed"
        if not bed_path.exists():
            log(f"downloading {row['file_accession']} -> {gz_path.name}")
            request = urllib.request.Request(row["href"], headers={"Accept": "*/*"})
            with urllib.request.urlopen(request, timeout=600) as response, open(gz_path, "wb") as out:
                shutil.copyfileobj(response, out)
            with gzip.open(gz_path, "rb") as src, open(bed_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
        row["local_bed"] = str(bed_path)
        row["local_sha256"] = sha256_file(bed_path)
        with open(bed_path, "rt", encoding="utf-8") as handle:
            row["peak_count"] = sum(1 for line in handle if line.strip() and not line.startswith("#"))
        log(f"  {bed_path.name}: {row['peak_count']} peaks")
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin/encode_hg38")
    parser.add_argument("--biosamples", default=";".join(DEFAULT_BIOSAMPLES))
    parser.add_argument("--targets", default=",".join(DEFAULT_TARGETS))
    parser.add_argument("--discover-only", action="store_true")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    biosamples = [b.strip() for b in args.biosamples.split(";") if b.strip()]
    targets = [t.strip() for t in args.targets.split(",") if t.strip()]

    rows = discover(biosamples, targets)
    log(f"{len(rows)} peak sets discovered")
    if not args.discover_only:
        rows = download(rows, outdir)

    manifest = outdir / "encode_hg38_peak_manifest.tsv"
    if rows:
        columns = list(rows[0].keys())
        with open(manifest, "wt", encoding="utf-8") as out:
            out.write("\t".join(columns) + "\n")
            for row in rows:
                out.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
        log(f"wrote {manifest}")
    else:
        log("no rows; manifest not written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
