#!/usr/bin/env python3
"""Collect per-sample sequencing QC from the QDUH cluster into Table S2 inputs.

Reads the `samtools flagstat` and `samtools markdup` reports produced by
`seq_qc_array.pbs` under /gpfs/data/gao/Covid-eccDNA/Revise/seq_qc, cross-checks
them against each other, and writes:

  revise/source_data/sequencing_qc_per_sample.tsv
  revise/source_data/sequencing_qc_summary.json

Definitions, which are also stated in the Table S2 legend:
  read pairs   - "paired in sequencing" / 2. The original FASTQs were not
                 retained, so this is the adapter- and quality-trimmed pair
                 count presented to BWA-MEM, not a pre-trimming raw count.
  mapping rate - primary mapped / primary alignments.
  duplicate rate - markdup DUPLICATE PRIMARY TOTAL / EXAMINED. Duplicates were
                 not marked in the archived alignments and were marked for this
                 revision.

Integrity check: markdup READ must equal flagstat "in total" for every sample.
A mismatch means the streaming pipeline was interrupted and the sample is
rejected rather than silently reported.
"""
from __future__ import annotations

import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path("/home/gao/eccDNA")
OUT_DIR = ROOT / "revise" / "source_data"
REMOTE = "/gpfs/data/gao/Covid-eccDNA/Revise/seq_qc"
HOST = "QDUH"


def fetch(kind: str) -> dict[str, str]:
    """Return {sample: file text} for every report of the given kind."""
    listing = subprocess.run(
        ["ssh", HOST, f"cd {REMOTE}/{kind} && for f in *.txt; do echo \"##$f\"; cat \"$f\"; done"],
        capture_output=True, text=True, timeout=600,
    )
    if listing.returncode != 0:
        raise SystemExit(f"ssh failed for {kind}: {listing.stderr[:400]}")
    out: dict[str, str] = {}
    name, buf = None, []
    for line in listing.stdout.splitlines():
        if line.startswith("##"):
            if name:
                out[name] = "\n".join(buf)
            name = line[2:].split(".")[0]
            buf = []
        else:
            buf.append(line)
    if name:
        out[name] = "\n".join(buf)
    return out


def parse_flagstat(text: str) -> dict[str, int]:
    def grab(pattern: str) -> int:
        m = re.search(rf"^(\d+) \+ \d+ {pattern}", text, re.M)
        if not m:
            raise ValueError(f"flagstat field not found: {pattern}")
        return int(m.group(1))

    return {
        "total": grab(r"in total"),
        "primary": grab(r"primary$"),
        "primary_mapped": grab(r"primary mapped"),
        "mapped_total": grab(r"mapped \("),
        "paired_in_sequencing": grab(r"paired in sequencing"),
    }


def parse_markdup(text: str) -> dict[str, int]:
    vals = {}
    for key in ("READ", "EXAMINED", "DUPLICATE PRIMARY TOTAL"):
        m = re.search(rf"^{re.escape(key)}: (\d+)$", text, re.M)
        if not m:
            raise ValueError(f"markdup field not found: {key}")
        vals[key] = int(m.group(1))
    return vals


def main() -> int:
    samples = [
        line.split("\t")[0]
        for line in (ROOT / "revise" / "source_data" / "seq_qc_samples.tsv").read_text().splitlines()
        if line.strip()
    ]
    fs, md = fetch("flagstat"), fetch("markdup")

    rows, bad, missing = [], [], []
    for s in samples:
        if s not in fs or s not in md:
            missing.append(s)
            continue
        f, m = parse_flagstat(fs[s]), parse_markdup(md[s])
        if m["READ"] != f["total"]:
            bad.append((s, m["READ"], f["total"]))
            continue
        rows.append({
            "sample": s,
            "read_pairs": f["paired_in_sequencing"] // 2,
            "mapping_rate": 100.0 * f["primary_mapped"] / f["primary"],
            "duplicate_rate": 100.0 * m["DUPLICATE PRIMARY TOTAL"] / m["EXAMINED"],
            # Carried so the EPM denominator already in Table S2 can be re-checked
            # against the alignments; the two must agree exactly.
            "mapped_total": f["mapped_total"],
        })

    if missing:
        print(f"still missing ({len(missing)}): {' '.join(missing)}")
    if bad:
        print("INTEGRITY FAILURE - markdup READ != flagstat total (rerun these):")
        for s, r, t in bad:
            print(f"  {s}: READ={r:,} total={t:,}")
    if missing or bad:
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tsv = OUT_DIR / "sequencing_qc_per_sample.tsv"
    with tsv.open("w") as fh:
        fh.write("sample\tread_pairs\tmapping_rate_pct\tduplicate_rate_pct\tmapped_total\n")
        for r in rows:
            fh.write(f"{r['sample']}\t{r['read_pairs']}\t{r['mapping_rate']:.4f}\t"
                     f"{r['duplicate_rate']:.4f}\t{r['mapped_total']}\n")

    pairs = [r["read_pairs"] / 1e6 for r in rows]
    mr = [r["mapping_rate"] for r in rows]
    dr = [r["duplicate_rate"] for r in rows]
    summary = {
        "n_samples": len(rows),
        "pairs_min_m": min(pairs), "pairs_max_m": max(pairs), "pairs_med_m": statistics.median(pairs),
        "maprate_min": min(mr), "maprate_max": max(mr), "maprate_med": statistics.median(mr),
        "duprate_min": min(dr), "duprate_max": max(dr), "duprate_med": statistics.median(dr),
    }
    (OUT_DIR / "sequencing_qc_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"wrote {tsv} ({len(rows)} samples)")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
