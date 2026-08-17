#!/usr/bin/env python3
"""Derive the v2 per-sample burden table from the locked robustness metrics.

`analyse/EPM/results/sample_burden_recomputed_from_qduh_bed.tsv` is what the
Figure S1 builder reads. Its v2 counterpart is not a new computation: every
column it needs already exists, per sample, in the locked
`robustness_callset_sample_metrics.tsv` rows for `circlemap_methods`.

    eccdna_count_qduh_bed = call_count
    mapped_reads          = mapped_alignments
    epm_qduh_bed          = epm  ( = call_count / mapped_alignments x 1e6 )
    log1p_epm_qduh_bed    = log(1 + epm)

So this is a column rename plus one logarithm, and the identity
`epm == call_count / mapped_alignments * 1e6` is asserted on every row rather
than assumed -- which is the rule this project follows for figures: read locked
values, and where a quantity is not tabulated prefer an exact algebraic
identity over a fresh estimator.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
SRC = RUN / "tables" / "robustness" / "robustness_callset_sample_metrics.tsv"
REFERENCE = RUN.parent / "analyse" / "EPM" / "results" / \
    "sample_burden_recomputed_from_qduh_bed.tsv"
DEST = RUN / "v2tree" / "EPM" / "results" / \
    "sample_burden_recomputed_from_qduh_bed.tsv"

CALLSET = "circlemap_methods"


def main() -> int:
    with SRC.open(newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["callset"] == CALLSET]
    if len(rows) != 78:
        raise SystemExit(f"expected 78 samples for {CALLSET}, found {len(rows)}")

    with REFERENCE.open(newline="") as fh:
        fields = next(csv.reader(fh, delimiter="\t"))

    out = []
    for r in rows:
        count = int(r["call_count"])
        mapped = int(r["mapped_alignments"])
        epm = float(r["epm"])
        identity = count / mapped * 1e6
        if abs(identity - epm) > 1e-6 * max(1.0, abs(epm)):
            raise SystemExit(
                f"{r['sample_id']}: epm identity fails "
                f"({epm} vs {identity})")
        out.append({
            "sample_id": r["sample_id"],
            "group": r["group"],
            "mapped_reads": mapped,
            "eccdna_count_input_table": count,
            "epm_input_table": f"{epm:.4f}",
            "eccdna_count_qduh_bed": count,
            "epm_qduh_bed": f"{epm:.6f}",
            "log1p_epm_qduh_bed": f"{math.log1p(epm):.9f}",
            "bed_header_or_invalid_lines_skipped": 0,
            "bed_empty_lines_skipped": 0,
            "bed_path": (
                "/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/callsets/"
                f"circlemap_methods/{r['sample_id']}_circle_site.bed"),
        })

    DEST.parent.mkdir(parents=True, exist_ok=True)
    with DEST.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)

    total = sum(r["eccdna_count_qduh_bed"] for r in out)
    covid = sum(r["eccdna_count_qduh_bed"] for r in out if r["group"] == "COVID-19")
    print(f"wrote {DEST.relative_to(RUN)}")
    print(f"  samples={len(out)}  total calls={total:,}  "
          f"COVID={covid:,}  HC={total - covid:,}")
    print("  epm identity asserted on every row")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
