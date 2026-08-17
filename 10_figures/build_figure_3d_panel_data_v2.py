#!/usr/bin/env python3
"""Emit the v2 panel data for Figure 3d, the COVID-19-specific exact intervals.

Figure 3d in the submitted artwork shows 27 intervals: the v1 result. Under the
methods-consistent call set the same rule -- identical chromosome, start and end
coordinates in at least 10 of the 39 COVID-19 samples and in no HC sample --
yields 16, which is what the Results text and the Figure 3 legend already state.

The v2 set is a strict subset of the v1 set, but the panel is a set of
concentric rings, one per interval, so eleven rings cannot simply be deleted:
the radial layout is a function of the number of intervals.  This writes the
membership matrix and the per-interval recurrence in the order the panel draws
them, so the ring layout can be rebuilt from it.

Row order matches the Figure 3 legend: decreasing recurrence, then genomic
coordinate.  Only COVID-19 samples are emitted, because no interval is detected
in a control.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
MATRIX = RUN / "tables" / "recurrent_intervals" / "covid_specific_exact.matrix.min10.tsv"
OUT = RUN / "figures" / "illustrator_panel_data"

MIN_SAMPLES = 10
N_COVID = 39


def chrom_key(name: str) -> tuple[int, int]:
    chrom, _, span = name.partition(":")
    start = int(span.split("_")[0])
    tag = chrom[3:]
    return (int(tag) if tag.isdigit() else 100 + ord(tag[0]), start)


def main() -> int:
    with MATRIX.open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    samples = rows[0][1:]
    if len(samples) != N_COVID:
        raise SystemExit(f"expected {N_COVID} COVID-19 columns, found {len(samples)}")

    intervals = []
    for row in rows[1:]:
        hits = [int(v) for v in row[1:]]
        intervals.append({"interval": row[0], "recurrence": sum(hits), "hits": hits})

    for item in intervals:
        if item["recurrence"] < MIN_SAMPLES:
            raise SystemExit(f"{item['interval']} has recurrence {item['recurrence']}")
    intervals.sort(key=lambda d: (-d["recurrence"], chrom_key(d["interval"])))
    print(f"{len(intervals)} intervals, recurrence "
          f"{intervals[-1]['recurrence']}-{intervals[0]['recurrence']}")
    if len(intervals) != 16:
        raise SystemExit(f"expected 16 v2 intervals, found {len(intervals)}")

    OUT.mkdir(parents=True, exist_ok=True)
    label = lambda s: re.sub(r"_", "–", s)          # chr2:79885063–79885444
    matrix_path = OUT / "Figure_3d_recurrent_intervals.tsv"
    with matrix_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["interval", "label", "recurrence_of_39_COVID"] + samples)
        for item in intervals:
            w.writerow([item["interval"], label(item["interval"]),
                        item["recurrence"]] + item["hits"])

    per_sample = {s: sum(item["hits"][i] for item in intervals)
                  for i, s in enumerate(samples)}
    summary = {
        "rule": ("identical chromosome, start and end coordinates in at least "
                 f"{MIN_SAMPLES} of {N_COVID} COVID-19 samples and in no HC sample"),
        "call_set": "circlemap_methods (v2)",
        "n_intervals_v2": len(intervals),
        "n_intervals_v1_in_submitted_artwork": 27,
        "row_order": "decreasing recurrence, then genomic coordinate",
        "recurrence_max": intervals[0]["recurrence"],
        "recurrence_min": intervals[-1]["recurrence"],
        "n_covid_samples_with_any": sum(1 for v in per_sample.values() if v),
        "n_hc_samples_with_any": 0,
        "intervals": [{"label": label(i["interval"]), "recurrence": i["recurrence"]}
                      for i in intervals],
        "per_sample_interval_count": per_sample,
    }
    (OUT / "Figure_3d_panel_summary.json").write_text(json.dumps(summary, indent=1))

    for item in intervals:
        print(f"  {label(item['interval']):34s} {item['recurrence']:2d}/39")
    print(f"wrote {matrix_path.name} and Figure_3d_panel_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
