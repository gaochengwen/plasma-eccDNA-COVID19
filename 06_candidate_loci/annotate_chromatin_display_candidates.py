#!/usr/bin/env python3
"""Derive the reproducible Figure 5 chromatin-display subset.

The primary eccGene eligibility and category rankings are not changed here.
This script applies a separate, explicitly recorded feasibility filter for a
locus-level *reference-chromatin* display: a strict-eligible gene must overlap
at least one H3K27ac peak in both of the two reference contexts already used
in Figure 5 (SARS-CoV-2-infected A549-ACE2 and CD14+ monocytes).  The test
window is the longest canonical hg38 refGene transcript plus 1.5 kb on each
side.  Within each category, eligible display loci retain the pre-existing
lexicographic order; track appearance and eccDNA coverage shape are never used
as selection criteria.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


A549_H3K27AC = "A549_GSE179184_H3K27ac_official"
CD14_H3K27AC = "CD14_H3K27ac"
CATEGORY_RANK_COLUMNS = {
    "Humoral immune response": "Humoral_rank",
    "Endothelial/coagulation": "Endothelial_coagulation_rank",
}


def read_tsv(path: Path) -> List[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def longest_refgene_records(path: Path, wanted: set[str]) -> Dict[str, Tuple[str, int, int]]:
    """Return the longest canonical-chromosome transcript per requested gene."""
    best: Dict[str, Tuple[str, int, int]] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16:
                continue
            gene, chrom = parts[12], parts[2]
            if gene not in wanted or chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}:
                continue
            start, end = int(parts[4]), int(parts[5])
            current = best.get(gene)
            if current is None or end - start > current[2] - current[1]:
                best[gene] = (chrom, start, end)
    missing = wanted - set(best)
    if missing:
        raise RuntimeError(f"canonical refGene record missing for: {', '.join(sorted(missing))}")
    return best


def read_bed3(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            intervals[chrom].append((int(start), int(end)))
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def overlap_count(intervals: Iterable[Tuple[int, int]], start: int, end: int) -> int:
    return sum(interval_end > start and interval_start < end for interval_start, interval_end in intervals)


def rank_value(text: str) -> float:
    return float(text) if text not in {"", "NA", "NaN", "nan"} else float("inf")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", required=True, type=Path)
    parser.add_argument("--refgene", required=True, type=Path)
    parser.add_argument("--peaks-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--criteria-out", required=True, type=Path)
    parser.add_argument("--flank-bp", type=int, default=1500)
    parser.add_argument(
        "--expect",
        default="",
        help="Optional pin, e.g. 'Humoral immune response=TFE3;Endothelial/coagulation=MCAM'.",
    )
    args = parser.parse_args()

    if args.flank_bp < 0:
        raise ValueError("--flank-bp must be non-negative")
    rows = read_tsv(args.ranking)
    strict_rows = [row for row in rows if row["Strict_selection_eligible"] == "Yes"]
    records = longest_refgene_records(args.refgene, {row["Gene"] for row in strict_rows})
    peaks = {
        A549_H3K27AC: read_bed3(
            args.peaks_dir / f"{A549_H3K27AC}.liftover.hg38.canonical.filtered.merged.bed"
        ),
        CD14_H3K27AC: read_bed3(
            args.peaks_dir / f"{CD14_H3K27AC}.liftover.hg38.canonical.filtered.merged.bed"
        ),
    }

    for row in rows:
        if row["Strict_selection_eligible"] != "Yes":
            row["Display_window_hg38"] = "NA"
            row["A549_H3K27ac_peak_count"] = "NA"
            row["CD14_H3K27ac_peak_count"] = "NA"
            row["Dual_reference_H3K27ac_eligible"] = "No"
            row["Chromatin_display_rank"] = "NA"
            row["Selected_for_Figure_5_display"] = "No"
            continue
        chrom, start, end = records[row["Gene"]]
        window_start = max(0, start - args.flank_bp)
        window_end = end + args.flank_bp
        a549_count = overlap_count(peaks[A549_H3K27AC].get(chrom, []), window_start, window_end)
        cd14_count = overlap_count(peaks[CD14_H3K27AC].get(chrom, []), window_start, window_end)
        row["Display_window_hg38"] = f"{chrom}:{window_start}-{window_end}"
        row["A549_H3K27ac_peak_count"] = str(a549_count)
        row["CD14_H3K27ac_peak_count"] = str(cd14_count)
        row["Dual_reference_H3K27ac_eligible"] = "Yes" if a549_count and cd14_count else "No"
        row["Chromatin_display_rank"] = "NA"
        row["Selected_for_Figure_5_display"] = "No"

    selected: dict[str, str] = {}
    for category, rank_column in CATEGORY_RANK_COLUMNS.items():
        category_rows = [
            row
            for row in rows
            if row["Category_for_row"] == category
            and row["Dual_reference_H3K27ac_eligible"] == "Yes"
        ]
        category_rows.sort(key=lambda row: (rank_value(row[rank_column]), row["Gene"]))
        if not category_rows:
            raise RuntimeError(f"no dual-reference display-eligible candidate in {category}")
        for rank, row in enumerate(category_rows, 1):
            row["Chromatin_display_rank"] = str(rank)
        category_rows[0]["Selected_for_Figure_5_display"] = "Yes"
        selected[category] = category_rows[0]["Gene"]

    # v1 hard-coded the display pair it happened to produce
    # ({"Humoral immune response": "HLA-E", "Endothelial/coagulation": "MCAM"})
    # and raised on anything else. That guard is why the v2 rerun stopped here:
    # on the methods call set the humoral winner is TFE3, not HLA-E.
    #
    # The real invariant is structural -- exactly one display locus per defined
    # category -- so that is what is enforced. The specific pair is reported,
    # and may be pinned with --expect when a run must reproduce a known result.
    if set(selected) != set(CATEGORY_RANK_COLUMNS):
        raise RuntimeError(
            f"expected one display locus per category {sorted(CATEGORY_RANK_COLUMNS)}, "
            f"got {selected}"
        )
    print(
        "display selections: "
        + ", ".join(f"{k}={v}" for k, v in sorted(selected.items())),
        flush=True,
    )
    if args.expect:
        pinned = dict(
            item.split("=", 1) for item in args.expect.split(";") if item.strip()
        )
        if selected != pinned:
            raise RuntimeError(
                f"display selections {selected} do not match --expect {pinned}"
            )

    extra_fields = [
        "Display_window_hg38",
        "A549_H3K27ac_peak_count",
        "CD14_H3K27ac_peak_count",
        "Dual_reference_H3K27ac_eligible",
        "Chromatin_display_rank",
        "Selected_for_Figure_5_display",
    ]
    write_tsv(args.out, rows, list(rows[0]) + extra_fields)
    criteria = [
        {
            "Stage": "8",
            "Criterion": "Reference-chromatin display feasibility",
            "Rule": (
                "Strict-eligible candidate overlaps at least one merged H3K27ac peak in both "
                "GSE179184 A549-ACE2 and CD14+ monocyte reference sets within the longest "
                f"canonical hg38 refGene transcript plus {args.flank_bp:,} bp on each side"
            ),
            "Used_for_ranking": "Display-subset eligibility only; not primary eccGene eligibility or ranking",
        },
        {
            "Stage": "9",
            "Criterion": "Display-subset selection",
            "Rule": (
                "Within each prespecified functional category, retain the existing lexicographic "
                "order and select the top display-eligible locus"
            ),
            "Used_for_ranking": "Figure 5d/e selection; coverage shape and visual inspection excluded",
        },
    ]
    write_tsv(args.criteria_out, criteria, ["Stage", "Criterion", "Rule", "Used_for_ranking"])
    print(
        "selected Figure 5 display loci: "
        + "; ".join(f"{category}={gene}" for category, gene in selected.items())
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
