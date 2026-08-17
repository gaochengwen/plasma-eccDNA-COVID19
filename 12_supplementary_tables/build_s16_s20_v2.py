#!/usr/bin/env python3
"""Rebuild Supplementary Tables S15, S16 and S20 on the v2 results.

These sheets hold several sub-tables stacked in one column. The target-specific
block in S15 and the validation-target blocks in S16 are restricted to CTNNA2,
SDK1 and TCF7L1, the three loci taken forward for wet-lab validation; S20c
reflects the updated candidate assignments. The sheets
are rebuilt rather than patched in place so the data regions have exact lengths.

Formatting is taken from the row being replaced so the rebuilt sheet keeps the
workbook's look.
"""

from __future__ import annotations

import csv
from copy import copy
from pathlib import Path

from openpyxl import load_workbook

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
T = RUN / "tables"
BOOK = T / "supplementary" / "Supplementary_Tables_v2.xlsx"

TARGET_ORDER = ["CTNNA2", "SDK1", "TCF7L1"]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(v):
    try:
        return float(v) if "." in str(v) or "e" in str(v).lower() else int(v)
    except (TypeError, ValueError):
        return v


def unmerge(ws, first_row: int) -> None:
    """Drop merged ranges below `first_row` so cells become writable."""
    for rng in list(ws.merged_cells.ranges):
        if rng.max_row >= first_row:
            ws.unmerge_cells(str(rng))


def blocks(ws, prefix: str) -> list[tuple[int, str]]:
    out = []
    for r in range(1, ws.max_row + 1):
        v = ws.cell(row=r, column=1).value
        if isinstance(v, str) and v.startswith(prefix) and "." in v[:5]:
            out.append((r, v))
    return out


def write_block(ws, start_row: int, header: list[str], rows: list[list],
                style_row: int) -> int:
    """Write a header + data block, copying cell style from `style_row`."""
    for c, name in enumerate(header, start=1):
        cell = ws.cell(row=start_row, column=c)
        cell.value = name
        cell._style = copy(ws.cell(row=style_row, column=c)._style)
    for i, data in enumerate(rows, start=1):
        for c, value in enumerate(data, start=1):
            cell = ws.cell(row=start_row + i, column=c)
            cell.value = value
            cell._style = copy(ws.cell(row=style_row + 1, column=c)._style)
    return start_row + len(rows) + 1


def rebuild_s15(wb) -> int:
    """Replace S15c with results for the three wet-lab targets only."""
    ws = wb["Table_S15"]
    section_start = next(
        row for row, value in blocks(ws, "S15") if value.startswith("S15c")
    )
    section_end = next(
        row for row, value in blocks(ws, "S15") if value.startswith("S15d")
    )
    rows = [
        row
        for row in read(T / "robustness" / "eccgene_candidate_stability.tsv")
        if row["gene"] in TARGET_ORDER
    ]
    callset_order = {
        "circlemap_methods": 0,
        "circlemap_high_support_masked": 1,
        "consensus_t10": 2,
    }
    definition_order = {"junction": 0, "interval": 1, "midpoint": 2}
    target_order = {target: index for index, target in enumerate(TARGET_ORDER)}
    rows.sort(
        key=lambda row: (
            callset_order.get(row["callset"], 99),
            definition_order.get(row["definition"], 99),
            target_order.get(row["gene"], 99),
        )
    )
    detection_cache: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    for row in rows:
        key = (row["callset"], row["definition"])
        if key not in detection_cache:
            detection_path = (
                RUN
                / "v2tree"
                / "Circle_finder"
                / "downstream"
                / row["callset"]
                / "eccGene"
                / "results"
                / f"{row['definition']}_detection_fisher.tsv"
            )
            detection_cache[key] = {
                item["Gene"]: item for item in read(detection_path)
            }
        detection = detection_cache[key][row["gene"]]
        row["detection_odds_ratio_COVID_vs_HC"] = detection[
            "odds_ratio_conditional_COVID_vs_HC"
        ]
        row["fisher_detection_p"] = detection["fisher_p"]
        row["fisher_detection_q_BH"] = detection["fisher_q_BH"]
        row["detection_direction"] = detection["direction"]
    header = list(rows[0])
    values = [[num(row[name]) for name in header] for row in rows]

    unmerge(ws, section_start)
    for row in range(section_start, section_end):
        for column in range(1, ws.max_column + 1):
            ws.cell(row=row, column=column).value = None
    ws.cell(row=section_start, column=1).value = (
        "S15c. Wet-lab-target eccGene stability for CTNNA2, SDK1 and TCF7L1"
    )
    write_block(
        ws,
        section_start + 1,
        header,
        values,
        style_row=section_start + 1,
    )
    return len(values)


def rebuild_s16(wb) -> tuple[int, int]:
    ws = wb["Table_S16"]
    caps = dict(blocks(ws, "S16"))
    cap_a = next(v for r, v in caps.items() if v.startswith("S16a"))
    cap_b = next(v for r, v in caps.items() if v.startswith("S16b"))

    support = [
        r for r in read(T / "robustness" /
                        "validation_target_support_summary.tsv")
        if r["target"] in TARGET_ORDER
    ]
    order = {t: i for i, t in enumerate(TARGET_ORDER)}
    support.sort(key=lambda r: (r["callset"], int(r["tolerance_bp"]),
                                order.get(r["target"], 99),
                                r["group"] != "COVID-19"))
    hdr_a = ["callset", "tolerance_bp", "target", "group", "n_samples",
             "supported_samples", "support_fraction", "matched_call_count"]
    rows_a = [[num(r[h]) for h in hdr_a] for r in support]

    mask = [
        r for r in read(T / "robustness" / "validation_target_mask_audit.tsv")
        if r["target"] in TARGET_ORDER
    ]
    mask.sort(key=lambda r: (order.get(r["target"], 99),
                             r["breakpoint_side"] != "start"))
    hdr_b = list(mask[0])
    rows_b = [[num(r[h]) for h in hdr_b] for r in mask]

    unmerge(ws, 3)
    for r in range(3, ws.max_row + 1):
        for c in range(1, ws.max_column + 1):
            ws.cell(row=r, column=c).value = None

    ws.cell(row=3, column=1).value = cap_a
    nxt = write_block(ws, 4, hdr_a, rows_a, style_row=4)
    ws.cell(row=nxt + 1, column=1).value = cap_b
    write_block(ws, nxt + 2, hdr_b, rows_b, style_row=4)
    return len(rows_a), len(rows_b)


def rebuild_s20(wb) -> int:
    ws = wb["Table_S20"]
    rank = read(T / "BCL3" / "chromatin_display_candidate_ranking.tsv")
    caps = blocks(ws, "S20")
    cap_c = next(v for _, v in caps if v.startswith("S20c"))
    start_c = next(r for r, v in caps if v.startswith("S20c"))
    start_d = next(r for r, v in caps if v.startswith("S20d"))

    hdr = [ws.cell(row=start_c + 1, column=c).value
           for c in range(1, ws.max_column + 1)]
    hdr = [h for h in hdr if h]
    rows = []
    for item in rank:
        rows.append([num(item.get(str(h), "")) for h in hdr])

    unmerge(ws, start_c)
    for r in range(start_c, start_d):
        for c in range(1, ws.max_column + 1):
            ws.cell(row=r, column=c).value = None
    ws.cell(row=start_c, column=1).value = cap_c
    write_block(ws, start_c + 1, hdr, rows, style_row=start_c + 1)

    title2 = ws.cell(row=2, column=1).value or ""
    ws.cell(row=2, column=1).value = str(title2).replace("281", str(len(rank)))
    return len(rows)


def main() -> int:
    wb = load_workbook(BOOK)
    s15 = rebuild_s15(wb)
    a, b = rebuild_s16(wb)
    n = rebuild_s20(wb)
    wb.save(BOOK)
    print(f"Table_S15: {s15} wet-lab-target eccGene rows")
    print(f"Table_S16: {a} support rows (3 wet-lab targets), {b} mask-audit rows")
    print(f"Table_S20: {n} candidate-category assignments")
    print(f"saved {BOOK.relative_to(RUN)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
