#!/usr/bin/env python3
"""Write the v2 results back into the submission workbook, format preserved.

The deliverable is a workbook in the same shape as
Submit/Supplementary_Tables_revised.xlsx: 21 sheets, each with its title row,
its header row and its column formatting. So the workbook is opened, only the
data regions of the affected sheets are rewritten, and everything else --
titles, column widths, number formats, untouched sheets -- is left alone.

Sheets rewritten here:

  Table_S3   cohort comparisons, all eight FDR families
  Table_S10  chromosome-normalized eccDNA per Mb, per sample
  Table_S12  eccDNA distribution across gene elements, per sample
  Table_S13  gene-feature observed/expected, per sample

Table_S11 (repeat classes) is deliberately NOT rewritten: its source analysis
runs on BAM reads, which are identical in v1 and v2, so the call-set migration
does not reach it.

Rows are written in the sheet's own column order and matched by header name.
Columns that no longer belong to the v2 data model are removed explicitly;
unmapped cells in retained columns are left blank rather than carrying stale
values from the previous revision.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

from openpyxl import load_workbook

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
T = RUN / "tables"
SRC = PROJECT / "Submit" / "Supplementary_Tables_revised.xlsx"
OUT = RUN / "tables" / "supplementary" / "Supplementary_Tables_v2.xlsx"

TITLE_ROW, HEADER_ROW, FIRST_DATA_ROW = 1, 2, 3
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ELEMENTS = ["3UTR", "5UTR", "CpG_islands", "Gene2KbD", "Gene2KbU", "exon", "intron"]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def rewrite(ws, rows: list[dict], mapping: dict[str, str]) -> int:
    """Replace the data region, matching sheet columns by header name."""
    headers = [ws.cell(row=HEADER_ROW, column=c).value
               for c in range(1, ws.max_column + 1)]
    for r in range(FIRST_DATA_ROW, ws.max_row + 1):
        for c in range(1, ws.max_column + 1):
            ws.cell(row=r, column=c).value = None
    for i, item in enumerate(rows):
        for c, head in enumerate(headers, start=1):
            key = mapping.get(str(head))
            if key is None or key not in item:
                continue
            ws.cell(row=FIRST_DATA_ROW + i, column=c).value = num(item[key])
    return len(rows)


def delete_columns_by_header(ws, names: set[str]) -> None:
    """Delete obsolete columns by their row-2 header, from right to left."""
    for column in range(ws.max_column, 0, -1):
        if str(ws.cell(row=HEADER_ROW, column=column).value) in names:
            ws.delete_cols(column)


def table_s3(wb) -> int:
    rows = read(T / "supplementary" / "Table_S3_v2.tsv")
    for row in rows:
        row["COVID-19 IQR"] = (
            float(row["COVID-19 Q3"]) - float(row["COVID-19 Q1"])
        )
        row["HC IQR"] = float(row["HC Q3"]) - float(row["HC Q1"])
    ws = wb["Table_S3"]
    delete_columns_by_header(
        ws,
        {"Rank-biserial correlation", "Mann-Whitney U", "Z statistic"},
    )
    mapping = {h: h for h in rows[0]}
    return rewrite(ws, rows, mapping)


def table_s10(wb) -> int:
    src = read(T / "genomic_distribution" / "chromosome_distribution.tsv")
    by: dict[tuple[str, str], dict[str, str]] = {}
    for r in src:
        by.setdefault((r["sample_id"], r["group"]), {})[r["chromosome"]] = \
            r["normalized_fraction_per_mb"]
    rows = []
    for (sample, group), values in sorted(by.items(),
                                          key=lambda kv: (kv[0][1] != "COVID-19",
                                                          kv[0][0])):
        item = {"Sample ID": sample, "Group": group}
        item.update({c: values.get(c, "") for c in CHROMS})
        rows.append(item)
    mapping = {"Sample ID": "Sample ID", "Group": "Group"}
    mapping.update({c: c for c in CHROMS})
    return rewrite(wb["Table_S10"], rows, mapping)


def gene_feature_rows() -> list[dict]:
    src = read(T / "genomic_distribution" / "gene_feature_enrichment_scores.tsv")
    order = {e: i for i, e in enumerate(ELEMENTS)}
    src.sort(key=lambda r: (r["group"] != "COVID-19", r["sample_id"],
                            order.get(r["element"], 99)))
    return src


def table_s12(wb) -> int:
    rows = gene_feature_rows()
    ws = wb["Table_S12"]
    delete_columns_by_header(
        ws,
        {"Coverage Ratio", "Nomalized EccDNA ratio"},
    )
    mapping = {"Sample ID": "sample_id", "Group": "group", "Element": "element",
               "Element eccDNA Count": "observed_overlap_o",
               "Total eccDNA": "total_eccdna_n",
               "EccDNA element ratio": "observed_fraction"}
    return rewrite(ws, rows, mapping)


def table_s13(wb) -> int:
    rows = gene_feature_rows()
    ws = wb["Table_S13"]
    delete_columns_by_header(
        ws,
        {"Background Total (B)", "Background Overlap (E)"},
    )
    for column in range(1, ws.max_column + 1):
        if ws.cell(row=HEADER_ROW, column=column).value == "Enrichment_score":
            ws.cell(row=HEADER_ROW, column=column).value = "Observed/expected ratio"
    mapping = {"Sample": "sample_id", "Group": "group", "Element": "element",
               "Total eccDNA (N)": "total_eccdna_n",
               "Observed Overlap (O)": "observed_overlap_o",
               "Observed Fraction": "observed_fraction",
               "Expected Fraction": "expected_fraction",
               "Observed/expected ratio": "oe_ratio"}
    return rewrite(ws, rows, mapping)


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SRC, OUT)
    wb = load_workbook(OUT)

    written = {"Table_S3": table_s3(wb), "Table_S10": table_s10(wb),
               "Table_S12": table_s12(wb), "Table_S13": table_s13(wb)}

    # Record the provenance of the sheet that is deliberately untouched.
    ws = wb["Table_S11"]
    note = ws.cell(row=TITLE_ROW, column=1).value or ""
    if "unchanged" not in str(note):
        ws.cell(row=TITLE_ROW, column=1).value = (
            f"{note} Values are unchanged from the previous revision: the "
            f"repeat-class analysis is computed from BAM reads, which are "
            f"identical under both call sets.")

    wb.save(OUT)
    print(f"wrote {OUT.relative_to(RUN)}")
    for sheet, n in written.items():
        print(f"  {sheet:10s} {n:>5d} data rows rewritten")
    print("  Table_S11      unchanged (BAM-level analysis), provenance noted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
