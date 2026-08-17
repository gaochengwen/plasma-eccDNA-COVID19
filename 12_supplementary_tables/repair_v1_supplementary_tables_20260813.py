#!/usr/bin/env python3
"""Rewrite the supplementary sheets that were never migrated to the v2 call set.

Background
----------
`build_supplementary_v2.py` created the v2 workbook by copying the v1
submission workbook and rewriting only Tables S3, S10, S12 and S13. Later
rounds added S15, S16 and S20. Everything else was carried over untouched, so a
group of sheets still holds numbers computed on the unfiltered archive call set
while the manuscript that cites them reports the methods-consistent v2 results.

An audit comparing every sheet cell-by-cell against the locked tables under
`tables/` found the following blocks still on v1:

  S2a          per-sample eccDNA count and EPM (sums to the v1 16,857,160)
  S5b-e        fragment-length peaks 196/365/571/760, spacings 169/206
  S14          the v1 3,402-gene significant list (v2 has 3,487)
  S17b-g       all chromatin eligibility, enrichment and abundance statistics
  S18a,c,e,f   downsampling, burden-window, HC3 and Spearman burden controls
  S19b,d       cross-cell-type and mappability-restricted enrichment

S4 is handled by `rebuild_age_matched_pairs_v2.py`, which must run after this
script because the age-matching inputs are read from Table S2a.

Sheets deliberately left alone: S11 (repeat classes, computed from BAM reads
that are identical under both call sets), S17a/S19a/S19c (peak-set provenance
and mappability, independent of the eccDNA call set), and every sheet already
verified against its locked source.

Method
------
Each block's data region is located from its header row, cleared, and rewritten
from the locked TSV, matching sheet columns to source columns by header name.
Titles, notes, headers, column widths and number formats are left untouched.
Only S14 changes row count, so rows are inserted there before writing.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import openpyxl

RUN = Path(__file__).resolve().parents[2]
WB = RUN / "tables" / "supplementary" / "Supplementary_Tables_v2.xlsx"
T = RUN / "tables"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(v):
    """Numbers as numbers, everything else verbatim; blanks stay blank."""
    if v is None or v == "":
        return None
    if v in ("NA", "nan", "NaN"):
        return "NA"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return v
    return int(f) if f.is_integer() and abs(f) < 2 ** 53 and "." not in v and "e" not in v.lower() else f


def block_rows(ws, header_row: int) -> int:
    """Number of data rows under `header_row`, stopping at a blank or a note."""
    n = 0
    r = header_row + 1
    while r <= ws.max_row:
        vals = [ws.cell(row=r, column=c).value for c in range(1, ws.max_column + 1)]
        filled = [v for v in vals if v is not None]
        if not filled:
            break
        if len(filled) == 1 and isinstance(vals[0], str):
            break
        n += 1
        r += 1
    return n


def rewrite(ws, header_row: int, rows: list[dict], *, label: str,
            rename: dict[str, str] | None = None,
            derive=None) -> None:
    """Clear the block's data region and write `rows` matched by header name."""
    rename = rename or {}
    ncols = ws.max_column
    headers = [ws.cell(row=header_row, column=c).value for c in range(1, ncols + 1)]
    keys = [rename.get(str(h), str(h)) for h in headers]

    have = block_rows(ws, header_row)
    if have != len(rows):
        raise SystemExit(
            f"{label}: sheet has {have} data rows but source has {len(rows)}; "
            "row-count changes must be handled explicitly")

    for i in range(have):
        for c in range(1, ncols + 1):
            ws.cell(row=header_row + 1 + i, column=c).value = None

    written = 0
    for i, item in enumerate(rows):
        if derive is not None:
            item = derive(item)
        for c, key in enumerate(keys, start=1):
            if key in item:
                ws.cell(row=header_row + 1 + i, column=c).value = num(item[key])
                written += 1
    missing = [k for k in keys if k not in (rows[0] if not derive else derive(rows[0]))]
    print(f"  {label:6s} rows={len(rows):5d} cells={written:6d}"
          + (f"  unmapped columns: {missing}" if missing else ""))


def main() -> None:
    wb = openpyxl.load_workbook(WB)
    print(f"repairing {WB.name}")

    # ---------- S2a: per-sample count and EPM ----------------------------
    rca = {r["sample_id"]: r for r in read(T / "RCA" / "RCA_sample_metrics.tsv")}
    ws = wb["Table_S2"]
    hdr = [ws.cell(row=3, column=c).value for c in range(1, ws.max_column + 1)]
    ci = {str(h): i + 1 for i, h in enumerate(hdr)}
    n = 0
    for r in range(4, 82):
        sid = ws.cell(row=r, column=1).value
        if not sid:
            continue
        s = rca[sid]
        ws.cell(row=r, column=ci["EccDNA Count"]).value = int(float(s["eccdna_count"]))
        ws.cell(row=r, column=ci["EPM"]).value = round(float(s["epm"]), 4)
        n += 1
    print(f"  S2a    rows={n:5d} (EccDNA Count, EPM; Mapped Reads/Read Pairs/rates already v2)")

    # ---------- S5: fragment-length peaks --------------------------------
    F = T / "fragment_length"
    ws = wb["Table_S5"]
    # the sheet calls the prominence column primary_peak_prominence; the source
    # calls it primary_prominence. It was left blank in the carried-over sheet.
    rewrite(ws, 27, read(F / "objective_peak_identification_and_stability.tsv"),
            label="S5b", rename={"primary_peak_prominence": "primary_prominence"})
    rewrite(ws, 36, read(F / "adjacent_peak_spacing.tsv"), label="S5c")
    rewrite(ws, 43, read(F / "peak_window_group_comparisons.tsv"), label="S5d")
    # S5e was never copied off QDUH with the other v2 fragment-length tables;
    # retrieved from the same run (the group table beside it is byte-identical
    # to the locked local copy).
    rewrite(ws, 51, read(F / "peak_window_sample_proportions.tsv"), label="S5e")

    # ---------- S17: chromatin -------------------------------------------
    C = T / "chromatin"
    ws = wb["Table_S17"]
    rewrite(ws, 19, read(C / "eccdna_sample_qc.tsv"), label="S17b")
    rewrite(ws, 102, read(C / "p0_within_group_enrichment.tsv"), label="S17c")
    rewrite(ws, 167, read(C / "relative_enrichment_group_stats.tsv"), label="S17d")
    rewrite(ws, 202, read(C / "absolute_burden_group_stats.tsv"), label="S17e")
    rewrite(ws, 236, read(C / "observed_fraction_group_stats.tsv"), label="S17f")
    rewrite(ws, 271, read(C / "sample_relative_enrichment.tsv"), label="S17g")

    # ---------- S18: burden controls -------------------------------------
    ws = wb["Table_S18"]
    rewrite(ws, 5, read(C / "p0_downsampled_group_stats.tsv"), label="S18a")
    rewrite(ws, 134, read(C / "p0_burden_matched_subset.tsv"), label="S18c")
    rewrite(ws, 251, read(C / "p0_covariate_adjusted_hc3.tsv"), label="S18e")
    rewrite(ws, 736, read(C / "p0_burden_association.tsv"), label="S18f")

    # ---------- S19: cross-cell-type and mappability ---------------------
    ws = wb["Table_S19"]
    rewrite(ws, 43, read(C / "ext_multicell_group_stats.tsv"), label="S19b")
    rewrite(ws, 129, read(C / "ext_mappability_group_stats.tsv"), label="S19d")

    wb.save(WB)
    print(f"saved {WB}")


if __name__ == "__main__":
    sys.exit(main())
