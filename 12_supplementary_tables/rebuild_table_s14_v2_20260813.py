#!/usr/bin/env python3
"""Rebuild Supplementary Table S14 on the v2 call set.

The carried-over sheet holds the v1 significant set: 3,402 genes, 3,094 of them
COVID-up. The manuscript, Figure 3a and the locked eccGene tables all report the
v2 set: 3,487 genes, 3,189 COVID-up. Only ~2,306 gene symbols are common to the
two lists, so this is a full rebuild rather than a value refresh.

Selection reproduces the Methods definition exactly: BH q < 0.05 on the
two-sided Wilcoxon rank-sum test across the 26,087 tested genes, together with
an absolute mean-EA log2 fold change >= 1. Per-sample values are the primary
junction EA matrix. Genes are listed alphabetically, as in the previous sheet.

Table_S14 is the only sheet whose row count changes, and nothing follows its
data region, so rows are simply written past the old extent and any surplus
rows from the shorter v1 list are cleared.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import openpyxl

RUN = Path(__file__).resolve().parents[2]
WB = RUN / "tables" / "supplementary" / "Supplementary_Tables_v2.xlsx"
ECC = RUN / "tables" / "eccGene"

TITLE_ROW, HEADER_ROW, FIRST_DATA_ROW = 1, 2, 3
Q_MAX, LFC_MIN = 0.05, 1.0


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    stats = {r["Gene"]: r for r in read(ECC / "junction_abundance_wilcoxon.tsv")}
    sig = sorted(
        g for g, r in stats.items()
        if float(r["wilcoxon_q_BH"]) < Q_MAX
        and abs(float(r["log2FC_mean_EA_COVID_vs_HC"])) >= LFC_MIN
    )
    up = sum(1 for g in sig if float(stats[g]["log2FC_mean_EA_COVID_vs_HC"]) > 0)
    print(f"tested genes {len(stats):,}; significant {len(sig):,}; COVID-up {up:,}")
    if (len(sig), up) != (3487, 3189):
        raise SystemExit(f"expected 3,487/3,189, got {len(sig)}/{up}")

    ea: dict[str, dict[str, str]] = {}
    with (ECC / "matrices" / "junction_EA.tsv").open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ea[row["Gene"]] = row

    wb = openpyxl.load_workbook(WB)
    ws = wb["Table_S14"]
    headers = [ws.cell(row=HEADER_ROW, column=c).value
               for c in range(1, ws.max_column + 1)]
    samples = [str(h) for h in headers[1:-4]]
    tail = [str(h) for h in headers[-4:]]
    if tail != ["log2FC", "Cliffs_delta", "p_wilcoxon", "adj_wilcoxon_BH"]:
        raise SystemExit(f"unexpected trailing columns: {tail}")
    missing = [s for s in samples if s not in next(iter(ea.values()))]
    if missing:
        raise SystemExit(f"samples absent from the EA matrix: {missing}")
    print(f"columns: 1 gene + {len(samples)} samples + {len(tail)} statistics")

    old_last = ws.max_row
    for i, gene in enumerate(sig):
        r = FIRST_DATA_ROW + i
        row = ea[gene]
        st = stats[gene]
        ws.cell(row=r, column=1).value = gene
        for c, s in enumerate(samples, start=2):
            ws.cell(row=r, column=c).value = round(float(row[s]), 4)
        n = len(samples) + 2
        ws.cell(row=r, column=n).value = float(st["log2FC_mean_EA_COVID_vs_HC"])
        ws.cell(row=r, column=n + 1).value = float(st["cliffs_delta_COVID_vs_HC"])
        ws.cell(row=r, column=n + 2).value = float(st["wilcoxon_p"])
        ws.cell(row=r, column=n + 3).value = float(st["wilcoxon_q_BH"])

    new_last = FIRST_DATA_ROW + len(sig) - 1
    cleared = 0
    for r in range(new_last + 1, old_last + 1):
        for c in range(1, ws.max_column + 1):
            if ws.cell(row=r, column=c).value is not None:
                ws.cell(row=r, column=c).value = None
                cleared += 1
    print(f"rows {FIRST_DATA_ROW}-{new_last} written; {cleared} surplus cells cleared")

    wb.save(WB)
    print(f"saved {WB}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
