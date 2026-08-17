#!/usr/bin/env python3
"""Add the objective Figure 5 locus-selection ranking as Supplementary Table S18.

The table is self-contained: it records the prespecified functional strata,
eligibility filters, lexicographic ranking order, and the category-specific
candidate rankings. It separately records the reference-chromatin feasibility
filter used only to select browser-display loci. Superseded illustrative loci
are intentionally omitted from the delivered supplementary table at the
authors' request.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


REVISE = Path(__file__).resolve().parent.parent
ANALYSIS = REVISE.parent / "analyse" / "BCL3" / "tables"
WORKBOOK = REVISE / "Supplementary_Tables_revised.xlsx"
BACKUP = REVISE / ".backup" / "Supplementary_Tables_pre_locus_selection.xlsx"
SHEET = "Table_S18"
OMIT_GENES = {"BCL3", "PROCR"}

TITLE = (
    "Supplementary Table S18. Category-based ranking and reference-chromatin "
    "display selection of candidate eccGene loci for Figure 5."
)

DISPLAY_COLUMNS = [
    ("Category_for_row", "Functional_category"),
    ("Category_rank", "Category_rank"),
    ("Gene", "Gene"),
    ("FDR", "BH_FDR"),
    ("log2FC", "log2FC_COVID_vs_HC"),
    ("Effect_size", "Cliffs_delta"),
    ("COVID_detected_n", "COVID_detected_n_of_39"),
    ("HC_detected_n", "HC_detected_n_of_39"),
    ("Prevalence_difference", "Detection_prevalence_difference"),
    ("Median_abundance_COVID", "Median_EA_COVID"),
    ("Median_abundance_HC", "Median_EA_HC"),
    ("Max_sample_contribution", "Maximum_COVID_sample_contribution"),
    ("quality_FDR", "Artifact_masked_BH_FDR"),
    ("quality_log2FC", "Artifact_masked_log2FC"),
    ("quality_effect_size", "Artifact_masked_Cliffs_delta"),
    ("quality_COVID_detected_n", "Artifact_masked_COVID_detected_n"),
    ("quality_HC_detected_n", "Artifact_masked_HC_detected_n"),
    ("artifact_recurrence_pass", "Independent_sample_support_pass"),
    ("blacklist_gap_pass", "Blacklist_and_gap_filter_pass"),
    ("Umap_k100_pass", "Umap_k100_filter_pass"),
    ("segmental_duplication_pass", "Segmental_duplication_filter_pass"),
    ("RepeatMasker_breakpoint_pass", "RepeatMasker_breakpoint_filter_pass"),
    ("quality_abundance_stable_all_definitions", "Artifact_masked_stable_all_assignments"),
    ("Strict_selection_eligible", "Eligible"),
    ("A549_H3K27ac_peak_count", "A549_H3K27ac_peaks_in_display_window"),
    ("CD14_H3K27ac_peak_count", "CD14_H3K27ac_peaks_in_display_window"),
    ("Dual_reference_H3K27ac_eligible", "Dual_reference_H3K27ac_display_eligible"),
    ("Chromatin_display_rank", "Display_subset_rank"),
    ("Selected_for_Figure_5_display", "Selected_for_Figure_5_d_or_e"),
    ("Exclusion_reason", "Exclusion_reason"),
]


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        values = line.split("\t")
        values += [""] * (len(header) - len(values))
        rows.append(dict(zip(header, values)))
    return header, rows


def value(text: str):
    if text in {"", "NA", "NaN", "nan"}:
        return None
    if text in {"Yes", "No", "Eligible"}:
        return text
    try:
        return float(text)
    except ValueError:
        return text


def style_section(ws, row: int, text: str) -> None:
    ws.cell(row=row, column=1, value=text)
    ws.cell(row=row, column=1).font = Font(bold=True)
    ws.cell(row=row, column=1).fill = PatternFill("solid", fgColor="E6E6ED")


def style_header(ws, row: int, columns: int) -> None:
    for cell in ws[row][:columns]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="253247")
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)


def main() -> int:
    criteria_header, criteria = read_tsv(ANALYSIS / "ranking_criteria.tsv")
    display_criteria_header, display_criteria = read_tsv(
        ANALYSIS / "chromatin_display_selection_criteria.tsv"
    )
    category_header, categories = read_tsv(ANALYSIS / "category_definitions.tsv")
    _, rankings = read_tsv(ANALYSIS / "chromatin_display_candidate_ranking.tsv")

    rankings = [row for row in rankings if row["Gene"] not in OMIT_GENES]
    for row in rankings:
        if row["Category_for_row"] == "Humoral immune response":
            row["Category_rank"] = row["Humoral_rank"]
        else:
            row["Category_rank"] = row["Endothelial_coagulation_rank"]

    category_order = {"Humoral immune response": 0, "Endothelial/coagulation": 1}

    def rank_value(text: str) -> float:
        return float(text) if text not in {"", "NA", "NaN", "nan"} else float("inf")

    rankings.sort(
        key=lambda row: (
            category_order[row["Category_for_row"]],
            rank_value(row["Category_rank"]),
            row["Gene"],
        )
    )

    if any(row["Gene"] in OMIT_GENES for row in rankings):
        raise RuntimeError("superseded locus escaped the export filter")
    selected = {
        row["Gene"]
        for row in rankings
        if row["Selected_for_Figure_5_display"] == "Yes"
        and row["Strict_selection_eligible"] == "Yes"
    }
    if selected != {"HLA-E", "MCAM"}:
        raise RuntimeError(f"unexpected selected loci: {sorted(selected)}")

    if not BACKUP.exists():
        BACKUP.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(WORKBOOK, BACKUP)

    wb = load_workbook(WORKBOOK)
    if SHEET in wb.sheetnames:
        del wb[SHEET]
    ws = wb.create_sheet(SHEET)

    ws.cell(row=1, column=1, value=TITLE)
    ws.cell(row=1, column=1).font = Font(bold=True, size=12)
    ws.cell(
        row=2,
        column=1,
        value=(
            "Primary eccGene eligibility and category ranking were completed independently of "
            "chromatin tracks. For Figure 5d/e only, strict-eligible candidates additionally "
            "had to overlap at least one H3K27ac peak in both predeclared reference contexts "
            "(GSE179184 A549-ACE2 and CD14+ monocytes), using the longest hg38 refGene "
            "transcript plus 1.5 kb on either side. The highest-ranked display-eligible locus "
            "in each category was selected before plotting. Visual signal intensity and coverage "
            "shape were not used. The endothelial/coagulation stratum was prespecified from "
            "Gene Ontology membership and was not itself a significantly enriched pathway in "
            "the locked enrichment analysis."
        ),
    )
    ws.cell(row=2, column=1).alignment = Alignment(wrap_text=True, vertical="top")

    row = 4
    style_section(ws, row, "S18a. Prespecified functional-category definitions")
    row += 1
    for col, name in enumerate(category_header, 1):
        ws.cell(row=row, column=col, value=name)
    style_header(ws, row, len(category_header))
    row += 1
    for record in categories:
        for col, name in enumerate(category_header, 1):
            ws.cell(row=row, column=col, value=value(record[name]))
        row += 1

    row += 1
    style_section(ws, row, "S18b. Eligibility filters and lexicographic ranking order")
    row += 1
    for col, name in enumerate(criteria_header, 1):
        ws.cell(row=row, column=col, value=name)
    style_header(ws, row, len(criteria_header))
    row += 1
    for record in criteria:
        for col, name in enumerate(criteria_header, 1):
            ws.cell(row=row, column=col, value=value(record[name]))
        row += 1
    for record in display_criteria:
        for col, name in enumerate(display_criteria_header, 1):
            ws.cell(row=row, column=col, value=value(record[name]))
        row += 1

    row += 1
    style_section(
        ws,
        row,
        "S18c. Category-specific candidate ranking and Figure 5 selection",
    )
    row += 1
    ranking_header_row = row
    for col, (_, label) in enumerate(DISPLAY_COLUMNS, 1):
        ws.cell(row=row, column=col, value=label)
    style_header(ws, row, len(DISPLAY_COLUMNS))
    row += 1
    for record in rankings:
        for col, (source, _) in enumerate(DISPLAY_COLUMNS, 1):
            ws.cell(row=row, column=col, value=value(record[source]))
        row += 1

    ws.freeze_panes = f"A{ranking_header_row + 1}"
    ws.auto_filter.ref = (
        f"A{ranking_header_row}:{get_column_letter(len(DISPLAY_COLUMNS))}{row - 1}"
    )
    ws.sheet_view.showGridLines = False
    ws.row_dimensions[1].height = 30
    ws.row_dimensions[2].height = 45
    for cells in ws.iter_rows():
        for cell in cells:
            cell.alignment = Alignment(
                vertical="top",
                wrap_text=cell.row <= ranking_header_row,
            )
    for col in range(1, len(DISPLAY_COLUMNS) + 1):
        label = ws.cell(ranking_header_row, col).value or ""
        width = min(max(13, len(str(label)) + 2), 34)
        ws.column_dimensions[get_column_letter(col)].width = width
    ws.column_dimensions["A"].width = 28
    ws.column_dimensions["C"].width = 18
    ws.column_dimensions[get_column_letter(len(DISPLAY_COLUMNS))].width = 42

    wb.save(WORKBOOK)
    print(
        f"wrote {SHEET}: {len(rankings)} category rows; selected "
        f"{', '.join(sorted(selected))}; omitted {', '.join(sorted(OMIT_GENES))}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
