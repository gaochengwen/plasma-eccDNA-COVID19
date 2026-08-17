#!/usr/bin/env python3
"""Add Supplementary Tables S11 and S12 to the supplementary workbook.

S11 collects everything needed to reproduce the corrected chromatin analysis:
liftOver QC, per-sample eccDNA eligibility, and the observed, expected, absolute
and relative quantities for all ten peak sets under all three localization
definitions, plus the within-cohort tests.

S12 collects the burden controls: downsampling and its agreement with the
full-sample estimates, the burden-matched subset, the HC3 regression models and
the burden-enrichment correlations.

Styling follows update_epm_statistics_revision.py so the new sheets match S1-S10.
"""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Sequence

from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

REVISE = Path(__file__).resolve().parent.parent
CHROMATIN = REVISE.parent / "analyse" / "chromatin"
TABLES_DIR = CHROMATIN / "tables"
WORKBOOK = REVISE / "Supplementary_Tables_revised.xlsx"

HEADER_FILL = PatternFill("solid", fgColor="D9EAF7")
SECTION_FILL = PatternFill("solid", fgColor="FCE4D6")


def read_tsv(name: str) -> List[Dict[str, str]]:
    path = TABLES_DIR / name
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_number(value: str):
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except ValueError:
        return value
    if number.is_integer() and abs(number) < 1e15 and "." not in value and "e" not in value.lower():
        return int(number)
    return number


def write_block(
    sheet,
    row: int,
    title: str,
    rows: Sequence[Dict[str, str]],
    columns: Sequence[str],
    note: str = "",
) -> int:
    """Write one titled block of rows; returns the next free row."""
    cell = sheet.cell(row=row, column=1, value=title)
    cell.font = Font(bold=True, size=11)
    cell.fill = SECTION_FILL
    cell.alignment = Alignment(horizontal="left", vertical="center")
    row += 1
    if note:
        note_cell = sheet.cell(row=row, column=1, value=note)
        note_cell.font = Font(italic=True, size=9)
        note_cell.alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)
        row += 1
    for col, name in enumerate(columns, 1):
        header = sheet.cell(row=row, column=col, value=name)
        header.font = Font(bold=True)
        header.fill = HEADER_FILL
        header.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    row += 1
    for record in rows:
        for col, name in enumerate(columns, 1):
            sheet.cell(row=row, column=col, value=as_number(record.get(name, "")))
        row += 1
    return row + 2


def finish_sheet(sheet, max_columns: int) -> None:
    for col in range(1, max_columns + 1):
        sheet.column_dimensions[get_column_letter(col)].width = 22 if col <= 2 else 16
    for row in sheet.iter_rows():
        for cell in row:
            if isinstance(cell.value, float):
                cell.number_format = (
                    "0.000000E+00" if cell.value != 0 and abs(cell.value) < 0.001 else "0.0000"
                )


def build_s11(workbook) -> None:
    if "Table_S11" in workbook.sheetnames:
        del workbook["Table_S11"]
    sheet = workbook.create_sheet("Table_S11")
    title = sheet.cell(
        row=1,
        column=1,
        value=(
            "Supplementary Table S11. Chromatin analysis in a single hg38 coordinate system: "
            "liftOver QC, eccDNA eligibility, and absolute and relative overlap with histone-mark "
            "peak sets under three localization definitions."
        ),
    )
    title.font = Font(bold=True, size=12)
    title.alignment = Alignment(horizontal="left", vertical="center")
    row = 3

    row = write_block(
        sheet,
        row,
        "S11a. liftOver QC for every hg19 peak set converted to hg38",
        read_tsv("chipseq_liftover_qc.tsv"),
        [
            "peak_set_id",
            "dataset",
            "mark",
            "analysis_tier",
            "original_peak_count",
            "liftover_success_count",
            "liftover_failure_count",
            "liftover_success_rate",
            "multiple_mapping_input_count",
            "abnormal_length_ratio_count_0.5_2",
            "analysis_unmerged_peak_count",
            "analysis_merged_peak_count",
            "analysis_merged_coverage_bp",
        ],
        note=(
            "Peaks were lifted individually with the UCSC hg19ToHg38 chain and retained only if "
            "canonical, within hg38 bounds, of positive length, and with a lifted/original length "
            "ratio between 0.5 and 2.0."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S11b. Per-sample eccDNA eligibility after canonical, bounds and blacklist/gap filtering",
        read_tsv("eccdna_sample_qc.tsv"),
        [
            "sample_id",
            "group",
            "mapped_alignments_idxstats",
            "original_eccdna_count",
            "valid_coordinate_count",
            "canonical_count",
            "within_chrom_bounds_count",
            "blacklist_gap_overlap_count",
            "eligible_eccdna_count",
        ],
    )

    row = write_block(
        sheet,
        row,
        "S11c. Within-cohort enrichment relative to the matched-placement expectation",
        read_tsv("p0_within_group_enrichment.tsv"),
        [
            "peak_set_id",
            "dataset",
            "mark",
            "analysis_tier",
            "mode",
            "group",
            "n_samples",
            "median_enrichment_ratio",
            "ratio_ci_low",
            "ratio_ci_high",
            "median_log2_enrichment",
            "n_samples_above_expected",
            "n_samples_below_expected",
            "wilcoxon_signed_rank_p_vs_0",
            "sign_test_p_vs_0",
            "bh_fdr_within_family",
        ],
        note=(
            "One-sample two-sided Wilcoxon signed-rank test of log2(observed/expected) against 0, "
            "computed separately in each cohort. Confidence intervals are percentile bootstrap "
            "intervals of the median over 10,000 resamples."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S11d. COVID-19 versus HC relative enrichment, all samples",
        read_tsv("relative_enrichment_group_stats.tsv"),
        [
            "peak_set_id",
            "dataset",
            "mark",
            "analysis_tier",
            "mode",
            "n_covid",
            "n_hc",
            "median_covid",
            "median_hc",
            "median_difference_covid_minus_hc",
            "mannwhitneyu_p",
            "cliffs_delta_covid_vs_hc",
            "bh_fdr_all_tests",
        ],
        note=(
            "Values are log2(observed/expected). These unadjusted comparisons are superseded by "
            "the burden-controlled analyses in Table S12 for interpretation."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S11e. COVID-19 versus HC absolute abundance within peak sets",
        read_tsv("absolute_burden_group_stats.tsv"),
        [
            "peak_set_id",
            "dataset",
            "mark",
            "analysis_tier",
            "mode",
            "n_covid",
            "n_hc",
            "median_covid",
            "median_hc",
            "median_difference_covid_minus_hc",
            "mannwhitneyu_p",
            "cliffs_delta_covid_vs_hc",
            "bh_fdr_all_tests",
        ],
        note=(
            "EPM(s,k) = 1e6 x (eccDNAs of sample s overlapping peak set k) / (total mapped "
            "alignments of sample s). This denominator differs from the peak-set-internal "
            "denominator used in the previous version of the analysis."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S11f. COVID-19 versus HC raw observed overlap fraction",
        read_tsv("observed_fraction_group_stats.tsv"),
        [
            "peak_set_id",
            "dataset",
            "mark",
            "analysis_tier",
            "mode",
            "median_covid",
            "median_hc",
            "median_difference_covid_minus_hc",
            "mannwhitneyu_p",
            "cliffs_delta_covid_vs_hc",
            "bh_fdr_all_tests",
        ],
    )

    row = write_block(
        sheet,
        row,
        "S11g. Per-sample observed, expected and relative enrichment for every peak set and mode",
        read_tsv("sample_relative_enrichment.tsv"),
        [
            "sample_id",
            "group",
            "peak_set_id",
            "mark",
            "analysis_tier",
            "mode",
            "eligible_eccdna_count",
            "denominator_units",
            "observed_overlap_units",
            "observed_fraction",
            "mean_shuffled_overlap_fraction",
            "enrichment_ratio",
            "log2_enrichment_ratio",
            "placement_bins",
            "matched_control_status",
        ],
        note=(
            "mean_shuffled_overlap_fraction is the exact expectation over all chromosome- and "
            "length-matched placements in the allowed space, not a Monte-Carlo estimate."
        ),
    )
    finish_sheet(sheet, 16)


def build_s12(workbook) -> None:
    if "Table_S12" in workbook.sheetnames:
        del workbook["Table_S12"]
    sheet = workbook.create_sheet("Table_S12")
    title = sheet.cell(
        row=1,
        column=1,
        value=(
            "Supplementary Table S12. Controls for the difference in detected eccDNA burden "
            "between cohorts: downsampling, burden-matched subset, covariate-adjusted regression, "
            "and the association between relative enrichment and eccDNA burden."
        ),
    )
    title.font = Font(bold=True, size=12)
    title.alignment = Alignment(horizontal="left", vertical="center")
    row = 3

    row = write_block(
        sheet,
        row,
        "S12a. COVID-19 versus HC after downsampling every sample to a common eccDNA count",
        read_tsv("p0_downsampled_group_stats.tsv"),
        [
            "peak_set_id",
            "mark",
            "analysis_tier",
            "mode",
            "downsample_target",
            "n_covid",
            "n_hc",
            "median_covid",
            "median_hc",
            "median_difference_covid_minus_hc",
            "mannwhitneyu_p",
            "cliffs_delta_covid_vs_hc",
            "bh_fdr_within_family",
        ],
        note=(
            "Each sample was reduced without replacement to the target count 100 times under a "
            "fixed seed, with the matched expectation recomputed from each subsample and the "
            "resulting log2 enrichment averaged."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S12b. Agreement between full-sample and downsampled estimates",
        read_tsv("p0_downsample_vs_full_agreement.tsv"),
        [
            "peak_set_id",
            "mark",
            "analysis_tier",
            "mode",
            "downsample_target",
            "n_samples",
            "pearson_r_full_vs_downsampled",
            "mean_absolute_difference",
            "max_absolute_difference",
        ],
        note=(
            "Because observed/expected is an intensive per-sample quantity, a random subsample is "
            "an unbiased estimator of the full-sample value. Downsampling therefore equalizes the "
            "precision of each sample's estimate rather than controlling for detection burden."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S12c. COVID-19 versus HC within the burden-matched eccDNA count window",
        read_tsv("p0_burden_matched_subset.tsv"),
        [
            "peak_set_id",
            "mark",
            "analysis_tier",
            "mode",
            "matched_window_low",
            "matched_window_high",
            "n_covid",
            "n_hc",
            "residual_burden_difference_p",
            "median_difference_covid_minus_hc",
            "mannwhitneyu_p",
            "cliffs_delta_covid_vs_hc",
            "bh_fdr_within_family",
        ],
    )

    row = write_block(
        sheet,
        row,
        "S12d. Samples included in the burden-matched window",
        read_tsv("p0_burden_matched_membership.tsv"),
        [
            "sample_id",
            "group",
            "eligible_eccdna_count",
            "in_matched_window",
            "matched_window_low",
            "matched_window_high",
        ],
    )

    row = write_block(
        sheet,
        row,
        "S12e. HC3 robust regression of log2 enrichment, nested models",
        read_tsv("p0_covariate_adjusted_hc3.tsv"),
        [
            "peak_set_id",
            "mark",
            "analysis_tier",
            "mode",
            "model",
            "term",
            "n_samples",
            "beta",
            "hc3_standard_error",
            "ci_low",
            "ci_high",
            "t_value",
            "p_value",
            "bh_fdr_within_family",
            "vif",
            "model_r_squared",
        ],
        note=(
            "Model conventions follow the RCA technical-bias analysis: log2 RCA concentration "
            "centred, age in units of 10 years from 70, male indicator, HC3 robust standard "
            "errors. VIF is reported because group and log10 eccDNA count are collinear."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S12f. Spearman association between relative enrichment and detected eccDNA count",
        read_tsv("p0_burden_association.tsv"),
        [
            "peak_set_id",
            "mark",
            "analysis_tier",
            "mode",
            "stratum",
            "n_samples",
            "spearman_rho",
            "spearman_p",
            "bh_fdr_within_family",
        ],
        note=(
            "A within-cohort association means the enrichment score tracks detection burden "
            "independently of group, which is why the group and burden contrasts cannot be "
            "separated in this design."
        ),
    )
    finish_sheet(sheet, 16)


def build_s13(workbook) -> None:
    if "Table_S13" in workbook.sheetnames:
        del workbook["Table_S13"]
    sheet = workbook.create_sheet("Table_S13")
    title = sheet.cell(
        row=1,
        column=1,
        value=(
            "Supplementary Table S13. Cross-cell-type replication using GRCh38-native ENCODE "
            "histone peak sets, and mappability controls."
        ),
    )
    title.font = Font(bold=True, size=12)
    title.alignment = Alignment(horizontal="left", vertical="center")
    row = 3

    row = write_block(
        sheet,
        row,
        "S13a. GRCh38-native ENCODE peak sets used, with provenance",
        read_tsv("ext_multicell_peak_qc.tsv"),
        [
            "peak_set_id", "biosample", "mark", "experiment_accession", "file_accession",
            "output_type", "assembly", "source_peak_count", "canonical_in_bounds_count",
            "merged_peak_count", "merged_coverage_bp", "local_sha256",
        ],
        note=(
            "Files were selected from ENCODE as released GRCh38 peak BEDs, preferring the "
            "experiment's default representation and replicated over pseudoreplicated peaks. "
            "Because they are already GRCh38 no liftOver was applied, so this comparison "
            "inherits no conversion artefact from the primary analysis."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S13b. Enrichment and COVID-19 versus HC comparison for every cell type and mark",
        read_tsv("ext_multicell_group_stats.tsv"),
        [
            "peak_set_id", "dataset", "mark", "mode", "n_covid", "n_hc",
            "median_ratio_covid", "ratio_ci_low_covid", "ratio_ci_high_covid",
            "median_ratio_hc", "ratio_ci_low_hc", "ratio_ci_high_hc",
            "n_above_expected_covid", "n_above_expected_hc",
            "wilcoxon_p_vs_0_covid", "bh_fdr_covid_vs_0",
            "cliffs_delta_covid_vs_hc", "mannwhitneyu_p", "bh_fdr_group_within_family",
        ],
        note=(
            "eccDNA lengths were rounded to 25 bp bins to keep the 33-peak-set computation "
            "tractable; S13e quantifies the effect of that approximation."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S13c. Mappability of each peak set relative to the genome background",
        read_tsv("ext_peak_set_mappability.tsv"),
        [
            "peak_set_id", "dataset", "mark", "analysis_tier", "peak_bp", "mappable_bp",
            "mappable_fraction", "genome_background_mappable_fraction",
            "mappability_ratio_vs_background",
        ],
        note="Umap k36 single-read mappability for GRCh38 (Bismap).",
    )

    row = write_block(
        sheet,
        row,
        "S13d. Enrichment recomputed with eccDNAs and placement space restricted to mappable sequence",
        read_tsv("ext_mappability_group_stats.tsv"),
        [
            "peak_set_id", "mark", "mode", "n_covid", "n_hc",
            "median_ratio_covid", "median_ratio_hc",
            "n_above_expected_covid", "n_above_expected_hc",
            "cliffs_delta_covid_vs_hc", "mannwhitneyu_p", "bh_fdr_group_within_family",
        ],
        note=(
            "The mappable mask was built by binning Umap k36 to 1 kb, keeping bins at least 80% "
            "mappable, merging, and dropping components under 5 kb; it retains 2,572,673,328 of "
            "2,831,289,217 allowed bp in 64,324 components."
        ),
    )

    row = write_block(
        sheet,
        row,
        "S13e. Effect of the 25 bp length binning, measured against the exact primary analysis",
        read_tsv("ext_binning_check_agreement.tsv"),
        [
            "peak_set_id", "mark", "mode", "length_bin_bp", "n_samples",
            "pearson_r_exact_vs_binned", "mean_absolute_difference", "max_absolute_difference",
        ],
        note=(
            "Across the six CD14+ reference marks binning changes per-sample log2 enrichment by "
            "at most 0.0054 (r >= 0.9999). The single large deviation is confined to the smallest "
            "legacy A549 peak set under the junction definition, which supports no conclusion."
        ),
    )
    finish_sheet(sheet, 19)


def main() -> int:
    if not WORKBOOK.exists():
        raise FileNotFoundError(WORKBOOK)
    backup = WORKBOOK.with_suffix(".pre_chromatin_bak.xlsx")
    if not backup.exists():
        shutil.copy2(WORKBOOK, backup)
        print(f"backed up {WORKBOOK.name} -> {backup.name}")

    workbook = load_workbook(WORKBOOK)
    before = list(workbook.sheetnames)
    build_s11(workbook)
    build_s12(workbook)
    build_s13(workbook)
    workbook.save(WORKBOOK)

    check = load_workbook(WORKBOOK, read_only=True)
    print(f"sheets before: {before}")
    print(f"sheets after:  {check.sheetnames}")
    for name in ("Table_S11", "Table_S12", "Table_S13"):
        sheet = check[name]
        print(f"  {name}: {sheet.max_row} rows x {sheet.max_column} columns")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
