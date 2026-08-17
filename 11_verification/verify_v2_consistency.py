#!/usr/bin/env python3
"""Stage S6 -- cross-check every v2 number across manuscript, letters and tables.

The failure this guards against is a partial migration: a number corrected in
the manuscript but left stale in a response letter, or a table that still
carries the v1 value the text now contradicts. Reading the three artefacts and
comparing them mechanically is the only way to be sure, so nothing here relies
on having read the documents by eye.

Three checks:

  1. no v1 value survives anywhere it should have been replaced
  2. every v2 value the manuscript prints is present
  3. the manuscript's headline numbers agree with the locked tables they cite
"""

from __future__ import annotations

import csv
import json
import re
import sys
import zipfile
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]

MANUSCRIPT = RUN / "manuscript" / "Maintext_revision_tracked_changes.docx"
LETTERS = [RUN / "response" / f"JARE-D-26-04883_Reviewer_{n}_response.docx"
           for n in (1, 2, 3)]
WORKBOOK = RUN / "tables" / "supplementary" / "Supplementary_Tables_v2.xlsx"
NUMBERS = RUN / "logs" / "v2_manuscript_numbers.json"

# v1 values that must not survive. The disclosure paragraph in the letters
# quotes the old total on purpose, so that one pair is exempted there.
STALE = [
    # round 1 (call-set migration)
    "15,006,670", "1,850,490", "384,786", "47,448", "291,263", "42,902",
    "3,143 [", "26,642", "3,402", "3,094", "3,271", "11,707",
    "281 candidate", "27 exact", "0.862", "196, 365, 571",
    "760 bp", "760-bp", "66.8%", "0.898 [", "0.0191 [",
    # round 2: values the hand-written round-1 list never reached. Each one was
    # found by sweeping every statistic-bearing sentence against a locked table,
    # which is why the list is now organised by the table it came from rather
    # than by the sentence that happens to quote it.
    #   Table S3, fragment-length bins
    "δ = 0.394", "0.00414", "−0.638", "3.67 × 10", "q ≥ 0.0788",
    #   matched pairs (regen_S1)
    "6.8×", "5.2×",
    #   tables/RCA
    #   ("64.0 [57.0" was listed here until 2026-08-13. RCA concentration does
    #   not depend on the call set; that entry tracked a precision change, and
    #   the author has since chosen to report the quartiles as integers, so the
    #   integer form is now the required value rather than a stale one.)
    "ρ = 0.400", "ρ = 0.322", "0.0071", "1.37 × 10", "ρ = 0.211",
    "1.42 × 10", "0.0070", "3.28-fold", "4.74-fold", "2.87 × 10", "2.46 × 10",
    "ρ = 0.511", "ρ = 0.095",
    #   tables/robustness
    "3.31–4.84", "3.31 to 4.84",
    #   tables/clinical
    "q = 0.94", "q = 1.00", "ρ = −0.575", "0.61; 0.03 to 0.93",
    #   tables/chromatin -- Figure 4 is drawn from full_interval, so any number
    #   that matches junction_start_end in this section is also an error
    "1.16–1.30", "0.77–0.80", "1.20 in COVID-19", "0.754 to 0.780",
    "0.740 to 0.769", "0.740–0.769", "δ = 0.675", "δ = −0.552", "δ = 0.324",
    "288,062", "41,480", "δ = 0.674", "42,527", "114,090", "18 HC",
    "q ≥ 0.97", "q ≥ 0.55", "P = 0.84", "ρ = 0.563", "−0.557", "r = 0.73",
    "1.24–1.38", "1.25–1.31", "ratio of 1.00", "0.766 to 0.756",
    "−0.05 in lung", "1.25 in COVID-19", "1.34 in HC", "δ = 0.716",
    "δ = 0.679", "δ = 0.040", "164,708", "46,555", "retaining all 78 samples",
    #   eccGene definition concordance (regen_S6)
    "0.865", "0.256", "0.966 vs", "12,793",
    #   pooled call burden
    "eight-fold", "15.0 versus 1.85",
    #   sentences whose antecedent no longer exists
    "five lower features",
]

STALE_EXEMPT_IN_LETTERS = {"16,857,160"}

REQUIRED = [
    "10,099,219", "8,830,916", "1,268,303", "226,434", "190,320",
    "2,053", "0.849", "0.776", "26,087", "3,487", "3,189", "9,910",
    "19,450", "196, 366, 571", "764", "278", "TFE3", "16 exact",
    "4.84", "0.910", "NAPSB", "three wet-lab validation targets", "CpG islands",
    "167–172", "202–207", "8,805,915", "1,262,947",
    "5,000-replicate", "TCF7L1",
    # round 2
    "0.215", "0.153", "−0.466", "0.00119", "6.1×", "4.4×", "0.0283",
    "64 [57–123]", "0.398", "0.283", "0.0209", "2.30 × 10", "0.0065",
    "3.20-fold", "4.63-fold", "3.31–5.08", "0.86", "0.97", "−0.568",
    "1.22–1.38", "0.75–0.78", "1.26 in COVID-19", "0.748 to 0.778",
    "0.546", "−0.393", "189,679", "25,542", "0.530", "−0.426",
    "25,887", "76,350", "19 HC", "0.52", "0.82", "0.430", "−0.349",
    "0.72", "1.33–1.52", "1.32–1.43", "ratio of 1.06", "0.747 to 0.736",
    "−0.18 in lung", "1.22 in COVID-19", "0.704", "0.677", "−0.174",
    "0.888", "0.869", "0.969", "0.957",
    "seven-fold", "8.83 versus 1.27",
]

# Values that appear only in the response letters, not in the manuscript body.
REQUIRED_IN_LETTERS = ["24,986", "5,342", "12,333", "0.72", "1.06"]


def accepted_text(path: Path) -> str:
    x = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    x = re.sub(r"<w:del\b[^>]*>.*?</w:del>", "", x, flags=re.S)
    return "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", x, re.S))


def letter_text(path: Path) -> str:
    from docx import Document
    doc = Document(str(path))
    parts = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            parts += [c.text for c in row.cells]
    return "\n".join(parts)


def main() -> int:
    failures: list[str] = []

    text = accepted_text(MANUSCRIPT)
    stale = [s for s in STALE if s in text]
    missing = [s for s in REQUIRED if s not in text]
    print(f"manuscript  stale={len(stale)}  missing={len(missing)}")
    if stale:
        failures.append(f"manuscript retains v1 values: {stale}")
    if missing:
        failures.append(f"manuscript missing v2 values: {missing}")
    if "PLACEHOLDER_" in text:
        failures.append("manuscript retains repository placeholders")
    if "most enriched at 5′UTRs and exons and least enriched at introns" in text:
        failures.append("manuscript retains the incorrect Summary gene-feature ordering")
    if (
        "chromosome-normalized and gene-feature analyses were reimplemented"
        " for this revision from the stated Methods because their original"
        " generating scripts were unavailable"
    ) in text:
        failures.append("manuscript misstates the recovered gene-feature script status")

    with (RUN / "tables" / "fragment_length" /
          "adjacent_peak_spacing.tsv").open(newline="") as handle:
        spacing_rows = list(csv.DictReader(handle, delimiter="\t"))
    peak_positions = [int(spacing_rows[0]["left_primary_position_bp"])] + [
        int(row["right_primary_position_bp"]) for row in spacing_rows
    ]
    peak_phrase = ", ".join(map(str, peak_positions[:-1])) + f" and {peak_positions[-1]}"
    if peak_phrase not in text:
        failures.append(f"manuscript does not report source-derived peaks: {peak_positions}")
    for row in spacing_rows:
        low = float(row["bootstrap_spacing_ci_low_bp"])
        high = float(row["bootstrap_spacing_ci_high_bp"])
        ci_phrase = f"{low:.0f}–{high:.0f}"
        if ci_phrase not in text:
            failures.append(f"manuscript lacks source-derived spacing CI {ci_phrase}")

    with (RUN / "tables" / "chromatin" /
          "eccdna_sample_qc.tsv").open(newline="") as handle:
        eligible_rows = list(csv.DictReader(handle, delimiter="\t"))
    eligible_totals = {}
    for row in eligible_rows:
        eligible_totals.setdefault(row["group"], 0)
        eligible_totals[row["group"]] += int(row["eligible_eccdna_count"])
    for group, value in eligible_totals.items():
        if f"{value:,}" not in text:
            failures.append(
                f"manuscript lacks source-derived eligible count for {group}: {value:,}"
            )

    for path in LETTERS:
        body = letter_text(path)
        bad = [s for s in STALE
               if s in body and s not in STALE_EXEMPT_IN_LETTERS]
        disclosed = "Note added in this revision" in body
        print(f"{path.stem[-20:]:22s} stale={len(bad)}  disclosure={'yes' if disclosed else 'NO'}")
        if bad:
            failures.append(f"{path.name} retains v1 values: {bad}")
        if not disclosed:
            failures.append(f"{path.name} has no disclosure paragraph")
        if "PLACEHOLDER_" in body:
            failures.append(f"{path.name} retains repository placeholders")
        if "10,099,219 of 16,857,160 calls (60.1%)" in body:
            failures.append(f"{path.name} retains the mixed-denominator percentage")
        if "10,099,219 entries (59.9%)" not in body:
            failures.append(f"{path.name} lacks the corrected archived-entry percentage")
        if "all numbers in the revised manuscript now come from it" in body:
            failures.append(f"{path.name} overstates the call-set correction scope")
        if (
            "gene-feature and recurrent exact-interval tables were reimplemented"
            " from the stated Methods because their original generating scripts"
            " were unavailable"
        ) in body:
            failures.append(f"{path.name} misstates the gene-feature script provenance")
    all_letters = "\n".join(letter_text(p) for p in LETTERS)
    absent_letters = [s for s in REQUIRED_IN_LETTERS if s not in all_letters]
    print(f"letter-only v2 values present: "
          f"{len(REQUIRED_IN_LETTERS) - len(absent_letters)}/{len(REQUIRED_IN_LETTERS)}")
    if absent_letters:
        failures.append(f"letters missing v2 values: {absent_letters}")

    # manuscript headline numbers vs the locked tables
    ref = json.loads(NUMBERS.read_text())
    burden = ref["burden"]["v2"]
    checks = {
        f"{burden['total_all']:,}": "total calls",
        f"{burden['COVID-19']['total']:,}": "COVID calls",
        f"{burden['HC']['total']:,}": "HC calls",
        f"{ref['eccGene']['tested']:,}": "tested genes",
        f"{ref['eccGene']['significant']:,}": "differential eccGenes",
        f"{ref['loci']['candidate_assignments']}": "candidate assignments",
        f"{ref['loci']['recurrent_min10_intervals']} exact": "recurrent intervals",
    }
    absent = [label for value, label in checks.items() if value not in text]
    print(f"table-to-text agreement: {len(checks) - len(absent)}/{len(checks)}")
    if absent:
        failures.append(f"manuscript does not print locked values for: {absent}")

    # Wet-lab target evidence. Read the locked tables so
    # the manuscript cannot silently drift from the actual per-call-set counts.
    selected = ["CTNNA2", "SDK1", "TCF7L1"]
    with (RUN / "tables" / "robustness" /
          "validation_target_support_summary.tsv").open(newline="") as handle:
        target_rows = list(csv.DictReader(handle, delimiter="\t"))

    def target_counts(callset: str, group: str = "COVID-19") -> list[int]:
        index = {
            row["target"]: int(row["supported_samples"])
            for row in target_rows
            if row["callset"] == callset
            and row["group"] == group
            and int(row["tolerance_bp"]) == 10
            and row["target"] in selected
        }
        return [index[target] for target in selected]

    methods_counts = target_counts("circlemap_methods")
    consensus_counts = target_counts("consensus_t10")
    hc_counts = target_counts("circlemap_methods", "HC")
    strict_sdk1 = target_counts("circlemap_strict")[1]
    def prose_counts(values: list[int]) -> str:
        return f"{values[0]}, {values[1]} and {values[2]}"

    target_text_checks = {
        prose_counts(methods_counts) + " of 39 COVID-19 samples":
            "methods target counts",
        prose_counts(consensus_counts) + " of 39":
            "consensus target counts",
        "one breakpoint window overlapped 62 bp of an L2a element":
            "SDK1 repeat overlap",
        "381, 380 and 380 bp": "matched target lengths",
    }
    absent_target_text = [
        label for value, label in target_text_checks.items() if value not in text
    ]
    if hc_counts != [0, 0, 0]:
        failures.append(f"selected-target HC counts are not all zero: {hc_counts}")
    if strict_sdk1 != 13:
        failures.append(f"SDK1 strict split/score support is {strict_sdk1}, expected 13")
    if absent_target_text:
        failures.append(f"manuscript target paragraph missing: {absent_target_text}")
    print(
        "selected-target agreement: methods="
        f"{methods_counts}, consensus={consensus_counts}, HC={hc_counts}, "
        f"SDK1_strict={strict_sdk1}"
    )

    from openpyxl import load_workbook
    wb = load_workbook(WORKBOOK, read_only=True)
    print(f"workbook sheets: {len(wb.sheetnames)}")
    if len(wb.sheetnames) != 21:
        failures.append(f"workbook has {len(wb.sheetnames)} sheets, expected 21")
    s16 = wb["Table_S16"]
    displayed_targets = {
        row[2]
        for row in s16.iter_rows(values_only=True)
        if len(row) > 2 and row[2] in {
            "CTNNA2", "SDK1", "TCF7L1", "CAB39", "chr22_target"
        }
    }
    if displayed_targets != set(selected):
        failures.append(
            f"Table S16 targets are {sorted(displayed_targets)}, expected {selected}"
        )
    print(f"Table S16 selected targets: {sorted(displayed_targets)}")

    s15 = wb["Table_S15"]
    s15_targets = {
        row[2]
        for row in s15.iter_rows(values_only=True)
        if len(row) > 2 and row[2] in {
            "CTNNA2", "SDK1", "TCF7L1", "CAB39", "chr22_target"
        }
    }
    if s15_targets != set(selected):
        failures.append(
            f"Table S15 target genes are {sorted(s15_targets)}, expected {selected}"
        )
    print(f"Table S15 wet-lab targets: {sorted(s15_targets)}")
    s15_headers = [cell.value for cell in s15[27]]
    required_s15_detection = {
        "detection_odds_ratio_COVID_vs_HC",
        "fisher_detection_p",
        "fisher_detection_q_BH",
        "detection_direction",
    }
    missing_s15_detection = required_s15_detection - set(s15_headers)
    if missing_s15_detection:
        failures.append(
            "Table S15 lacks target detection fields: "
            f"{sorted(missing_s15_detection)}"
        )
    else:
        header_index = {value: index for index, value in enumerate(s15_headers)}
        tcf7l1_rows = [
            row
            for row in s15.iter_rows(min_row=30, max_row=56, values_only=True)
            if row[2] == "TCF7L1"
        ]
        abundance_max_q = max(float(row[8]) for row in tcf7l1_rows)
        detection_max_q = max(
            float(row[header_index["fisher_detection_q_BH"]])
            for row in tcf7l1_rows
        )
        if len(tcf7l1_rows) != 9:
            failures.append(
                f"Table S15 has {len(tcf7l1_rows)} TCF7L1 stability rows, expected 9"
            )
        if abs(abundance_max_q - 0.0306909131) > 1e-12:
            failures.append(
                f"TCF7L1 abundance max q is {abundance_max_q}, expected 0.0306909131"
            )
        if abs(detection_max_q - 0.004463907019) > 1e-12:
            failures.append(
                f"TCF7L1 detection max q is {detection_max_q}, expected 0.004463907019"
            )
        print(
            "TCF7L1 endpoint maxima: abundance q="
            f"{abundance_max_q:.10g}, detection q={detection_max_q:.10g}"
        )

    s12 = wb["Table_S12"]
    s12_headers = [cell.value for cell in s12[2]]
    obsolete_s12 = {"Coverage Ratio", "Nomalized EccDNA ratio"} & set(s12_headers)
    if obsolete_s12:
        failures.append(f"Table S12 retains obsolete columns: {sorted(obsolete_s12)}")

    s3 = wb["Table_S3"]
    s3_headers = [cell.value for cell in s3[2]]
    s3_index = {value: index for index, value in enumerate(s3_headers)}
    for label in ("COVID-19", "HC"):
        required = {f"{label} Q1", f"{label} Q3", f"{label} IQR"}
        if not required <= set(s3_headers):
            failures.append(f"Table S3 lacks {label} quartile/IQR columns")
            continue
        bad_iqr = 0
        for row in s3.iter_rows(min_row=3, values_only=True):
            q1 = row[s3_index[f"{label} Q1"]]
            q3 = row[s3_index[f"{label} Q3"]]
            iqr = row[s3_index[f"{label} IQR"]]
            if q1 is None and q3 is None and iqr is None:
                continue
            if None in (q1, q3, iqr) or abs(float(iqr) - (float(q3) - float(q1))) > 1e-10:
                bad_iqr += 1
        if bad_iqr:
            failures.append(f"Table S3 has {bad_iqr} invalid {label} IQR values")

    s13 = wb["Table_S13"]
    s13_headers = [cell.value for cell in s13[2]]
    if "Observed/expected ratio" not in s13_headers:
        failures.append("Table S13 lacks the observed/expected-ratio column")
    else:
        oe_column = s13_headers.index("Observed/expected ratio")
        workbook_oe = {}
        for row in s13.iter_rows(min_row=3, values_only=True):
            if len(row) <= oe_column or row[oe_column] is None:
                continue
            workbook_oe[(str(row[0]), str(row[1]), str(row[2]))] = float(
                row[oe_column]
            )
        oe_count = len(workbook_oe)
        if oe_count != 546:
            failures.append(
                f"Table S13 has {oe_count} observed/expected values, expected 546"
            )
        with (RUN / "tables" / "genomic_distribution" /
              "gene_feature_enrichment_scores.tsv").open(newline="") as handle:
            source_oe = {
                (row["sample_id"], row["group"], row["element"]): float(row["oe_ratio"])
                for row in csv.DictReader(handle, delimiter="\t")
            }
        if set(workbook_oe) != set(source_oe):
            failures.append("Table S13 keys do not match the locked gene-feature source")
        else:
            bad_oe = sum(
                abs(workbook_oe[key] - source_oe[key]) > 1e-12 for key in source_oe
            )
            if bad_oe:
                failures.append(
                    f"Table S13 has {bad_oe} O/E values differing from locked source"
                )
        print(f"Table S13 observed/expected values: {oe_count}/546")

    print()
    if failures:
        print(f"FAIL -- {len(failures)} problem(s):")
        for f in failures:
            print(f"  {f}")
        return 1
    print("PASS -- manuscript, letters and tables agree on every checked value")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
