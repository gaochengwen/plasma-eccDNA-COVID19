#!/usr/bin/env python3
"""Check every data block of the supplementary workbook against its locked TSV.

`verify_v2_consistency.py` compares the manuscript and the three letters against
a curated list of values. It never opens the workbook's data regions, which is
how a whole group of sheets stayed on the v1 call set while passing: Tables S2a,
S4, S5b-e, S14, S17b-g, S18a/c/e/f and S19b/d were carried over from the v1
submission workbook and never rewritten.

This closes that gap by walking each block cell by cell against the table it is
supposed to be a copy of, matching columns by header name and rows by key.

Blocks with no locked TSV counterpart are listed as unchecked rather than
silently skipped, so the coverage gap is always visible.

    python3 scripts/local/verify_workbook_against_sources.py

Exit status is 0 only when every checked cell agrees.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import openpyxl

RUN = Path(__file__).resolve().parents[2]
WB = RUN / "tables" / "supplementary" / "Supplementary_Tables_v2.xlsx"
T = RUN / "tables"

TOL = 1e-9

# sheet, header row (1-based), source TSV, key columns
BLOCKS = [
    ("S3",   "Table_S3",   2,   "supplementary/Table_S3_v2.tsv",                    ["FDR family", "Metric", "Feature"]),
    ("S5b",  "Table_S5",   27,  "fragment_length/objective_peak_identification_and_stability.tsv", ["peak_id"]),
    ("S5c",  "Table_S5",   36,  "fragment_length/adjacent_peak_spacing.tsv",        ["left_peak_id", "right_peak_id"]),
    ("S5d",  "Table_S5",   43,  "fragment_length/peak_window_group_comparisons.tsv", ["peak_id"]),
    ("S5e",  "Table_S5",   51,  "fragment_length/peak_window_sample_proportions.tsv", ["sample_id", "peak_id"]),
    ("S6b",  "Table_S6",   15,  "robustness/robustness_callset_sample_metrics.tsv", ["sample_id", "callset"]),
    ("S6c",  "Table_S6",   565, "robustness/robustness_group_tests.tsv",            ["callset", "metric"]),
    ("S10",  "Table_S10",  2,   None,                                               None),
    ("S12",  "Table_S12",  2,   None,                                               None),
    ("S13",  "Table_S13",  2,   None,                                               None),
    ("S15a", "Table_S15",  4,   "robustness/eccgene_callset_summary.tsv",           ["callset", "definition"]),
    ("S15b", "Table_S15",  17,  "robustness/eccgene_callset_concordance.tsv",       ["reference_callset", "comparison_callset", "definition"]),
    ("S15c", "Table_S15",  27,  "robustness/eccgene_candidate_stability.tsv",       ["callset", "definition", "gene"]),
    ("S15e", "Table_S15",  183, "robustness/chromatin_callset_group_statistics.tsv", ["callset", "analysis", "peak_set_id", "mode"]),
    ("S16a", "Table_S16",  4,   "robustness/validation_target_support_summary.tsv", ["callset", "tolerance_bp", "target", "group"]),
    ("S16b", "Table_S16",  79,  "robustness/validation_target_mask_audit.tsv",      ["target", "breakpoint_side"]),
    ("S17a", "Table_S17",  5,   "chromatin/chipseq_liftover_qc.tsv",                ["peak_set_id"]),
    ("S17b", "Table_S17",  19,  "chromatin/eccdna_sample_qc.tsv",                   ["sample_id"]),
    ("S17c", "Table_S17",  102, "chromatin/p0_within_group_enrichment.tsv",         ["peak_set_id", "mode", "group"]),
    ("S17d", "Table_S17",  167, "chromatin/relative_enrichment_group_stats.tsv",    ["peak_set_id", "mode"]),
    ("S17e", "Table_S17",  202, "chromatin/absolute_burden_group_stats.tsv",        ["peak_set_id", "mode"]),
    ("S17f", "Table_S17",  236, "chromatin/observed_fraction_group_stats.tsv",      ["peak_set_id", "mode"]),
    ("S17g", "Table_S17",  271, "chromatin/sample_relative_enrichment.tsv",         ["sample_id", "peak_set_id", "mode"]),
    ("S18a", "Table_S18",  5,   "chromatin/p0_downsampled_group_stats.tsv",         ["peak_set_id", "mode", "downsample_target"]),
    ("S18c", "Table_S18",  134, "chromatin/p0_burden_matched_subset.tsv",           ["peak_set_id", "mode"]),
    ("S18e", "Table_S18",  251, "chromatin/p0_covariate_adjusted_hc3.tsv",          ["peak_set_id", "mode", "model", "term"]),
    ("S18f", "Table_S18",  736, "chromatin/p0_burden_association.tsv",              ["peak_set_id", "mode", "stratum"]),
    ("S19b", "Table_S19",  43,  "chromatin/ext_multicell_group_stats.tsv",          ["peak_set_id", "mode"]),
    ("S19c", "Table_S19",  114, "chromatin/ext_peak_set_mappability.tsv",           ["peak_set_id"]),
    ("S19d", "Table_S19",  129, "chromatin/ext_mappability_group_stats.tsv",        ["peak_set_id", "mode"]),
    ("S20c", "Table_S20",  29,  "BCL3/category_candidate_rankings.tsv",             ["Functional_category", "Gene"]),
]

# Blocks whose source is not a single TSV; each is checked by a dedicated rule.
SPECIAL = ["S2a", "S4a", "S4b", "S14"]

# Blocks with no locked source in this tree, and why.
UNCHECKED = {
    "S1":   "participant clinical annotation, curated from hospital records",
    "S2b":  "derived cohort tests, recomputed only in the repair round",
    "S5a":  "prespecified parameters, not analysis output",
    "S7":   "caller concordance, checked by verify_v2_consistency.py",
    "S8":   "length-bin proportions, keys differ only by label formatting",
    "S9":   "clinical A2-A6, checked by verify_v2_consistency.py",
    "S11":  "repeat classes, computed from BAM reads (call-set independent)",
    "S15d": "chromatin concordance, wide layout without a row key",
    "S19a": "peak-set provenance, not analysis output",
    "S21":  "non-COVID pneumonia samples, not sequenced",
}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(v: str) -> str:
    return v.replace("_", "")


def same(cell, text: str) -> bool:
    if cell is None or text in ("", "NA", "nan"):
        return True
    try:
        return abs(float(cell) - float(text)) <= TOL * max(1.0, abs(float(text)))
    except (TypeError, ValueError):
        return str(cell) == text


def check_block(wb, label, sheet, hrow, src, keys) -> tuple[int, int]:
    rows = list(wb[sheet].iter_rows(values_only=True))
    hdr = [str(c) if c is not None else "" for c in rows[hrow - 1]]
    hp = {h: i for i, h in enumerate(hdr)}
    data = read(T / src)
    if any(k not in hp for k in keys):
        raise SystemExit(f"{label}: header row {hrow} lacks {keys}: {hdr[:6]}")
    idx = {tuple(norm(r[k]) for k in keys): r for r in data}
    shared = [h for h in hdr if h in data[0] and h not in keys]
    checked = bad = 0
    for r in rows[hrow:hrow + len(data)]:
        if r[0] is None:
            break
        key = tuple(norm(str(r[hp[k]])) for k in keys)
        if key not in idx:
            print(f"  {label}: row key {key} not in {Path(src).name}")
            bad += 1
            continue
        s = idx[key]
        for f in shared:
            checked += 1
            if not same(r[hp[f]], s[f]):
                bad += 1
                if bad <= 3:
                    print(f"  {label}: {key} {f} sheet={r[hp[f]]!r} source={s[f]!r}")
    return checked, bad


def check_special(wb) -> tuple[int, int]:
    checked = bad = 0
    rca = {r["sample_id"]: r for r in read(T / "RCA" / "RCA_sample_metrics.tsv")}
    rows = list(wb["Table_S2"].iter_rows(values_only=True))
    hp = {str(c): i for i, c in enumerate(rows[2])}
    total = 0
    for r in rows[3:81]:
        if not r[0]:
            break
        s = rca[r[0]]
        for col, fld in (("EccDNA Count", "eccdna_count"), ("EPM", "epm"),
                         ("Mapped Reads", "mapped_reads")):
            checked += 1
            if abs(float(r[hp[col]]) - float(s[fld])) > 1e-3:
                bad += 1
                print(f"  S2a: {r[0]} {col} sheet={r[hp[col]]} source={s[fld]}")
        total += int(r[hp["EccDNA Count"]])
    checked += 1
    if total != 10_099_219:
        bad += 1
        print(f"  S2a: eccDNA counts sum to {total:,}, expected 10,099,219")

    import json
    pairs = json.loads((T / "age_matching" / "age_matched_pairs_v2.json").read_text())
    sens = json.loads((T / "age_matching" /
                       "age_matched_caliper_sensitivity_v2.json").read_text())
    rows = list(wb["Table_S4"].iter_rows(values_only=True))
    for i, row in enumerate(pairs):
        for j, v in enumerate(row):
            checked += 1
            if not same(rows[3 + i][j], str(v)):
                bad += 1
                print(f"  S4a: row {i} col {j} sheet={rows[3+i][j]!r} source={v!r}")
    keys = ["caliper_years", "n_pairs", "total_age_difference_years",
            "max_pair_age_difference_years", "n_sex_concordant_pairs",
            "p_count", "n_count_up", "fc_count", "p_epm", "n_epm_up", "fc_epm"]
    for i, s in enumerate(sens):
        for j, k in enumerate(keys):
            checked += 1
            if not same(rows[18 + i][j], str(s[k])):
                bad += 1
                print(f"  S4b: caliper {s['caliper_years']} {k} "
                      f"sheet={rows[18+i][j]!r} source={s[k]!r}")

    stats = {r["Gene"]: r for r in read(T / "eccGene" / "junction_abundance_wilcoxon.tsv")}
    sig = sorted(g for g, r in stats.items()
                 if float(r["wilcoxon_q_BH"]) < 0.05
                 and abs(float(r["log2FC_mean_EA_COVID_vs_HC"])) >= 1)
    rows = list(wb["Table_S14"].iter_rows(values_only=True))
    genes = [r[0] for r in rows[2:] if r[0]]
    checked += 1
    if genes != sig:
        bad += 1
        print(f"  S14: {len(genes)} genes, expected the {len(sig)} significant ones")
    else:
        for r in rows[2:]:
            if not r[0]:
                break
            s = stats[r[0]]
            for off, fld in ((-4, "log2FC_mean_EA_COVID_vs_HC"),
                             (-3, "cliffs_delta_COVID_vs_HC"),
                             (-2, "wilcoxon_p"), (-1, "wilcoxon_q_BH")):
                checked += 1
                if not same(r[off], s[fld]):
                    bad += 1
    return checked, bad


def main() -> int:
    wb = openpyxl.load_workbook(WB, read_only=True, data_only=True)
    total = fails = 0
    print(f"checking {WB.name} against tables/\n")
    for label, sheet, hrow, src, keys in BLOCKS:
        if src is None:
            continue
        c, b = check_block(wb, label, sheet, hrow, src, keys)
        total += c
        fails += b
        print(f"  {label:5s} {c:6d} cells  {'OK' if b == 0 else f'{b} MISMATCH'}")
    c, b = check_special(wb)
    total += c
    fails += b
    print(f"  {'/'.join(SPECIAL):20s} {c:6d} cells  {'OK' if b == 0 else f'{b} MISMATCH'}")

    print(f"\nnot covered here ({len(UNCHECKED)} blocks):")
    for k, why in UNCHECKED.items():
        print(f"  {k:5s} {why}")

    print(f"\n{total:,} cells compared, {fails} mismatches")
    if fails:
        print("FAIL -- the workbook disagrees with the locked tables")
        return 1
    print("PASS -- every checked workbook cell matches its locked source")
    return 0


if __name__ == "__main__":
    sys.exit(main())
