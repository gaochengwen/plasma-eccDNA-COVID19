#!/usr/bin/env python3
"""Extract every manuscript-facing number from the v2 tables, beside its v1 value.

This is the backbone of stages S5 and S6: S5 edits the manuscript from this
report, S6 re-runs it to prove the manuscript and the tables agree. Nothing is
recomputed from raw data -- every v2 value is read out of a locked table, and
every v1 value is the figure currently printed in the submitted manuscript.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "tables"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(text: str) -> float:
    return float(text)


def eccgene_summary() -> dict[str, object]:
    ab = {r["Gene"]: r for r in rows(T / "eccGene" / "junction_abundance_wilcoxon.tsv")}
    tested = len(ab)

    def is_sig(r: dict[str, str]) -> bool:
        q, lfc = r["wilcoxon_q_BH"], r["log2FC_mean_EA_COVID_vs_HC"]
        if q in ("", "NA") or lfc in ("", "NA", "inf", "-inf"):
            return False
        return fnum(q) < 0.05 and abs(fnum(lfc)) >= 1

    sig = {g for g, r in ab.items() if is_sig(r)}
    up = {g for g in sig if fnum(ab[g]["log2FC_mean_EA_COVID_vs_HC"]) > 0}

    det = {r["Gene"]: r for r in rows(T / "eccGene" / "junction_detection_fisher.tsv")}
    hc_zero = [r for r in det.values() if r["hc_detected"] == "0"]
    hc_zero.sort(key=lambda r: (fnum(r["fisher_q_BH"]), -int(r["covid_detected"])))

    return {
        "tested_genes": tested,
        "differential_eccgenes": len(sig),
        "differential_up_in_covid": len(up),
        "fig3b_top30": [r["Gene"] for r in hc_zero[:30]],
    }


def fragment_peaks() -> dict[str, object]:
    peaks = rows(T / "fragment_length" / "objective_peak_identification_and_stability.tsv")
    spacing = rows(T / "fragment_length" / "adjacent_peak_spacing.tsv")
    return {
        "peaks_bp": [int(fnum(p["primary_peak_position_bp"])) for p in peaks],
        "hc_bootstrap_detection": [
            round(fnum(p["hc_bootstrap_detection_frequency"]), 3) for p in peaks
        ],
        "adjacent_spacing_bp": [int(fnum(s["primary_spacing_bp"])) for s in spacing],
        "adjacent_spacing_bootstrap_median_bp": [
            int(fnum(s["bootstrap_spacing_median_bp"])) for s in spacing
        ],
    }


def chromatin_flips() -> dict[str, object]:
    out: dict[str, object] = {}
    for name in (
        "relative_enrichment_group_stats",
        "absolute_burden_group_stats",
        "observed_fraction_group_stats",
    ):
        rs = rows(T / "chromatin" / f"{name}.tsv")
        sig = sum(1 for r in rs if r["bh_fdr_all_tests"] not in ("", "NA")
                  and fnum(r["bh_fdr_all_tests"]) < 0.05)
        out[name] = {"tests": len(rs), "significant_q_lt_0.05": sig}
    return out


def locus_selection() -> dict[str, object]:
    rank = rows(T / "BCL3" / "chromatin_display_candidate_ranking.tsv")
    selected = {
        r["Category_for_row"]: r["Gene"]
        for r in rank
        if r["Selected_for_Figure_5_display"] == "Yes"
        and r["Strict_selection_eligible"] == "Yes"
    }
    return {
        "candidate_category_assignments": len(rank),
        "display_selected": selected,
    }


def validation_targets() -> dict[str, object]:
    path = T / "chromatin" / "validation_target_support_summary.tsv"
    if not path.is_file():
        path = T / "validation_target_support_summary.tsv"
    if not path.is_file():
        return {"note": "pending S3 core robustness run"}
    rs = [r for r in rows(path)
          if r["callset"] == "circlemap_methods" and r["tolerance_bp"] == "10"]
    return {
        r["target"] + "_" + r["group"]: f"{r['supported_samples']}/{r['n_samples']}"
        for r in rs
    }


# Values as currently printed in Submit/Maintext_revised_clean.docx.
V1 = {
    "tested_genes": 26642,
    "differential_eccgenes": 3402,
    "differential_up_in_covid": 3094,
    "peaks_bp": [196, 365, 571, 760],
    "adjacent_spacing_bp": [169, 206],
    "candidate_category_assignments": 281,
    "display_selected": {
        "Humoral immune response": "HLA-E",
        "Endothelial/coagulation": "MCAM",
    },
}


def main() -> int:
    report = {
        "eccGene": eccgene_summary(),
        "fragment_length": fragment_peaks(),
        "chromatin": chromatin_flips(),
        "locus_selection": locus_selection(),
        "validation_targets": validation_targets(),
    }
    out = ROOT / "logs" / "v2_key_numbers.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    flat = {**report["eccGene"], **report["fragment_length"], **report["locus_selection"]}
    print(f"{'quantity':40s} {'v1 (submitted)':>26s}   {'v2 (methods)':>26s}")
    print("-" * 98)
    for key, v1 in V1.items():
        v2 = flat.get(key)
        flag = "" if v1 == v2 else "   <-- CHANGES"
        print(f"{key:40s} {str(v1):>26s}   {str(v2):>26s}{flag}")

    print()
    print("Fig 3b top-30 (HC-zero, smallest q):")
    print("  " + ", ".join(flat["fig3b_top30"][:10]) + " ...")
    print()
    print(f"Report written: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
