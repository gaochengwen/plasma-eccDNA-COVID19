#!/usr/bin/env python3
"""Every Results number the manuscript prints, computed from the v2 locked tables.

This is the source for the tracked manuscript edits: each entry pairs the value
currently printed in Submit/Maintext_Covid-eccDNA.docx with the v2 value, so the
edit script replaces strings rather than re-deriving statistics inside the
document builder.

Nothing here recomputes an analysis. Group medians, IQRs, means and ranges are
order statistics of the locked per-sample columns; Wilcoxon/Cliff's delta values
are read from the locked group-test tables where they exist.
"""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
T = RUN / "tables"

LENGTH_BINS = ["lt200", "200_399", "400_599", "600_999", "1000_1999", "ge2000"]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def quartiles(values: list[float]) -> tuple[float, float, float]:
    s = sorted(values)
    n = len(s)
    return (statistics.median(s[: n // 2]), statistics.median(s),
            statistics.median(s[(n + 1) // 2:]))


def burden() -> dict:
    rows = read(T / "robustness" / "robustness_callset_sample_metrics.tsv")
    out = {}
    for callset, tag in (("circlemap_current", "v1"), ("circlemap_methods", "v2")):
        sub = [r for r in rows if r["callset"] == callset]
        block = {}
        for group in ("COVID-19", "HC"):
            g = [r for r in sub if r["group"] == group]
            counts = [float(r["call_count"]) for r in g]
            epms = [float(r["epm"]) for r in g]
            clo, cmed, chi = quartiles(counts)
            elo, emed, ehi = quartiles(epms)
            fractions = {}
            for k in LENGTH_BINS:
                vals = [int(r[f"length_{k}_count"]) for r in g
                        if r[f"length_{k}_count"] not in ("", "nan")]
                totals = [sum(int(r[f"length_{b}_count"]) for b in LENGTH_BINS)
                          for r in g if r[f"length_{k}_count"] not in ("", "nan")]
                if vals:
                    per = [v / t for v, t in zip(vals, totals) if t]
                    lo, med, hi = quartiles(per)
                    fractions[k] = {"median": med, "q1": lo, "q3": hi}
            block[group] = {
                "n": len(g), "total": int(sum(counts)),
                "mean_count": sum(counts) / len(counts),
                "min_count": int(min(counts)), "max_count": int(max(counts)),
                "count_median": cmed, "count_q1": clo, "count_q3": chi,
                "epm_median": emed, "epm_q1": elo, "epm_q3": ehi,
                "length_fractions": fractions,
            }
        block["total_all"] = block["COVID-19"]["total"] + block["HC"]["total"]
        out[tag] = block
    return out


def group_tests() -> dict:
    rows = read(T / "robustness" / "robustness_group_tests.tsv")
    out = {}
    for r in rows:
        out[f"{r['callset']}|{r['metric']}"] = {
            "cliffs_delta": float(r["cliffs_delta_COVID_vs_HC"]),
            "p": float(r["p_value_two_sided"]),
            "q": float(r["q_value_bh_within_metric"]),
            "median_ratio": float(r["median_ratio_COVID_over_HC"]),
        }
    return out


def regression() -> dict:
    rows = [r for r in read(T / "robustness" / "robustness_regression_HC3.tsv")
            if r["term"] == "COVID-19_vs_HC"]
    return {r["callset"]: {
        "fold": float(r["multiplicative_effect"]),
        "ci_lo": float(r["multiplicative_ci_low"]),
        "ci_hi": float(r["multiplicative_ci_high"]),
    } for r in rows}


def fragment() -> dict:
    peaks = read(T / "fragment_length" /
                 "objective_peak_identification_and_stability.tsv")
    spacing = read(T / "fragment_length" / "adjacent_peak_spacing.tsv")
    return {
        "peaks_bp": [int(float(p["primary_peak_position_bp"])) for p in peaks],
        "hc_bootstrap": [round(float(p["hc_bootstrap_detection_frequency"]), 3)
                         for p in peaks],
        "spacing_bp": [int(float(s["primary_spacing_bp"])) for s in spacing],
        "spacing_bootstrap_median": [int(float(s["bootstrap_spacing_median_bp"]))
                                     for s in spacing],
    }


def eccgene() -> dict:
    ab = read(T / "eccGene" / "junction_abundance_wilcoxon.tsv")

    def sig(r):
        q, lfc = r["wilcoxon_q_BH"], r["log2FC_mean_EA_COVID_vs_HC"]
        if q in ("", "NA") or lfc in ("", "NA", "inf", "-inf"):
            return False
        return float(q) < 0.05 and abs(float(lfc)) >= 1

    s = [r for r in ab if sig(r)]
    up = [r for r in s if float(r["log2FC_mean_EA_COVID_vs_HC"]) > 0]
    det = read(T / "eccGene" / "junction_detection_fisher.tsv")
    hc0 = [r for r in det if r["hc_detected"] == "0"]
    hc0.sort(key=lambda r: (float(r["fisher_q_BH"]), -int(r["covid_detected"])))
    return {"tested": len(ab), "significant": len(s), "up": len(up),
            "top30": [r["Gene"] for r in hc0[:30]]}


def loci() -> dict:
    rank = read(T / "BCL3" / "chromatin_display_candidate_ranking.tsv")
    selected = {r["Category_for_row"]: r["Gene"] for r in rank
                if r["Selected_for_Figure_5_display"] == "Yes"}
    matrix = T / "recurrent_intervals" / "covid_specific_exact.matrix.min10.tsv"
    with matrix.open() as fh:
        n_intervals = sum(1 for _ in fh) - 1
    return {"candidate_assignments": len(rank), "display_selected": selected,
            "recurrent_min10_intervals": n_intervals}


def targets() -> dict:
    rows = [r for r in read(T / "robustness" /
                            "validation_target_support_summary.tsv")
            if r["tolerance_bp"] == "10" and r["group"] == "COVID-19"]
    out = {}
    for r in rows:
        out.setdefault(r["target"], {})[r["callset"]] = \
            f"{r['supported_samples']}/{r['n_samples']}"
    return out


def main() -> int:
    report = {
        "burden": burden(), "group_tests": group_tests(),
        "regression": regression(), "fragment_length": fragment(),
        "eccGene": eccgene(), "loci": loci(), "validation_targets": targets(),
    }
    (RUN / "logs" / "v2_manuscript_numbers.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    b = report["burden"]
    print("Results paragraph 1 — burden\n")
    print(f"{'quantity':34s} {'v1 (printed)':>22s} {'v2':>22s}")
    print("-" * 80)
    rows = [
        ("total eccDNA calls", f"{b['v1']['total_all']:,}", f"{b['v2']['total_all']:,}"),
        ("COVID-19 calls", f"{b['v1']['COVID-19']['total']:,}", f"{b['v2']['COVID-19']['total']:,}"),
        ("HC calls", f"{b['v1']['HC']['total']:,}", f"{b['v2']['HC']['total']:,}"),
        ("COVID mean per sample", f"{b['v1']['COVID-19']['mean_count']:,.0f}", f"{b['v2']['COVID-19']['mean_count']:,.0f}"),
        ("COVID range", f"{b['v1']['COVID-19']['min_count']:,}-{b['v1']['COVID-19']['max_count']:,}",
                        f"{b['v2']['COVID-19']['min_count']:,}-{b['v2']['COVID-19']['max_count']:,}"),
        ("HC mean per sample", f"{b['v1']['HC']['mean_count']:,.0f}", f"{b['v2']['HC']['mean_count']:,.0f}"),
        ("HC range", f"{b['v1']['HC']['min_count']:,}-{b['v1']['HC']['max_count']:,}",
                     f"{b['v2']['HC']['min_count']:,}-{b['v2']['HC']['max_count']:,}"),
        ("COVID count median [IQR]", f"{b['v1']['COVID-19']['count_median']:,.0f}", f"{b['v2']['COVID-19']['count_median']:,.0f}"),
        ("HC count median [IQR]", f"{b['v1']['HC']['count_median']:,.0f}", f"{b['v2']['HC']['count_median']:,.0f}"),
        ("COVID EPM median", f"{b['v1']['COVID-19']['epm_median']:,.0f}", f"{b['v2']['COVID-19']['epm_median']:,.0f}"),
        ("HC EPM median", f"{b['v1']['HC']['epm_median']:,.0f}", f"{b['v2']['HC']['epm_median']:,.0f}"),
    ]
    for name, v1, v2 in rows:
        print(f"{name:34s} {v1:>22s} {v2:>22s}")

    gt, rg = report["group_tests"], report["regression"]
    print()
    for metric in ("call_count", "epm"):
        a = gt[f"circlemap_current|{metric}"]
        c = gt[f"circlemap_methods|{metric}"]
        print(f"{metric:12s} delta {a['cliffs_delta']:.3f} -> {c['cliffs_delta']:.3f}   "
              f"q {a['q']:.2e} -> {c['q']:.2e}")
    print(f"{'adj fold':12s} {rg['circlemap_current']['fold']:.2f} "
          f"({rg['circlemap_current']['ci_lo']:.2f}-{rg['circlemap_current']['ci_hi']:.2f})"
          f" -> {rg['circlemap_methods']['fold']:.2f} "
          f"({rg['circlemap_methods']['ci_lo']:.2f}-{rg['circlemap_methods']['ci_hi']:.2f})")

    e, l, f = report["eccGene"], report["loci"], report["fragment_length"]
    print(f"\neccGene tested {e['tested']:,}  significant {e['significant']:,}  up {e['up']:,}")
    print(f"peaks {f['peaks_bp']}  spacing {f['spacing_bp']}")
    print(f"candidate assignments {l['candidate_assignments']}  "
          f"display {l['display_selected']}  recurrent min10 {l['recurrent_min10_intervals']}")
    print(f"\nwritten: logs/v2_manuscript_numbers.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
