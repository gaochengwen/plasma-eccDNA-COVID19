#!/usr/bin/env python3
"""Emit v1-to-v2 source data for the data panels of the Illustrator main figures.

Figures 1, 2 and 3 were assembled in Adobe Illustrator and have no generator in
this repository, so their data panels cannot be re-rendered by script. What can
be delivered exactly is the numbers behind each panel, v1 beside v2, so whoever
holds the .ai files can replace the panels without re-deriving anything.

Covered here (all inputs present in the v2 tree):

  Figure 1b  per-sample total eccDNA counts
  Figure 1c  per-sample EPM
  Figure 1e  per-sample length-bin composition
  Figure 3a  eccGenes passing BH q < 0.05 and |log2FC| >= 1
  Figure 3b  top-30 HC-zero detection genes

Not covered, because the code that produced their inputs is not in this
repository (see docs/S4_figure_regeneration_QA.md):

  Figure 1a  workflow schematic, carries no data
  Figure 2   chromosome / repeat-class / gene-element distributions
  Figure 3c  Metascape terms (external service)
  Figure 3d  recurrent exact intervals
  Figure 3e  conceptual model, carries no data
"""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
OUT = RUN / "figures" / "illustrator_panel_data"

LENGTH_BINS = ["lt200", "200_399", "400_599", "600_999", "1000_1999", "ge2000"]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def iqr(values: list[float]) -> tuple[float, float, float]:
    s = sorted(values)
    n = len(s)
    med = statistics.median(s)
    lo = statistics.median(s[: n // 2])
    hi = statistics.median(s[(n + 1) // 2:])
    return lo, med, hi


def figure1(out: Path) -> dict[str, object]:
    v2 = [r for r in read(RUN / "tables" / "robustness" /
                          "robustness_callset_sample_metrics.tsv")
          if r["callset"] == "circlemap_methods"]
    v1 = [r for r in read(RUN / "tables" / "robustness" /
                          "robustness_callset_sample_metrics.tsv")
          if r["callset"] == "circlemap_current"]
    by = {"v1": {r["sample_id"]: r for r in v1},
          "v2": {r["sample_id"]: r for r in v2}}

    rows = []
    for sample in sorted(by["v2"]):
        a, b = by["v1"][sample], by["v2"][sample]
        row = {
            "sample_id": sample,
            "group": b["group"],
            "count_v1": a["call_count"], "count_v2": b["call_count"],
            "epm_v1": f"{float(a['epm']):.1f}", "epm_v2": f"{float(b['epm']):.1f}",
        }
        for tag, src in (("v1", a), ("v2", b)):
            total = sum(int(src[f"length_{k}_count"]) for k in LENGTH_BINS
                        if src[f"length_{k}_count"] not in ("", "nan"))
            for k in LENGTH_BINS:
                raw = src[f"length_{k}_count"]
                row[f"frac_{k}_{tag}"] = (
                    f"{int(raw) / total:.6f}" if raw not in ("", "nan") and total else "NA")
        rows.append(row)

    path = out / "Figure_1bce_per_sample.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    summary = {}
    for tag in ("v1", "v2"):
        for group in ("COVID-19", "HC"):
            counts = [float(r[f"count_{tag}"]) for r in rows if r["group"] == group]
            epms = [float(r[f"epm_{tag}"]) for r in rows if r["group"] == group]
            lo, med, hi = iqr(epms)
            summary[f"{group}_{tag}"] = {
                "n": len(counts),
                "total_calls": int(sum(counts)),
                "count_median": int(statistics.median(counts)),
                "epm_median": round(med, 1),
                "epm_iqr": [round(lo, 1), round(hi, 1)],
            }
    return {"per_sample_table": path.name, "summary": summary}


def figure3(out: Path) -> dict[str, object]:
    ab = {r["Gene"]: r for r in read(RUN / "tables" / "eccGene" /
                                     "junction_abundance_wilcoxon.tsv")}

    def sig(r: dict[str, str]) -> bool:
        q, lfc = r["wilcoxon_q_BH"], r["log2FC_mean_EA_COVID_vs_HC"]
        if q in ("", "NA") or lfc in ("", "NA", "inf", "-inf"):
            return False
        return float(q) < 0.05 and abs(float(lfc)) >= 1

    genes = sorted(g for g, r in ab.items() if sig(r))
    path_a = out / "Figure_3a_significant_eccgenes.tsv"
    with path_a.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Gene", "log2FC", "cliffs_delta", "wilcoxon_q_BH",
                    "covid_detected", "hc_detected"])
        for g in genes:
            r = ab[g]
            w.writerow([g, r["log2FC_mean_EA_COVID_vs_HC"],
                        r["cliffs_delta_COVID_vs_HC"], r["wilcoxon_q_BH"],
                        r["covid_detected"], r["hc_detected"]])

    det = read(RUN / "tables" / "eccGene" / "junction_detection_fisher.tsv")
    hc_zero = [r for r in det if r["hc_detected"] == "0"]
    hc_zero.sort(key=lambda r: (float(r["fisher_q_BH"]), -int(r["covid_detected"])))
    top = hc_zero[:30]
    # Figure 3b orders rows by decreasing COVID-19 recurrence, then q, then symbol.
    top.sort(key=lambda r: (-int(r["covid_detected"]), float(r["fisher_q_BH"]),
                            r["Gene"]))
    path_b = out / "Figure_3b_top30_detection.tsv"
    with path_b.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rank", "Gene", "covid_detected_n", "hc_detected_n",
                    "fisher_q_BH"])
        for i, r in enumerate(top, 1):
            w.writerow([i, r["Gene"], r["covid_detected"], r["hc_detected"],
                        r["fisher_q_BH"]])

    return {
        "significant_eccgenes_v2": len(genes),
        "significant_eccgenes_v1": 3402,
        "tables": [path_a.name, path_b.name],
        "top30_v2": [r["Gene"] for r in top],
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"figure_1": figure1(OUT), "figure_3": figure3(OUT)}
    (OUT / "panel_data_summary.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n")

    f1 = report["figure_1"]["summary"]
    print("Figure 1b/1c — cohort summary")
    print(f"  {'':10s} {'total calls':>14s} {'median count':>13s} {'median EPM':>11s}")
    for group in ("COVID-19", "HC"):
        for tag in ("v1", "v2"):
            s = f1[f"{group}_{tag}"]
            print(f"  {group + ' ' + tag:10s} {s['total_calls']:>14,d} "
                  f"{s['count_median']:>13,d} {s['epm_median']:>11.1f}")
    print()
    f3 = report["figure_3"]
    print(f"Figure 3a — significant eccGenes: v1 {f3['significant_eccgenes_v1']:,} "
          f"-> v2 {f3['significant_eccgenes_v2']:,}")
    print(f"Figure 3b — top-30 (v2): {', '.join(f3['top30_v2'][:8])} ...")
    print()
    print(f"written to figures/illustrator_panel_data/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
