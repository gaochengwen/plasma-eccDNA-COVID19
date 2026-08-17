#!/usr/bin/env python3
"""Rebuild Supplementary Table S3 (cohort comparisons) on the v2 call set.

Table S3 is the master statistics table: eight named FDR families, each
BH-corrected within itself. The test is the two-sided Wilcoxon rank-sum with
the tie-corrected normal approximation and no continuity correction -- an
implementation confirmed by reproducing the v1 sample_burden p values exactly
(5.7283e-11 and 5.1847e-09) before any v2 value was computed.

Family provenance under v2:

  sample_burden, sample_burden_log1p   locked per-sample metrics
  chromosome_distribution              recomputed to Methods (amendment A4)
  repeat_class_enrichment              CARRIED OVER FROM v1 UNCHANGED -- the
                                       original script works from BAM reads,
                                       which are identical in v1 and v2, so the
                                       call-set migration does not touch it
  gene_feature_*                       recomputed to Methods (amendment A4)
  fragment_length_bin_proportion       locked per-sample proportions

The three gene_feature families are three views of one set of counts
(observed fraction, O/E ratio, log2 O/E). Because log2 is monotone, Cliff's
delta and the p value are identical across them by construction; v1 reported
slightly different values for the enrichment-score family, which cannot arise
from the same counts and is noted in the QA document.
"""

from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path

import numpy as np
from openpyxl import load_workbook
from scipy import stats

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
T = RUN / "tables"
V1_BOOK = PROJECT / "Submit" / "Supplementary_Tables_revised.xlsx"
OUT = RUN / "tables" / "supplementary" / "Table_S3_v2.tsv"

TEST = ("two-sided Wilcoxon rank-sum (Mann-Whitney U, tie-corrected normal "
        "approximation)")
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ELEMENTS = ["3UTR", "5UTR", "CpG_islands", "Gene2KbD", "Gene2KbU", "exon", "intron"]
BINS = [("< 1 kb", ("lt200", "200_399", "400_599", "600_999")),
        ("1-2 kb", ("1000_1999",)), ("> 2 kb", ("ge2000",))]

HEADER = ["Family", "FDR family", "Source", "Metric", "Feature", "Unit",
          "Transform", "Test", "COVID-19 n", "HC n", "COVID-19 median",
          "COVID-19 Q1", "COVID-19 Q3", "HC median", "HC Q1", "HC Q3",
          "Median difference", "Cliff's delta", "Raw P value",
          "BH-adjusted q value", "Direction", "Notes"]


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def quart(v):
    s = sorted(v)
    n = len(s)
    return statistics.median(s[: n // 2]), statistics.median(s), statistics.median(s[(n + 1) // 2:])


def cliffs(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(((a[:, None] > b[None, :]).sum() - (a[:, None] < b[None, :]).sum())
                 / (len(a) * len(b)))


def bh(ps):
    ps = np.asarray(ps, float)
    n = len(ps)
    order = np.argsort(ps)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        prev = min(prev, ps[i] * n / (n - rank + 1))
        q[i] = prev
    return q


def row(family, fdr_family, source, metric, feature, unit, transform,
        covid, hc, note="") -> dict:
    c1, cm, c3 = quart(covid)
    h1, hm, h3 = quart(hc)
    d = cliffs(covid, hc)
    p = float(stats.ranksums(covid, hc).pvalue)
    return {"Family": family, "FDR family": fdr_family, "Source": source,
            "Metric": metric, "Feature": feature, "Unit": unit,
            "Transform": transform, "Test": TEST,
            "COVID-19 n": len(covid), "HC n": len(hc),
            "COVID-19 median": cm, "COVID-19 Q1": c1, "COVID-19 Q3": c3,
            "HC median": hm, "HC Q1": h1, "HC Q3": h3,
            "Median difference": cm - hm, "Cliff's delta": d,
            "Raw P value": p, "BH-adjusted q value": None,
            "Direction": "COVID-19 higher" if d > 0 else "HC higher",
            "Notes": note}


def main() -> int:
    rows: list[dict] = []

    # --- sample burden -----------------------------------------------------
    m = [r for r in read(T / "robustness" / "robustness_callset_sample_metrics.tsv")
         if r["callset"] == "circlemap_methods"]
    cov = lambda k, g: [float(r[k]) for r in m if r["group"] == g]  # noqa: E731
    rows.append(row("sample_burden", "sample_burden", "v2 methods call set",
                    "total_eccDNA_count", "all_samples", "count", "none",
                    cov("call_count", "COVID-19"), cov("call_count", "HC")))
    rows.append(row("sample_burden", "sample_burden", "v2 methods call set",
                    "total_EPM", "all_samples", "eccDNA per million mapped reads",
                    "none", cov("epm", "COVID-19"), cov("epm", "HC")))
    rows.append(row("sample_burden", "sample_burden_log1p", "v2 methods call set",
                    "total_EPM", "all_samples", "log1p(EPM)", "log1p",
                    [math.log1p(x) for x in cov("epm", "COVID-19")],
                    [math.log1p(x) for x in cov("epm", "HC")]))

    # --- chromosome distribution ------------------------------------------
    chrom = read(T / "genomic_distribution" / "chromosome_distribution.tsv")
    for c in CHROMS:
        sub = [r for r in chrom if r["chromosome"] == c]
        rows.append(row("chromosome_distribution", "chromosome_distribution",
                        "v2 methods call set, reimplemented (A4)",
                        "chromosome_normalized_fraction_per_mb", c,
                        "fraction per Mb", "none",
                        [float(r["normalized_fraction_per_mb"]) for r in sub
                         if r["group"] == "COVID-19"],
                        [float(r["normalized_fraction_per_mb"]) for r in sub
                         if r["group"] == "HC"]))

    # --- repeat classes: unchanged, carried over from v1 -------------------
    wb = load_workbook(V1_BOOK, read_only=True)
    ws = wb["Table_S3"]
    v1 = list(ws.iter_rows(values_only=True))
    hdr = [str(c) if c is not None else "" for c in v1[1]]
    carried = 0
    for r in v1[2:]:
        if r[1] != "repeat_class_enrichment":
            continue
        rows.append({h: r[i] for i, h in enumerate(hdr) if h in HEADER}
                    | {"Notes": "unchanged: computed from BAM reads, which are "
                                "identical in v1 and v2"})
        carried += 1

    # --- gene features: three views of one count set -----------------------
    gf = read(T / "genomic_distribution" / "gene_feature_enrichment_scores.tsv")
    views = [("gene_feature_element_ratio", "gene_feature_eccDNA_element_ratio",
              "fraction", "none", lambda r: float(r["observed_fraction"])),
             ("gene_feature_normalized_ratio", "gene_feature_normalized_eccDNA_ratio",
              "observed/expected", "none", lambda r: float(r["oe_ratio"])),
             ("gene_feature_enrichment_score", "gene_feature_enrichment_score",
              "log2(observed/expected)", "log2",
              lambda r: math.log2(float(r["oe_ratio"])))]
    for fam, metric, unit, transform, get in views:
        for el in ELEMENTS:
            sub = [r for r in gf if r["element"] == el]
            rows.append(row("gene_feature", fam,
                            "v2 methods call set, reimplemented (A4)",
                            metric, el, unit, transform,
                            [get(r) for r in sub if r["group"] == "COVID-19"],
                            [get(r) for r in sub if r["group"] == "HC"]))

    # --- fragment-length bins ---------------------------------------------
    fs = [r for r in read(T / "robustness" /
                          "robustness_fragment_size_sample_proportions.tsv")
          if r["callset"] == "circlemap_methods"]
    by: dict[tuple[str, str], dict[str, float]] = {}
    for r in fs:
        by.setdefault((r["sample_id"], r["group"]), {})[r["length_bin"]] = \
            float(r["proportion"])
    for label, keys in BINS:
        rows.append(row("fragment_length_bin_proportion",
                        "fragment_length_bin_proportion", "v2 methods call set",
                        "fragment_length_bin_proportion", label, "proportion",
                        "none",
                        [sum(v[k] for k in keys) for (s, g), v in by.items()
                         if g == "COVID-19"],
                        [sum(v[k] for k in keys) for (s, g), v in by.items()
                         if g == "HC"]))

    # --- BH within each named family --------------------------------------
    for fam in {r["FDR family"] for r in rows}:
        members = [r for r in rows if r["FDR family"] == fam
                   and r["BH-adjusted q value"] is None]
        if not members:
            continue
        for r, q in zip(members, bh([r["Raw P value"] for r in members])):
            r["BH-adjusted q value"] = float(q)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT.relative_to(RUN)}: {len(rows)} rows "
          f"({carried} carried over unchanged)")
    print(f"\n{'FDR family':34s} {'n':>3s} {'q<0.05':>7s}")
    for fam in sorted({r["FDR family"] for r in rows}):
        members = [r for r in rows if r["FDR family"] == fam]
        sig = sum(1 for r in members
                  if r["BH-adjusted q value"] is not None
                  and float(r["BH-adjusted q value"]) < 0.05)
        print(f"{fam:34s} {len(members):>3d} {sig:>7d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
