#!/usr/bin/env python3
"""Generate the final human-readable robustness report and checksum manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


CALLSET_LABELS = {
    "circlemap_current": "Current Circle-Map",
    "circlemap_methods": "Methods-consistent Circle-Map",
    "circlemap_strict": "Circle-Map split≥10, score>1,000",
    "circlemap_very_strict": "Circle-Map split≥20, score>2,000",
    "circlemap_artifact_masked": "Circle-Map artifact-masked",
    "circlemap_high_support_masked": "Circle-Map high-support + artifact-masked",
    "consensus_t10": "Circle-Map/Circle_finder consensus (10 bp)",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "NA"
    if value != 0 and (abs(value) < 0.001 or abs(value) >= 10_000):
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def validate_required(root: Path) -> None:
    required = (
        root / "circlefinder_validation_all.txt",
        root / "consensus_validation_all.txt",
        root / "core_robustness_validation.txt",
        root / "artifact_filter_validation.txt",
        root / "downstream_summary_validation.txt",
        root / "downstream" / "circlemap_methods" / "downstream_validation.txt",
        root
        / "downstream"
        / "circlemap_high_support_masked"
        / "downstream_validation.txt",
        root / "downstream" / "consensus_t10" / "downstream_validation.txt",
    )
    failures = []
    for path in required:
        if not path.is_file() or path.read_text(errors="replace").splitlines()[0] != "PASS":
            failures.append(str(path))
    if failures:
        raise RuntimeError(f"Required validation files absent or not PASS: {failures}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    validate_required(root)

    group_tests = read_tsv(root / "tables" / "robustness_group_tests.tsv")
    regressions = read_tsv(root / "tables" / "robustness_regression_HC3.tsv")
    concordance = read_tsv(root / "tables" / "caller_concordance_all.tsv")
    target_summary = read_tsv(root / "tables" / "validation_target_support_summary.tsv")
    target_masks = read_tsv(root / "tables" / "validation_target_mask_audit.tsv")
    eccgene_candidates = read_tsv(root / "tables" / "eccgene_candidate_stability.tsv")
    eccgene_concordance = read_tsv(root / "tables" / "eccgene_callset_concordance.tsv")
    chromatin_concordance = read_tsv(root / "tables" / "chromatin_callset_concordance.tsv")

    group_effects = {
        row["callset"]: row
        for row in regressions
        if row["term"] == "COVID-19_vs_HC"
    }
    epm_tests = {
        row["callset"]: row for row in group_tests if row["metric"] == "epm"
    }
    callset_rows = []
    for callset in CALLSET_LABELS:
        if callset not in epm_tests or callset not in group_effects:
            continue
        test = epm_tests[callset]
        regression = group_effects[callset]
        callset_rows.append(
            [
                CALLSET_LABELS[callset],
                fmt(float(test["covid_median"]), 1),
                fmt(float(test["hc_median"]), 1),
                fmt(float(test["median_ratio_COVID_over_HC"]), 2),
                fmt(float(test["p_value_two_sided"]), 2),
                fmt(float(test["cliffs_delta_COVID_vs_HC"]), 2),
                (
                    f"{fmt(float(regression['multiplicative_effect']), 2)} "
                    f"({fmt(float(regression['multiplicative_ci_low']), 2)}–"
                    f"{fmt(float(regression['multiplicative_ci_high']), 2)})"
                ),
                fmt(float(regression["p_value_two_sided"]), 2),
            ]
        )

    concordance_rows = []
    for tolerance in (2, 5, 10, 20):
        selected = [row for row in concordance if int(row["tolerance_bp"]) == tolerance]
        concordance_rows.append(
            [
                str(tolerance),
                fmt(
                    statistics.median(
                        float(row["circlemap_retention_fraction"]) for row in selected
                    ),
                    3,
                ),
                fmt(
                    statistics.median(
                        float(row["circlefinder_retention_fraction"]) for row in selected
                    ),
                    3,
                ),
                fmt(
                    statistics.median(float(row["jaccard_index"]) for row in selected),
                    3,
                ),
                f"{sum(int(row['consensus_count']) for row in selected):,}",
            ]
        )

    target_rows = []
    for target in ("CTNNA2", "CAB39", "chr22_target", "SDK1", "TCF7L1"):
        values = []
        for callset in (
            "circlemap_methods",
            "circlemap_high_support_masked",
            "consensus_t10",
        ):
            selected = [
                row
                for row in target_summary
                if row["target"] == target
                and row["callset"] == callset
                and row["group"] == "COVID-19"
                and int(row["tolerance_bp"]) == 10
            ]
            values.append(
                f"{selected[0]['supported_samples']}/39" if selected else "NA"
            )
        target_rows.append([target, *values])

    mask_rows = []
    for target in ("CTNNA2", "CAB39", "chr22_target", "SDK1", "TCF7L1"):
        selected = [row for row in target_masks if row["target"] == target]
        minimum_umap = min(
            float(row["umap_k100_unique_overlap_fraction"]) for row in selected
        )
        repeats = sorted(
            {
                annotation.split(":")[0]
                for row in selected
                for annotation in row["repeatmasker_annotations"].split(";")
                if annotation
            }
        )
        mask_rows.append(
            [
                target,
                fmt(minimum_umap, 2),
                "; ".join(repeats) if repeats else "None",
                (
                    "Pass"
                    if minimum_umap >= 0.90 and not repeats
                    else "Fail strict breakpoint mask"
                ),
            ]
        )

    candidate_rows = []
    for gene in ("CTNNA2", "CAB39", "SDK1", "TCF7L1"):
        for callset in (
            "circlemap_methods",
            "circlemap_high_support_masked",
            "consensus_t10",
        ):
            selected = [
                row
                for row in eccgene_candidates
                if row["gene"] == gene
                and row["callset"] == callset
                and row["definition"] == "junction"
            ]
            if not selected:
                continue
            row = selected[0]
            candidate_rows.append(
                [
                    gene,
                    CALLSET_LABELS[callset],
                    fmt(float(row["log2FC_mean_EA_COVID_vs_HC"]), 2),
                    f"{row['covid_detected']}/39",
                    f"{row['hc_detected']}/39",
                    fmt(float(row["wilcoxon_q_BH"]), 2),
                    row["direction"],
                ]
            )

    eccgene_consensus = [
        row
        for row in eccgene_concordance
        if row["comparison_callset"] == "consensus_t10"
    ]
    chromatin_consensus = [
        row
        for row in chromatin_concordance
        if row["comparison_callset"] == "consensus_t10"
    ]
    chromatin_direction = (
        sum(int(row["direction_concordant"]) for row in chromatin_consensus)
        / len(chromatin_consensus)
        if chromatin_consensus
        else math.nan
    )
    chromatin_significance = (
        sum(int(row["significance_concordant"]) for row in chromatin_consensus)
        / len(chromatin_consensus)
        if chromatin_consensus
        else math.nan
    )

    report = f"""# Circle-Map robustness analysis

Generated: {datetime.now(timezone.utc).isoformat()}

## Executive conclusion

The COVID-19 versus HC eccDNA-burden conclusion was evaluated under the
current Circle-Map callset, the manuscript methods filter, progressively
stricter split-read/score thresholds, a stringent breakpoint artifact mask,
and an independent two-caller consensus. The primary normalized outcome is
circles per million mapped alignments (EPM). All effect estimates below are
computed at sample level (39 COVID-19 and 39 HC); adjusted models include
log2 RCA concentration, age, and sex with HC3 standard errors.

{markdown_table(
    ["Callset", "COVID median EPM", "HC median EPM", "Median ratio", "MWU P", "Cliff's δ", "Adjusted fold (95% CI)", "Adjusted P"],
    callset_rows,
)}

## Caller concordance

Circle_finder was run independently from paired FASTQ records reconstructed
from the complete BWA BAM records. Consensus requires the same chromosome and
both breakpoints within the stated tolerance; calls are matched one-to-one.
The prespecified primary tolerance is 10 bp.

{markdown_table(
    ["Tolerance (bp)", "Median CM retention", "Median CF retention", "Median Jaccard", "Consensus calls"],
    concordance_rows,
)}

## Four validation targets

Support is the number of COVID-19 samples with both breakpoints within 10 bp
of the fixed target. No HC sample supported any of the four targets in the
methods-consistent Circle-Map callset.

{markdown_table(
    ["Target", "Methods CM", "High-support + masked", "Two-caller consensus"],
    target_rows,
)}

The strict mask audit explains whether target loss is caused by generic
breakpoint artifacts. Both 200-bp breakpoint windows must have at least 90%
Umap k=100 unique coverage and no RepeatMasker/segmental-duplication/
blacklist-gap overlap.

{markdown_table(
    ["Target", "Minimum Umap fraction", "Breakpoint repeats", "Strict-mask status"],
    mask_rows,
)}

## eccGene stability

Candidate results below use the prespecified junction gene assignment. The
complete genome-wide abundance/detection tables for junction, interval, and
midpoint definitions are under `downstream/<callset>/eccGene/`.

{markdown_table(
    ["Gene", "Callset", "log2FC", "COVID detected", "HC detected", "BH q", "Direction"],
    candidate_rows,
)}

Across the three eccGene definitions, methods-Circle-Map versus consensus
log2FC Spearman correlations were:
{", ".join(f"{row['definition']}={fmt(float(row['log2FC_spearman_rho']), 3)}" for row in eccgene_consensus)}.

## Chromatin stability

Across {len(chromatin_consensus)} matched chromatin tests, methods-Circle-Map
and consensus agreed in direction for {fmt(chromatin_direction * 100, 1)}%
and in BH-significance status for {fmt(chromatin_significance * 100, 1)}%.
Full absolute and relative enrichment results are under
`downstream/<callset>/chromatin/tables/`.

## Reproducibility and interpretation

- Circle_finder upstream commit:
  `3eb333db2ea6277dde36cbf640be9afeb710c717`.
- Circle-Map methods filter: split reads ≥3, score >200, mean coverage > SD,
  start/end coverage increase ≥0.3, uncovered fraction <0.1.
- Artifact mask: no full-circle or breakpoint-window blacklist/gap overlap;
  each ±100-bp breakpoint window has Umap k=100 unique fraction ≥0.90 and no
  segmental-duplication or RepeatMasker overlap. Internal repeats are retained.
- High-support masked: artifact mask plus split reads ≥10 and score >1,000.
- The benchmarked caller and implementation provenance are documented in the
  [Circle_finder repository](https://github.com/pk7zuva/Circle_finder) and a
  [2024 eccDNA caller benchmark](https://www.nature.com/articles/s41467-024-53496-8).
- Consensus improves specificity but is not a gold standard: both callers use
  short-read evidence and BWA-family alignment, and Circle_finder FASTQ was
  reconstructed from BAM rather than retained as the original FASTQ files.
"""
    (root / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

    candidates = []
    for directory in (root / "tables", root / "figures", root / "downstream"):
        for path in directory.rglob("*"):
            if path.is_file() and (
                path.suffix in {".tsv", ".txt", ".json", ".pdf", ".svg", ".png"}
                or path.name.endswith(".tsv.gz")
            ):
                candidates.append(path)
    for path in root.glob("*"):
        if path.is_file() and path.suffix in {".md", ".txt", ".json"}:
            candidates.append(path)
    unique = sorted(set(candidates))
    manifest = root / "final_output_manifest.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "sha256"])
        for path in unique:
            writer.writerow([path.relative_to(root), path.stat().st_size, sha256(path)])

    with (root / "final_analysis_validation.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write(f"manifested_files={len(unique)}\n")
        handle.write(f"completed_utc={datetime.now(timezone.utc).isoformat()}\n")


if __name__ == "__main__":
    main()
