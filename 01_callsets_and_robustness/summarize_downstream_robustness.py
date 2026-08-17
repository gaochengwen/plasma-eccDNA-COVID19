#!/usr/bin/env python3
"""Summarize eccGene and chromatin stability across robustness callsets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats


CALLSETS = ("circlemap_methods", "circlemap_high_support_masked", "consensus_t10")
DEFINITIONS = ("junction", "interval", "midpoint")
CANDIDATE_GENES = ("CTNNA2", "CAB39", "SDK1", "TCF7L1")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def eccgene_summary(root: Path) -> tuple[list[dict], list[dict], list[dict]]:
    gene_tables: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    summary_rows = []
    candidate_rows = []
    for callset in CALLSETS:
        for definition in DEFINITIONS:
            path = (
                root
                / "downstream"
                / callset
                / "eccGene"
                / "results"
                / f"{definition}_abundance_wilcoxon.tsv"
            )
            rows = read_tsv(path)
            by_gene = {row["Gene"]: row for row in rows}
            gene_tables[(callset, definition)] = by_gene
            significant = [row for row in rows if float(row["wilcoxon_q_BH"]) < 0.05]
            summary_rows.append(
                {
                    "callset": callset,
                    "definition": definition,
                    "tested_genes": len(rows),
                    "significant_q_lt_0_05": len(significant),
                    "significant_COVID_up": sum(row["direction"] == "COVID_up" for row in significant),
                    "significant_HC_up": sum(row["direction"] == "HC_up" for row in significant),
                }
            )
            for gene in CANDIDATE_GENES:
                row = by_gene[gene]
                candidate_rows.append(
                    {
                        "callset": callset,
                        "definition": definition,
                        "gene": gene,
                        "log2FC_mean_EA_COVID_vs_HC": float(
                            row["log2FC_mean_EA_COVID_vs_HC"]
                        ),
                        "cliffs_delta_COVID_vs_HC": float(
                            row["cliffs_delta_COVID_vs_HC"]
                        ),
                        "covid_detected": int(row["covid_detected"]),
                        "hc_detected": int(row["hc_detected"]),
                        "wilcoxon_p": float(row["wilcoxon_p"]),
                        "wilcoxon_q_BH": float(row["wilcoxon_q_BH"]),
                        "direction": row["direction"],
                    }
                )

    concordance_rows = []
    reference = "circlemap_methods"
    for callset in CALLSETS:
        if callset == reference:
            continue
        for definition in DEFINITIONS:
            ref = gene_tables[(reference, definition)]
            comparison = gene_tables[(callset, definition)]
            common = sorted(set(ref) & set(comparison))
            ref_significant = {
                gene for gene in common if float(ref[gene]["wilcoxon_q_BH"]) < 0.05
            }
            comparison_significant = {
                gene
                for gene in common
                if float(comparison[gene]["wilcoxon_q_BH"]) < 0.05
            }
            union = ref_significant | comparison_significant
            rho = stats.spearmanr(
                [float(ref[gene]["log2FC_mean_EA_COVID_vs_HC"]) for gene in common],
                [
                    float(comparison[gene]["log2FC_mean_EA_COVID_vs_HC"])
                    for gene in common
                ],
            )
            concordance_rows.append(
                {
                    "reference_callset": reference,
                    "comparison_callset": callset,
                    "definition": definition,
                    "common_tested_genes": len(common),
                    "reference_significant_genes": len(ref_significant),
                    "comparison_significant_genes": len(comparison_significant),
                    "significant_gene_intersection": len(
                        ref_significant & comparison_significant
                    ),
                    "significant_gene_union": len(union),
                    "significant_gene_jaccard": (
                        len(ref_significant & comparison_significant) / len(union)
                        if union
                        else math.nan
                    ),
                    "log2FC_spearman_rho": float(rho.statistic),
                    "log2FC_spearman_p": float(rho.pvalue),
                }
            )
    return summary_rows, candidate_rows, concordance_rows


def chromatin_summary(root: Path) -> tuple[list[dict], list[dict]]:
    combined_rows = []
    for callset in CALLSETS:
        for analysis, filename in (
            ("relative", "relative_enrichment_group_stats.tsv"),
            ("absolute", "absolute_burden_group_stats.tsv"),
        ):
            path = root / "downstream" / callset / "chromatin" / "tables" / filename
            for row in read_tsv(path):
                combined_rows.append(
                    {
                        "callset": callset,
                        "analysis": analysis,
                        "peak_set_id": row["peak_set_id"],
                        "dataset": row["dataset"],
                        "mark": row["mark"],
                        "analysis_tier": row["analysis_tier"],
                        "mode": row["mode"],
                        "metric": row["metric"],
                        "n_covid": int(row["n_covid"]),
                        "n_hc": int(row["n_hc"]),
                        "median_covid": float(row["median_covid"]),
                        "median_hc": float(row["median_hc"]),
                        "median_difference_covid_minus_hc": float(
                            row["median_difference_covid_minus_hc"]
                        ),
                        "mannwhitneyu_p": float(row["mannwhitneyu_p"]),
                        "cliffs_delta_covid_vs_hc": float(
                            row["cliffs_delta_covid_vs_hc"]
                        ),
                        "bh_fdr_all_tests": float(row["bh_fdr_all_tests"]),
                    }
                )

    concordance_rows = []
    key_fields = ("analysis", "peak_set_id", "mode", "metric")
    reference = {
        tuple(row[field] for field in key_fields): row
        for row in combined_rows
        if row["callset"] == "circlemap_methods"
    }
    for callset in CALLSETS[1:]:
        comparison = {
            tuple(row[field] for field in key_fields): row
            for row in combined_rows
            if row["callset"] == callset
        }
        for key in sorted(set(reference) & set(comparison)):
            ref = reference[key]
            comp = comparison[key]
            concordance_rows.append(
                {
                    "reference_callset": "circlemap_methods",
                    "comparison_callset": callset,
                    **{field: value for field, value in zip(key_fields, key)},
                    "reference_median_difference": ref[
                        "median_difference_covid_minus_hc"
                    ],
                    "comparison_median_difference": comp[
                        "median_difference_covid_minus_hc"
                    ],
                    "direction_concordant": int(
                        np.sign(ref["median_difference_covid_minus_hc"])
                        == np.sign(comp["median_difference_covid_minus_hc"])
                    ),
                    "reference_q_lt_0_05": int(ref["bh_fdr_all_tests"] < 0.05),
                    "comparison_q_lt_0_05": int(comp["bh_fdr_all_tests"] < 0.05),
                    "significance_concordant": int(
                        (ref["bh_fdr_all_tests"] < 0.05)
                        == (comp["bh_fdr_all_tests"] < 0.05)
                    ),
                }
            )
    return combined_rows, concordance_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    tables = args.root / "tables"
    eccgene, candidates, eccgene_concordance = eccgene_summary(args.root)
    chromatin, chromatin_concordance = chromatin_summary(args.root)
    write_tsv(tables / "eccgene_callset_summary.tsv", eccgene)
    write_tsv(tables / "eccgene_candidate_stability.tsv", candidates)
    write_tsv(tables / "eccgene_callset_concordance.tsv", eccgene_concordance)
    write_tsv(tables / "chromatin_callset_group_statistics.tsv", chromatin)
    write_tsv(tables / "chromatin_callset_concordance.tsv", chromatin_concordance)

    parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "callsets": CALLSETS,
        "eccgene_definitions": DEFINITIONS,
        "candidate_genes": CANDIDATE_GENES,
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.root / "downstream_summary_parameters.json").open("w") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")
    with (args.root / "downstream_summary_validation.txt").open("w") as handle:
        handle.write("PASS\n")
        handle.write(f"eccgene_summary_rows={len(eccgene)}\n")
        handle.write(f"chromatin_group_statistic_rows={len(chromatin)}\n")


if __name__ == "__main__":
    main()
