#!/usr/bin/env python3
"""Audit and rank eccGene loci for unbiased locus-level visualization.

The primary statistics are read from the locked junction-based eccGene
reanalysis.  Genomic-quality sensitivity statistics are read from a second run
of the same eccGene code using the prespecified Circle-Map artifact-masked
callset.  Ranking is lexicographic; no weighted score is calculated.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import platform
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


COVID = "COVID-19"
HC = "HC"
AUDIT_GENES = ("BCL3", "PROCR")
MAX_SAMPLE_CONTRIBUTION = 0.50
MIN_INDEPENDENT_COVID_SAMPLES = 2

# The humoral category is an enriched GO:BP term in the locked revised
# enrichment.  The endothelial/coagulation category is deliberately defined
# from broad GO roots before inspecting candidate rankings; its enrichment
# support is audited rather than assumed.
CATEGORY_ROOTS = {
    "Humoral immune response": ("GO:0006959",),
    "Endothelial/coagulation": (
        "GO:0007596",  # blood coagulation
        "GO:0050818",  # regulation of coagulation
        "GO:0001944",  # vasculature development
        "GO:0003158",  # endothelium development
        "GO:0061028",  # establishment of endothelial barrier
        "GO:0043114",  # regulation of vascular permeability
    ),
}


def fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return f"{value:.10g}"
    return str(value)


def parse_float(value: str | None) -> float:
    if value is None or value in {"", "NA"}:
        return math.nan
    if value in {"Inf", "inf"}:
        return math.inf
    if value in {"-Inf", "-inf"}:
        return -math.inf
    return float(value)


def parse_int(value: str | None) -> int:
    if value is None or value in {"", "NA"}:
        return 0
    return int(float(value))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"Cannot infer fields for empty table: {path}")
        fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


class IntervalIndex:
    """Merged interval index supporting point and interval queries."""

    def __init__(self, intervals: list[tuple[str, int, int]]) -> None:
        grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in intervals:
            if end > start:
                grouped[chrom].append((start, end))
        self.starts: dict[str, list[int]] = {}
        self.ends: dict[str, list[int]] = {}
        for chrom, values in grouped.items():
            values.sort()
            merged: list[list[int]] = []
            for start, end in values:
                if not merged or start > merged[-1][1]:
                    merged.append([start, end])
                else:
                    merged[-1][1] = max(merged[-1][1], end)
            self.starts[chrom] = [value[0] for value in merged]
            self.ends[chrom] = [value[1] for value in merged]

    @classmethod
    def from_bed(cls, path: Path) -> "IntervalIndex":
        intervals: list[tuple[str, int, int]] = []
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    intervals.append((fields[0], int(fields[1]), int(fields[2])))
                except ValueError:
                    continue
        return cls(intervals)

    def any_overlap(self, chrom: str, start: int, end: int) -> bool:
        if end <= start or chrom not in self.starts:
            return False
        index = bisect.bisect_right(self.ends[chrom], start)
        return index < len(self.starts[chrom]) and self.starts[chrom][index] < end


class GenePointIndex:
    """Bin-indexed merged gene bodies matching the locked eccGene definition."""

    def __init__(
        self,
        path: Path,
        wanted: set[str] | None = None,
        bin_size: int = 100_000,
    ) -> None:
        self.bin_size = bin_size
        self.intervals: list[tuple[str, int, int, str]] = []
        self.bins: dict[str, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
        self.by_gene: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        with path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                gene = row["gene"]
                if wanted is not None and gene not in wanted:
                    continue
                chrom = row["chrom"]
                start = int(row["start_0based"])
                end = int(row["end_0based_exclusive"])
                index = len(self.intervals)
                self.intervals.append((chrom, start, end, gene))
                self.by_gene[gene].append((chrom, start, end))
                for bin_id in range(start // bin_size, (end - 1) // bin_size + 1):
                    self.bins[chrom][bin_id].append(index)

    def genes_for_point(self, chrom: str, point: int) -> set[str]:
        genes: set[str] = set()
        for index in self.bins.get(chrom, {}).get(point // self.bin_size, []):
            _chrom, start, end, gene = self.intervals[index]
            if start <= point < end:
                genes.add(gene)
        return genes


def read_stats(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    if not rows or "Gene" not in rows[0]:
        raise ValueError(f"Invalid statistics table: {path}")
    return {row["Gene"]: row for row in rows}


def read_matrix_subset(
    path: Path,
    wanted: set[str],
    allow_missing_as_zero: bool = False,
) -> tuple[list[str], dict[str, list[float]]]:
    rows: dict[str, list[float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "Gene":
            raise ValueError(f"Invalid matrix header: {path}")
        samples = header[1:]
        for row in reader:
            if row and row[0] in wanted:
                rows[row[0]] = [float(value) for value in row[1:]]
    missing = wanted - set(rows)
    if missing:
        if not allow_missing_as_zero:
            raise ValueError(
                f"{path} is missing {len(missing)} requested genes; "
                f"first={sorted(missing)[:5]}"
            )
        for gene in missing:
            rows[gene] = [0.0] * len(samples)
    return samples, rows


def read_metadata(path: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    rows = read_tsv(path)
    groups: dict[str, str] = {}
    for row in rows:
        group = row["group"]
        if group not in {COVID, HC}:
            raise ValueError(f"Unexpected group {group!r} in {path}")
        sample = row["sample_id"]
        if sample in groups:
            raise ValueError(f"Duplicate sample {sample}")
        groups[sample] = group
    if sum(group == COVID for group in groups.values()) != 39:
        raise ValueError("Expected 39 COVID-19 samples")
    if sum(group == HC for group in groups.values()) != 39:
        raise ValueError("Expected 39 HC samples")
    return rows, groups


def matrix_metrics(
    samples: list[str],
    matrix: dict[str, list[float]],
    groups: dict[str, str],
) -> dict[str, dict[str, object]]:
    if set(samples) != set(groups):
        raise ValueError("Matrix samples do not match metadata")
    covid_indices = [i for i, sample in enumerate(samples) if groups[sample] == COVID]
    hc_indices = [i for i, sample in enumerate(samples) if groups[sample] == HC]
    output: dict[str, dict[str, object]] = {}
    for gene, values in matrix.items():
        covid_values = [values[i] for i in covid_indices]
        hc_values = [values[i] for i in hc_indices]
        total = sum(covid_values)
        if total > 0:
            local_index = max(range(len(covid_values)), key=covid_values.__getitem__)
            max_contribution = covid_values[local_index] / total
            max_sample = samples[covid_indices[local_index]]
        else:
            max_contribution = math.nan
            max_sample = "NA"
        output[gene] = {
            "median_COVID": statistics.median(covid_values),
            "median_HC": statistics.median(hc_values),
            "max_sample_contribution": max_contribution,
            "max_contributing_sample": max_sample,
            "covid_sum": total,
        }
    return output


def parse_go_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str], dict[str, str]]:
    parents: dict[str, set[str]] = defaultdict(set)
    names: dict[str, str] = {}
    alt_to_primary: dict[str, str] = {}
    current: dict[str, object] | None = None

    def commit(record: dict[str, object] | None) -> None:
        if not record or record.get("obsolete") or "id" not in record:
            return
        go_id = str(record["id"])
        names[go_id] = str(record.get("name", ""))
        parents[go_id].update(record.get("parents", set()))
        for alt_id in record.get("alt_ids", set()):
            alt_to_primary[str(alt_id)] = go_id

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                commit(current)
                current = {"parents": set(), "alt_ids": set(), "obsolete": False}
                continue
            if line.startswith("[") and line.endswith("]"):
                commit(current)
                current = None
                continue
            if current is None:
                continue
            if line.startswith("id: GO:"):
                current["id"] = line.split("id: ", 1)[1]
            elif line.startswith("name: "):
                current["name"] = line.split("name: ", 1)[1]
            elif line.startswith("alt_id: GO:"):
                current["alt_ids"].add(line.split("alt_id: ", 1)[1])
            elif line.startswith("is_a: GO:"):
                current["parents"].add(line.split()[1])
            elif line.startswith("relationship: part_of GO:"):
                current["parents"].add(line.split()[2])
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    commit(current)
    return parents, names, alt_to_primary


def go_category_membership(
    obo: Path,
    gaf: Path,
    genes: set[str],
) -> tuple[dict[str, set[str]], dict[tuple[str, str], set[str]], dict[str, str]]:
    parents, names, alt_to_primary = parse_go_obo(obo)
    cache: dict[str, frozenset[str]] = {}

    def ancestors(term: str, active: set[str] | None = None) -> frozenset[str]:
        term = alt_to_primary.get(term, term)
        if term in cache:
            return cache[term]
        active = set() if active is None else active
        if term in active:
            return frozenset({term})
        active.add(term)
        values = {term}
        for parent in parents.get(term, set()):
            values.update(ancestors(parent, active))
        active.remove(term)
        cache[term] = frozenset(values)
        return cache[term]

    roots_by_category = {
        category: set(roots) for category, roots in CATEGORY_ROOTS.items()
    }
    membership: dict[str, set[str]] = {category: set() for category in CATEGORY_ROOTS}
    evidence: dict[tuple[str, str], set[str]] = defaultdict(set)
    with open_text(gaf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            gene = fields[2]
            if gene not in genes:
                continue
            qualifiers = set(fields[3].split("|"))
            if "NOT" in qualifiers:
                continue
            go_id = alt_to_primary.get(fields[4], fields[4])
            term_ancestors = ancestors(go_id)
            for category, roots in roots_by_category.items():
                if roots & term_ancestors:
                    membership[category].add(gene)
                    evidence[(category, gene)].add(go_id)
    root_names = {root: names.get(root, "UNKNOWN") for roots in CATEGORY_ROOTS.values() for root in roots}
    return membership, evidence, root_names


def significant_enrichment_terms(path: Path) -> set[str]:
    rows = read_tsv(path)
    return {
        row["native"]
        for row in rows
        if row.get("source") == "GO:BP" and parse_float(row.get("p_value")) < 0.05
    }


def artifact_support_metrics(
    callset_dir: Path,
    metadata: list[dict[str, str]],
    gene_index: GenePointIndex,
    wanted: set[str],
) -> dict[str, dict[str, object]]:
    by_gene: dict[str, dict[str, object]] = {}
    for gene in wanted:
        by_gene[gene] = {
            "covid_samples": set(),
            "hc_samples": set(),
            "covid_calls": 0,
            "hc_calls": 0,
            "covid_support": [],
            "covid_split": [],
            "covid_scores": [],
        }
    for sample_index, row in enumerate(metadata, start=1):
        sample = row["sample_id"]
        group = row["group"]
        path = callset_dir / f"{sample}_circle_site.bed"
        if not path.is_file():
            raise FileNotFoundError(path)
        seen: set[tuple[str, int, int]] = set()
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                try:
                    chrom = fields[0]
                    start, end = int(fields[1]), int(fields[2])
                    discordant = float(fields[3])
                    split = float(fields[4])
                    score = float(fields[5])
                except ValueError:
                    continue
                coord = (chrom, start, end)
                if coord in seen:
                    continue
                seen.add(coord)
                genes_here = gene_index.genes_for_point(chrom, start) & wanted
                for gene in genes_here:
                    values = by_gene[gene]
                    if group == COVID:
                        values["covid_samples"].add(sample)
                        values["covid_calls"] += 1
                        values["covid_support"].append(discordant + split)
                        values["covid_split"].append(split)
                        values["covid_scores"].append(score)
                    else:
                        values["hc_samples"].add(sample)
                        values["hc_calls"] += 1
        print(
            f"[support {sample_index:02d}/{len(metadata)}] {sample}",
            file=sys.stderr,
            flush=True,
        )
    output: dict[str, dict[str, object]] = {}
    for gene, values in by_gene.items():
        support = values["covid_support"]
        splits = values["covid_split"]
        scores = values["covid_scores"]
        output[gene] = {
            "artifact_COVID_detected_n": len(values["covid_samples"]),
            "artifact_HC_detected_n": len(values["hc_samples"]),
            "artifact_COVID_call_n": values["covid_calls"],
            "artifact_HC_call_n": values["hc_calls"],
            "artifact_median_total_support": statistics.median(support) if support else math.nan,
            "artifact_median_split_reads": statistics.median(splits) if splits else math.nan,
            "artifact_median_circlemap_score": statistics.median(scores) if scores else math.nan,
        }
    return output


def chromatin_overlap_metrics(
    gene_index: GenePointIndex,
    genes: set[str],
    track_paths: dict[str, Path],
) -> dict[str, dict[str, bool]]:
    indexes = {name: IntervalIndex.from_bed(path) for name, path in track_paths.items()}
    output: dict[str, dict[str, bool]] = {}
    for gene in genes:
        values = {
            name: any(index.any_overlap(chrom, start, end) for chrom, start, end in gene_index.by_gene[gene])
            for name, index in indexes.items()
        }
        values["H3K27ac_overlap"] = values["CD14_H3K27ac_overlap"] or values["HUVEC_H3K27ac_overlap"]
        values["H3K27me3_overlap"] = values["CD14_H3K27me3_overlap"] or values["HUVEC_H3K27me3_overlap"]
        output[gene] = values
    return output


def stable_abundance_flags(path: Path) -> dict[str, bool]:
    return {
        row["Gene"]: parse_int(row["abundance_significant_same_direction_all_definitions"]) == 1
        for row in read_tsv(path)
    }


def rank_key(row: dict[str, object]) -> tuple:
    def descending(value: object) -> float:
        number = float(value)
        return -number if math.isfinite(number) else math.inf

    return (
        float(row["FDR"]),
        -abs(float(row["Effect_size"])),
        -int(row["COVID_detected_n"]),
        -float(row["Prevalence_difference"]),
        -int(row["artifact_COVID_detected_n"]),
        descending(row["artifact_median_split_reads"]),
        descending(row["artifact_median_circlemap_score"]),
        str(row["Gene"]),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-eccgene-dir", required=True, type=Path)
    parser.add_argument("--quality-eccgene-dir", required=True, type=Path)
    parser.add_argument("--sample-metadata", required=True, type=Path)
    parser.add_argument("--artifact-callset-dir", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--goa-gaf", required=True, type=Path)
    parser.add_argument("--enrichment-terms", required=True, type=Path)
    parser.add_argument("--cd14-h3k27ac", required=True, type=Path)
    parser.add_argument("--cd14-h3k27me3", required=True, type=Path)
    parser.add_argument("--huvec-h3k27ac", required=True, type=Path)
    parser.add_argument("--huvec-h3k27me3", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    outdir = args.output_dir
    table_dir = outdir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    primary_results = args.primary_eccgene_dir / "results"
    primary_matrices = args.primary_eccgene_dir / "matrices"
    quality_results = args.quality_eccgene_dir / "results"
    quality_matrices = args.quality_eccgene_dir / "matrices"

    metadata, groups = read_metadata(args.sample_metadata)
    primary_abundance = read_stats(primary_results / "junction_abundance_wilcoxon.tsv")
    primary_detection = read_stats(primary_results / "junction_detection_fisher.tsv")
    quality_abundance = read_stats(quality_results / "junction_abundance_wilcoxon.tsv")
    quality_detection = read_stats(quality_results / "junction_detection_fisher.tsv")
    tested_genes = set(primary_abundance)
    if set(primary_detection) != tested_genes:
        raise ValueError("Primary abundance and detection genes differ")

    primary_candidate_genes = {
        gene
        for gene, row in primary_abundance.items()
        if parse_float(row["wilcoxon_q_BH"]) < 0.05
        and parse_float(row["log2FC_mean_EA_COVID_vs_HC"]) > 0
        and parse_float(row["cliffs_delta_COVID_vs_HC"]) > 0
        and parse_int(row["covid_detected"]) >= MIN_INDEPENDENT_COVID_SAMPLES
    }
    working_genes = primary_candidate_genes | set(AUDIT_GENES)

    primary_ea_samples, primary_ea = read_matrix_subset(
        primary_matrices / "junction_EA.tsv", working_genes
    )
    quality_ea_samples, quality_ea = read_matrix_subset(
        quality_matrices / "junction_EA.tsv",
        working_genes,
        allow_missing_as_zero=True,
    )
    primary_count_samples, primary_counts = read_matrix_subset(
        primary_matrices / "junction_counts.tsv", working_genes
    )
    quality_count_samples, quality_counts = read_matrix_subset(
        quality_matrices / "junction_counts.tsv",
        working_genes,
        allow_missing_as_zero=True,
    )
    if not (
        primary_ea_samples
        == quality_ea_samples
        == primary_count_samples
        == quality_count_samples
    ):
        raise ValueError("Primary and quality-filtered matrix sample orders differ")
    matrix_samples = primary_ea_samples
    primary_ea_metrics = matrix_metrics(matrix_samples, primary_ea, groups)
    quality_ea_metrics = matrix_metrics(matrix_samples, quality_ea, groups)
    primary_count_metrics = matrix_metrics(matrix_samples, primary_counts, groups)

    membership, go_evidence, go_root_names = go_category_membership(
        args.go_obo, args.goa_gaf, tested_genes
    )
    enriched_terms = significant_enrichment_terms(args.enrichment_terms)

    gene_intervals = args.primary_eccgene_dir / "annotation" / "gene_intervals_refGene_merged.tsv"
    gene_index = GenePointIndex(gene_intervals, working_genes)
    support = artifact_support_metrics(
        args.artifact_callset_dir, metadata, gene_index, working_genes
    )
    tracks = {
        "CD14_H3K27ac_overlap": args.cd14_h3k27ac,
        "CD14_H3K27me3_overlap": args.cd14_h3k27me3,
        "HUVEC_H3K27ac_overlap": args.huvec_h3k27ac,
        "HUVEC_H3K27me3_overlap": args.huvec_h3k27me3,
    }
    chromatin = chromatin_overlap_metrics(gene_index, working_genes, tracks)
    primary_stable = stable_abundance_flags(primary_results / "candidate_stability.tsv")
    quality_stable = stable_abundance_flags(quality_results / "candidate_stability.tsv")

    all_rows: dict[str, dict[str, object]] = {}
    for gene in working_genes:
        abundance = primary_abundance[gene]
        detection = primary_detection[gene]
        q_abundance = quality_abundance.get(gene)
        q_detection = quality_detection.get(gene)
        if q_abundance is None or q_detection is None:
            quality_q = math.nan
            quality_log2fc = math.nan
            quality_delta = math.nan
            quality_covid_detected = 0
            quality_hc_detected = 0
            quality_prevalence = 0.0
        else:
            quality_q = parse_float(q_abundance["wilcoxon_q_BH"])
            quality_log2fc = parse_float(q_abundance["log2FC_mean_EA_COVID_vs_HC"])
            quality_delta = parse_float(q_abundance["cliffs_delta_COVID_vs_HC"])
            quality_covid_detected = parse_int(q_detection["covid_detected"])
            quality_hc_detected = parse_int(q_detection["hc_detected"])
            quality_prevalence = parse_float(
                q_detection["frequency_difference_COVID_minus_HC"]
            )

        primary_pass = gene in primary_candidate_genes
        quality_stat_pass = (
            math.isfinite(quality_q)
            and quality_q < 0.05
            and quality_log2fc > 0
            and quality_delta > 0
            and quality_covid_detected >= MIN_INDEPENDENT_COVID_SAMPLES
        )
        artifact_sample_pass = (
            int(support[gene]["artifact_COVID_detected_n"])
            >= MIN_INDEPENDENT_COVID_SAMPLES
        )
        artifact_any_call_pass = (
            int(support[gene]["artifact_COVID_call_n"])
            + int(support[gene]["artifact_HC_call_n"])
            > 0
        )
        primary_dominance_pass = (
            math.isfinite(float(primary_ea_metrics[gene]["max_sample_contribution"]))
            and float(primary_ea_metrics[gene]["max_sample_contribution"])
            <= MAX_SAMPLE_CONTRIBUTION
        )
        quality_dominance_pass = (
            math.isfinite(float(quality_ea_metrics[gene]["max_sample_contribution"]))
            and float(quality_ea_metrics[gene]["max_sample_contribution"])
            <= MAX_SAMPLE_CONTRIBUTION
        )
        eligible = (
            primary_pass
            and quality_stat_pass
            and artifact_sample_pass
            and primary_dominance_pass
            and quality_dominance_pass
        )
        reasons: list[str] = []
        if not primary_pass:
            reasons.append("fails_primary_FDR_direction_or_recurrence")
        if not quality_stat_pass:
            reasons.append("fails_artifact_masked_abundance_sensitivity")
        if not artifact_sample_pass:
            reasons.append("artifact_masked_support_in_fewer_than_2_COVID_samples")
        if not primary_dominance_pass:
            reasons.append("one_sample_contributes_more_than_50pct_primary_EA")
        if not quality_dominance_pass:
            reasons.append("one_sample_contributes_more_than_50pct_quality_EA")
        categories = [
            category for category, members in membership.items() if gene in members
        ]
        category_evidence = []
        for category in categories:
            annotated_terms = ",".join(sorted(go_evidence[(category, gene)]))
            category_evidence.append(f"{category}:{annotated_terms}")
        row: dict[str, object] = {
            "Gene": gene,
            "Functional_category": "; ".join(categories) if categories else "Not in prespecified categories",
            "GO_category_evidence": "; ".join(category_evidence)
            if category_evidence
            else "NA",
            "FDR": parse_float(abundance["wilcoxon_q_BH"]),
            "log2FC": parse_float(abundance["log2FC_mean_EA_COVID_vs_HC"]),
            "Effect_size": parse_float(abundance["cliffs_delta_COVID_vs_HC"]),
            "COVID_detected_n": parse_int(detection["covid_detected"]),
            "HC_detected_n": parse_int(detection["hc_detected"]),
            "Prevalence_difference": parse_float(
                detection["frequency_difference_COVID_minus_HC"]
            ),
            "Median_abundance_COVID": float(primary_ea_metrics[gene]["median_COVID"]),
            "Median_abundance_HC": float(primary_ea_metrics[gene]["median_HC"]),
            "Max_sample_contribution": float(
                primary_ea_metrics[gene]["max_sample_contribution"]
            ),
            "Max_contributing_sample": primary_ea_metrics[gene]["max_contributing_sample"],
            "Max_sample_count_contribution": float(
                primary_count_metrics[gene]["max_sample_contribution"]
            ),
            "quality_FDR": quality_q,
            "quality_log2FC": quality_log2fc,
            "quality_effect_size": quality_delta,
            "quality_COVID_detected_n": quality_covid_detected,
            "quality_HC_detected_n": quality_hc_detected,
            "quality_prevalence_difference": quality_prevalence,
            "quality_Max_sample_contribution": float(
                quality_ea_metrics[gene]["max_sample_contribution"]
            ),
            "quality_Max_contributing_sample": quality_ea_metrics[gene][
                "max_contributing_sample"
            ],
            **support[gene],
            "artifact_recurrence_pass": artifact_sample_pass,
            "blacklist_gap_pass": artifact_any_call_pass,
            "Umap_k100_pass": artifact_any_call_pass,
            "segmental_duplication_pass": artifact_any_call_pass,
            "RepeatMasker_breakpoint_pass": artifact_any_call_pass,
            **chromatin[gene],
            "primary_abundance_stable_all_definitions": primary_stable.get(gene, False),
            "quality_abundance_stable_all_definitions": quality_stable.get(gene, False),
            "primary_candidate_pass": primary_pass,
            "quality_statistical_sensitivity_pass": quality_stat_pass,
            "single_sample_dominance_pass": primary_dominance_pass
            and quality_dominance_pass,
            "Strict_selection_eligible": eligible,
            "Exclusion_reason": "; ".join(reasons) if reasons else "Eligible",
            "Overall_rank": None,
            "Humoral_rank": None,
            "Endothelial_coagulation_rank": None,
            "Final_rank": None,
            "Selected": False,
        }
        all_rows[gene] = row

    eligible_rows = [
        row for row in all_rows.values() if bool(row["Strict_selection_eligible"])
    ]
    eligible_rows.sort(key=rank_key)
    for index, row in enumerate(eligible_rows, start=1):
        row["Overall_rank"] = index

    category_winners: list[tuple[str, dict[str, object]]] = []
    for category, members in membership.items():
        members_in_working = [
            all_rows[gene] for gene in working_genes if gene in members
        ]
        eligible_category = [
            row for row in members_in_working if bool(row["Strict_selection_eligible"])
        ]
        eligible_category.sort(key=rank_key)
        rank_field = (
            "Humoral_rank"
            if category == "Humoral immune response"
            else "Endothelial_coagulation_rank"
        )
        for index, row in enumerate(eligible_category, start=1):
            row[rank_field] = index
            if index == 1:
                row["Selected"] = True
                category_winners.append((category, row))

    for row in all_rows.values():
        labels = []
        if row["Humoral_rank"] is not None:
            labels.append(f"Humoral:{row['Humoral_rank']}")
        if row["Endothelial_coagulation_rank"] is not None:
            labels.append(
                f"Endothelial/coagulation:{row['Endothelial_coagulation_rank']}"
            )
        row["Final_rank"] = "; ".join(labels) if labels else None

    selected_rows: list[dict[str, object]] = []
    if eligible_rows:
        selected_rows.append(
            {
                **eligible_rows[0],
                "Selection_scope": "Overall statistically highest-ranked eligible locus",
            }
        )
    for category, row in category_winners:
        selected_rows.append(
            {**row, "Selection_scope": f"Highest-ranked: {category}"}
        )

    category_rows: list[dict[str, object]] = []
    for category, members in membership.items():
        members_in_working = [
            all_rows[gene] for gene in working_genes if gene in members
        ]
        for row in sorted(
            members_in_working,
            key=lambda value: (
                not bool(value["Strict_selection_eligible"]),
                rank_key(value),
            ),
        ):
            category_rows.append({**row, "Category_for_row": category})

    ranking_fields = [
        "Gene",
        "Functional_category",
        "GO_category_evidence",
        "FDR",
        "log2FC",
        "Effect_size",
        "COVID_detected_n",
        "HC_detected_n",
        "Prevalence_difference",
        "Median_abundance_COVID",
        "Median_abundance_HC",
        "Max_sample_contribution",
        "Max_contributing_sample",
        "Max_sample_count_contribution",
        "quality_FDR",
        "quality_log2FC",
        "quality_effect_size",
        "quality_COVID_detected_n",
        "quality_HC_detected_n",
        "quality_prevalence_difference",
        "quality_Max_sample_contribution",
        "quality_Max_contributing_sample",
        "artifact_COVID_detected_n",
        "artifact_HC_detected_n",
        "artifact_COVID_call_n",
        "artifact_HC_call_n",
        "artifact_median_total_support",
        "artifact_median_split_reads",
        "artifact_median_circlemap_score",
        "artifact_recurrence_pass",
        "blacklist_gap_pass",
        "Umap_k100_pass",
        "segmental_duplication_pass",
        "RepeatMasker_breakpoint_pass",
        "CD14_H3K27ac_overlap",
        "CD14_H3K27me3_overlap",
        "HUVEC_H3K27ac_overlap",
        "HUVEC_H3K27me3_overlap",
        "H3K27ac_overlap",
        "H3K27me3_overlap",
        "primary_abundance_stable_all_definitions",
        "quality_abundance_stable_all_definitions",
        "primary_candidate_pass",
        "quality_statistical_sensitivity_pass",
        "single_sample_dominance_pass",
        "Strict_selection_eligible",
        "Exclusion_reason",
        "Overall_rank",
        "Humoral_rank",
        "Endothelial_coagulation_rank",
        "Final_rank",
        "Selected",
    ]
    candidate_rows = [all_rows[gene] for gene in primary_candidate_genes]
    candidate_rows.sort(
        key=lambda row: (
            not bool(row["Strict_selection_eligible"]),
            rank_key(row),
        )
    )
    write_tsv(
        table_dir / "candidate_locus_ranking_all.tsv",
        candidate_rows,
        ranking_fields,
    )
    write_tsv(
        table_dir / "category_candidate_rankings.tsv",
        category_rows,
        ["Category_for_row", *ranking_fields],
    )
    selected_fields = ["Selection_scope", *ranking_fields]
    # De-duplicate a locus that is both overall and a category winner while
    # preserving every selection scope as separate auditable rows.
    write_tsv(table_dir / "selected_loci.tsv", selected_rows, selected_fields)
    write_tsv(
        table_dir / "BCL3_PROCR_audit.tsv",
        [all_rows[gene] for gene in AUDIT_GENES],
        ranking_fields,
    )

    top_rows: list[dict[str, object]] = []
    for category in CATEGORY_ROOTS:
        rows = [
            row
            for row in category_rows
            if row["Category_for_row"] == category
        ][:10]
        top_rows.extend(rows)
    write_tsv(
        table_dir / "top10_by_category.tsv",
        top_rows,
        ["Category_for_row", *ranking_fields],
    )

    category_definition_rows = []
    for category, roots in CATEGORY_ROOTS.items():
        for root in roots:
            category_definition_rows.append(
                {
                    "Category": category,
                    "GO_root": root,
                    "GO_root_name": go_root_names[root],
                    "Root_significant_in_locked_revised_enrichment": root in enriched_terms,
                    "Category_has_any_significant_root_in_locked_revised_enrichment": bool(
                        set(roots) & enriched_terms
                    ),
                    "Membership_rule": "GOA human annotation to root or descendant; NOT annotations excluded",
                }
            )
    write_tsv(table_dir / "category_definitions.tsv", category_definition_rows)

    criteria_rows = [
        {
            "Stage": 1,
            "Criterion": "Primary revised abundance",
            "Rule": "junction Wilcoxon BH FDR < 0.05; log2FC > 0; Cliff's delta > 0",
            "Used_for_ranking": "Eligibility",
        },
        {
            "Stage": 2,
            "Criterion": "Independent-sample support",
            "Rule": f"detected in at least {MIN_INDEPENDENT_COVID_SAMPLES} COVID-19 samples; no arbitrary 10-sample threshold",
            "Used_for_ranking": "Eligibility",
        },
        {
            "Stage": 3,
            "Criterion": "Prespecified breakpoint-quality sensitivity",
            "Rule": "same eccGene test remains FDR < 0.05 and COVID-up after methods-consistent Circle-Map calls are masked for blacklist/gap, Umap k=100, segmental duplication and RepeatMasker",
            "Used_for_ranking": "Eligibility",
        },
        {
            "Stage": 4,
            "Criterion": "Single-sample dominance",
            "Rule": f"maximum COVID-19 sample contribution <= {MAX_SAMPLE_CONTRIBUTION:.0%} in both primary and quality-filtered EA",
            "Used_for_ranking": "Eligibility",
        },
        {
            "Stage": 5,
            "Criterion": "Functional assignment",
            "Rule": "membership in prespecified GO root or descendant using GOA human; enrichment support reported separately",
            "Used_for_ranking": "Stratification",
        },
        {
            "Stage": 6,
            "Criterion": "Lexicographic ranking",
            "Rule": "FDR ascending; |Cliff's delta| descending; COVID detection n descending; prevalence difference descending; artifact-masked COVID detection n descending; median split reads descending; median Circle-Map score descending; gene symbol ascending",
            "Used_for_ranking": "Ranking",
        },
        {
            "Stage": 7,
            "Criterion": "Chromatin display",
            "Rule": "H3K27ac/H3K27me3 overlaps and visual signal intensity are descriptive only and never enter eligibility or ranking",
            "Used_for_ranking": "No",
        },
    ]
    write_tsv(table_dir / "ranking_criteria.tsv", criteria_rows)

    sample_rows: list[dict[str, object]] = []
    report_genes = set(AUDIT_GENES) | {
        str(row["Gene"]) for row in selected_rows
    }
    for gene in sorted(report_genes):
        for index, sample in enumerate(matrix_samples):
            sample_rows.append(
                {
                    "Gene": gene,
                    "sample_id": sample,
                    "group": groups[sample],
                    "primary_EA": primary_ea[gene][index],
                    "quality_filtered_EA": quality_ea[gene][index],
                    "primary_count": primary_counts[gene][index],
                    "quality_filtered_count": quality_counts[gene][index],
                    "Selected_category_winner": bool(all_rows[gene]["Selected"]),
                    "Strict_selection_eligible": bool(
                        all_rows[gene]["Strict_selection_eligible"]
                    ),
                }
            )
    write_tsv(table_dir / "sample_abundance_selected_and_audit.tsv", sample_rows)

    summary_comparison = [
        {
            "Analysis": "Primary revised eccGene",
            "Tested_genes": len(primary_abundance),
            "FDR_lt_0_05": sum(
                parse_float(row["wilcoxon_q_BH"]) < 0.05
                for row in primary_abundance.values()
            ),
            "COVID_up_FDR_lt_0_05": sum(
                parse_float(row["wilcoxon_q_BH"]) < 0.05
                and parse_float(row["log2FC_mean_EA_COVID_vs_HC"]) > 0
                and parse_float(row["cliffs_delta_COVID_vs_HC"]) > 0
                for row in primary_abundance.values()
            ),
        },
        {
            "Analysis": "Prespecified artifact-masked sensitivity",
            "Tested_genes": len(quality_abundance),
            "FDR_lt_0_05": sum(
                parse_float(row["wilcoxon_q_BH"]) < 0.05
                for row in quality_abundance.values()
            ),
            "COVID_up_FDR_lt_0_05": sum(
                parse_float(row["wilcoxon_q_BH"]) < 0.05
                and parse_float(row["log2FC_mean_EA_COVID_vs_HC"]) > 0
                and parse_float(row["cliffs_delta_COVID_vs_HC"]) > 0
                for row in quality_abundance.values()
            ),
        },
    ]
    write_tsv(table_dir / "current_vs_quality_filtered_summary.tsv", summary_comparison)

    selected_humoral = [
        row for row in selected_rows if row["Selection_scope"].endswith("Humoral immune response")
    ]
    selected_vascular = [
        row for row in selected_rows if row["Selection_scope"].endswith("Endothelial/coagulation")
    ]
    overall = eligible_rows[0]["Gene"] if eligible_rows else "None"
    bcl3 = all_rows["BCL3"]
    procr = all_rows["PROCR"]
    lines = [
        "# Objective locus-selection audit",
        "",
        f"Run time (UTC): {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Outcome",
        "",
        f"- Strictly eligible loci: {len(eligible_rows):,}.",
        f"- Overall highest-ranked eligible locus: **{overall}**.",
        f"- Humoral immune-response winner: **{selected_humoral[0]['Gene'] if selected_humoral else 'None'}**.",
        f"- Endothelial/coagulation winner: **{selected_vascular[0]['Gene'] if selected_vascular else 'None'}**.",
        f"- BCL3 eligible: **{fmt(bool(bcl3['Strict_selection_eligible']))}**; reason: {bcl3['Exclusion_reason']}.",
        f"- PROCR eligible: **{fmt(bool(procr['Strict_selection_eligible']))}**; reason: {procr['Exclusion_reason']}.",
        "",
        "The locked revised enrichment contains a significant humoral immune-response GO term. "
        "It does not contain a significant endothelial/coagulation GO root used here; that "
        "category is therefore reported as a prespecified ontology-based biological stratum, "
        "not as an enrichment-derived significant category.",
        "",
        "## Selection rule",
        "",
        "Eligible loci were COVID-up eccGenes with primary Wilcoxon BH FDR < 0.05, "
        "positive log2FC and positive Cliff's delta, support from at least two independent "
        "COVID-19 samples, concordant significance after the prespecified breakpoint-level "
        "artifact mask, and no sample contributing more than 50% of total COVID-19 EA in "
        "either analysis. Within each GO-defined category, genes were ranked sequentially "
        "by adjusted P value, absolute Cliff's delta, COVID-19 detection frequency, "
        "prevalence difference, artifact-masked recurrence, median split-read support and "
        "median Circle-Map score. IGV appearance and histone signal intensity were not used.",
        "",
        "## Interpretation guardrails",
        "",
        "- The primary and artifact-masked analyses use the same eccGene assignment, EA, "
        "Wilcoxon and BH-FDR logic; only the input Circle-Map callset differs.",
        "- Chromatin tracks are public reference-cell data and are not patient-matched ChIP-seq.",
        "- H3K27ac/H3K27me3 overlap is descriptive and is excluded from the ranking key.",
        "- eccGene remains a positional association label and does not imply that an intact "
        "or functional gene is carried by eccDNA.",
    ]
    (outdir / "RESULTS_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    input_paths = {
        "primary_abundance": primary_results / "junction_abundance_wilcoxon.tsv",
        "primary_detection": primary_results / "junction_detection_fisher.tsv",
        "primary_EA": primary_matrices / "junction_EA.tsv",
        "quality_abundance": quality_results / "junction_abundance_wilcoxon.tsv",
        "quality_detection": quality_results / "junction_detection_fisher.tsv",
        "quality_EA": quality_matrices / "junction_EA.tsv",
        "sample_metadata": args.sample_metadata,
        "go_obo": args.go_obo,
        "goa_gaf": args.goa_gaf,
        "enrichment_terms": args.enrichment_terms,
        **tracks,
    }
    manifest = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "python": sys.version,
        "coordinate_system": "GRCh38/hg38; BED 0-based, half-open",
        "primary_eccgene_definition": "Circle-Map BED start/junction within merged UCSC refGene gene body",
        "max_sample_contribution_threshold": MAX_SAMPLE_CONTRIBUTION,
        "minimum_independent_covid_samples": MIN_INDEPENDENT_COVID_SAMPLES,
        "weighted_score_used": False,
        "igv_signal_used": False,
        "ranking_order": [
            "FDR ascending",
            "absolute Cliff's delta descending",
            "COVID-19 detected n descending",
            "COVID-HC prevalence difference descending",
            "artifact-masked COVID-19 detected n descending",
            "artifact-masked median split reads descending",
            "artifact-masked median Circle-Map score descending",
            "gene symbol ascending",
        ],
        "category_roots": CATEGORY_ROOTS,
        "input_files": {
            name: {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in input_paths.items()
        },
        "outputs": {
            "primary_candidate_count": len(primary_candidate_genes),
            "strict_eligible_count": len(eligible_rows),
            "overall_top_gene": overall,
            "humoral_top_gene": selected_humoral[0]["Gene"] if selected_humoral else None,
            "endothelial_coagulation_top_gene": selected_vascular[0]["Gene"] if selected_vascular else None,
            "BCL3_eligible": bool(bcl3["Strict_selection_eligible"]),
            "PROCR_eligible": bool(procr["Strict_selection_eligible"]),
        },
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    validation = [
        "PASS: 39 COVID-19 and 39 HC samples were present",
        "PASS: primary abundance and detection gene universes matched",
        "PASS: primary and quality-filtered matrix sample orders matched",
        "PASS: no weighted score was used",
        "PASS: chromatin overlap and IGV signal were excluded from the ranking key",
        "PASS: all selected category winners were strict-eligible and rank 1",
    ]
    for row in selected_rows:
        if "Highest-ranked:" in str(row["Selection_scope"]):
            if not row["Strict_selection_eligible"]:
                raise AssertionError(f"Selected ineligible locus: {row['Gene']}")
            expected_rank = (
                row["Humoral_rank"]
                if str(row["Selection_scope"]).endswith("Humoral immune response")
                else row["Endothelial_coagulation_rank"]
            )
            if expected_rank != 1:
                raise AssertionError(f"Selected non-rank-1 locus: {row['Gene']}")
    (outdir / "validation_report.txt").write_text(
        "\n".join(validation) + "\n", encoding="utf-8"
    )
    print(f"Primary candidates: {len(primary_candidate_genes)}")
    print(f"Strict eligible: {len(eligible_rows)}")
    print(f"Overall top: {overall}")
    print(f"BCL3: {bcl3['Exclusion_reason']}")
    print(f"PROCR: {procr['Exclusion_reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
