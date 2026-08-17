#!/usr/bin/env python3
"""Per-sample chromosome distribution and gene-element O/E, per the Methods.

REIMPLEMENTED under pre-declaration amendments A3/A4. Where the original code
and the manuscript Methods disagree, the Methods wins -- that is the user's
explicit instruction, and it has the side effect of making the Methods text
true for the first time.

Two quantities, both defined verbatim by the manuscript.

Chromosome distribution (Figure 2b, Table S10)

    normalized_fraction_per_mb
        = (eccDNAs on the chromosome / eccDNAs in the sample x 100)
          / chromosome length in Mb

    "After normalizing by chromosomal length (percentage of eccDNAs per Mb)".

Gene-element enrichment (Figure 2d, Tables S12/S13)

    O/E = ratio of unique eccDNAs falling in a given element
          / ratio of the genome occupied by that element        (Methods eq. 1)

    with element membership decided by the eccDNA START coordinate:
    "we counted unique eccDNA molecules by their start coordinates and mapped
    them to seven genomic elements ... using bedtools" (Methods).

Two deliberate departures from the v1 artefacts, both recorded in A4:

  * v1 intersected the FULL eccDNA interval (`bedtools intersect -u`), not the
    start coordinate.
  * v1's expected fraction was a pooled cohort eccDNA background (a constant
    16,717,103 across all samples), not the genome-occupied fraction the
    Methods equation specifies.

So v2 element O/E values are not comparable value-for-value with v1's; only
directions and between-cohort contrasts are. Chromosome values are comparable,
since that definition is unchanged.

Element membership is evaluated against the merged element BEDs, so an eccDNA
counts at most once per element. Elements are not mutually exclusive, by
design -- an exon start also lies inside a gene body.
"""

from __future__ import annotations

import bisect
import csv
import sys
from pathlib import Path

ELEMENTS = ["exon", "intron", "5UTR", "3UTR", "CpG_islands", "Gene2KbU", "Gene2KbD"]
CANON = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CANON_SET = set(CANON)


def load_merged(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    """Merged, sorted intervals per chromosome as parallel start/end lists."""
    starts: dict[str, list[int]] = {}
    ends: dict[str, list[int]] = {}
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in CANON_SET:
                continue
            try:
                a, b = int(f[1]), int(f[2])
            except ValueError:
                continue
            if b <= a:
                continue
            starts.setdefault(f[0], []).append(a)
            ends.setdefault(f[0], []).append(b)
    out = {}
    for chrom in starts:
        pairs = sorted(zip(starts[chrom], ends[chrom]))
        merged: list[tuple[int, int]] = []
        for a, b in pairs:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        out[chrom] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def covered_bp(index: dict[str, tuple[list[int], list[int]]]) -> int:
    return sum(e - s for chrom in index
               for s, e in zip(index[chrom][0], index[chrom][1]))


def contains(index: dict[str, tuple[list[int], list[int]]], chrom: str,
             pos: int) -> bool:
    """Is the 0-based point `pos` inside any merged interval on `chrom`?"""
    entry = index.get(chrom)
    if entry is None:
        return False
    starts, ends = entry
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < ends[i]


def main() -> int:
    adapter, annot, outdir = (Path(sys.argv[1]), Path(sys.argv[2]),
                              Path(sys.argv[3]))
    outdir.mkdir(parents=True, exist_ok=True)

    sizes = {}
    with (annot / "hg38.chrom.sizes").open() as fh:
        for line in fh:
            f = line.split()
            if f and f[0] in CANON_SET:
                sizes[f[0]] = int(f[1])
    genome_bp = sum(sizes.values())
    print(f"canonical genome: {genome_bp:,} bp over {len(sizes)} chromosomes",
          flush=True)

    index = {}
    genome_fraction = {}
    for element in ELEMENTS:
        path = annot / f"{element}_merged.bed"
        if not path.is_file():
            raise SystemExit(f"missing element annotation: {path}")
        index[element] = load_merged(path)
        bp = covered_bp(index[element])
        genome_fraction[element] = bp / genome_bp
        print(f"  {element:12s} {bp:>13,} bp  genome fraction "
              f"{genome_fraction[element]:.10f}", flush=True)

    with (outdir / "element_genome_fraction.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["element", "covered_bp", "genome_bp", "genome_fraction"])
        for element in ELEMENTS:
            w.writerow([element, covered_bp(index[element]), genome_bp,
                        f"{genome_fraction[element]:.10f}"])

    samples = [(p.name.replace("_circle_site.bed", ""), "COVID-19", p)
               for p in sorted((adapter / "covid").glob("*_circle_site.bed"))]
    samples += [(p.name.replace("_circle_site.bed", ""), "HC", p)
                for p in sorted((adapter / "normal").glob("*_circle_site.bed"))]
    if len(samples) != 78:
        raise SystemExit(f"expected 78 samples, found {len(samples)}")

    chrom_rows, element_rows = [], []
    for n, (sample, group, bed) in enumerate(samples, 1):
        per_chrom = dict.fromkeys(CANON, 0)
        per_element = dict.fromkeys(ELEMENTS, 0)
        total = 0
        seen: set[tuple[str, int, int]] = set()
        with bed.open() as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or f[0] not in CANON_SET:
                    continue
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    continue
                if end <= start:
                    continue
                key = (f[0], start, end)
                if key in seen:          # "unique eccDNA molecules"
                    continue
                seen.add(key)
                total += 1
                per_chrom[f[0]] += 1
                for element in ELEMENTS:
                    if contains(index[element], f[0], start):
                        per_element[element] += 1

        for chrom in CANON:
            chrom_rows.append({
                "sample_id": sample, "group": group, "chromosome": chrom,
                "eccdna_count": per_chrom[chrom], "sample_total": total,
                "chromosome_length_bp": sizes[chrom],
                "normalized_fraction_per_mb":
                    f"{(per_chrom[chrom] / total * 100) / (sizes[chrom] / 1e6):.6f}"
                    if total else "NA",
            })
        for element in ELEMENTS:
            observed = per_element[element] / total if total else float("nan")
            expected = genome_fraction[element]
            element_rows.append({
                "sample_id": sample, "group": group, "element": element,
                "total_eccdna_n": total,
                "observed_overlap_o": per_element[element],
                "observed_fraction": f"{observed:.6f}",
                "expected_fraction": f"{expected:.6f}",
                "oe_ratio": f"{observed / expected:.6f}",
            })
        print(f"[{n}/78] {sample} ({group}): {total:,} unique eccDNAs", flush=True)

    for name, rows in (("chromosome_distribution.tsv", chrom_rows),
                       ("gene_feature_enrichment_scores.tsv", element_rows)):
        path = outdir / name
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {name}: {len(rows)} rows", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
