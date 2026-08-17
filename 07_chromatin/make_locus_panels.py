#!/usr/bin/env python3
"""Regenerate the BCL3 and PROCR locus panels (Figure 5c, d) in hg38.

The published panels were IGV screenshots annotated with hg19 coordinates, which
contradicted the rest of the revised analysis after everything was harmonized to
hg38 (Reviewer 1 comment 8, Reviewer 3 major comment 4). These panels are drawn
directly from the hg38 data instead: eccDNA coverage aggregated per cohort from
the same eligible eccDNA calls used elsewhere, the liftOver-converted histone
peak sets, and hg38 refGene exon structure.

Coverage is expressed on the same scale as the EPM used throughout the study, so
the two cohorts are directly comparable, and both cohorts share one y axis per
locus as in the original figure.
"""

from __future__ import annotations

import argparse
import json
import socket
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402

BIN_BP = 25
FLANK_BP = 1500

# Peak sets shown beneath the coverage tracks. The GSE179184 A549-ACE2 sets are the
# infected-cell context for Figure 5; the CD14+ sets are the reference immune-cell
# context used in Figure 4, included so the locus can be read against both.
PEAK_TRACKS = [
    ("A549_GSE179184_H3K27ac_official", "A549 H3K27ac", "#C93E3F"),
    ("A549_GSE179184_H3K4me3_official", "A549 H3K4me3", "#EC6F00"),
    ("CD14_H3K27ac", "CD14$^+$ H3K27ac", "#459434"),
    ("CD14_H3K27me3", "CD14$^+$ H3K27me3", "#6F7B91"),
]

HC_COLOR, COVID_COLOR = "#0272B2", "#EC6F00"
GENE_COLOR = "#253247"


@dataclass
class Locus:
    gene: str
    chrom: str
    start: int
    end: int
    panel: str


def log(message: str) -> None:
    base.log(message)


def read_refgene(path: Path, genes: Sequence[str]) -> Dict[str, dict]:
    """Longest hg38 refGene transcript per requested gene, with exon structure."""
    wanted = set(genes)
    best: Dict[str, dict] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16:
                continue
            symbol = parts[12]
            if symbol not in wanted:
                continue
            chrom = parts[2]
            if chrom not in base.CANON_SET:
                continue
            start, end = int(parts[4]), int(parts[5])
            exon_starts = [int(v) for v in parts[9].rstrip(",").split(",") if v]
            exon_ends = [int(v) for v in parts[10].rstrip(",").split(",") if v]
            record = {
                "gene": symbol,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": parts[3],
                "transcript": parts[1],
                "exons": list(zip(exon_starts, exon_ends)),
            }
            if symbol not in best or (end - start) > (best[symbol]["end"] - best[symbol]["start"]):
                best[symbol] = record
    missing = wanted - set(best)
    if missing:
        raise RuntimeError(f"genes not found in refGene: {sorted(missing)}")
    return best


def locus_coverage(
    sample_manifest: pd.DataFrame,
    loci: Sequence[Locus],
    chrom_sizes: Dict[str, int],
    forbidden_arrays,
) -> Dict[Tuple[str, str], np.ndarray]:
    """Per-cohort mean eccDNA coverage per million mapped alignments, per bin.

    Only eligible eccDNAs are counted, using the same canonical-chromosome,
    in-bounds and blacklist/gap filters as every other analysis, so the panels
    cannot show signal that the statistics excluded.
    """
    grids = {
        locus.gene: np.arange(locus.start, locus.end + BIN_BP, BIN_BP) for locus in loci
    }
    accum: Dict[Tuple[str, str], List[np.ndarray]] = {
        (locus.gene, group): [] for locus in loci for group in ("COVID", "HC")
    }

    samples = sample_manifest.sort_values(["group", "sample_id"]).to_dict(orient="records")
    for index, sample in enumerate(samples, 1):
        ecc, _ = base.load_sample_eccdna(
            Path(sample["eccdna_bed"]), chrom_sizes, forbidden_arrays
        )
        mapped = int(sample["mapped_alignments_idxstats"])
        scale = 1e6 / mapped if mapped > 0 else 0.0
        chrom_values = ecc["chrom"].to_numpy()
        starts = ecc["start"].to_numpy(dtype=np.int64)
        ends = ecc["end"].to_numpy(dtype=np.int64)
        for locus in loci:
            grid = grids[locus.gene]
            counts = np.zeros(len(grid) - 1, dtype=np.float64)
            keep = (
                (chrom_values == locus.chrom)
                & (ends > locus.start)
                & (starts < locus.end)
            )
            if keep.any():
                for s, e in zip(starts[keep].tolist(), ends[keep].tolist()):
                    lo = max(s, locus.start)
                    hi = min(e, locus.end)
                    first = (lo - locus.start) // BIN_BP
                    last = (hi - 1 - locus.start) // BIN_BP
                    counts[first : last + 1] += 1.0
            accum[(locus.gene, sample["group"])].append(counts * scale)
        if index % 10 == 0 or index == len(samples):
            log(f"coverage {index}/{len(samples)} samples")

    out = {}
    for key, stack in accum.items():
        out[key] = np.mean(np.vstack(stack), axis=0) if stack else np.zeros(0)
    return out, grids


def peaks_in_window(bed: Path, locus: Locus) -> List[Tuple[int, int]]:
    intervals = []
    with open(bed, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            parts = line.split("\t")
            if parts[0] != locus.chrom:
                continue
            start, end = int(parts[1]), int(parts[2])
            if end > locus.start and start < locus.end:
                intervals.append((max(start, locus.start), min(end, locus.end)))
    return intervals


def make_figure(
    loci: Sequence[Locus],
    coverage: Dict[Tuple[str, str], np.ndarray],
    grids: Dict[str, np.ndarray],
    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]],
    genes: Dict[str, dict],
    dirs: base.Dirs,
    stem: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    base.setup_matplotlib()
    body_font = mpl.rcParams["font.family"]
    body_font = body_font[0] if isinstance(body_font, (list, tuple)) else body_font
    mpl.rcParams.update(
        {
            "mathtext.fontset": "custom",
            "mathtext.rm": body_font,
            "mathtext.it": f"{body_font}:italic",
            "mathtext.bf": f"{body_font}:bold",
            "mathtext.cal": f"{body_font}:italic",
            "mathtext.tt": body_font,
            "mathtext.default": "regular",
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
        }
    )
    LABEL_SIZE, TICK_SIZE, TRACK_SIZE, PANEL_SIZE = 6.0, 5.5, 5.0, 8.0

    fig = plt.figure(figsize=(183 / 25.4, 78 / 25.4))
    # Two loci side by side; within each, two coverage tracks, the peak tracks and
    # the gene model, with height ratios matching their information content.
    heights = [1.5, 1.5] + [0.32] * len(PEAK_TRACKS) + [0.5]
    outer = GridSpec(
        1, 2, figure=fig, wspace=0.30, left=0.115, right=0.99, top=0.90, bottom=0.16
    )

    for locus in loci:
        column = 0 if locus.panel == "c" else 1
        inner = outer[0, column].subgridspec(len(heights), 1, height_ratios=heights, hspace=0.18)
        grid = grids[locus.gene]
        centres = (grid[:-1] + grid[1:]) / 2.0 / 1000.0
        vmax = max(
            float(np.max(coverage[(locus.gene, "COVID")])),
            float(np.max(coverage[(locus.gene, "HC")])),
            1e-9,
        )

        axes_list = []
        for row, (group, color, label) in enumerate(
            (("COVID", COVID_COLOR, "COVID-19"), ("HC", HC_COLOR, "HC"))
        ):
            ax = fig.add_subplot(inner[row, 0])
            axes_list.append(ax)
            values = coverage[(locus.gene, group)]
            ax.fill_between(centres, 0, values, color=color, linewidth=0, alpha=0.9, step="mid")
            ax.set_ylim(0, vmax * 1.12)
            ax.set_xlim(locus.start / 1000.0, locus.end / 1000.0)
            ax.set_yticks([0, vmax])
            ax.set_yticklabels(["0", f"{vmax:.2g}"])
            if row == 0:
                ax.set_ylabel(
                    "eccDNA per 10$^6$\nmapped alignments", fontsize=TRACK_SIZE, labelpad=2
                )
                ax.yaxis.set_label_coords(-0.105, -0.05)
            ax.tick_params(labelsize=TICK_SIZE, length=1.8, pad=1.2)
            ax.set_xticks([])
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.text(
                0.985, 0.88, label, transform=ax.transAxes, fontsize=TRACK_SIZE,
                ha="right", va="top", color=color, fontweight="bold",
            )
            if row == 0:
                ax.set_title(
                    f"{locus.gene}   {locus.chrom}:{locus.start:,}-{locus.end:,} (hg38)",
                    fontsize=LABEL_SIZE, pad=3,
                )
                ax.text(
                    -0.13, 1.20, locus.panel, transform=ax.transAxes, fontsize=PANEL_SIZE,
                    fontweight="bold", va="bottom", ha="left",
                )

        for offset, (peak_id, peak_label, peak_color) in enumerate(PEAK_TRACKS):
            ax = fig.add_subplot(inner[2 + offset, 0])
            axes_list.append(ax)
            for start, end in peaks[(locus.gene, peak_id)]:
                ax.axvspan(start / 1000.0, end / 1000.0, color=peak_color, linewidth=0)
            ax.set_xlim(locus.start / 1000.0, locus.end / 1000.0)
            ax.set_ylim(0, 1)
            ax.set_xticks([])
            ax.set_yticks([])
            for side in ("top", "right", "left", "bottom"):
                ax.spines[side].set_visible(False)
            ax.text(
                -0.012, 0.5, peak_label, transform=ax.transAxes, fontsize=TRACK_SIZE,
                ha="right", va="center",
            )

        ax = fig.add_subplot(inner[len(heights) - 1, 0])
        axes_list.append(ax)
        record = genes[locus.gene]
        ax.plot(
            [max(record["start"], locus.start) / 1000.0, min(record["end"], locus.end) / 1000.0],
            [0.5, 0.5], "-", color=GENE_COLOR, linewidth=0.6,
        )
        for exon_start, exon_end in record["exons"]:
            if exon_end <= locus.start or exon_start >= locus.end:
                continue
            ax.add_patch(
                plt.Rectangle(
                    (max(exon_start, locus.start) / 1000.0, 0.18),
                    (min(exon_end, locus.end) - max(exon_start, locus.start)) / 1000.0,
                    0.64, facecolor=GENE_COLOR, edgecolor="none",
                )
            )
        # Strand arrows along the transcript body.
        span = (min(record["end"], locus.end) - max(record["start"], locus.start)) / 1000.0
        if span > 0:
            for frac in np.linspace(0.08, 0.92, 7):
                xpos = max(record["start"], locus.start) / 1000.0 + frac * span
                ax.annotate(
                    "", xy=(xpos + (0.012 if record["strand"] == "+" else -0.012) * span, 0.5),
                    xytext=(xpos, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="white", linewidth=0.4,
                                    mutation_scale=3, shrinkA=0, shrinkB=0),
                )
        ax.set_xlim(locus.start / 1000.0, locus.end / 1000.0)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.tick_params(labelsize=TICK_SIZE, length=1.8, pad=1.2)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
        ax.xaxis.set_major_formatter(
            mpl.ticker.FuncFormatter(lambda v, _pos: f"{v:,.0f}")
        )
        ax.set_xlabel(f"{locus.chrom} position (kb, hg38)", fontsize=LABEL_SIZE)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.text(
            -0.012, 0.5, f"{record['gene']} ({record['strand']})", transform=ax.transAxes,
            fontsize=TRACK_SIZE, ha="right", va="center", style="italic",
        )

    for ext in ("pdf", "svg", "png"):
        fig.savefig(dirs.figures / f"{stem}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    log(f"wrote figures/{stem}.{{pdf,svg,png}}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default="/gpfs/data/gao/Covid-eccDNA")
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument(
        "--refgene", default="/gpfs/data/gao/Covid-eccDNA/eccDNA_abundance3/refGene.txt"
    )
    parser.add_argument("--stem", default="figure5_cd_loci_hg38")
    args = parser.parse_args(argv)

    dirs = base.ensure_dirs(Path(args.outdir))
    chrom_sizes = base.read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    base.build_forbidden_and_allowed(dirs, chrom_sizes)
    forbidden_df = base.read_bed3(
        dirs.reference / "hg38_forbidden.blacklist_gaps.canonical.merged.bed", chrom_sizes
    )
    forbidden_arrays = base.interval_dict_to_arrays(base.merge_intervals_from_df(forbidden_df))

    genes = read_refgene(Path(args.refgene), ["BCL3", "PROCR"])
    loci = []
    for gene, panel in (("BCL3", "c"), ("PROCR", "d")):
        record = genes[gene]
        loci.append(
            Locus(
                gene=gene,
                chrom=record["chrom"],
                start=max(0, record["start"] - FLANK_BP),
                end=min(chrom_sizes[record["chrom"]], record["end"] + FLANK_BP),
                panel=panel,
            )
        )
        log(f"{gene}: {record['chrom']}:{loci[-1].start:,}-{loci[-1].end:,} hg38 ({record['strand']})")

    peaks: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for locus in loci:
        for peak_id, _, _ in PEAK_TRACKS:
            bed = dirs.peaks / f"{base.sanitize_id(peak_id)}.liftover.hg38.canonical.filtered.merged.bed"
            if not bed.exists():
                raise RuntimeError(f"missing {bed}")
            peaks[(locus.gene, peak_id)] = peaks_in_window(bed, locus)
            log(f"{locus.gene} / {peak_id}: {len(peaks[(locus.gene, peak_id)])} peaks in window")

    sample_manifest = base.make_sample_manifest(Path(args.project_dir), dirs, args.samtools)
    coverage, grids = locus_coverage(sample_manifest, loci, chrom_sizes, forbidden_arrays)

    rows = []
    for locus in loci:
        grid = grids[locus.gene]
        for group in ("COVID", "HC"):
            for i, value in enumerate(coverage[(locus.gene, group)]):
                rows.append(
                    {
                        "gene": locus.gene,
                        "panel": locus.panel,
                        "chrom": locus.chrom,
                        "bin_start": int(grid[i]),
                        "bin_end": int(grid[i + 1]),
                        "group": group,
                        "mean_eccdna_per_million_mapped": float(value),
                    }
                )
    frame = pd.DataFrame(rows)
    frame.to_csv(dirs.tables / "figure5_cd_locus_coverage.tsv", sep="\t", index=False)
    log(f"wrote figure5_cd_locus_coverage.tsv ({len(frame)} rows)")

    make_figure(loci, coverage, grids, peaks, genes, dirs, args.stem)

    (dirs.outdir / "run_manifest_figure5cd.json").write_text(
        json.dumps(
            {
                "loci": [
                    {"gene": l.gene, "chrom": l.chrom, "start": l.start, "end": l.end, "panel": l.panel}
                    for l in loci
                ],
                "bin_bp": BIN_BP,
                "flank_bp": FLANK_BP,
                "peak_tracks": [p[0] for p in PEAK_TRACKS],
                "refgene": args.refgene,
                "assembly": "hg38",
                "python": sys.version,
                "hostname": socket.gethostname(),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    log("Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
