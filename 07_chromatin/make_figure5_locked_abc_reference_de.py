#!/usr/bin/env python3
"""Render Figure 5 with locked panels a--c and reference-style d/e tracks.

Panels a--c are intentionally *not redrawn*.  They are loaded verbatim from
the archived 600-dpi manuscript figure that was approved before the locus
selection revision.  This preserves their geometry, palette, data treatment,
and typography exactly.  Only the right-hand d/e region is replaced.

The new d/e panels replicate the supplied original locus-browser hierarchy:
chromosome-position overview, COVID-19 and HC eccDNA tracks, two chromatin
rows, and a RefSeq gene model.  The locus data are HLA-E and MCAM selected by
the preregistered display-eligibility procedure.  The two chromatin rows are
reference H3K27ac peak annotations (A549-ACE2 and CD14+ monocytes), not
patient ChIP-seq signal.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import run_chromatin_hg38_analysis as base  # noqa: E402


# The historical image is used solely as the immutable a--c layer.  It was
# rasterized at 600 dpi directly from mei:/home/gao/eccDNA/revise/Figure_5.pdf.
# Its hash makes any accidental change to those panels a hard error rather than
# a visual drift that could be missed during figure review.
LOCKED_PNG_SHA256 = "3bc52d33d016d44c0f879f63188f696ee988e92a3879be1d408c0511f40c28fe"
LOCKED_SIZE_PX = (4323, 3845)
FINAL_SIZE_INCH = (LOCKED_SIZE_PX[0] / 600.0, LOCKED_SIZE_PX[1] / 600.0)

LOCUS_PANELS = (
    ("HLA-E", "d", "Immune-response stratum"),
    ("MCAM", "e", "Endothelial/coagulation stratum"),
)
PEAK_TRACKS = (
    ("A549_GSE179184_H3K27ac_official", "A549-ACE2 H3K27ac", "#009999"),
    ("CD14_H3K27ac", "CD14$^+$ H3K27ac", "#006666"),
)

# These are the legacy Figure 5 red/blue roles.  They deliberately match the
# locked panels a--c and the supplied D/E reference, rather than introducing a
# new palette into one figure.
COVID_COLOR = "#C84F50"
HC_COLOR = "#36617B"
GENE_COLOR = "#0000B2"
TRACK_BORDER = "#111111"


def format_q(value: float) -> str:
    """Short BH-q representation suitable for the compact locus header."""
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    return f"{mantissa:.2g}×10$^{{{exponent}}}$"


def load_locked_background(path: Path) -> np.ndarray:
    """Read and validate the exact approved a--c layer from the mei PDF export."""
    from PIL import Image

    if not path.exists():
        raise FileNotFoundError(f"locked a--c PNG is missing: {path}")
    payload = path.read_bytes()
    observed_hash = hashlib.sha256(payload).hexdigest()
    if observed_hash != LOCKED_PNG_SHA256:
        raise RuntimeError(
            "locked a--c image hash differs from the approved source: "
            f"{observed_hash} != {LOCKED_PNG_SHA256}"
        )
    with Image.open(path) as image:
        rgba = np.asarray(image.convert("RGBA"))
    if (rgba.shape[1], rgba.shape[0]) != LOCKED_SIZE_PX:
        raise RuntimeError(
        "locked a--c image dimensions differ from the approved mei export: "
            f"{(rgba.shape[1], rgba.shape[0])} != {LOCKED_SIZE_PX}"
        )
    return rgba


def read_refgene(path: Path, genes: Sequence[str]) -> Dict[str, dict]:
    """Obtain the longest canonical hg38 refGene transcript for each locus."""
    wanted = set(genes)
    canonical = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    best: Dict[str, dict] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16 or parts[12] not in wanted or parts[2] not in canonical:
                continue
            start, end = int(parts[4]), int(parts[5])
            record = {
                "gene": parts[12],
                "chrom": parts[2],
                "start": start,
                "end": end,
                "strand": parts[3],
                "exons": list(
                    zip(
                        [int(value) for value in parts[9].rstrip(",").split(",") if value],
                        [int(value) for value in parts[10].rstrip(",").split(",") if value],
                    )
                ),
            }
            if record["gene"] not in best or end - start > best[record["gene"]]["end"] - best[record["gene"]]["start"]:
                best[record["gene"]] = record
    missing = wanted - set(best)
    if missing:
        raise RuntimeError(f"refGene records missing for: {', '.join(sorted(missing))}")
    return best


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    sizes = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def peaks_in_window(bed: Path, chrom: str, start: int, end: int) -> List[Tuple[int, int]]:
    intervals: List[Tuple[int, int]] = []
    with bed.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != chrom:
                continue
            left, right = int(parts[1]), int(parts[2])
            if right > start and left < end:
                intervals.append((max(left, start), min(right, end)))
    return intervals


def locus_arrays(
    coverage: pd.DataFrame, gene: str, record: dict
) -> Tuple[int, int, np.ndarray, Dict[str, np.ndarray]]:
    """Recover one common 25-bp grid and cohort-mean coverage vectors."""
    subset = coverage[coverage["gene"] == gene].copy()
    if set(subset["group"]) != {"COVID", "HC"}:
        raise RuntimeError(f"{gene}: coverage input lacks a cohort")
    starts = ends = None
    values: Dict[str, np.ndarray] = {}
    for group in ("COVID", "HC"):
        group_frame = subset[subset["group"] == group].sort_values("bin_start")
        local_starts = group_frame["bin_start"].to_numpy(int)
        local_ends = group_frame["bin_end"].to_numpy(int)
        if starts is None:
            starts, ends = local_starts, local_ends
        elif not np.array_equal(starts, local_starts) or not np.array_equal(ends, local_ends):
            raise RuntimeError(f"{gene}: cohort coverage grids differ")
        values[group] = group_frame["mean_eccdna_per_million_mapped"].to_numpy(float)
    assert starts is not None and ends is not None
    display_start, display_end = record["start"] - 1500, record["end"] + 1500
    bin_width = int(starts[1] - starts[0])
    if starts[0] != display_start or ends[-1] < display_end or ends[-1] - display_end >= bin_width:
        raise RuntimeError(f"{gene}: coverage does not match the predeclared refGene ±1.5 kb window")
    grid = np.concatenate([starts, [ends[-1]]])
    return display_start, display_end, grid, values


def prepare_row(ax) -> None:
    """Use the thin black framed rows of the supplied original browser panel."""
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(TRACK_BORDER)
        spine.set_linewidth(0.65)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.tick_params(axis="x", length=0, pad=1.5, labelsize=4.8)


def draw_chromosome_overview(
    ax,
    chrom: str,
    start: int,
    end: int,
    chrom_size: int,
    selection: pd.Series,
) -> None:
    """Chromosome-position strip modeled on the original cytoband header.

    The alternating blocks are an orientation strip, rather than invented
    cytoband calls.  The red tick is placed from the actual locus coordinate and
    therefore remains informative for the hg38 locus being shown.
    """
    import matplotlib.patches as patches
    prepare_row(ax)
    ax.set_ylim(0, 1)
    left, width = 0.02, 0.96
    y, height = 0.62, 0.18
    for index in range(12):
        shade = "#E6E6E6" if index % 2 == 0 else "#BFBFBF"
        ax.add_patch(
            patches.Rectangle(
                (left + index * width / 12.0, y),
                width / 12.0,
                height,
                transform=ax.transAxes,
                facecolor=shade,
                edgecolor="black",
                linewidth=0.35,
            )
        )
    midpoint_fraction = ((start + end) / 2.0) / chrom_size
    marker_width = max(0.012, (end - start) / chrom_size)
    marker_left = left + width * min(max(midpoint_fraction - marker_width / 2.0, 0.0), 1.0 - marker_width)
    ax.add_patch(
        patches.Rectangle(
            (marker_left, y - 0.025),
            width * marker_width,
            height + 0.05,
            transform=ax.transAxes,
            facecolor=COVID_COLOR,
            edgecolor="none",
            zorder=3,
        )
    )
    ax.text(
        0.01,
        0.97,
        f"{chrom}:{start:,}–{end:,} (hg38)",
        transform=ax.transAxes,
        fontsize=5.9,
        ha="left",
        va="top",
    )
    ax.text(
        0.99,
        0.42,
        f"BH q={format_q(float(selection['FDR']))}; δ={float(selection['Effect_size']):.2f}",
        transform=ax.transAxes,
        fontsize=4.55,
        ha="right",
        va="center",
        color="#333333",
    )
    # Put coordinate ticks *inside* the header.  Default Matplotlib tick labels
    # extend beyond an axes boundary and would otherwise overlap the next row.
    tick_positions = np.linspace(start, end, 4)
    for coordinate in tick_positions:
        fraction = (coordinate - start) / (end - start)
        ax.plot(
            [fraction, fraction], [0.13, 0.22], transform=ax.transAxes,
            color=TRACK_BORDER, linewidth=0.55, clip_on=False,
        )
        ax.text(
            fraction,
            0.025,
            f"{coordinate / 1000.0:,.0f} kb",
            transform=ax.transAxes,
            fontsize=4.35,
            ha="center",
            va="bottom",
        )
    ax.set_xlim(start / 1000.0, end / 1000.0)


def draw_coverage_row(
    ax,
    grid: np.ndarray,
    values: np.ndarray,
    start: int,
    end: int,
    vmax: float,
    group: str,
) -> None:
    """Render the cohort-mean eccDNA profile within one original-style boxed row."""
    prepare_row(ax)
    centres = (grid[:-1] + grid[1:]) / 2000.0
    color = COVID_COLOR if group == "COVID" else HC_COLOR
    label = "COVID-19 eccDNA" if group == "COVID" else "HC eccDNA"
    # Reserve a white header band for the track label and a white footer for
    # the scale.  This keeps every annotation clear of the coverage profile.
    baseline = 0.20
    height = baseline + 0.48 * values / vmax
    ax.fill_between(centres, baseline, height, step="mid", color=color, linewidth=0, alpha=0.98)
    ax.set_xlim(start / 1000.0, end / 1000.0)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.text(0.97, 0.90, label, transform=ax.transAxes, fontsize=5.35, ha="right", va="top")
    ax.text(0.012, 0.055, f"0–{vmax:.3g} EPM", transform=ax.transAxes, fontsize=3.75, ha="left", va="bottom")


def draw_peak_row(
    ax,
    intervals: Sequence[Tuple[int, int]],
    start: int,
    end: int,
    label: str,
    color: str,
) -> None:
    """Render a fine reference-H3K27ac peak track with a separate label band."""
    import matplotlib.patches as patches

    prepare_row(ax)
    ax.set_xlim(start / 1000.0, end / 1000.0)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    # A thin horizontal reference line makes sparse peaks legible without
    # turning the annotation rows into large solid blocks.
    ax.axhline(0.48, color="#D3D3D3", linewidth=0.35, zorder=1)
    for left, right in intervals:
        ax.add_patch(
            patches.Rectangle(
                (left / 1000.0, 0.42),
                (right - left) / 1000.0,
                0.12,
                facecolor=color,
                edgecolor="none",
                zorder=2,
            )
        )
    ax.text(0.97, 0.88, label, transform=ax.transAxes, fontsize=5.0, ha="right", va="top")


def draw_gene_row(ax, record: dict, start: int, end: int) -> None:
    """Draw a RefSeq-style gene model in the last framed row."""
    import matplotlib.pyplot as plt

    prepare_row(ax)
    ax.set_xlim(start / 1000.0, end / 1000.0)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    gene_start = max(record["start"], start) / 1000.0
    gene_end = min(record["end"], end) / 1000.0
    ax.plot([gene_start, gene_end], [0.48, 0.48], color=GENE_COLOR, linewidth=0.72)
    for exon_start, exon_end in record["exons"]:
        if exon_end <= start or exon_start >= end:
            continue
        left = max(exon_start, start) / 1000.0
        width = (min(exon_end, end) - max(exon_start, start)) / 1000.0
        ax.add_patch(plt.Rectangle((left, 0.36), width, 0.25, facecolor=GENE_COLOR, edgecolor="none"))
    span = gene_end - gene_start
    for fraction in np.linspace(0.08, 0.92, 11):
        x = gene_start + fraction * span
        delta = (0.013 if record["strand"] == "+" else -0.013) * span
        ax.annotate(
            "",
            xy=(x + delta, 0.48),
            xytext=(x, 0.48),
            arrowprops=dict(
                arrowstyle="-|>", color=GENE_COLOR, linewidth=0.45, mutation_scale=4.0,
                shrinkA=0, shrinkB=0,
            ),
        )
    ax.text(0.50, 0.10, record["gene"], transform=ax.transAxes, fontsize=6.0, style="italic", ha="center", va="bottom")
    ax.text(0.97, 0.88, "RefSeq", transform=ax.transAxes, fontsize=4.95, ha="right", va="top")


def add_locus_panel(
    fig,
    rect: Tuple[float, float, float, float],
    gene: str,
    panel_label: str,
    record: dict,
    selection: pd.Series,
    coverage: pd.DataFrame,
    dirs: base.Dirs,
    chrom_sizes: Dict[str, int],
) -> None:
    """Place one d/e locus browser in the original right-column geometry."""
    x0, y0, width, height = rect
    start, end, grid, signals = locus_arrays(coverage, gene, record)
    chrom = record["chrom"]
    if chrom not in chrom_sizes:
        raise RuntimeError(f"missing chromosome size for {chrom}")
    peak_map = {}
    expected = {
        "A549_GSE179184_H3K27ac_official": int(selection["A549_H3K27ac_peak_count"]),
        "CD14_H3K27ac": int(selection["CD14_H3K27ac_peak_count"]),
    }
    for peak_id, _label, _color in PEAK_TRACKS:
        path = dirs.peaks / f"{base.sanitize_id(peak_id)}.liftover.hg38.canonical.filtered.merged.bed"
        found = peaks_in_window(path, chrom, start, end)
        if len(found) != expected[peak_id]:
            raise RuntimeError(f"{gene}/{peak_id}: observed {len(found)} peaks, expected {expected[peak_id]}")
        peak_map[peak_id] = found

    # Fixed row proportions reproduce the visual rhythm of the supplied d/e
    # browser: coordinate overview, two cohort tracks, two chromatin rows, and
    # a RefSeq model.  Rows touch so their thin frames form one compact panel.
    fractions = (0.205, 0.168, 0.168, 0.152, 0.152, 0.155)
    top = y0 + height
    axes = []
    for fraction in fractions:
        row_height = height * fraction
        top -= row_height
        axes.append(fig.add_axes([x0, top, width, row_height], zorder=3))

    axes[0].text(-0.034, 1.04, panel_label, transform=axes[0].transAxes, fontsize=8.0, fontweight="bold", ha="left", va="bottom")
    draw_chromosome_overview(axes[0], chrom, start, end, chrom_sizes[chrom], selection)
    vmax = max(float(np.max(signals["COVID"])), float(np.max(signals["HC"])), 1e-12)
    draw_coverage_row(axes[1], grid, signals["COVID"], start, end, vmax, "COVID")
    draw_coverage_row(axes[2], grid, signals["HC"], start, end, vmax, "HC")
    for axis, (peak_id, label, color) in zip(axes[3:5], PEAK_TRACKS):
        draw_peak_row(axis, peak_map[peak_id], start, end, label, color)
    draw_gene_row(axes[5], record, start, end)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(SCRIPT_DIR.parent))
    parser.add_argument("--refgene", default=str(SCRIPT_DIR.parent / "inputs" / "figure5_refGene_hg38.txt"))
    parser.add_argument(
        "--coverage",
        default=str(SCRIPT_DIR.parent / "tables" / "figure5_screen_hlae_mcam_coverage.tsv"),
    )
    parser.add_argument(
        "--display-ranking",
        default=str(SCRIPT_DIR.parent.parent / "BCL3" / "tables" / "chromatin_display_candidate_ranking.tsv"),
    )
    parser.add_argument(
        "--locked-png",
        required=True,
        help="600-dpi rasterization of mei:/home/gao/eccDNA/revise/Figure_5.pdf.",
    )
    parser.add_argument("--stem", default="figure5_hg38")
    args = parser.parse_args(argv)

    dirs = base.ensure_dirs(Path(args.outdir))
    coverage = pd.read_csv(args.coverage, sep="\t")
    ranking = pd.read_csv(args.display_ranking, sep="\t")
    genes = [gene for gene, _label, _category in LOCUS_PANELS]
    selection = ranking[
        ranking["Selected_for_Figure_5_display"].eq("Yes") & ranking["Gene"].isin(genes)
    ].copy()
    if set(selection["Gene"]) != set(genes) or len(selection) != len(genes):
        raise RuntimeError("display-ranking table does not contain exactly HLA-E and MCAM")
    if not pd.to_numeric(selection["Chromatin_display_rank"], errors="coerce").eq(1).all():
        raise RuntimeError("a selected locus is not rank 1 in the display-eligible subset")
    if not selection["Dual_reference_H3K27ac_eligible"].eq("Yes").all():
        raise RuntimeError("a selected locus lacks dual-reference H3K27ac eligibility")
    if set(coverage["gene"]) != set(genes):
        raise RuntimeError("coverage input does not contain exactly HLA-E and MCAM")

    records = read_refgene(Path(args.refgene), genes)
    chrom_sizes = read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    background = load_locked_background(Path(args.locked_png))

    # Deposit exact display inputs alongside the final figure.
    coverage.to_csv(dirs.tables / "figure5_display_locus_coverage.tsv", sep="\t", index=False)
    selection.to_csv(dirs.tables / "figure5_display_locus_selection.tsv", sep="\t", index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base.setup_matplotlib()
    body_font = matplotlib.rcParams["font.family"]
    body_font = body_font[0] if isinstance(body_font, (list, tuple)) else body_font
    matplotlib.rcParams.update(
        {
            "font.family": body_font,
            "mathtext.fontset": "custom",
            "mathtext.rm": body_font,
            "mathtext.it": f"{body_font}:italic",
            "mathtext.bf": f"{body_font}:bold",
            "mathtext.default": "regular",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=FINAL_SIZE_INCH, dpi=600)
    fig.patch.set_facecolor("white")
    fig.figimage(background, xo=0, yo=0, zorder=0)

    # Mask every historical d/e pixel before laying down replacement panels.
    # Its left boundary lies safely right of locked panel c and its label.
    # A real axes (rather than a figure-level patch) is used because figimage
    # is painted above figure artists by some Matplotlib backends.  This blank
    # axes therefore reliably masks every old d/e pixel, including tick labels
    # that extend beyond their old axes boundaries.
    mask_ax = fig.add_axes([0.406, 0.025, 0.586, 0.955], zorder=1)
    mask_ax.set_facecolor("white")
    mask_ax.set_xticks([])
    mask_ax.set_yticks([])
    for spine in mask_ax.spines.values():
        spine.set_visible(False)
    positions = (
        (0.447, 0.555, 0.526, 0.400),
        (0.447, 0.060, 0.526, 0.400),
    )
    for (gene, panel_label, _category), rect in zip(LOCUS_PANELS, positions):
        row = selection[selection["Gene"] == gene]
        if len(row) != 1:
            raise RuntimeError(f"expected one selected row for {gene}, found {len(row)}")
        add_locus_panel(
            fig,
            rect,
            gene,
            panel_label,
            records[gene],
            row.iloc[0],
            coverage,
            dirs,
            chrom_sizes,
        )

    for extension in ("pdf", "svg", "png"):
        fig.savefig(dirs.figures / f"{args.stem}.{extension}", dpi=600 if extension == "png" else None)
    plt.close(fig)
    base.log(f"wrote figures/{args.stem}.{{pdf,svg,png}} with locked a--c and updated d/e")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
