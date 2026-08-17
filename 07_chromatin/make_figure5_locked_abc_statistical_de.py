#!/usr/bin/env python3
"""Render an evidence-focused Figure 5 while preserving mei panels a--c.

The approved 600-dpi mei Figure 5 export is used as an immutable background,
and every pixel in its left a--c region is retained.  The historical right
side is replaced with:

  d. the complete candidate-category landscape used for objective locus
     selection; and
  e. sample-level primary junction-EA distributions for HLA-E and MCAM,
     accompanied by recurrence, artifact-masked sensitivity, dominance and
     reference-H3K27ac annotations.

No locus-browser coverage shape or visual signal intensity enters selection.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import run_chromatin_hg38_analysis as base  # noqa: E402


LOCKED_PNG_SHA256 = "3bc52d33d016d44c0f879f63188f696ee988e92a3879be1d408c0511f40c28fe"
LOCKED_SIZE_PX = (4323, 3845)
FINAL_SIZE_INCH = (LOCKED_SIZE_PX[0] / 600.0, LOCKED_SIZE_PX[1] / 600.0)

SELECTED_GENES = ("HLA-E", "MCAM")
CATEGORIES = (
    "Humoral immune response",
    "Endothelial/coagulation",
)
CATEGORY_MARKERS = {
    "Humoral immune response": "o",
    "Endothelial/coagulation": "^",
}

# COVID/HC retain the exact semantic colours already locked into panels a--c.
COVID_COLOR = "#C84F50"
HC_COLOR = "#36617B"

# Additional roles use the exact Nature scientific-illustration palette.
GREY_LIGHT = "#C8CEDA"
GREY_MID = "#99A3B4"
GREY_DARK = "#49566D"
TEAL_LIGHT = "#CCE7EE"
TEAL = "#019AA3"
PURPLE = "#A84E94"
PURPLE_DARK = "#792C74"
GRID = "#E6E6ED"
TEXT = "#253247"


def load_locked_background(path: Path) -> np.ndarray:
    """Read and hash-check the approved mei a--c raster layer."""
    from PIL import Image

    payload = path.read_bytes()
    observed = hashlib.sha256(payload).hexdigest()
    if observed != LOCKED_PNG_SHA256:
        raise RuntimeError(
            "locked a--c image hash differs from the approved mei export: "
            f"{observed} != {LOCKED_PNG_SHA256}"
        )
    with Image.open(path) as image:
        rgba = np.asarray(image.convert("RGBA"))
    if (rgba.shape[1], rgba.shape[0]) != LOCKED_SIZE_PX:
        raise RuntimeError(
            "locked a--c image dimensions differ from the approved mei export: "
            f"{(rgba.shape[1], rgba.shape[0])} != {LOCKED_SIZE_PX}"
        )
    return rgba


def format_q(value: float) -> str:
    """Compact BH-q notation for the locus headers."""
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    return f"{mantissa:.2g}×10$^{{{exponent}}}$"


def point_area(detected_n: pd.Series | np.ndarray | float) -> np.ndarray:
    """Map COVID detection recurrence (5--39 of 39) to readable point area."""
    values = np.asarray(detected_n, dtype=float)
    return 8.0 + 0.58 * values


def read_inputs(
    ranking_path: Path,
    matrix_path: Path,
    metadata_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read and cross-check the candidate and sample-level inputs."""
    ranking = pd.read_csv(ranking_path, sep="\t")
    required = {
        "Category_for_row",
        "Gene",
        "FDR",
        "Effect_size",
        "COVID_detected_n",
        "HC_detected_n",
        "Prevalence_difference",
        "Max_sample_contribution",
        "quality_COVID_detected_n",
        "quality_HC_detected_n",
        "Strict_selection_eligible",
        "Dual_reference_H3K27ac_eligible",
        "Selected_for_Figure_5_display",
    }
    missing = required - set(ranking)
    if missing:
        raise RuntimeError(f"candidate table lacks columns: {', '.join(sorted(missing))}")
    if set(ranking["Category_for_row"]) != set(CATEGORIES):
        raise RuntimeError("candidate table contains unexpected functional categories")
    if len(ranking) != 281:
        raise RuntimeError(f"expected 281 candidate-category rows, found {len(ranking)}")
    if int(ranking["Strict_selection_eligible"].eq("Yes").sum()) != 149:
        raise RuntimeError("strict-eligible candidate count no longer equals 149")
    if int(ranking["Dual_reference_H3K27ac_eligible"].eq("Yes").sum()) != 54:
        raise RuntimeError("dual-reference H3K27ac candidate count no longer equals 54")
    selected = ranking.loc[
        ranking["Selected_for_Figure_5_display"].eq("Yes"), "Gene"
    ].tolist()
    if selected != list(SELECTED_GENES):
        raise RuntimeError(f"unexpected selected loci or order: {selected}")

    metadata = pd.read_csv(metadata_path, sep="\t")
    if not {"sample_id", "group"}.issubset(metadata):
        raise RuntimeError("sample metadata lacks sample_id/group")
    if metadata["sample_id"].duplicated().any():
        raise RuntimeError("duplicate sample IDs in metadata")
    group_counts = metadata["group"].value_counts().to_dict()
    if group_counts != {"COVID-19": 39, "HC": 39}:
        raise RuntimeError(f"expected 39 samples per cohort, found {group_counts}")

    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    missing_genes = set(SELECTED_GENES) - set(matrix.index)
    if missing_genes:
        raise RuntimeError(f"junction-EA matrix lacks: {', '.join(sorted(missing_genes))}")
    matrix_samples = list(matrix.columns)
    metadata_samples = metadata["sample_id"].tolist()
    if matrix_samples != metadata_samples:
        raise RuntimeError("junction-EA matrix and metadata sample orders differ")

    records = []
    group_map = metadata.set_index("sample_id")["group"].to_dict()
    for gene in SELECTED_GENES:
        values = pd.to_numeric(matrix.loc[gene], errors="raise")
        for sample, value in values.items():
            records.append(
                {
                    "Gene": gene,
                    "sample_id": sample,
                    "group": group_map[sample],
                    "primary_junction_EA": float(value),
                    "detected": int(float(value) > 0),
                }
            )
    samples = pd.DataFrame.from_records(records)

    # The displayed dots must reproduce every primary ranking-table quantity
    # that can be recomputed from the matrix.
    for gene in SELECTED_GENES:
        row = ranking.loc[
            ranking["Gene"].eq(gene)
            & ranking["Selected_for_Figure_5_display"].eq("Yes")
        ]
        if len(row) != 1:
            raise RuntimeError(f"expected one selected candidate row for {gene}")
        row = row.iloc[0]
        gene_samples = samples[samples["Gene"].eq(gene)]
        detected = (
            gene_samples.groupby("group")["detected"].sum().astype(int).to_dict()
        )
        expected = {
            "COVID-19": int(row["COVID_detected_n"]),
            "HC": int(row["HC_detected_n"]),
        }
        if detected != expected:
            raise RuntimeError(
                f"{gene}: primary sample detections {detected} != ranking table {expected}"
            )
        covid = gene_samples.loc[
            gene_samples["group"].eq("COVID-19"), "primary_junction_EA"
        ].to_numpy(float)
        contribution = float(np.max(covid) / np.sum(covid))
        if not np.isclose(
            contribution, float(row["Max_sample_contribution"]), rtol=1e-8, atol=1e-10
        ):
            raise RuntimeError(f"{gene}: maximum-sample contribution does not reproduce")
    return ranking, samples


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.55)
    ax.spines["bottom"].set_linewidth(0.55)
    ax.tick_params(labelsize=5.0, width=0.5, length=2.2, pad=1.5)
    ax.xaxis.label.set_size(5.6)
    ax.yaxis.label.set_size(5.6)


def add_panel_label(fig, rect: tuple[float, float, float, float], label: str) -> None:
    x0, y0, _width, height = rect
    fig.text(
        x0 - 0.022,
        y0 + height + 0.004,
        label,
        fontsize=8.0,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="black",
    )


def draw_candidate_panel(
    fig,
    rect: tuple[float, float, float, float],
    ranking: pd.DataFrame,
) -> None:
    """Panel d: all candidate-category assignments and display selection."""
    x0, y0, width, height = rect
    add_panel_label(fig, rect, "d")
    ax = fig.add_axes(
        [x0 + 0.085 * width, y0 + 0.17 * height, 0.66 * width, 0.75 * height],
        zorder=3,
    )
    frame = ranking.copy()
    frame["minus_log10_fdr"] = -np.log10(pd.to_numeric(frame["FDR"], errors="raise"))
    frame["Effect_size"] = pd.to_numeric(frame["Effect_size"], errors="raise")
    frame["COVID_detected_n"] = pd.to_numeric(
        frame["COVID_detected_n"], errors="raise"
    )

    states = (
        ("Not strict-eligible", ~frame["Strict_selection_eligible"].eq("Yes")),
        (
            "Strict-eligible",
            frame["Strict_selection_eligible"].eq("Yes")
            & ~frame["Dual_reference_H3K27ac_eligible"].eq("Yes"),
        ),
        (
            "Dual-reference H3K27ac",
            frame["Dual_reference_H3K27ac_eligible"].eq("Yes")
            & ~frame["Selected_for_Figure_5_display"].eq("Yes"),
        ),
        ("Displayed", frame["Selected_for_Figure_5_display"].eq("Yes")),
    )
    style = {
        "Not strict-eligible": dict(facecolors="white", edgecolors=GREY_LIGHT, linewidths=0.45, alpha=0.80),
        "Strict-eligible": dict(facecolors=GREY_MID, edgecolors="white", linewidths=0.35, alpha=0.72),
        "Dual-reference H3K27ac": dict(facecolors=TEAL_LIGHT, edgecolors=TEAL, linewidths=0.65, alpha=0.95),
        "Displayed": dict(facecolors=PURPLE, edgecolors=PURPLE_DARK, linewidths=0.75, alpha=1.0),
    }
    for state, state_mask in states:
        for category, marker in CATEGORY_MARKERS.items():
            subset = frame[state_mask & frame["Category_for_row"].eq(category)]
            if subset.empty:
                continue
            ax.scatter(
                subset["Effect_size"],
                subset["minus_log10_fdr"],
                s=point_area(subset["COVID_detected_n"]),
                marker=marker,
                zorder=5 if state == "Displayed" else 2,
                **style[state],
            )

    label_offsets = {
        "HLA-E": (-29, 13),
        "MCAM": (9, 8),
    }
    for gene in SELECTED_GENES:
        row = frame.loc[
            frame["Gene"].eq(gene)
            & frame["Selected_for_Figure_5_display"].eq("Yes")
        ].iloc[0]
        ax.annotate(
            gene,
            xy=(float(row["Effect_size"]), float(row["minus_log10_fdr"])),
            xytext=label_offsets[gene],
            textcoords="offset points",
            fontsize=5.5,
            fontstyle="italic",
            fontweight="bold",
            color=PURPLE_DARK,
            ha="left",
            va="center",
            arrowprops=dict(
                arrowstyle="-",
                color=PURPLE_DARK,
                linewidth=0.55,
                shrinkA=2,
                shrinkB=2,
            ),
            bbox=dict(
                boxstyle="round,pad=0.12",
                facecolor="white",
                edgecolor="none",
                alpha=0.92,
            ),
            zorder=8,
        )

    ax.set_xlim(0.10, 0.66)
    ax.set_ylim(1.25, 4.48)
    ax.set_xlabel("Cliff's δ (COVID-19 vs HC)")
    ax.set_ylabel("−log$_{10}$(BH FDR)")
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    ax.set_yticks([1.3, 2.0, 3.0, 4.0])
    ax.grid(axis="y", color=GRID, linewidth=0.45, zorder=0)
    style_axis(ax)
    ax.text(
        0.0,
        1.035,
        "All candidate–category assignments",
        transform=ax.transAxes,
        fontsize=5.6,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=TEXT,
    )
    legend_ax = fig.add_axes(
        [x0 + 0.77 * width, y0 + 0.17 * height, 0.22 * width, 0.75 * height],
        zorder=3,
    )
    legend_ax.axis("off")
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.text(
        0.0,
        1.03,
        "Key",
        fontsize=5.4,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=TEXT,
    )
    legend_ax.text(
        0.0,
        0.95,
        "n=281 assignments\nstrict=149; dual=54",
        fontsize=5.0,
        ha="left",
        va="top",
        color=GREY_DARK,
    )
    key_rows = [
        (0.78, "o", GREY_MID, "white", "Humoral\nimmune"),
        (0.63, "^", GREY_MID, "white", "Endothelial/\ncoagulation"),
        (0.48, "o", "white", GREY_LIGHT, "Not strict-\neligible"),
        (0.35, "o", GREY_MID, "white", "Strict-\neligible"),
        (0.22, "o", TEAL_LIGHT, TEAL, "Dual-reference\nH3K27ac"),
        (0.08, "o", PURPLE, PURPLE_DARK, "Displayed\nlocus"),
    ]
    for y, marker, face, edge, label in key_rows:
        legend_ax.scatter(
            [0.12],
            [y],
            s=20,
            marker=marker,
            facecolors=face,
            edgecolors=edge,
            linewidths=0.6,
        )
        legend_ax.text(
            0.27,
            y,
            label,
            fontsize=5.0,
            ha="left",
            va="center",
            color="black",
            linespacing=0.90,
        )
    fig.text(
        x0 + 0.50 * width,
        y0 + 0.018 * height,
        "Point size indicates COVID-19 detection recurrence (5–39 of 39); visual track appearance was not used.",
        fontsize=5.0,
        ha="center",
        va="bottom",
        color=GREY_DARK,
    )


def draw_header(
    ax,
    gene: str,
    row: pd.Series,
) -> None:
    """A separate annotation strip so no text overlaps sample points."""
    ax.axis("off")
    ax.text(
        0.5,
        0.91,
        gene,
        fontsize=6.5,
        fontstyle="italic",
        fontweight="bold",
        ha="center",
        va="top",
        color=TEXT,
    )
    ax.text(
        0.5,
        0.68,
        f"BH q={format_q(float(row['FDR']))}; δ={float(row['Effect_size']):.2f}",
        fontsize=5.0,
        ha="center",
        va="top",
        color="black",
    )
    ax.text(
        0.5,
        0.46,
        (
            f"Detected: COVID-19 {int(row['COVID_detected_n'])}/39  |  "
            f"HC {int(row['HC_detected_n'])}/39"
        ),
        fontsize=5.0,
        ha="center",
        va="top",
        color="black",
    )
    ax.text(
        0.5,
        0.25,
        (
            f"Masked: {int(row['quality_COVID_detected_n'])}/39  |  "
            f"{int(row['quality_HC_detected_n'])}/39;  "
            f"max sample {100 * float(row['Max_sample_contribution']):.1f}%"
        ),
        fontsize=5.0,
        ha="center",
        va="top",
        color=GREY_DARK,
    )
    ax.text(
        0.5,
        0.04,
        "●  dual-reference H3K27ac overlap",
        fontsize=5.0,
        ha="center",
        va="bottom",
        color=TEAL,
    )


def draw_sample_axis(
    ax,
    frame: pd.DataFrame,
    gene: str,
    show_ylabel: bool,
) -> None:
    """Primary junction-EA distribution with every zero and non-zero sample."""
    groups = ("COVID-19", "HC")
    transformed = []
    for group in groups:
        values = frame.loc[
            frame["Gene"].eq(gene) & frame["group"].eq(group),
            "primary_junction_EA",
        ].to_numpy(float)
        if len(values) != 39:
            raise RuntimeError(f"{gene}/{group}: expected 39 samples, found {len(values)}")
        transformed.append(np.log10(values + 1.0))

    bp = ax.boxplot(
        transformed,
        positions=[0, 1],
        widths=0.50,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 0.75},
        boxprops={"color": "black", "linewidth": 0.55},
        whiskerprops={"color": "black", "linewidth": 0.55},
        capprops={"color": "black", "linewidth": 0.55},
    )
    for patch, color in zip(bp["boxes"], (COVID_COLOR, HC_COLOR)):
        patch.set_facecolor(color)
        patch.set_alpha(0.24)

    for index, (group, values, color) in enumerate(
        zip(groups, transformed, (COVID_COLOR, HC_COLOR))
    ):
        seed = 5100 + index + (0 if gene == "HLA-E" else 20)
        rng = np.random.default_rng(seed)
        ax.scatter(
            np.full(len(values), index) + rng.uniform(-0.13, 0.13, len(values)),
            values,
            s=5.2,
            facecolors=color,
            edgecolors="white",
            linewidths=0.18,
            alpha=0.90,
            zorder=3,
        )

    tick_values = np.array([0, 10, 100, 1000], dtype=float)
    ax.set_yticks(np.log10(tick_values + 1.0))
    ax.set_yticklabels(["0", "10", "100", "1,000"])
    ax.set_ylim(-0.08, 3.27)
    ax.set_xlim(-0.52, 1.52)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(groups)
    ax.grid(axis="y", color=GRID, linewidth=0.45, zorder=0)
    if show_ylabel:
        ax.set_ylabel("Junction-level EA\n(log$_{10}$[EA + 1])")
    else:
        ax.tick_params(labelleft=False)
        ax.spines["left"].set_visible(False)
    style_axis(ax)


def draw_sample_panel(
    fig,
    rect: tuple[float, float, float, float],
    ranking: pd.DataFrame,
    samples: pd.DataFrame,
) -> None:
    """Panel e: selected-locus sample distributions and audit annotations."""
    x0, y0, width, height = rect
    add_panel_label(fig, rect, "e")
    fig.text(
        x0,
        y0 + height + 0.005,
        "Selected loci: sample-level primary abundance",
        fontsize=5.6,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=TEXT,
    )

    facet_width = 0.435 * width
    facet_x = (x0 + 0.065 * width, x0 + 0.555 * width)
    for index, (gene, left) in enumerate(zip(SELECTED_GENES, facet_x)):
        row = ranking.loc[
            ranking["Gene"].eq(gene)
            & ranking["Selected_for_Figure_5_display"].eq("Yes")
        ].iloc[0]
        header = fig.add_axes(
            [left, y0 + 0.73 * height, facet_width, 0.25 * height],
            zorder=3,
        )
        draw_header(header, gene, row)
        plot = fig.add_axes(
            [left, y0 + 0.20 * height, facet_width, 0.48 * height],
            zorder=3,
        )
        draw_sample_axis(plot, samples, gene, show_ylabel=index == 0)

    fig.text(
        x0 + 0.50 * width,
        y0 + 0.105 * height,
        "Points are individual samples; q and δ are from the primary junction-level analysis.",
        fontsize=5.0,
        ha="center",
        va="bottom",
        color=GREY_DARK,
    )
    fig.text(
        x0 + 0.50 * width,
        y0 + 0.040 * height,
        "Masked = prespecified breakpoint-quality sensitivity; H3K27ac annotations are reference-cell data.",
        fontsize=5.0,
        ha="center",
        va="bottom",
        color=GREY_DARK,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(SCRIPT_DIR.parent))
    parser.add_argument("--ranking", required=True, type=Path)
    parser.add_argument("--primary-ea", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--locked-png", required=True, type=Path)
    parser.add_argument(
        "--locked-pdf",
        required=True,
        type=Path,
        help="Vector source for the immutable a--c panels.",
    )
    parser.add_argument("--stem", default="figure5_statistical_de_v1")
    args = parser.parse_args(argv)

    dirs = base.ensure_dirs(Path(args.outdir))
    background = load_locked_background(args.locked_png)
    ranking, samples = read_inputs(args.ranking, args.primary_ea, args.metadata)

    candidate_columns = [
        "Category_for_row",
        "Gene",
        "FDR",
        "Effect_size",
        "COVID_detected_n",
        "HC_detected_n",
        "Prevalence_difference",
        "Max_sample_contribution",
        "quality_FDR",
        "quality_effect_size",
        "quality_COVID_detected_n",
        "quality_HC_detected_n",
        "Strict_selection_eligible",
        "Dual_reference_H3K27ac_eligible",
        "Selected_for_Figure_5_display",
    ]
    ranking[candidate_columns].to_csv(
        dirs.tables / f"{args.stem}_candidates.tsv", sep="\t", index=False
    )
    samples.to_csv(
        dirs.tables / f"{args.stem}_sample_EA.tsv", sep="\t", index=False
    )

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
            "svg.fonttype": "none",
        }
    )
    def draw_replacement(fig, with_white_mask: bool) -> None:
        if with_white_mask:
            # The historical panel label begins just left of the former 40.6%
            # lock boundary.  Starting at 38.8% removes it while retaining the
            # complete approved a--c plotting region.
            mask_ax = fig.add_axes([0.388, 0.020, 0.607, 0.965], zorder=1)
            mask_ax.set_facecolor("white")
            mask_ax.set_xticks([])
            mask_ax.set_yticks([])
            for spine in mask_ax.spines.values():
                spine.set_visible(False)
        draw_candidate_panel(fig, (0.447, 0.545, 0.526, 0.395), ranking)
        draw_sample_panel(fig, (0.447, 0.060, 0.526, 0.405), ranking, samples)

    # A 600-dpi whole-figure preview is retained for Word embedding and visual
    # QA.  The submission PDF is assembled separately from vector layers.
    preview = plt.figure(figsize=FINAL_SIZE_INCH, dpi=600)
    preview.patch.set_facecolor("white")
    preview.figimage(background, xo=0, yo=0, zorder=0)
    draw_replacement(preview, with_white_mask=True)
    preview.savefig(dirs.figures / f"{args.stem}.png", dpi=600)
    plt.close(preview)

    # Build the publication PDF from the locked vector a--c source plus a
    # vector statistical overlay.  Cropping the source prevents obsolete
    # right-panel text from remaining invisibly embedded beneath a white box.
    from pypdf import PdfReader, PdfWriter
    from pypdf._page import PageObject

    with tempfile.TemporaryDirectory(prefix="figure5_overlay_") as tmp:
        tmpdir = Path(tmp)
        overlay_path = Path(tmp) / "overlay.pdf"
        overlay = plt.figure(figsize=FINAL_SIZE_INCH)
        overlay.patch.set_alpha(0.0)
        draw_replacement(overlay, with_white_mask=False)
        overlay.savefig(overlay_path, transparent=True)
        plt.close(overlay)

        source_reader = PdfReader(str(args.locked_pdf))
        overlay_reader = PdfReader(str(overlay_path))
        if len(source_reader.pages) != 1 or len(overlay_reader.pages) != 1:
            raise RuntimeError("Figure 5 vector inputs must each contain one page")
        source_page = source_reader.pages[0]
        width = float(source_page.mediabox.width)
        height = float(source_page.mediabox.height)
        crop_width = 0.388 * width
        bounded_crop_path = tmpdir / "left_bounded_crop.pdf"
        # A PDF CropBox is only a viewing hint: text extractors can still see
        # clipped historical d/e content.  Rendering the source through a
        # fixed-size device establishes the intended a--c boundary.
        subprocess.run(
            [
                "gs",
                "-q",
                "-dSAFER",
                "-dBATCH",
                "-dNOPAUSE",
                "-sDEVICE=pdfwrite",
                "-dAutoRotatePages=/None",
                "-dFIXEDMEDIA",
                f"-dDEVICEWIDTHPOINTS={crop_width:.6f}",
                f"-dDEVICEHEIGHTPOINTS={height:.6f}",
                f"-sOutputFile={bounded_crop_path}",
                str(args.locked_pdf),
            ],
            check=True,
        )
        # pdfwrite may retain out-of-bounds drawing operators even on a
        # physically smaller page.  The vector PDF -> EPS -> vector PDF
        # round-trip below resolves the page boundary and removes those
        # operators, rather than merely hiding them from view.  This prevents
        # obsolete historical d/e labels from leaking into PDF text exports.
        eps_crop_path = tmpdir / "left_clean_crop.eps"
        subprocess.run(
            [
                "pdftops",
                "-eps",
                "-rasterize",
                "never",
                str(bounded_crop_path),
                str(eps_crop_path),
            ],
            check=True,
        )
        clean_crop_path = tmpdir / "left_clean_crop.pdf"
        subprocess.run(
            [
                "gs",
                "-q",
                "-dSAFER",
                "-dBATCH",
                "-dNOPAUSE",
                "-sDEVICE=pdfwrite",
                "-dEPSCrop",
                "-dAutoRotatePages=/None",
                f"-sOutputFile={clean_crop_path}",
                str(eps_crop_path),
            ],
            check=True,
        )
        clean_source_page = PdfReader(str(clean_crop_path)).pages[0]
        output_page = PageObject.create_blank_page(width=width, height=height)
        output_page.merge_page(clean_source_page)
        output_page.merge_page(overlay_reader.pages[0])
        writer = PdfWriter()
        writer.add_page(output_page)
        pdf_path = dirs.figures / f"{args.stem}.pdf"
        with pdf_path.open("wb") as handle:
            writer.write(handle)

    svg_path = dirs.figures / f"{args.stem}.svg"
    subprocess.run(
        ["pdftocairo", "-svg", str(pdf_path), str(svg_path)],
        check=True,
    )
    base.log(
        f"wrote figures/{args.stem}.{{pdf,svg,png}} with vector a--c and statistical d/e"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
