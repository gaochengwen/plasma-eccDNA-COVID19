#!/usr/bin/env python3
"""Create the revised, evidence-bounded graphical abstract.

The artwork summarizes only analyses present in the revision package.  It
deliberately distinguishes observed plasma-eccDNA associations from
mechanistic or clinical claims that were not tested in this study.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


MM = 1 / 25.4

# Exact Nature-branded colours used for newly created artwork.
RED = "#C93E3F"
BLUE = "#0272B2"
TEAL = "#019AA3"
TEAL_LIGHT = "#CCE7EE"
PURPLE = "#A84E94"
PURPLE_LIGHT = "#E5C4D9"
NAVY = "#253247"
GREY_DARK = "#49566D"
GREY_MID = "#99A3B4"
GREY_LIGHT = "#E6E6ED"
PAPER = "#F7F8FB"
WHITE = "#FFFFFF"


def rounded_box(
    ax,
    xy: tuple[float, float],
    width: float,
    height: float,
    *,
    facecolor: str,
    edgecolor: str,
    linewidth: float = 0.7,
    radius: float = 0.018,
    zorder: int = 1,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.010,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def panel_label(ax, x: float, y: float, label: str) -> None:
    ax.text(
        x,
        y,
        label,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="top",
        color="black",
    )


def draw_people(ax, centre_x: float, y: float, color: str) -> None:
    """Three simple cohort glyphs, constructed as vectors."""
    for offset, scale in ((-0.034, 0.90), (0, 1.08), (0.034, 0.90)):
        x = centre_x + offset
        ax.add_patch(
            Circle(
                (x, y + 0.028 * scale),
                0.010 * scale,
                facecolor=color,
                edgecolor="none",
            )
        )
        ax.add_patch(
            FancyBboxPatch(
                (x - 0.014 * scale, y - 0.021 * scale),
                0.028 * scale,
                0.039 * scale,
                boxstyle=f"round,pad=0.001,rounding_size={0.009 * scale}",
                facecolor=color,
                edgecolor="none",
            )
        )


def draw_plasma_tube(ax, x: float, y: float) -> None:
    rounded_box(
        ax,
        (x - 0.025, y - 0.052),
        0.050,
        0.100,
        facecolor=WHITE,
        edgecolor=GREY_DARK,
        linewidth=0.8,
        radius=0.006,
    )
    ax.add_patch(
        FancyBboxPatch(
            (x - 0.021, y - 0.045),
            0.042,
            0.049,
            boxstyle="round,pad=0.001,rounding_size=0.004",
            facecolor="#F1C96A",
            edgecolor="none",
        )
    )
    ax.plot([x - 0.027, x + 0.027], [y + 0.050, y + 0.050], color=GREY_DARK, lw=1.4)


def arrow(ax, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.8,
            color=GREY_MID,
            shrinkA=1,
            shrinkB=1,
        )
    )


def observation_card(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    title: str,
    color: str,
    lines: tuple[str, ...],
) -> None:
    rounded_box(
        ax,
        (x, y),
        width,
        height,
        facecolor=WHITE,
        edgecolor=GREY_LIGHT,
        linewidth=0.65,
        radius=0.012,
    )
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            0.010,
            height,
            boxstyle="round,pad=0,rounding_size=0.006",
            facecolor=color,
            edgecolor="none",
        )
    )
    ax.text(
        x + 0.025,
        y + height - 0.025,
        title,
        fontsize=6.2,
        fontweight="bold",
        color=NAVY,
        ha="left",
        va="top",
    )
    ax.text(
        x + 0.025,
        y + height - 0.061,
        "\n".join(lines),
        fontsize=5.2,
        color=GREY_DARK,
        ha="left",
        va="top",
        linespacing=1.25,
    )


def create_figure(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(183 * MM, 120 * MM), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.035,
        0.955,
        "Plasma eccDNA profiling in COVID-19",
        fontsize=11,
        fontweight="bold",
        color=NAVY,
        ha="left",
        va="top",
    )
    ax.text(
        0.035,
        0.900,
        "Study design, reproducible observations and interpretation boundary",
        fontsize=6.4,
        color=GREY_DARK,
        ha="left",
        va="top",
    )
    ax.plot([0.035, 0.965], [0.865, 0.865], color=GREY_LIGHT, lw=0.8)

    # a: study design
    panel_label(ax, 0.035, 0.835, "a")
    ax.text(
        0.070,
        0.833,
        "Study and analysis",
        fontsize=7.2,
        fontweight="bold",
        color=NAVY,
        ha="left",
        va="top",
    )
    rounded_box(
        ax,
        (0.035, 0.190),
        0.275,
        0.595,
        facecolor=PAPER,
        edgecolor=GREY_LIGHT,
        linewidth=0.7,
        radius=0.018,
    )
    draw_people(ax, 0.104, 0.675, RED)
    draw_people(ax, 0.239, 0.675, BLUE)
    ax.text(0.104, 0.610, "COVID-19\nn = 39", fontsize=5.7, ha="center", va="top", color=NAVY)
    ax.text(0.239, 0.610, "Healthy controls\nn = 39", fontsize=5.7, ha="center", va="top", color=NAVY)
    ax.text(0.172, 0.535, "plasma", fontsize=5.5, ha="center", va="center", color=GREY_DARK)
    arrow(ax, (0.104, 0.595), (0.154, 0.528))
    arrow(ax, (0.239, 0.595), (0.190, 0.528))
    draw_plasma_tube(ax, 0.172, 0.462)
    arrow(ax, (0.172, 0.398), (0.172, 0.355))
    rounded_box(
        ax,
        (0.067, 0.255),
        0.210,
        0.085,
        facecolor=WHITE,
        edgecolor=TEAL,
        linewidth=0.8,
        radius=0.012,
    )
    ax.text(
        0.172,
        0.312,
        "Circle-seq → hg38 eccDNA calls",
        fontsize=5.7,
        fontweight="bold",
        ha="center",
        va="center",
        color=NAVY,
    )
    ax.text(
        0.172,
        0.278,
        "7 prespecified call sets · age-adjusted checks",
        fontsize=4.9,
        ha="center",
        va="center",
        color=GREY_DARK,
    )

    arrow(ax, (0.315, 0.487), (0.347, 0.487))

    # b: results
    panel_label(ax, 0.345, 0.835, "b")
    ax.text(
        0.380,
        0.833,
        "Reproducible observations",
        fontsize=7.2,
        fontweight="bold",
        color=NAVY,
        ha="left",
        va="top",
    )
    rounded_box(
        ax,
        (0.345, 0.190),
        0.340,
        0.595,
        facecolor=PAPER,
        edgecolor=GREY_LIGHT,
        linewidth=0.7,
        radius=0.018,
    )
    observation_card(
        ax,
        0.372,
        0.605,
        0.286,
        0.135,
        title="Higher plasma eccDNA abundance",
        color=RED,
        lines=(
            "COVID-19/HC median ratio 3.02–3.86",
            "Age-adjusted ratio 3.31–5.08 across 7 call sets",
        ),
    )
    observation_card(
        ax,
        0.372,
        0.425,
        0.286,
        0.135,
        title="Stable fragment-length modes",
        color=TEAL,
        lines=(
            "Objective consensus peaks: 196, 366 and 571 bp",
            "Peak positions remained stable in sensitivity analyses",
        ),
    )
    observation_card(
        ax,
        0.372,
        0.245,
        0.286,
        0.135,
        title="Genomic-context associations",
        color=PURPLE,
        lines=(
            "Positional eccGene and pathway associations",
            "Overlap with reference-cell chromatin annotations",
        ),
    )

    arrow(ax, (0.690, 0.487), (0.722, 0.487))

    # c: interpretation
    panel_label(ax, 0.720, 0.835, "c")
    ax.text(
        0.755,
        0.833,
        "Interpretation boundary",
        fontsize=7.2,
        fontweight="bold",
        color=NAVY,
        ha="left",
        va="top",
    )
    rounded_box(
        ax,
        (0.720, 0.190),
        0.245,
        0.595,
        facecolor=PAPER,
        edgecolor=GREY_LIGHT,
        linewidth=0.7,
        radius=0.018,
    )
    ax.text(
        0.842,
        0.716,
        "Supported",
        fontsize=6.4,
        fontweight="bold",
        ha="center",
        va="center",
        color=TEAL,
    )
    ax.text(
        0.752,
        0.667,
        "●  Cross-sectional plasma-eccDNA associations\n"
        "●  Robustness across prespecified definitions\n"
        "●  Reference-annotation overlap",
        fontsize=5.2,
        ha="left",
        va="top",
        color=NAVY,
        linespacing=1.45,
    )
    ax.plot([0.750, 0.935], [0.548, 0.548], color=GREY_LIGHT, lw=0.8)
    ax.text(
        0.842,
        0.510,
        "Not tested",
        fontsize=6.4,
        fontweight="bold",
        ha="center",
        va="center",
        color=RED,
    )
    ax.text(
        0.752,
        0.462,
        "×  Patient-matched chromatin state\n"
        "×  Molecular mechanism or tissue origin\n"
        "×  Disease specificity or clinical utility",
        fontsize=5.2,
        ha="left",
        va="top",
        color=NAVY,
        linespacing=1.45,
    )
    ax.text(
        0.842,
        0.285,
        "Requires molecular, external,\n"
        "disease-control and clinical validation",
        fontsize=5.3,
        fontweight="bold",
        ha="center",
        va="center",
        color=GREY_DARK,
        linespacing=1.25,
    )

    rounded_box(
        ax,
        (0.125, 0.070),
        0.750,
        0.072,
        facecolor=TEAL_LIGHT,
        edgecolor=TEAL,
        linewidth=0.8,
        radius=0.014,
    )
    ax.text(
        0.500,
        0.106,
        "A descriptive plasma-eccDNA map for hypothesis generation and independent validation",
        fontsize=6.4,
        fontweight="bold",
        color=NAVY,
        ha="center",
        va="center",
    )

    stem = output_dir / "Graphical_Abstract"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches=None)
    fig.savefig(stem.with_suffix(".jpeg"), dpi=300, bbox_inches=None, pil_kwargs={"quality": 95})
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    create_figure(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
