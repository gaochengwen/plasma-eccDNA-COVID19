#!/usr/bin/env python3
"""Create revised eccGene figures and their machine-readable source data."""

from __future__ import annotations

import csv
import json
import math
import textwrap
import urllib.request
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analyse" / "eccGene"
REVISE = ROOT / "revise"
RESULTS = ANALYSIS / "results"
MATRICES = ANALYSIS / "matrices"
SOURCE_DATA = REVISE / "source_data"

COVID = "#C93E3F"
HC = "#0272B2"
ORANGE = "#EC6F00"
GREEN = "#459434"
PURPLE = "#A84E94"
TEAL = "#019AA3"
YELLOW = "#CCA02C"
GREY_1 = "#E6E6ED"
GREY_2 = "#C8CEDA"
GREY_3 = "#99A3B4"
GREY_4 = "#6F7B91"
GREY_5 = "#49566D"
STONE_1 = "#F7F5EF"

MM = 1 / 25.4

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,
        "axes.linewidth": 0.55,
        "axes.labelsize": 6,
        "axes.titlesize": 6.5,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial",
        "mathtext.default": "regular",
        "savefig.transparent": False,
    }
)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix, kwargs in (
        (".pdf", {}),
        (".svg", {}),
        (".png", {"dpi": 300}),
    ):
        fig.savefig(REVISE / f"{stem}{suffix}", **kwargs)


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.05) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def prepare_abundance_panel() -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    stability = read_tsv(RESULTS / "candidate_stability.tsv")
    abundance = read_tsv(RESULTS / "junction_abundance_wilcoxon.tsv")
    abundance_lookup = abundance.set_index("Gene")

    selected = stability[
        (stability["abundance_significant_same_direction_all_definitions"] == 1)
        & (stability["junction_abundance_log2FC"] >= 1)
    ].copy()
    selected["cliffs_delta"] = selected["Gene"].map(
        abundance_lookup["cliffs_delta_COVID_vs_HC"]
    )
    selected = selected[selected["cliffs_delta"] > 0]
    selected = selected.sort_values(
        ["junction_abundance_q", "cliffs_delta", "Gene"],
        ascending=[True, False, True],
    ).head(30)

    metadata = read_tsv(ANALYSIS / "sample_metadata.tsv")
    sample_order = metadata["sample_id"].tolist()
    groups = metadata.set_index("sample_id")["group"].to_dict()
    covid_samples = [sample for sample in sample_order if groups[sample] == "COVID-19"]
    hc_samples = [sample for sample in sample_order if groups[sample] == "HC"]
    ordered_samples = covid_samples + hc_samples

    ea = read_tsv(MATRICES / "junction_EA.tsv").set_index("Gene")
    values = ea.loc[selected["Gene"], ordered_samples].astype(float)
    transformed = np.log1p(values.to_numpy())
    means = transformed.mean(axis=1, keepdims=True)
    stds = transformed.std(axis=1, ddof=0, keepdims=True)
    stds[stds == 0] = 1
    z = (transformed - means) / stds
    z_frame = pd.DataFrame(z, index=values.index, columns=ordered_samples)

    export = selected[
        [
            "Gene",
            "junction_abundance_q",
            "junction_abundance_log2FC",
            "cliffs_delta",
            "interval_abundance_q",
            "interval_abundance_log2FC",
            "midpoint_abundance_q",
            "midpoint_abundance_log2FC",
        ]
    ].copy()
    export = export.rename(
        columns={
            "cliffs_delta": "junction_cliffs_delta",
            "junction_abundance_q": "junction_wilcoxon_q_BH",
            "interval_abundance_q": "interval_wilcoxon_q_BH",
            "midpoint_abundance_q": "midpoint_wilcoxon_q_BH",
        }
    )
    write_tsv(SOURCE_DATA / "Figure_3A_abundance_genes.tsv", export)
    z_frame.reset_index().to_csv(
        SOURCE_DATA / "Figure_3A_log1p_EA_row_zscore.tsv", sep="\t", index=False
    )
    return z_frame, export, covid_samples, hc_samples


def prepare_detection_panel() -> pd.DataFrame:
    stability = read_tsv(RESULTS / "candidate_stability.tsv")
    detection = read_tsv(RESULTS / "junction_detection_fisher.tsv")
    stable_genes = set(
        stability.loc[
            stability["detection_significant_same_direction_all_definitions"] == 1,
            "Gene",
        ]
    )
    selected = detection[
        (detection["Gene"].isin(stable_genes))
        & (detection["fisher_q_BH"] < 0.05)
        & (detection["frequency_difference_COVID_minus_HC"] > 0)
    ].copy()
    selected = selected.sort_values(
        ["fisher_q_BH", "frequency_difference_COVID_minus_HC", "Gene"],
        ascending=[True, False, True],
    ).head(15)
    write_tsv(SOURCE_DATA / "Figure_3B_detection_genes.tsv", selected)
    return selected


EXPECT_QUERY = None
EXPECT_BACKGROUND = None


def gprofiler_enrichment() -> tuple[pd.DataFrame, dict]:
    abundance = read_tsv(RESULTS / "junction_abundance_wilcoxon.tsv")
    query = abundance.loc[
        (abundance["wilcoxon_q_BH"] < 0.05)
        & (abundance["log2FC_mean_EA_COVID_vs_HC"] >= 1)
        & (abundance["cliffs_delta_COVID_vs_HC"] > 0),
        "Gene",
    ].tolist()
    background = abundance["Gene"].tolist()
    payload = {
        "organism": "hsapiens",
        "query": query,
        "sources": ["GO:BP", "REAC", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": background,
        "no_evidences": True,
        "highlight": True,
    }
    request_path = SOURCE_DATA / "eccGene_gProfiler_request.json"
    response_path = SOURCE_DATA / "eccGene_gProfiler_response.json"
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    request_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    (SOURCE_DATA / "eccGene_gProfiler_query_genes.txt").write_text(
        "\n".join(query) + "\n", encoding="utf-8"
    )
    (SOURCE_DATA / "eccGene_gProfiler_background_genes.txt").write_text(
        "\n".join(background) + "\n", encoding="utf-8"
    )

    if response_path.exists():
        response = json.loads(response_path.read_text(encoding="utf-8"))
    else:
        req = urllib.request.Request(
            "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            data=json.dumps(payload).encode("utf-8"),
            headers={
                "Content-Type": "application/json",
                "User-Agent": "eccDNA-revision/1.0",
            },
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=180) as handle:
            response = json.loads(handle.read().decode("utf-8"))
        response_path.write_text(
            json.dumps(response, indent=2) + "\n", encoding="utf-8"
        )

    results = pd.DataFrame(response.get("result", []))
    if results.empty:
        raise RuntimeError("g:Profiler returned no significant enrichment terms")
    # v1 pinned the set sizes it was built with (query=3094, background=26642)
    # and aborted on anything else. Under the methods call set the differential
    # gene list legitimately differs, so the pin is recorded rather than
    # enforced; EXPECT_QUERY/EXPECT_BACKGROUND can re-pin it once the v2 sizes
    # are themselves locked.
    print(f"enrichment set sizes: query={len(query)}, background={len(background)}",
          flush=True)
    if EXPECT_QUERY is not None and (
        len(query) != EXPECT_QUERY or len(background) != EXPECT_BACKGROUND
    ):
        raise RuntimeError(
            f"Unexpected enrichment set sizes: query={len(query)}, "
            f"background={len(background)}"
        )
    results["minus_log10_fdr"] = -np.log10(results["p_value"])
    results["submitted_query_gene_count"] = len(query)
    results["submitted_background_gene_count"] = len(background)
    all_term_columns = [
        "source",
        "native",
        "name",
        "p_value",
        "minus_log10_fdr",
        "intersection_size",
        "term_size",
        "query_size",
        "effective_domain_size",
        "precision",
        "recall",
        "group_id",
        "highlighted",
        "submitted_query_gene_count",
        "submitted_background_gene_count",
    ]
    write_tsv(
        SOURCE_DATA / "eccGene_gProfiler_all_significant_terms.tsv",
        results.sort_values(["p_value", "source", "native"])[all_term_columns],
    )

    representatives = (
        results.sort_values(["p_value", "source", "native"])
        .groupby("group_id", as_index=False, sort=False)
        .first()
        .sort_values(["p_value", "source", "native"])
        .head(12)
        .copy()
    )
    representatives["minus_log10_fdr"] = -np.log10(representatives["p_value"])
    representatives["query_gene_count"] = len(query)
    representatives["background_gene_count"] = len(background)
    keep = [
        "source",
        "native",
        "name",
        "p_value",
        "minus_log10_fdr",
        "intersection_size",
        "term_size",
        "query_size",
        "effective_domain_size",
        "group_id",
        "query_gene_count",
        "background_gene_count",
    ]
    write_tsv(
        SOURCE_DATA / "Figure_3C_gProfiler_representative_terms.tsv",
        representatives[keep],
    )
    return representatives, response


def prepare_circle_panel() -> pd.DataFrame:
    path = SOURCE_DATA / "covid_specific_exact.matrix.min10.tsv"
    circles = pd.read_csv(path, sep="\t", index_col=0)
    circles = circles.astype(int)
    recurrence = circles.sum(axis=1)
    column_total = circles.sum(axis=0)
    circles = circles.loc[
        recurrence.sort_values(ascending=False, kind="stable").index,
        column_total.sort_values(ascending=False, kind="stable").index,
    ]
    export = circles.copy()
    export.insert(0, "recurrence", export.sum(axis=1))
    export.index.name = "eccDNA"
    export.reset_index().to_csv(
        SOURCE_DATA / "Figure_3D_sorted_recurrent_circles.tsv",
        sep="\t",
        index=False,
    )
    return circles


def plot_figure_3() -> None:
    z_frame, _, covid_samples, hc_samples = prepare_abundance_panel()
    detection = prepare_detection_panel()
    enrichment, _ = gprofiler_enrichment()
    circles = prepare_circle_panel()

    fig = plt.figure(figsize=(183 * MM, 220 * MM), facecolor="white")
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.32, 1.0],
        height_ratios=[1.02, 1.0],
        left=0.205,
        right=0.865,
        bottom=0.045,
        top=0.965,
        hspace=0.34,
        wspace=0.43,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    heat_cmap = LinearSegmentedColormap.from_list(
        "nature_diverging", [HC, "#FFFFFF", COVID]
    )
    image = ax_a.imshow(
        z_frame.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        cmap=heat_cmap,
        vmin=-2.5,
        vmax=2.5,
        rasterized=True,
    )
    ax_a.set_yticks(np.arange(len(z_frame)))
    ax_a.set_yticklabels(z_frame.index, fontsize=5)
    ax_a.set_xticks([])
    ax_a.set_xlabel("Individual plasma samples", labelpad=30)
    ax_a.set_title(
        "Top stable abundance genes", loc="left", pad=5
    )
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    n_covid = len(covid_samples)
    ax_a.add_patch(
        Rectangle(
            (0, -0.045),
            n_covid / len(z_frame.columns),
            0.015,
            color=COVID,
            clip_on=False,
            linewidth=0,
            transform=ax_a.transAxes,
        )
    )
    ax_a.add_patch(
        Rectangle(
            (n_covid / len(z_frame.columns), -0.045),
            len(hc_samples) / len(z_frame.columns),
            0.015,
            color=HC,
            clip_on=False,
            linewidth=0,
            transform=ax_a.transAxes,
        )
    )
    ax_a.text(
        n_covid / (2 * len(z_frame.columns)),
        -0.068,
        "COVID-19 (n=39)",
        ha="center",
        va="top",
        fontsize=5,
        transform=ax_a.transAxes,
    )
    ax_a.text(
        (n_covid + len(hc_samples) / 2) / len(z_frame.columns),
        -0.068,
        "HC (n=39)",
        ha="center",
        va="top",
        fontsize=5,
        transform=ax_a.transAxes,
    )
    cbar = fig.colorbar(image, ax=ax_a, orientation="horizontal", fraction=0.045, pad=0.16)
    cbar.set_label("Row z-score of log1p(EA)", labelpad=1)
    cbar.set_ticks([-2, 0, 2])
    cbar.outline.set_linewidth(0.4)
    panel_label(ax_a, "a", x=-0.13, y=1.04)

    ax_b = fig.add_subplot(gs[0, 1])
    detection_plot = detection.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(detection_plot))
    for i, row in detection_plot.iterrows():
        ax_b.plot(
            [row["hc_detection_frequency"], row["covid_detection_frequency"]],
            [i, i],
            color=GREY_3,
            linewidth=0.7,
            zorder=1,
        )
    ax_b.scatter(
        detection_plot["hc_detection_frequency"],
        y,
        s=17,
        color=HC,
        edgecolor="white",
        linewidth=0.35,
        label="HC",
        zorder=3,
    )
    ax_b.scatter(
        detection_plot["covid_detection_frequency"],
        y,
        s=17,
        color=COVID,
        edgecolor="white",
        linewidth=0.35,
        label="COVID-19",
        zorder=3,
    )
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(detection_plot["Gene"])
    ax_b.set_xlim(-0.02, 1.02)
    ax_b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_b.set_xticklabels(["0", "25", "50", "75", "100"])
    ax_b.set_xlabel("Detection frequency (%)")
    ax_b.set_title(
        "Top assignment-stable eccGene detection differences", loc="left", pad=5
    )
    ax_b.grid(axis="x", color=GREY_1, linewidth=0.5)
    ax_b.spines[["top", "right", "left"]].set_visible(False)
    ax_b.tick_params(axis="y", length=0)
    ax_b.legend(frameon=False, loc="lower right", handletextpad=0.4)
    panel_label(ax_b, "b", x=-0.20, y=1.04)

    ax_c = fig.add_subplot(gs[1, 0])
    enrich_plot = enrichment.sort_values("minus_log10_fdr").reset_index(drop=True)
    y = np.arange(len(enrich_plot))
    source_colors = {"GO:BP": GREEN, "KEGG": ORANGE, "REAC": PURPLE}
    colors = [source_colors[source] for source in enrich_plot["source"]]
    sizes = 12 + 1.6 * np.sqrt(enrich_plot["intersection_size"]) ** 2
    ax_c.scatter(
        enrich_plot["minus_log10_fdr"],
        y,
        s=sizes,
        c=colors,
        edgecolor="white",
        linewidth=0.4,
        zorder=3,
    )
    ax_c.set_yticks(y)
    ax_c.set_yticklabels(
        [textwrap.fill(name, width=39) for name in enrich_plot["name"]], fontsize=5
    )
    ax_c.set_xlabel("−log10(BH-FDR)")
    ax_c.set_title(
        "Functional enrichment of COVID-up abundance genes", loc="left", pad=5
    )
    ax_c.grid(axis="x", color=GREY_1, linewidth=0.5)
    ax_c.spines[["top", "right", "left"]].set_visible(False)
    ax_c.tick_params(axis="y", length=0)
    legend = [
        Patch(facecolor=source_colors[source], edgecolor="none", label=source)
        for source in ("GO:BP", "KEGG", "REAC")
        if source in set(enrich_plot["source"])
    ]
    ax_c.legend(handles=legend, frameon=False, loc="lower right", ncol=3)
    panel_label(ax_c, "c", x=-0.13, y=1.04)

    ax_d = fig.add_subplot(gs[1, 1])
    circle_cmap = ListedColormap([GREY_1, COVID])
    ax_d.imshow(
        circles.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        cmap=circle_cmap,
        vmin=0,
        vmax=1,
        rasterized=True,
    )
    labels = []
    for circle, recurrence in circles.sum(axis=1).items():
        labels.append(f"{circle.replace('_', '–')}  ({int(recurrence)})")
    ax_d.set_yticks(np.arange(len(circles)))
    ax_d.set_yticklabels(labels, fontsize=4.1)
    ax_d.yaxis.tick_right()
    ax_d.tick_params(axis="y", length=0, pad=1.2)
    ax_d.set_xticks([])
    ax_d.set_xlabel("COVID-19 samples, ordered by recurrent-circle count")
    ax_d.set_title(
        "Recurrent exact circles absent from HC", loc="left", pad=5
    )
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    ax_d.legend(
        handles=[
            Patch(facecolor=COVID, edgecolor="none", label="Detected"),
            Patch(facecolor=GREY_1, edgecolor="none", label="Not detected"),
        ],
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0, -0.07),
        ncol=2,
    )
    panel_label(ax_d, "d", x=-0.20, y=1.04)

    save_figure(fig, "Figure_3")
    plt.close(fig)


def plot_figure_s3() -> None:
    summary = read_tsv(RESULTS / "definition_summary.tsv")
    overlap = read_tsv(RESULTS / "definition_overlap_summary.tsv")
    stability = read_tsv(RESULTS / "candidate_stability.tsv")

    # Taller canvas at unchanged width. The margins and the gap between the two
    # left-column panels are pinned to the absolute sizes they had at 92 mm, so the
    # added height is spent on the panel boxes rather than on white space; everything
    # set in points (type, line weights) is untouched.
    fig_h_mm = 115.0
    top_mm, bottom_mm, gap_mm = 9.2, 14.72, 18.0
    rows_mm = fig_h_mm - top_mm - bottom_mm - gap_mm

    fig = plt.figure(figsize=(183 * MM, fig_h_mm * MM), facecolor="white")
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 1.58],
        height_ratios=[0.76, 1.0],
        left=0.075,
        right=0.93,
        bottom=bottom_mm / fig_h_mm,
        top=1.0 - top_mm / fig_h_mm,
        wspace=0.43,
        hspace=gap_mm / (rows_mm / 2.0),
    )

    stable_counts = pd.DataFrame(
        {
            "endpoint": ["Abundance", "Detection"],
            "n_genes_q005_same_direction_all_definitions": [
                int(
                    stability[
                        "abundance_significant_same_direction_all_definitions"
                    ].sum()
                ),
                int(
                    stability[
                        "detection_significant_same_direction_all_definitions"
                    ].sum()
                ),
            ],
            "criterion": [
                "BH q < 0.05 with the same direction under junction, full-interval, and midpoint definitions"
            ]
            * 2,
        }
    )

    ax_a = fig.add_subplot(gs[0, 0])
    y = np.arange(len(stable_counts))
    counts = stable_counts[
        "n_genes_q005_same_direction_all_definitions"
    ].to_numpy()
    # Bar thickness is set in axis units, so it is trimmed to keep the two bars at
    # roughly their previous physical weight on the taller panel.
    bars = ax_a.barh(y, counts, height=0.38, color=[TEAL, PURPLE])
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(stable_counts["endpoint"])
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, counts.max() * 1.23)
    ax_a.set_xlabel("Stable genes (BH $q$ < 0.05, same direction)")
    ax_a.set_title("Stable under all definitions", loc="left")
    for bar, value in zip(bars, counts):
        ax_a.text(
            bar.get_width() + counts.max() * 0.025,
            bar.get_y() + bar.get_height() / 2,
            f"{int(value):,}",
            va="center",
            ha="left",
            fontsize=6,
            fontweight="bold",
        )
    ax_a.spines[["top", "right"]].set_visible(False)
    panel_label(ax_a, "a", x=-0.20, y=1.08)

    ax_b = fig.add_subplot(gs[1, 0])
    labels = ["J–I", "J–M", "I–M"]
    abundance_rows = overlap[
        (overlap["endpoint"] == "abundance_wilcoxon")
        & (overlap["set_filter"] == "all_significant")
    ]
    detection_rows = overlap[
        (overlap["endpoint"] == "detection_fisher")
        & (overlap["set_filter"] == "all_significant")
    ]
    values = np.vstack(
        [abundance_rows["jaccard"].to_numpy(), detection_rows["jaccard"].to_numpy()]
    )
    image = ax_b.pcolormesh(
        np.arange(values.shape[1] + 1) - 0.5,
        np.arange(values.shape[0] + 1) - 0.5,
        values,
        cmap=LinearSegmentedColormap.from_list("blue_scale", [GREY_1, HC]),
        vmin=0,
        vmax=1,
        shading="flat",
    )
    ax_b.set_xlim(-0.5, values.shape[1] - 0.5)
    ax_b.set_ylim(values.shape[0] - 0.5, -0.5)
    ax_b.set_xticks(np.arange(3))
    ax_b.set_xticklabels(labels)
    ax_b.set_yticks([0, 1])
    ax_b.set_yticklabels(["Abundance", "Detection"])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax_b.text(
                j,
                i,
                f"{values[i, j]:.2f}",
                ha="center",
                va="center",
                color="white" if values[i, j] > 0.55 else "black",
                fontsize=6,
            )
    ax_b.add_patch(
        Rectangle(
            (0.5, -0.5),
            1.0,
            2.0,
            fill=False,
            edgecolor=HC,
            linewidth=0.9,
        )
    )
    ax_b.set_title("Significant-set concordance", loc="left")
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax_b, fraction=0.05, pad=0.03)
    cbar.solids.set_rasterized(False)
    cbar.set_label("Jaccard index")
    panel_label(ax_b, "b", x=-0.20, y=1.08)

    ax_c = fig.add_subplot(gs[:, 1])
    stable = stability[
        (stability["abundance_significant_same_direction_all_definitions"] == 1)
        & (stability["junction_abundance_log2FC"] >= 1)
    ].copy()
    stable = stable.sort_values(
        ["junction_abundance_q", "Gene"], ascending=[True, True]
    ).head(15)
    values = stable[
        [
            "junction_abundance_log2FC",
            "interval_abundance_log2FC",
            "midpoint_abundance_log2FC",
        ]
    ].to_numpy()
    if np.nanmin(values) <= 0:
        raise ValueError("Figure S3c selection must be COVID-19-up under every definition")
    cap = 5.0
    image = ax_c.pcolormesh(
        np.arange(values.shape[1] + 1) - 0.5,
        np.arange(values.shape[0] + 1) - 0.5,
        np.clip(values, 0, cap),
        cmap=LinearSegmentedColormap.from_list(
            "covid_up_scale", ["#FFF7F5", "#F2C6C5", COVID]
        ),
        vmin=0,
        vmax=cap,
        shading="flat",
    )
    ax_c.set_xlim(-0.5, values.shape[1] - 0.5)
    ax_c.set_ylim(values.shape[0] - 0.5, -0.5)
    ax_c.set_xticks([0, 1, 2])
    ax_c.set_xticklabels(["Junction", "Interval", "Midpoint"])
    ax_c.set_yticks(np.arange(len(stable)))
    ax_c.set_yticklabels(stable["Gene"])
    for label in ax_c.get_yticklabels():
        label.set_fontstyle("italic")
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            actual = float(values[row_index, column_index])
            display = min(actual, cap)
            ax_c.text(
                column_index,
                row_index,
                "≥5" if actual >= cap else f"{actual:.1f}",
                ha="center",
                va="center",
                color="white" if display >= 3.2 else "black",
                fontsize=4.7,
            )
    ax_c.set_title("Effect sizes for top stable genes", loc="left")
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax_c, fraction=0.05, pad=0.03)
    cbar.solids.set_rasterized(False)
    cbar.set_ticks([0, 1, 2, 3, 4, 5])
    cbar.set_ticklabels(["0", "1", "2", "3", "4", "≥5"])
    cbar.set_label("Mean-EA log2 fold change", fontsize=5.5, labelpad=2)
    panel_label(ax_c, "c", x=-0.15, y=1.04)

    write_tsv(SOURCE_DATA / "Figure_S3_definition_summary.tsv", summary)
    write_tsv(SOURCE_DATA / "Figure_S3_definition_overlap.tsv", overlap)
    write_tsv(SOURCE_DATA / "Figure_S3_stable_counts.tsv", stable_counts)
    stable_export = stable[
        [
            "Gene",
            "junction_abundance_q",
            "junction_abundance_log2FC",
            "interval_abundance_q",
            "interval_abundance_log2FC",
            "midpoint_abundance_q",
            "midpoint_abundance_log2FC",
        ]
    ].copy()
    stable_export.insert(0, "selection_rank", np.arange(1, len(stable_export) + 1))
    stable_export["display_cap_log2FC"] = cap
    write_tsv(
        SOURCE_DATA / "Figure_S3_top_candidate_effects.tsv",
        stable_export,
    )
    save_figure(fig, "Figure_S3")
    plt.close(fig)


def add_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    lines: list[str],
    edge: str,
    fill: str = "white",
) -> None:
    x, y = xy
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.01,rounding_size=0.012",
        linewidth=0.8,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(box)
    ax.text(
        x + 0.02,
        y + height - 0.05,
        title,
        fontsize=7,
        fontweight="bold",
        va="top",
    )
    ax.text(
        x + 0.02,
        y + height - 0.14,
        "\n".join(lines),
        fontsize=5.7,
        va="top",
        linespacing=1.35,
    )


def plot_graphical_abstract() -> None:
    fig, ax = plt.subplots(figsize=(183 * MM, 88 * MM), facecolor="white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_box(
        ax,
        (0.02, 0.18),
        0.19,
        0.66,
        "Plasma cohorts",
        [
            "COVID-19 (n=39)",
            "Healthy controls (n=39)",
            "",
            "Circle-seq profiling",
            "16.86 million calls",
        ],
        edge=GREY_4,
        fill=STONE_1,
    )
    add_box(
        ax,
        (0.27, 0.18),
        0.21,
        0.66,
        "Gene assignment",
        [
            "Main: BED start/junction",
            "within merged gene body",
            "",
            "Sensitivity:",
            "full-interval intersection",
            "midpoint within gene body",
        ],
        edge=TEAL,
    )
    ax.plot([0.31, 0.44], [0.34, 0.34], color=GREY_5, lw=1.1)
    ax.add_patch(
        Rectangle(
            (0.345, 0.315),
            0.07,
            0.05,
            facecolor=GREY_1,
            edgecolor=GREY_4,
            lw=0.5,
        )
    )
    ax.scatter([0.365], [0.34], s=18, color=TEAL, zorder=3)
    ax.text(0.38, 0.285, "junction-associated gene body", ha="center", fontsize=4.8)

    add_box(
        ax,
        (0.54, 0.18),
        0.20,
        0.66,
        "Separate endpoints",
        [
            "Relative abundance (EA)",
            "two-sided Wilcoxon",
            "BH-FDR, log2FC, Cliff's δ",
            "",
            "Detection frequency",
            "Fisher's exact test",
            "OR and exact 95% CI",
        ],
        edge=PURPLE,
    )

    add_box(
        ax,
        (0.80, 0.18),
        0.18,
        0.66,
        "Main observations",
        [
            "Higher Circle-seq burden",
            "in COVID-19",
            "",
            "Altered gene-associated",
            "abundance and detection",
            "",
            "Association with reference",
            "chromatin annotations",
        ],
        edge=ORANGE,
        fill=STONE_1,
    )

    for start, end in ((0.215, 0.265), (0.485, 0.535), (0.745, 0.795)):
        ax.add_patch(
            FancyArrowPatch(
                (start, 0.51),
                (end, 0.51),
                arrowstyle="-|>",
                mutation_scale=8,
                linewidth=0.8,
                color=GREY_5,
            )
        )

    ax.text(
        0.5,
        0.93,
        "Plasma eccDNA profiling identifies COVID-19-associated genomic patterns",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.075,
        "Associations are descriptive and do not establish gene function, pathway activation, tissue of origin, or clinical utility.",
        ha="center",
        va="center",
        fontsize=5.5,
        color=GREY_5,
    )

    save_figure(fig, "Graphical_abstract")
    plt.close(fig)


def write_run_metadata(response: dict) -> None:
    meta = response.get("meta", {})
    metadata = {
        "figure_script": str(Path(__file__).resolve()),
        "main_definition": "Circle-Map BED start/junction within merged UCSC refGene gene body",
        "abundance_display_selection": (
            "Top 30 by junction Wilcoxon BH q among genes significant in the same "
            "direction under junction, interval, and midpoint definitions, with "
            "junction mean-EA log2FC >= 1 and positive Cliff's delta"
        ),
        "detection_display_selection": (
            "Top 15 by junction Fisher BH q among genes significant in the same "
            "direction under all three assignment definitions"
        ),
        "enrichment_query": (
            "junction abundance BH q < 0.05, mean-EA log2FC >= 1, "
            "positive Cliff's delta"
        ),
        "enrichment_background": "all 26,642 genes tested in the junction abundance analysis",
        "enrichment_sources": ["GO:BP", "REAC", "KEGG"],
        "enrichment_correction": "BH-FDR",
        "gprofiler_meta": meta,
    }
    (SOURCE_DATA / "eccGene_figure_run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    _, response = gprofiler_enrichment()
    plot_figure_3()
    plot_figure_s3()
    plot_graphical_abstract()
    write_run_metadata(response)
    print("Created revised Figure 3, Figure S3, and graphical abstract.")


if __name__ == "__main__":
    main()
