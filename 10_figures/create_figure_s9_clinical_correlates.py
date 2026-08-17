#!/usr/bin/env python3
"""Create Supplementary Figure S9 (Submit numbering): association of plasma eccDNA
burden with clinical and laboratory variables in the COVID-19 cohort.

Core conclusion defended by the figure
--------------------------------------
Within the COVID-19 cohort plasma eccDNA burden is not associated with inflammatory,
coagulation or lymphocyte markers, nor with severity, ICU admission, respiratory
support or in-hospital mortality; and the association seen when the two cohorts are
pooled is a between-cohort effect that disappears within each cohort.

Panel d does NOT show a clinically adjusted cohort contrast, and must not be read
as one. The laboratory markers, comorbidity count and days from diagnosis were
recorded only in the COVID-19 participants, so none of them can enter a
between-cohort model. The panel therefore shows the within-COVID-19 effect of each
clinical variable, with the unadjusted cohort contrast plotted above them for
scale and labelled as the reference model.

Every plotted estimate is read from the locked tables under
``deposit/tables/08_clinical_correlates``; nothing is recomputed here. The only
derived quantity is the per-interquartile-range rescaling of the regression
coefficients in panel d, which is the exact algebraic identity
``ratio_per_IQR = 2 ** (coef_log2 * IQR)``; the script asserts this against the
per-unit ratio stored in the locked table. Panel b draws Theil-Sen slopes, which
are rank based and therefore consistent with the Spearman statistics quoted on
the panel; the slopes are drawn for visual guidance only and no new test is run.

Inks follow the locked project palette used by Figure 1, Figure 4 and
Submit/Figure_S2.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "deposit" / "tables" / "08_clinical_correlates"
OUT_STEM = ROOT / "revise" / "Figure_S9_clinical_correlates"
COPIES = (ROOT / "Submit" / "Figure_S9.pdf", ROOT / "deposit" / "figures" / "Figure_S9.pdf")

MM = 1.0 / 25.4
FIG_W_MM, FIG_H_MM = 183.0, 138.0

POINT = {"COVID-19": "#D70000", "HC": "#034E61"}
FILL = {"COVID-19": "#FF8080", "HC": "#8BABD3"}
NEUTRAL_DARK = "#3C4C5A"
NEUTRAL_MID = "#6E7A87"
NEUTRAL_LIGHT = "#AFBAC4"
BAND = "#F1F3F6"
RULE = "#49566D"
NULL_LINE = "#9AA5B1"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6.4,
    "ytick.labelsize": 6.4,
    "legend.fontsize": 6.2,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

LAB_LABEL = {"CRP": "CRP", "IL6": "IL-6", "Ferritin": "Ferritin",
             "D_dimer": "D-dimer", "Lymphocyte_Count": "Lymphocyte count"}
METRICS = ("epm", "eccdna_count")
METRIC_STYLE = {"epm": dict(marker="o", filled=True, label="EPM"),
                "eccdna_count": dict(marker="s", filled=False, label="eccDNA count")}


def panel_tag(ax, letter, dx=-0.155, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="top", ha="left")


def forest(ax, rows, xlabel, null_at, xlim, show_n=True):
    """rows: list of (label, {metric: (est, lo, hi, n, q)})."""
    y = np.arange(len(rows))[::-1]
    for yi, (label, per_metric) in zip(y, rows):
        if yi % 2 == 0:
            ax.axhspan(yi - 0.5, yi + 0.5, color=BAND, lw=0, zorder=0)
        for k, metric in enumerate(METRICS):
            if metric not in per_metric:
                continue
            est, lo, hi, n, q = per_metric[metric]
            off = 0.19 if k == 0 else -0.19
            st = METRIC_STYLE[metric]
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo, hi], [yi + off] * 2, color=NEUTRAL_DARK, lw=0.8,
                        solid_capstyle="butt", zorder=2)
            ax.plot([est], [yi + off], marker=st["marker"], ms=3.4,
                    mfc=NEUTRAL_DARK if st["filled"] else "white",
                    mec=NEUTRAL_DARK, mew=0.8, ls="none", zorder=3)
    ax.axvline(null_at, color=NULL_LINE, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([lb for lb, _ in rows])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    if show_n:
        for yi, (_, per_metric) in zip(y, rows):
            n = next(iter(per_metric.values()))[3]
            ax.text(1.008, yi, f"{int(n)}", transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=5.8, color=NEUTRAL_MID)
        ax.text(1.008, len(rows) - 0.35, "n", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=5.8, color=NEUTRAL_MID, style="italic")


def main() -> None:
    a2 = pd.read_csv(TABLES / "A2_burden_laboratory_correlation_covid.tsv", sep="\t")
    a3 = pd.read_csv(TABLES / "A3_lymphocyte_cross_cohort.tsv", sep="\t")
    a4 = pd.read_csv(TABLES / "A4_burden_by_clinical_strata_covid.tsv", sep="\t")
    a5 = pd.read_csv(TABLES / "A5_regression_clinical_covariates.tsv", sep="\t")
    merged = pd.read_csv(TABLES / "analysis_input_merged.tsv", sep="\t")

    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    gs = fig.add_gridspec(2, 2, left=0.195, right=0.955, bottom=0.085, top=0.935,
                          wspace=0.62, hspace=0.52,
                          width_ratios=[1.0, 0.94], height_ratios=[1.0, 1.0])

    # ---------------------------------------------------------------- a
    ax = fig.add_subplot(gs[0, 0])
    rows = []
    for lab in ("CRP", "IL6", "Ferritin", "D_dimer", "Lymphocyte_Count"):
        per = {}
        for m in METRICS:
            r = a2[(a2.lab == lab) & (a2.burden_metric == m)].iloc[0]
            per[m] = (r.spearman_rho, r.ci_lo, r.ci_hi, r.n, r.q_value_BH_within_metric)
        rows.append((LAB_LABEL[lab], per))
    forest(ax, rows, "Spearman \u03c1 (95% CI)", 0.0, (-0.95, 0.99))
    ax.set_title("Burden vs laboratory markers, within COVID-19", pad=5, loc="left",
                 fontsize=7, color=NEUTRAL_DARK)
    panel_tag(ax, "a")
    # No marker survives correction; quote the smallest q rather than singling
    # out one laboratory marker.
    qmin = float(a2.q_value_BH_within_metric.min())
    ax.text(-0.92, -0.45,
            f"no marker significant after BH correction (smallest q = {qmin:.3f})",
            fontsize=5.6, color=NEUTRAL_MID, ha="left", va="center")


    # ---------------------------------------------------------------- b (hero)
    ax = fig.add_subplot(gs[0, 1])
    sub = merged.dropna(subset=["Lymphocyte_Count", "epm"])
    for g in ("HC", "COVID-19"):
        s = sub[sub.group == g]
        ax.plot(s.Lymphocyte_Count, s.epm, ls="none", marker="o", ms=2.9,
                mfc=FILL[g], mec=POINT[g], mew=0.5, alpha=0.95, zorder=3, label=g)

    def theil(s):
        x = s.Lymphocyte_Count.to_numpy(float)
        y = np.log2(s.epm.to_numpy(float))
        sl, ic, _, _ = stats.theilslopes(y, x)
        xs = np.linspace(x.min(), x.max(), 50)
        return xs, 2 ** (ic + sl * xs)

    xs, ys = theil(sub)
    ax.plot(xs, ys, color=NEUTRAL_DARK, lw=1.3, zorder=4)
    for g in ("COVID-19", "HC"):
        xs, ys = theil(sub[sub.group == g])
        ax.plot(xs, ys, color=POINT[g], lw=1.1, ls=(0, (3.5, 1.8)), zorder=4)

    rho_all = a3[(a3.stratum == "all 78") & (a3.burden_metric == "epm")].iloc[0]
    rho_cov = a3[(a3.stratum == "COVID-19") & (a3.burden_metric == "epm")].iloc[0]
    rho_hc = a3[(a3.stratum == "HC") & (a3.burden_metric == "epm")].iloc[0]
    ax.set_yscale("log")
    ax.set_xlabel("Lymphocyte count ($\\mathregular{10^{9}}$ l$\\mathregular{^{-1}}$)")
    ax.set_ylabel("eccDNA burden (EPM)")
    ax.set_title("Pooled association is a between-cohort effect", pad=5, loc="left",
                 fontsize=7, color=NEUTRAL_DARK)
    panel_tag(ax, "b", dx=-0.20)
    ax.text(0.975, 0.965,
            f"pooled  \u03c1 = {rho_all.spearman_rho:.2f},  q = {rho_all.q_value_BH:.1e}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.9, color=NEUTRAL_DARK)
    ax.text(0.975, 0.885,
            f"COVID-19  \u03c1 = {rho_cov.spearman_rho:.2f},  "
            f"q = {rho_cov.q_value_BH:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.9, color=POINT["COVID-19"])
    ax.text(0.975, 0.808,
            f"HC  \u03c1 = {rho_hc.spearman_rho:.2f},  q = {rho_hc.q_value_BH:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.9, color=POINT["HC"])
    ax.legend(loc="lower left", handletextpad=0.35, borderpad=0.1, labelspacing=0.25,
              bbox_to_anchor=(-0.02, -0.02))

    # ---------------------------------------------------------------- c
    ax = fig.add_subplot(gs[1, 0])
    order = [("ICU vs ward (at sampling)", "ICU vs ward"),
             ("Invasive ventilation/ECMO vs lower support", "IMV or ECMO vs lower"),
             ("Critical vs Severe/Moderate", "Critical vs rest"),
             ("Critical/Severe vs Moderate", "Critical/severe vs moderate"),
             ("Died vs survived (in-hospital)", "Died vs survived")]
    rows = []
    for key, label in order:
        per = {}
        for m in METRICS:
            r = a4[(a4.comparison == key) & (a4.burden_metric == m)].iloc[0]
            per[m] = (r.cliffs_delta, np.nan, np.nan,
                      r.n_group1 + r.n_group2, r.q_value_BH)
        rows.append((label, per))
    forest(ax, rows, "Cliff's \u03b4", 0.0, (-1.02, 1.02))
    ax.set_title("Burden by clinical stratum, within COVID-19", pad=5, loc="left",
                 fontsize=7, color=NEUTRAL_DARK)
    panel_tag(ax, "c")
    qvals = a4.q_value_BH.to_numpy(float)
    if np.allclose(qvals, qvals[0]):
        q_note = f"all q = {qvals[0]:.2f}"
    else:
        q_note = f"q = {qvals.min():.2f}\u2013{qvals.max():.2f}"
    ax.set_xlabel(f"Cliff's \u03b4   ({q_note})")


    # ---------------------------------------------------------------- d
    ax = fig.add_subplot(gs[1, 1])
    ref = a5[(a5.model.str.startswith("Full cohort")) & (a5.term == "group_covid")].iloc[0]
    entries = [("COVID-19 vs HC\n(reference model)", ref.ratio, ref.ci_lo, ref.ci_hi,
                ref.n, POINT["COVID-19"])]
    cov_terms = [("Comorbidity count", "n_comorbidity", "Within COVID-19 + comorbidity count"),
                 ("Days from diagnosis", "Days_Diagnosis_to_Sampling",
                  "Within COVID-19 + days from diagnosis"),
                 ("CRP", "CRP", "Within COVID-19 + CRP"),
                 ("IL-6", "IL6", "Within COVID-19 + IL6"),
                 ("Ferritin", "Ferritin", "Within COVID-19 + Ferritin"),
                 ("D-dimer", "D_dimer", "Within COVID-19 + D_dimer"),
                 ("Lymphocyte count", "Lymphocyte_Count", "Within COVID-19 + Lymphocyte_Count")]
    covid = merged[merged.group == "COVID-19"]
    for label, term, model in cov_terms:
        r = a5[(a5.model == model) & (a5.term == term)].iloc[0]
        v = covid[term].dropna()
        iqr = float(v.quantile(0.75) - v.quantile(0.25))
        # exact algebraic rescaling of the locked per-unit coefficient
        assert abs(2 ** r.coef_log2 - r.ratio) < 1e-9, f"per-unit ratio mismatch for {term}"
        entries.append((label, 2 ** (r.coef_log2 * iqr),
                        2 ** (np.log2(r.ci_lo) * iqr), 2 ** (np.log2(r.ci_hi) * iqr),
                        r.n, NEUTRAL_DARK))
    y = np.arange(len(entries))[::-1]
    for yi, (label, est, lo, hi, n, ink) in zip(y, entries):
        if yi % 2 == 0:
            ax.axhspan(yi - 0.5, yi + 0.5, color=BAND, lw=0, zorder=0)
        ax.plot([lo, hi], [yi, yi], color=ink, lw=0.9, solid_capstyle="butt", zorder=2)
        ax.plot([est], [yi], marker="o", ms=3.6, mfc=ink, mec=ink, ls="none", zorder=3)
        ax.text(1.008, yi, f"{int(n)}", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=5.8, color=NEUTRAL_MID)
    ax.axvline(1.0, color=NULL_LINE, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.28, 14)
    ax.set_xticks([0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.5", "1", "2", "5", "10"])
    ax.set_yticks(y)
    ax.set_yticklabels([e[0] for e in entries])
    ax.set_ylim(-0.7, len(entries) - 0.3)
    ax.set_xlabel("Fold change in EPM (95% CI)\nclinical terms within COVID-19, per IQR")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.text(1.008, len(entries) - 0.35, "n", transform=ax.get_yaxis_transform(),
            va="center", ha="left", fontsize=5.8, color=NEUTRAL_MID, style="italic")
    ax.set_title("No clinical covariate predicts burden", pad=5,
                 loc="left", fontsize=7, color=NEUTRAL_DARK)
    panel_tag(ax, "d", dx=-0.29)

    handles = [Line2D([], [], marker=METRIC_STYLE[m]["marker"], ls="none", ms=3.4,
                      mfc=NEUTRAL_DARK if METRIC_STYLE[m]["filled"] else "white",
                      mec=NEUTRAL_DARK, mew=0.8, label=METRIC_STYLE[m]["label"])
               for m in METRICS]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.195, 0.995),
               ncol=2, handletextpad=0.35, columnspacing=1.1)

    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{OUT_STEM}.pdf")
    fig.savefig(f"{OUT_STEM}.svg")
    plt.close(fig)
    for dest in COPIES:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(f"{OUT_STEM}.pdf", dest)
    print("written:", OUT_STEM.with_suffix(".pdf"))
    for dest in COPIES:
        print("  copied ->", dest)


if __name__ == "__main__":
    main()
