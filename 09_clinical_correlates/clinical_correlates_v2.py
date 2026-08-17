#!/usr/bin/env python3
"""Clinical-correlate analysis on the v2 burden (Figure S5, Tables S9/A1-A6).

REIMPLEMENTED under pre-declaration amendments A3/A4. The original code exists
neither in the repository nor on the cluster; the definitions below are taken
from the manuscript Methods and from the column structure of the v1 artefacts.

Confidence intervals use 5,000 percentile bootstrap resamples under a fixed
seed, not a Fisher z transform. That was established empirically: the v1 CIs do
not reproduce under Fisher z (for the pooled lymphocyte correlation, rho =
-0.5748 with n = 78, Fisher z gives [-0.706, -0.403] against the recorded
[-0.685, -0.427]), and the Methods specify percentile bootstrap for the
parallel RCA correlations. Amendment A3 item 5 is corrected accordingly.

Only the eccDNA burden columns differ between v1 and v2; every clinical
variable is carried through unchanged, and `mapped_reads` is asserted identical
because the BAMs are shared.

This version is the single authoritative generator for A2--A6 and the clinical
block of Supplementary Table S9.  Earlier versions generated only A2--A4,
leaving A5/A6 and Table S9 on a different call-set version.  All bootstrap
intervals now use the same recorded 5,000-resample convention.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from openpyxl import load_workbook
from scipy import stats

RUN = Path(__file__).resolve().parents[2]
MERGED = RUN / "v2tree" / "clinical" / "analysis_input_merged.tsv"
OUT = RUN / "tables" / "clinical"
MIRROR = RUN / "v2tree" / "clinical_tables"
WORKBOOK = RUN / "tables" / "supplementary" / "Supplementary_Tables_v2.xlsx"

LABS = ["CRP", "IL6", "Ferritin", "D_dimer", "Lymphocyte_Count"]
METRICS = ["epm", "eccdna_count"]
BOOTSTRAP = 5000
SEED = 20260809


def load() -> list[dict[str, str]]:
    with MERGED.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(rows, key):
    out = []
    for r in rows:
        v = r.get(key, "")
        out.append(float(v) if v not in ("", "NA", "None") else math.nan)
    return np.asarray(out, dtype=float)


def bh(p: list[float]) -> list[float]:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        prev = min(prev, p[i] * n / (n - rank + 1))
        q[i] = prev
    return list(q)


def spearman_ci(x, y, rng):
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4:
        return n, math.nan, math.nan, math.nan, math.nan
    rho, p = stats.spearmanr(x, y)
    boots = np.empty(BOOTSTRAP)
    for b in range(BOOTSTRAP):
        idx = rng.integers(0, n, n)
        if len(np.unique(x[idx])) < 2 or len(np.unique(y[idx])) < 2:
            boots[b] = np.nan
            continue
        boots[b] = stats.spearmanr(x[idx], y[idx]).statistic
    finite = boots[np.isfinite(boots)]
    lo, hi = np.percentile(finite, [2.5, 97.5])
    return n, float(rho), float(lo), float(hi), float(p)


def cliffs_delta(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    gt = sum((a[:, None] > b[None, :]).sum(axis=1))
    lt = sum((a[:, None] < b[None, :]).sum(axis=1))
    return float((gt - lt) / (len(a) * len(b)))


def write(name: str, rows: list[dict], fields: list[str]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    MIRROR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT / name, MIRROR / name)
    print(f"  wrote {name}: {len(rows)} rows")


def stable_rng(*parts: str) -> np.random.Generator:
    """Return a deterministic per-analysis RNG independent of loop ordering."""
    token = "|".join((str(SEED),) + tuple(parts)).encode()
    seed = int.from_bytes(hashlib.sha256(token).digest()[:8], "little")
    return np.random.default_rng(seed)


def fit_hc3(rows: list[dict[str, str]], predictors: list[tuple[str, str, float]],
            model: str) -> list[dict]:
    """Fit log2(EPM+1) OLS with HC3 SEs and t-based inference.

    ``predictors`` contains (output term, input column, scale divisor).  The
    divisor makes the full-cohort age term explicitly per 10 years while the
    within-COVID models retain their established per-year coefficient.
    """
    keep = ["epm"] + [source for _, source, _ in predictors]
    frame = []
    for row in rows:
        item = {}
        try:
            for key in keep:
                value = row.get(key, "")
                if value in ("", "NA", "None"):
                    raise ValueError
                item[key] = float(value)
        except (TypeError, ValueError):
            continue
        frame.append(item)
    if not frame:
        raise RuntimeError(f"no complete cases for model: {model}")

    y = np.log2(np.asarray([r["epm"] for r in frame], float) + 1.0)
    names = [term for term, _, _ in predictors]
    x = np.column_stack([
        np.asarray([r[source] for r in frame], float) / divisor
        for _, source, divisor in predictors
    ])
    x = sm.add_constant(x, has_constant="add")
    fit = sm.OLS(y, x).fit(cov_type="HC3")
    term_names = ["intercept"] + names
    tq = float(stats.t.ppf(0.975, fit.df_resid))
    out = []
    for i, term in enumerate(term_names):
        coef = float(fit.params[i])
        se = float(fit.bse[i])
        p = float(2.0 * stats.t.sf(abs(coef / se), fit.df_resid))
        out.append({
            "model": model,
            "n": int(fit.nobs),
            "term": term,
            "coef_log2": coef,
            "ratio": float(2.0 ** coef),
            "ci_lo": float(2.0 ** (coef - tq * se)),
            "ci_hi": float(2.0 ** (coef + tq * se)),
            "p_value": p,
        })
    return out


def rebuild_table_s9(a2, a3, a4, a5, a6) -> None:
    """Rewrite every clinical result in Table S9 from the v2 locked tables."""
    wb = load_workbook(WORKBOOK)
    ws = wb["Table_S9"]
    for r in range(3, ws.max_row + 1):
        for c in range(1, ws.max_column + 1):
            ws.cell(r, c).value = None

    rows = []
    for item in a2:
        rows.append([
            "A. Burden–laboratory correlation (COVID-19 only)",
            f"{item['burden_metric']} vs {item['lab']}", item["n"],
            "Spearman rho", round(item["spearman_rho"], 3),
            round(item["ci_lo"], 3), round(item["ci_hi"], 3),
            item["p_value"], item["q_value_BH_within_metric"],
        ])
    for item in a3:
        label = f"{item['stratum']} — {item['burden_metric']}"
        if item["burden_metric"] == "-":
            rows.append([
                "B. Lymphocyte count, pooled vs within-cohort", label,
                item["n"], "Mann–Whitney U", None, None, None,
                item["p_value"], None,
            ])
        else:
            rows.append([
                "B. Lymphocyte count, pooled vs within-cohort", label,
                item["n"], "Spearman rho", round(item["spearman_rho"], 3),
                round(item["ci_lo"], 3), round(item["ci_hi"], 3),
                item["p_value"], item["q_value_BH"],
            ])
    for item in a4:
        rows.append([
            "C. Burden by clinical stratum (COVID-19 only)",
            f"{item['comparison']} — {item['burden_metric']} "
            f"(n={item['n_group1']} vs {item['n_group2']})",
            item["n_group1"] + item["n_group2"], "Cliff's delta",
            round(item["cliffs_delta"], 3), None, None,
            item["p_value"], item["q_value_BH"],
        ])

    a5_index = {(r["model"], r["term"]): r for r in a5}
    d_specs = [
        ("Full cohort: group + RCA + age + sex (reference model, no clinical "
         "covariate; identical to robustness_regression_HC3.tsv, callset "
         "circlemap_methods)", "group_covid",
         "Full cohort: group + RCA + age + sex (reference model — no clinical covariate)"),
        ("Within COVID-19 + comorbidity count", "n_comorbidity",
         "Within COVID-19 + comorbidity count — term: n_comorbidity"),
        ("Within COVID-19 + days from diagnosis", "Days_Diagnosis_to_Sampling",
         "Within COVID-19 + days from diagnosis — term: Days_Diagnosis_to_Sampling"),
        ("Within COVID-19 + comorbidity + days", "n_comorbidity",
         "Within COVID-19 + comorbidity + days — term: n_comorbidity"),
        ("Within COVID-19 + comorbidity + days", "Days_Diagnosis_to_Sampling",
         "Within COVID-19 + comorbidity + days — term: Days_Diagnosis_to_Sampling"),
        ("Within COVID-19 + CRP", "CRP", "Within COVID-19 + CRP — term: CRP"),
        ("Within COVID-19 + IL6", "IL6", "Within COVID-19 + IL6 — term: IL6"),
        ("Within COVID-19 + Ferritin", "Ferritin",
         "Within COVID-19 + Ferritin — term: Ferritin"),
        ("Within COVID-19 + D_dimer", "D_dimer",
         "Within COVID-19 + D_dimer — term: D_dimer"),
        ("Within COVID-19 + Lymphocyte_Count", "Lymphocyte_Count",
         "Within COVID-19 + Lymphocyte_Count — term: Lymphocyte_Count"),
    ]
    for model, term, label in d_specs:
        item = a5_index[(model, term)]
        rows.append([
            "D. HC3 robust regression of log2(EPM+1)", label, item["n"],
            "Fold change per unit", round(item["ratio"], 3),
            round(item["ci_lo"], 3), round(item["ci_hi"], 3),
            item["p_value"], None,
        ])

    for item in a6:
        rows.append([
            "E. Sensitivity analyses (EPM)",
            f"{item['subset']} — {item['lab']}", item["n"],
            "Spearman rho / Cliff's delta", round(item["spearman_rho"], 3),
            None if item["ci_lo"] is None else round(item["ci_lo"], 3),
            None if item["ci_hi"] is None else round(item["ci_hi"], 3),
            item["p_value"], None,
        ])

    if len(rows) != 61:
        raise RuntimeError(f"Table S9 expected 61 data rows, found {len(rows)}")
    for r, values in enumerate(rows, start=3):
        for c, value in enumerate(values, start=1):
            ws.cell(r, c).value = value
    ws.cell(65, 1).value = (
        "Q values are Benjamini–Hochberg values computed within the family "
        "stated for each panel (panel A, within each burden metric; panels B "
        "and C, within the panel). Spearman confidence intervals in panels A, "
        "B and E are percentile bootstrap intervals from 5,000 resamples with "
        f"base seed {SEED}. Panel D reports the coefficient of the named term "
        "using HC3 robust standard errors and t-based confidence intervals; "
        "the laboratory markers, comorbidity count and days from diagnosis "
        "were recorded only in the COVID-19 participants, so models containing "
        "them are fitted within that cohort."
    )
    wb.save(WORKBOOK)
    print("  rebuilt Table_S9 from A2–A6")


def main() -> int:
    rows = load()
    covid = [r for r in rows if r["group"] == "COVID-19"]
    hc = [r for r in rows if r["group"] == "HC"]
    print(f"loaded {len(rows)} samples ({len(covid)} COVID-19, {len(hc)} HC)")
    rng = np.random.default_rng(SEED)

    # ---- A2: burden vs laboratory values, COVID-19 only -------------------
    a2, pv = [], {m: [] for m in METRICS}
    for metric in METRICS:
        for lab in LABS:
            n, rho, lo, hi, p = spearman_ci(num(covid, metric), num(covid, lab), rng)
            a2.append({"cohort": "COVID-19", "burden_metric": metric, "lab": lab,
                       "n": n, "spearman_rho": rho, "ci_lo": lo, "ci_hi": hi,
                       "p_value": p, "q_value_BH_within_metric": None,
                       "ci_method": f"percentile bootstrap, {BOOTSTRAP} resamples"})
            pv[metric].append(p)
    for metric in METRICS:
        qs = bh(pv[metric])
        for row, q in zip([r for r in a2 if r["burden_metric"] == metric], qs):
            row["q_value_BH_within_metric"] = q
    write("A2_burden_laboratory_correlation_covid.tsv", a2, list(a2[0]))

    # ---- A3: lymphocyte count across strata -------------------------------
    a3, praw = [], []
    for label, subset in (("all 78", rows), ("COVID-19", covid), ("HC", hc)):
        for metric in METRICS:
            n, rho, lo, hi, p = spearman_ci(num(subset, metric),
                                            num(subset, "Lymphocyte_Count"), rng)
            a3.append({"stratum": label, "burden_metric": metric, "n": n,
                       "spearman_rho": rho, "ci_lo": lo, "ci_hi": hi,
                       "p_value": p, "q_value_BH": None,
                       "ci_method": f"percentile bootstrap, {BOOTSTRAP} resamples"})
            praw.append(p)
    for row, q in zip(a3, bh(praw)):
        row["q_value_BH"] = q
    lc = num(covid, "Lymphocyte_Count"), num(hc, "Lymphocyte_Count")
    ok = [v[~np.isnan(v)] for v in lc]
    u = stats.mannwhitneyu(ok[0], ok[1], alternative="two-sided")
    a3.append({"stratum": "COVID-19 vs HC lymphocyte count", "burden_metric": "-",
               "n": len(ok[0]) + len(ok[1]), "spearman_rho": None,
               "ci_lo": None, "ci_hi": None, "p_value": float(u.pvalue),
               "q_value_BH": None,
               "ci_method": f"two-sided Mann-Whitney U; Cliff's delta "
                            f"{cliffs_delta(ok[0], ok[1]):.6f}"})
    write("A3_lymphocyte_cross_cohort.tsv", a3, list(a3[0]))

    # ---- A4: burden by clinical strata, COVID-19 only ---------------------
    # Field values are taken from the merged table as recorded, not guessed:
    #   Ward_ICU              Yes / No           (ICU at sampling)
    #   Mortality_in_hospital Yes / No
    #   Severity              Critical / Severe / Moderate / Not specified
    #   Respiratory_Support   free text; the invasive arm is invasive
    #                         mechanical ventilation, and "No record" is
    #                         excluded from the comparison arm while
    #                         "None documented" / "None (room air)" are lower
    #                         support and stay in it.
    NO_SUPPORT_RECORD = {"No record", "", "NA"}
    INVASIVE = {"Invasive mechanical ventilation", "ECMO"}
    strata = [
        ("ICU vs ward (at sampling)",
         lambda r: r["Ward_ICU"] == "Yes", lambda r: r["Ward_ICU"] == "No"),
        ("Died vs survived (in-hospital)",
         lambda r: r["Mortality_in_hospital"] == "Yes",
         lambda r: r["Mortality_in_hospital"] == "No"),
        ("Critical vs Severe/Moderate",
         lambda r: r["Severity"] == "Critical",
         lambda r: r["Severity"] in ("Severe", "Moderate")),
        ("Critical/Severe vs Moderate",
         lambda r: r["Severity"] in ("Critical", "Severe"),
         lambda r: r["Severity"] == "Moderate"),
        ("Invasive ventilation/ECMO vs lower support",
         lambda r: r["Respiratory_Support"] in INVASIVE,
         lambda r: r["Respiratory_Support"] not in INVASIVE | NO_SUPPORT_RECORD),
    ]
    a4, praw = [], []
    for label, f1, f2 in strata:
        g1, g2 = [r for r in covid if f1(r)], [r for r in covid if f2(r)]
        for metric in METRICS:
            x, y = num(g1, metric), num(g2, metric)
            x, y = x[~np.isnan(x)], y[~np.isnan(y)]
            if len(x) < 2 or len(y) < 2:
                continue
            u = stats.mannwhitneyu(x, y, alternative="two-sided")
            a4.append({"comparison": label, "burden_metric": metric,
                       "n_group1": len(x), "n_group2": len(y),
                       "median1": float(np.median(x)), "median2": float(np.median(y)),
                       "cliffs_delta": cliffs_delta(x, y),
                       "p_value": float(u.pvalue), "q_value_BH": None})
            praw.append(float(u.pvalue))
    for row, q in zip(a4, bh(praw)):
        row["q_value_BH"] = q
    write("A4_burden_by_clinical_strata_covid.tsv", a4, list(a4[0]))

    # ---- A5: HC3 clinical-covariate models --------------------------------
    reference_model = (
        "Full cohort: group + RCA + age + sex (reference model, no clinical "
        "covariate; identical to robustness_regression_HC3.tsv, callset "
        "circlemap_methods)"
    )
    # The merged clinical input stores EPM rounded to four decimals.  Use the
    # authoritative full-precision RCA table for the identical reference model
    # so it remains exactly the 4.8382-fold fit reported in the main Results.
    with (RUN / "tables" / "RCA" / "RCA_regression_HC3.tsv").open(
            newline="") as fh:
        rca = list(csv.DictReader(fh, delimiter="\t"))
    rca = [r for r in rca if r["outcome"] == "epm" and
           r["model"] == "full_covariate_adjusted"]
    term_map = {"intercept": "intercept", "COVID-19_vs_HC": "group_covid",
                "log2_RCA_per_doubling": "log2_rca",
                "age_per_10_years": "age_per_10_years",
                "male_vs_female": "male"}
    a5 = [{"model": reference_model, "n": int(r["n"]),
           "term": term_map[r["term"]],
           "coef_log2": float(r["estimate_log2_scale"]),
           "ratio": float(r["multiplicative_effect"]),
           "ci_lo": float(r["multiplicative_ci_low"]),
           "ci_hi": float(r["multiplicative_ci_high"]),
           "p_value": float(r["p_value_two_sided"])} for r in rca]
    model_specs = [
        ("Within COVID-19 + CRP", [("CRP", "CRP", 1.0)]),
        ("Within COVID-19 + IL6", [("IL6", "IL6", 1.0)]),
        ("Within COVID-19 + Ferritin", [("Ferritin", "Ferritin", 1.0)]),
        ("Within COVID-19 + D_dimer", [("D_dimer", "D_dimer", 1.0)]),
        ("Within COVID-19 + Lymphocyte_Count",
         [("Lymphocyte_Count", "Lymphocyte_Count", 1.0)]),
        ("Within COVID-19 + comorbidity count",
         [("n_comorbidity", "n_comorbidity", 1.0)]),
        ("Within COVID-19 + days from diagnosis",
         [("Days_Diagnosis_to_Sampling", "Days_Diagnosis_to_Sampling", 1.0)]),
        ("Within COVID-19 + comorbidity + days",
         [("n_comorbidity", "n_comorbidity", 1.0),
          ("Days_Diagnosis_to_Sampling", "Days_Diagnosis_to_Sampling", 1.0)]),
    ]
    within_base = [("log2_rca", "log2_rca", 1.0),
                   ("age", "age", 1.0), ("male", "male", 1.0)]
    for model, extra in model_specs:
        a5.extend(fit_hc3(covid, within_base + extra, model))
    write("A5_regression_clinical_covariates.tsv", a5, list(a5[0]))

    ref = next(r for r in a5 if r["model"] == reference_model
               and r["term"] == "group_covid")
    if not (abs(ref["ratio"] - 4.838182184) < 1e-9
            and abs(ref["p_value"] - 2.299234665e-06) < 1e-15):
        raise RuntimeError("A5 reference HC3 model did not reproduce the locked fit")

    # ---- A6: prespecified clinical-record sensitivity subsets -------------
    subsets = [
        ("all COVID (n=39)", covid),
        ("excluding manually reconstructed records (XS120, NS46)",
         [r for r in covid if r["sample_id"] not in {"XS120", "NS46"}]),
        ("documented severity label only",
         [r for r in covid if r.get("Severity_documented", "") not in
          ("", "NA", "None")]),
        ("excluding low-confidence index date (XS57, NS167, NS271)",
         [r for r in covid if r["sample_id"] not in {"XS57", "NS167", "NS271"}]),
    ]
    a6 = []
    for subset_name, subset in subsets:
        for lab in LABS:
            n, rho, lo, hi, p = spearman_ci(
                num(subset, "epm"), num(subset, lab),
                stable_rng("A6", subset_name, lab),
            )
            a6.append({"subset": subset_name, "n_total": len(subset),
                       "lab": lab, "n": n, "spearman_rho": rho,
                       "ci_lo": lo, "ci_hi": hi, "p_value": p})
        icu = num([r for r in subset if r["Ward_ICU"] == "Yes"], "epm")
        ward = num([r for r in subset if r["Ward_ICU"] == "No"], "epm")
        u = stats.mannwhitneyu(icu, ward, alternative="two-sided")
        a6.append({"subset": subset_name, "n_total": len(subset),
                   "lab": "ICU vs ward (EPM)", "n": len(icu) + len(ward),
                   "spearman_rho": cliffs_delta(icu, ward),
                   "ci_lo": None, "ci_hi": None, "p_value": float(u.pvalue)})
    write("A6_sensitivity_analyses.tsv", a6, list(a6[0]))

    invalid_ci = [r for r in a6 if r["ci_lo"] is not None and
                  not (r["ci_lo"] <= r["spearman_rho"] <= r["ci_hi"])]
    if invalid_ci:
        raise RuntimeError(f"A6 contains {len(invalid_ci)} CIs excluding rho")

    rebuild_table_s9(a2, a3, a4, a5, a6)

    (OUT / "clinical_run_parameters.json").write_text(json.dumps({
        "reimplemented": True,
        "amendments": ["A3", "A4"],
        "input": str(MERGED.relative_to(RUN)),
        "bootstrap_resamples": BOOTSTRAP,
        "seed": SEED,
        "ci_method": "percentile bootstrap",
        "generated_tables": ["A2", "A3", "A4", "A5", "A6", "Table_S9"],
        "hc3_inference": "HC3 robust standard errors with t-based intervals and p values",
        "note": ("No v1 code exists to gate against. Clinical variables are "
                 "unchanged; only the eccDNA burden columns come from the v2 "
                 "methods call set."),
    }, indent=2) + "\n")
    shutil.copy2(OUT / "clinical_run_parameters.json",
                 MIRROR / "clinical_run_parameters.json")
    print(f"\nwritten to {OUT.relative_to(RUN)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
