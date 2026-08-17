#!/usr/bin/env python3
"""Re-run module 02 and Supplementary Figure S3 without the BAMs.

deposit/scripts/02_rca_technical_bias/run_rca_bias_analysis.py needs the BAM
indices only inside collect_sample_metrics(), which derives the per-sample read
and call counts. Those counts do not depend on the RCA concentration, so they
are taken unchanged from the already-deposited RCA_sample_metrics.tsv and the
original module's own group_statistics / correlation_statistics /
regression_statistics / figure functions are driven directly. The statistics and
the figure are therefore produced by the published code, not a reimplementation.

The run is validated first: with the UNROUNDED concentrations every regenerated
table must reproduce the currently deposited one, otherwise the script stops
before writing anything.

    python3 revise/scripts/fixes/rerun_module02_locally.py --check   # validate only
    python3 revise/scripts/fixes/rerun_module02_locally.py           # validate, then write

Run from the repository root.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MOD = ROOT / "deposit" / "scripts" / "02_rca_technical_bias" / "run_rca_bias_analysis.py"
T02 = ROOT / "deposit" / "tables" / "02_rca_technical_bias"
METRICS = T02 / "RCA_sample_metrics.tsv"

INT_COLS = ("eccdna_count", "discordant_read_pairs", "split_reads",
            "total_support_events", "header_lines_skipped", "mapped_reads")
FLOAT_COLS = ("rca_concentration_ng_ul", "epm", "support_events_per_million_mapped_reads",
              "age_years")

FIELDS = {
    "RCA_group_summary.tsv": ("metric", "metric_label", "group", "n", "mean", "sd",
                              "median", "q1", "q3", "minimum", "maximum"),
    "RCA_group_comparisons.tsv": ("metric", "metric_label", "comparison", "n_covid", "n_hc",
                                  "mann_whitney_u", "p_value_two_sided",
                                  "q_value_bh_across_group_tests",
                                  "cliffs_delta_covid_minus_hc", "median_covid", "median_hc"),
    "RCA_correlations.tsv": ("scope", "n", "x_metric", "y_metric", "y_metric_label",
                             "spearman_rho", "bootstrap_ci_low", "bootstrap_ci_high",
                             "p_value_two_sided", "q_value_bh_within_scope",
                             "q_value_bh_across_all_correlations", "bootstrap_iterations"),
    "RCA_regression_HC3.tsv": ("outcome", "outcome_label", "outcome_transform", "model", "term",
                               "estimate_log2_scale", "hc3_standard_error", "ci_low_log2_scale",
                               "ci_high_log2_scale", "multiplicative_effect",
                               "multiplicative_ci_low", "multiplicative_ci_high", "t_value",
                               "p_value_two_sided", "q_value_bh_within_term", "n",
                               "degrees_of_freedom", "r_squared"),
}


def load_module():
    spec = importlib.util.spec_from_file_location("rca_mod", MOD)
    m = importlib.util.module_from_spec(spec)
    sys.modules["rca_mod"] = m
    spec.loader.exec_module(m)
    return m


def sample_rows_from(df: pd.DataFrame) -> list[dict]:
    rows = []
    for _, r in df.iterrows():
        d = {k: r[k] for k in df.columns}
        for c in INT_COLS:
            if c in d and pd.notna(d[c]):
                d[c] = int(d[c])
        for c in FLOAT_COLS:
            if c in d and pd.notna(d[c]):
                d[c] = float(d[c])
        rows.append(d)
    return rows


def compute(m, rows, boot=5000):
    gs, cmp_ = m.group_statistics(rows)
    corr = m.correlation_statistics(rows, boot)
    reg = m.regression_statistics(rows)
    return {"RCA_group_summary.tsv": gs, "RCA_group_comparisons.tsv": cmp_,
            "RCA_correlations.tsv": corr, "RCA_regression_HC3.tsv": reg}


def to_frame(rows, fields):
    return pd.DataFrame([{k: r.get(k) for k in fields} for r in rows])


def compare(new_rows, path: Path, fields, tol=1e-9) -> list[str]:
    if not path.exists():
        return [f"{path.name}: reference missing"]
    a = to_frame(new_rows, fields)
    b = pd.read_csv(path, sep="\t")[list(fields)]
    if a.shape != b.shape:
        return [f"{path.name}: shape {a.shape} vs {b.shape}"]
    bad = []
    for c in fields:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(y) and pd.api.types.is_numeric_dtype(pd.to_numeric(x, errors="coerce")):
            xv = pd.to_numeric(x, errors="coerce").to_numpy(float)
            yv = y.to_numpy(float)
            d = np.abs(xv - yv)
            rel = d / np.maximum(np.abs(yv), 1e-12)
            m_ = np.isfinite(xv) & np.isfinite(yv)
            if ((d > tol) & (rel > 1e-6) & m_).any():
                bad.append(f"{path.name}:{c} max rel diff {np.nanmax(rel[m_]):.3g}")
        else:
            if (x.astype(str).values != y.astype(str).values).any():
                bad.append(f"{path.name}:{c} text differs")
    return bad


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="validate only, write nothing")
    args = ap.parse_args()

    m = load_module()
    df = pd.read_csv(METRICS, sep="\t")

    # ---- validate: unrounded concentrations must reproduce the deposited tables
    pre = T02.parent.parent.parent / "revise" / "backups" / "rca_rounding" / \
        "deposit_tables_02_rca_technical_bias_RCA_sample_metrics.tsv.preround"
    unrounded = pd.read_csv(pre, sep="\t") if pre.exists() else df
    ref_dir = ROOT / "revise" / "backups" / "rca_rounding"
    got = compute(m, sample_rows_from(unrounded))
    problems = []
    for name, fields in FIELDS.items():
        ref = ref_dir / f"deposit_tables_02_rca_technical_bias_{name}.preround"
        problems += compare(got[name], ref, fields)
    if problems:
        print("VALIDATION FAILED - nothing written:")
        for p in problems:
            print("   ", p)
        sys.exit(1)
    print("  validation: the original module code reproduces all four deposited "
          "tables exactly from the unrounded concentrations")

    if args.check:
        return

    # ---- recompute from the rounded concentrations --------------------------
    out = compute(m, sample_rows_from(df))
    for name, fields in FIELDS.items():
        m.write_tsv(T02 / name, out[name], fields)
        print(f"  wrote {name}")

    # ---- figures, using the module's own plotting code ----------------------
    fig_dir = ROOT / "revise" / "module02_rerun"
    (fig_dir / "figures").mkdir(parents=True, exist_ok=True)
    m.configure_matplotlib()
    m.make_correlation_figure(sample_rows_from(df), out["RCA_group_comparisons.tsv"],
                              out["RCA_correlations.tsv"], fig_dir)
    m.make_adjusted_effect_figure(out["RCA_regression_HC3.tsv"], fig_dir)
    src = fig_dir / "figures" / "Figure_RCA_technical_bias.pdf"
    for dest in (ROOT / "Submit" / "Figure_S3.pdf",
                 ROOT / "deposit" / "figures" / "Figure_S3.pdf"):
        shutil.copy2(src, dest)
        print(f"  copied -> {dest.relative_to(ROOT)}")

    marker = T02 / "STALE_AWAITING_HPC_RERUN.txt"
    if marker.exists():
        marker.unlink()
        print("  removed the STALE marker (module 02 is current again)")


if __name__ == "__main__":
    main()
