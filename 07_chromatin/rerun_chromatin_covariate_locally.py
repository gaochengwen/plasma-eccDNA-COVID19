#!/usr/bin/env python3
"""Regenerate deposit/tables/06_chromatin/p0_covariate_adjusted_hc3.tsv with the
rounded RCA concentrations.

The covariate stage of run_p0_burden_control.py needs nothing but
`sample_relative_enrichment.tsv` (already deposited) and the RCA metadata, so it
can be driven here without the BAMs, peak BEDs or the reference genome. The
module's own `load_covariates`, `build_models`, `ols_hc3`, `bh_within_family` and
`analysis_covariate` are called directly -- this is the published code, not a
reimplementation.

The run is validated first: driven with the UNROUNDED concentrations it must
reproduce the currently deposited table to machine precision, otherwise nothing
is written. That guard is what caught an earlier bad refit of this same file.

    python3 revise/scripts/fixes/rerun_chromatin_covariate_locally.py --check
    python3 revise/scripts/fixes/rerun_chromatin_covariate_locally.py

Run from the repository root.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MODDIR = ROOT / "deposit" / "scripts" / "06_chromatin"
T06 = ROOT / "deposit" / "tables" / "06_chromatin"
T02 = ROOT / "deposit" / "tables" / "02_rca_technical_bias"
BACKUP = ROOT / "revise" / "backups" / "rca_rounding"

TARGET = T06 / "p0_covariate_adjusted_hc3.tsv"
PREROUND = BACKUP / "deposit_tables_06_chromatin_p0_covariate_adjusted_hc3.tsv.preround"
RELATIVE = T06 / "sample_relative_enrichment.tsv"

# the module's own input metadata: rounded copy and the pre-rounding original
META_ROUNDED = T02 / "RCA_sample_metrics.tsv"
META_UNROUNDED = BACKUP / "deposit_tables_02_rca_technical_bias_RCA_sample_metrics.tsv.preround"

KEY = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode", "model", "term"]
NUMERIC = ["beta", "hc3_standard_error", "t_value", "p_value", "ci_low", "ci_high",
           "vif", "model_r_squared", "bh_fdr_within_family"]
EXACT = ["n_samples", "n_samples_dropped_nonfinite", "degrees_of_freedom", "fdr_family"]


def load_module():
    sys.path.insert(0, str(MODDIR))
    spec = importlib.util.spec_from_file_location("p0mod", MODDIR / "run_p0_burden_control.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["p0mod"] = m
    spec.loader.exec_module(m)
    return m


def run_covariate(m, meta_path: Path) -> pd.DataFrame:
    """Drive the module's covariate stage into a throwaway output tree."""
    with tempfile.TemporaryDirectory(prefix="p0cov_") as tmp:
        dirs = m.base.ensure_dirs(Path(tmp))
        shutil.copy2(RELATIVE, dirs.tables / "sample_relative_enrichment.tsv")
        relative = pd.read_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t")
        covariates = m.load_covariates(meta_path)
        return m.analysis_covariate(relative, covariates, dirs)


def align(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(KEY).reset_index(drop=True)


def compare(new: pd.DataFrame, ref_path: Path, tol=1e-9) -> list[str]:
    ref = align(pd.read_csv(ref_path, sep="\t"))
    got = align(new[ref.columns.tolist()])
    if got.shape != ref.shape:
        return [f"shape {got.shape} vs {ref.shape}"]
    bad = []
    for c in EXACT:
        if (got[c].astype(str).values != ref[c].astype(str).values).any():
            bad.append(f"{c}: differs")
    for c in NUMERIC:
        x = got[c].to_numpy(float)
        y = ref[c].to_numpy(float)
        fin = np.isfinite(x) & np.isfinite(y)
        if (np.isfinite(x) != np.isfinite(y)).any():
            bad.append(f"{c}: finiteness differs")
            continue
        d = np.abs(x[fin] - y[fin])
        rel = d / np.maximum(np.abs(y[fin]), 1e-12)
        if ((d > tol) & (rel > 1e-6)).any():
            bad.append(f"{c}: max rel diff {rel.max():.3g}")
    return bad


def headline(df: pd.DataFrame, label: str) -> None:
    prim = df[(df.analysis_tier == "primary") & (df["mode"] == "full_interval")
              & (df.term == "COVID_vs_HC")]
    print(f"  {label}")
    for mod in ("m1_group_only", "m2_group_rca_age_sex", "m3_plus_eccdna_burden",
                "m4_plus_burden_and_depth"):
        x = prim[prim.model == mod]
        if x.empty:
            continue
        print(f"    {mod:<26} q {x.bh_fdr_within_family.min():.4g} - "
              f"{x.bh_fdr_within_family.max():.4g}   min P {x.p_value.min():.4g}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="validate only, write nothing")
    args = ap.parse_args()

    m = load_module()

    ref_path = PREROUND if PREROUND.exists() else TARGET
    got = run_covariate(m, META_UNROUNDED if META_UNROUNDED.exists() else META_ROUNDED)
    problems = compare(got, ref_path)
    if problems:
        print("VALIDATION FAILED - nothing written:")
        for p in problems:
            print("   ", p)
        sys.exit(1)
    print("  validation: the original module code reproduces the deposited table "
          "exactly from the unrounded concentrations")
    headline(got, "before rounding:")

    if args.check:
        return

    out = run_covariate(m, META_ROUNDED)
    headline(out, "after rounding:")

    ref = align(pd.read_csv(ref_path, sep="\t"))
    new = align(out[ref.columns.tolist()])
    moved = {c: int((np.abs(new[c].to_numpy(float) - ref[c].to_numpy(float))
                     > 1e-12).sum()) for c in NUMERIC}
    print("  rows changed per column: " + ", ".join(f"{k}={v}" for k, v in moved.items()))

    out[ref.columns.tolist()].to_csv(TARGET, sep="\t", index=False)
    print(f"  wrote {TARGET.relative_to(ROOT)} ({len(out)} rows)")

    # the module's own deposited input file must carry the rounded values too
    meta = pd.read_csv(T02 / "rca_metadata.tsv", sep="\t")
    if (meta.rca_concentration_ng_ul % 1 != 0).any():
        shutil.copy2(T02 / "rca_metadata.tsv",
                     BACKUP / "deposit_tables_02_rca_technical_bias_rca_metadata.tsv.preround")
        meta["rca_concentration_ng_ul"] = meta.rca_concentration_ng_ul.round().astype(float)
        meta.to_csv(T02 / "rca_metadata.tsv", sep="\t", index=False)
        print("  rounded deposit/tables/02_rca_technical_bias/rca_metadata.tsv "
              "(the module input, previously still unrounded)")


if __name__ == "__main__":
    main()
