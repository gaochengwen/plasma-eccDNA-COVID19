#!/usr/bin/env python3
"""P0 revision analyses: prove (or refute) burden-independence of chromatin enrichment.

Reviewer 2 major point 2 and Reviewer 3 major point 4 (final paragraph) both ask
whether the higher eccDNA signal in histone-marked regions in COVID-19 merely
reflects the much larger number of detected eccDNAs. The hg38 redo answered this
with an observed/expected ratio, but that ratio is itself correlated with the
number of detected eccDNAs, so the question is not yet settled.

This script adds the four analyses needed before the chromatin section can be
rewritten:

  within_group   One-sample tests that log2(observed/expected) differs from 0
                 inside each cohort separately. This is the formal support for
                 "plasma eccDNA preferentially derives from marked chromatin",
                 which so far rested only on the Fig. 4A meta-profiles.

  matched_subset COVID-19 vs HC restricted to samples whose eligible eccDNA
                 count falls in the range shared by both cohorts.

  covariate      HC3 robust OLS of log2 enrichment on group while adjusting for
                 detected eccDNA count, RCA yield, age and sex. Model conventions
                 (centred log2 RCA, age per 10 years from 70, male dummy, HC3)
                 follow Revise/RCA/run_rca_bias_analysis.py so the two analyses
                 are directly comparable.

  downsample     Every sample is randomly reduced to a common number of eligible
                 eccDNAs, and both the observed and the expected fraction are
                 recomputed from that subsample. This breaks the group/burden
                 collinearity by construction and is the decisive test.

The peak tracks are rebuilt from the merged hg38 BEDs already written by
run_chromatin_hg38_analysis.py, so liftOver is not repeated and the peak sets are
byte-identical to the main analysis.
"""

from __future__ import annotations

import argparse
import json
import math
import socket
import sys
import time
import zlib
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402


P0_SEED = 20260728
DEFAULT_TARGETS = (4000, 20000)
DEFAULT_REPEATS = 100
PRIMARY_MARK_ORDER = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3"]
GROUPS = ("COVID", "HC")


def log(message: str) -> None:
    base.log(message)


# --------------------------------------------------------------------------
# shared inputs
# --------------------------------------------------------------------------


def load_covariates(path: Path) -> pd.DataFrame:
    """Read the RCA metadata table used by the manuscript's RCA analysis."""
    meta = pd.read_csv(path, sep="\t")
    required = {"sample_id", "rca_concentration_ng_ul", "sex", "age_years", "group"}
    missing = required - set(meta.columns)
    if missing:
        raise RuntimeError(f"{path} is missing required columns: {sorted(missing)}")
    meta = meta.copy()
    meta["group"] = meta["group"].replace({"COVID-19": "COVID", "Healthy control": "HC"})
    meta["rca_concentration_ng_ul"] = meta["rca_concentration_ng_ul"].astype(float)
    meta["age_years"] = meta["age_years"].astype(float)
    meta["male"] = (meta["sex"].astype(str).str.lower() == "male").astype(float)
    return meta[["sample_id", "group", "rca_concentration_ng_ul", "age_years", "male"]]


def build_tracks_from_merged_beds(
    peak_manifest: pd.DataFrame,
    dirs: base.Dirs,
    chrom_sizes: Mapping[str, int],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
) -> Dict[str, base.Track]:
    """Reconstruct Track objects from the merged hg38 peak BEDs on disk.

    `allowed` is required, not optional: p_full_interval is derived from the
    peak-free components of the placement space, so a Track whose
    free_components are left empty silently reports an expected full-interval
    overlap fraction of 1.0.
    """
    tracks: Dict[str, base.Track] = {}
    for source in peak_manifest.to_dict(orient="records"):
        peak_id = source["peak_set_id"]
        stem = base.sanitize_id(peak_id)
        merged_bed = dirs.peaks / f"{stem}.liftover.hg38.canonical.filtered.merged.bed"
        if not merged_bed.exists():
            raise RuntimeError(
                f"missing {merged_bed}; run run_chromatin_hg38_analysis.py --stage analyze first"
            )
        df = base.read_bed3(merged_bed, chrom_sizes)
        merged = base.merge_intervals_from_df(df)
        arrays = base.interval_dict_to_arrays(merged)
        tracks[peak_id] = base.Track(
            peak_set_id=peak_id,
            dataset=source["dataset"],
            mark=source["mark"],
            analysis_tier=source["analysis_tier"],
            intervals=arrays,
            free_components={},
            prefix=base.prefix_tracks(arrays),
            merged_interval_count=int(sum(len(v) for v in merged.values())),
            coverage_bp=int(sum(e - s for vals in merged.values() for s, e in vals)),
        )
        log(f"track {peak_id}: {tracks[peak_id].merged_interval_count} merged hg38 intervals")
    base.finalize_track_free_components(tracks, allowed)
    return tracks


def bh_within_family(df: pd.DataFrame, family_cols: Sequence[str], p_col: str, out_col: str) -> pd.DataFrame:
    """Benjamini-Hochberg correction applied separately inside each family."""
    df = df.copy()
    df[out_col] = np.nan
    for _, idx in df.groupby(list(family_cols), dropna=False).groups.items():
        sub = df.loc[idx, p_col].tolist()
        df.loc[idx, out_col] = base.bh_fdr(sub)
    return df


def median_ci(values: np.ndarray, rng: np.random.Generator, n_boot: int = 10000) -> Tuple[float, float]:
    if len(values) == 0:
        return math.nan, math.nan
    draws = rng.integers(0, len(values), size=(n_boot, len(values)))
    medians = np.median(values[draws], axis=1)
    return float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


# --------------------------------------------------------------------------
# analysis 1: within-group enrichment vs the matched expectation
# --------------------------------------------------------------------------


def analysis_within_group(relative: pd.DataFrame, dirs: base.Dirs) -> pd.DataFrame:
    rng = np.random.default_rng(P0_SEED)
    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode"]
    for key_vals, sub in relative.groupby(keys, dropna=False):
        for group in GROUPS:
            vals = (
                sub.loc[sub["group"] == group, "log2_enrichment_ratio"]
                .replace([np.inf, -np.inf], np.nan)
                .dropna()
                .to_numpy(float)
            )
            n = len(vals)
            n_dropped = int((sub["group"] == group).sum() - n)
            if n == 0:
                continue
            wilcoxon_p = float(stats.wilcoxon(vals, alternative="two-sided").pvalue) if n >= 6 else math.nan
            n_pos = int((vals > 0).sum())
            n_neg = int((vals < 0).sum())
            sign_p = (
                float(stats.binomtest(n_pos, n_pos + n_neg, 0.5, alternative="two-sided").pvalue)
                if (n_pos + n_neg) > 0
                else math.nan
            )
            lo, hi = median_ci(vals, rng)
            rows.append(
                {
                    **dict(zip(keys, key_vals)),
                    "group": group,
                    "n_samples": n,
                    "n_samples_dropped_nonfinite": n_dropped,
                    "median_log2_enrichment": float(np.median(vals)),
                    "median_log2_ci_low": lo,
                    "median_log2_ci_high": hi,
                    "median_enrichment_ratio": float(2 ** np.median(vals)),
                    "ratio_ci_low": float(2 ** lo),
                    "ratio_ci_high": float(2 ** hi),
                    "n_samples_above_expected": n_pos,
                    "n_samples_below_expected": n_neg,
                    "wilcoxon_signed_rank_p_vs_0": wilcoxon_p,
                    "sign_test_p_vs_0": sign_p,
                }
            )
    out = pd.DataFrame(rows)
    out["fdr_family"] = out["analysis_tier"] + "|" + out["mode"] + "|" + out["group"]
    out = bh_within_family(out, ["fdr_family"], "wilcoxon_signed_rank_p_vs_0", "bh_fdr_within_family")
    out = out.sort_values(["analysis_tier", "peak_set_id", "mode", "group"])
    out.to_csv(dirs.tables / "p0_within_group_enrichment.tsv", sep="\t", index=False)
    log(f"wrote p0_within_group_enrichment.tsv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------
# analysis 2: burden-matched subset
# --------------------------------------------------------------------------


def analysis_matched_subset(relative: pd.DataFrame, dirs: base.Dirs) -> Tuple[pd.DataFrame, pd.DataFrame]:
    counts = (
        relative[["sample_id", "group", "eligible_eccdna_count"]]
        .drop_duplicates("sample_id")
        .reset_index(drop=True)
    )
    covid_n = counts.loc[counts["group"] == "COVID", "eligible_eccdna_count"]
    hc_n = counts.loc[counts["group"] == "HC", "eligible_eccdna_count"]
    low = int(covid_n.min())
    high = int(hc_n.max())
    counts["in_matched_window"] = counts["eligible_eccdna_count"].between(low, high)
    counts["matched_window_low"] = low
    counts["matched_window_high"] = high
    counts.to_csv(dirs.tables / "p0_burden_matched_membership.tsv", sep="\t", index=False)

    keep = set(counts.loc[counts["in_matched_window"], "sample_id"])
    sub_all = relative[relative["sample_id"].isin(keep)]
    kept = counts[counts["in_matched_window"]]
    n_covid = int((kept["group"] == "COVID").sum())
    n_hc = int((kept["group"] == "HC").sum())
    burden_p = float(
        stats.mannwhitneyu(
            kept.loc[kept["group"] == "COVID", "eligible_eccdna_count"].to_numpy(float),
            kept.loc[kept["group"] == "HC", "eligible_eccdna_count"].to_numpy(float),
            alternative="two-sided",
        ).pvalue
    )
    log(
        f"burden-matched window [{low}, {high}]: COVID n={n_covid}, HC n={n_hc}; "
        f"residual burden difference P={burden_p:.3g}"
    )

    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode"]
    for key_vals, sub in sub_all.groupby(keys, dropna=False):
        covid = (
            sub.loc[sub["group"] == "COVID", "log2_enrichment_ratio"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
        )
        hc = (
            sub.loc[sub["group"] == "HC", "log2_enrichment_ratio"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .to_numpy(float)
        )
        if len(covid) == 0 or len(hc) == 0:
            continue
        test = stats.mannwhitneyu(covid, hc, alternative="two-sided")
        rows.append(
            {
                **dict(zip(keys, key_vals)),
                "matched_window_low": low,
                "matched_window_high": high,
                "n_covid": len(covid),
                "n_hc": len(hc),
                "residual_burden_difference_p": burden_p,
                "median_covid": float(np.median(covid)),
                "median_hc": float(np.median(hc)),
                "median_difference_covid_minus_hc": float(np.median(covid) - np.median(hc)),
                "mannwhitneyu_p": float(test.pvalue),
                "cliffs_delta_covid_vs_hc": base.cliffs_delta(covid, hc),
            }
        )
    out = pd.DataFrame(rows)
    out["fdr_family"] = out["analysis_tier"] + "|" + out["mode"]
    out = bh_within_family(out, ["fdr_family"], "mannwhitneyu_p", "bh_fdr_within_family")
    out = out.sort_values(["analysis_tier", "peak_set_id", "mode"])
    out.to_csv(dirs.tables / "p0_burden_matched_subset.tsv", sep="\t", index=False)
    log(f"wrote p0_burden_matched_subset.tsv ({len(out)} rows)")
    return out, counts


# --------------------------------------------------------------------------
# analysis 3: covariate-adjusted HC3 robust regression
# --------------------------------------------------------------------------


def ols_hc3(y: np.ndarray, x: np.ndarray) -> dict:
    """OLS with HC3 robust standard errors (same implementation as the RCA analysis)."""
    n, p = x.shape
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residuals = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, xtx_inv, x)
    denominator = np.maximum(1.0 - leverage, np.finfo(float).eps)
    scaled_residual_sq = (residuals / denominator) ** 2
    meat = x.T @ (x * scaled_residual_sq[:, None])
    covariance = xtx_inv @ meat @ xtx_inv
    standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    dof = n - p
    t_values = beta / standard_errors
    p_values = 2.0 * stats.t.sf(np.abs(t_values), dof)
    critical = stats.t.ppf(0.975, dof)
    centered = y - np.mean(y)
    return {
        "beta": beta,
        "standard_error": standard_errors,
        "t_value": t_values,
        "p_value": p_values,
        "ci_low": beta - critical * standard_errors,
        "ci_high": beta + critical * standard_errors,
        "degrees_of_freedom": dof,
        "r_squared": 1.0 - float(np.sum(residuals**2)) / float(np.sum(centered**2)),
        "n": n,
    }


def variance_inflation_factors(x: np.ndarray, names: Sequence[str]) -> Dict[str, float]:
    """VIF for every non-intercept column, so collinearity is reported not hidden."""
    out: Dict[str, float] = {}
    for j, name in enumerate(names):
        if name == "intercept":
            continue
        others = [k for k in range(x.shape[1]) if k != j]
        fit = ols_hc3(x[:, j], x[:, others])
        r2 = fit["r_squared"]
        out[name] = float("inf") if r2 >= 1.0 else float(1.0 / (1.0 - r2))
    return out


def build_models(frame: pd.DataFrame) -> List[Tuple[str, List[str], np.ndarray]]:
    n = len(frame)
    ones = np.ones(n)
    covid = (frame["group"].to_numpy() == "COVID").astype(float)
    log2_rca = np.log2(frame["rca_concentration_ng_ul"].to_numpy(float))
    log2_rca_c = log2_rca - log2_rca.mean()
    age_per_10 = (frame["age_years"].to_numpy(float) - 70.0) / 10.0
    male = frame["male"].to_numpy(float)
    log10_n = np.log10(frame["eligible_eccdna_count"].to_numpy(float))
    log10_n_c = log10_n - log10_n.mean()
    log2_depth = np.log2(frame["mapped_alignments_idxstats"].to_numpy(float))
    log2_depth_c = log2_depth - log2_depth.mean()
    return [
        ("m1_group_only", ["intercept", "COVID_vs_HC"], np.column_stack([ones, covid])),
        (
            "m2_group_rca_age_sex",
            ["intercept", "COVID_vs_HC", "log2_RCA_per_doubling", "age_per_10_years", "male_vs_female"],
            np.column_stack([ones, covid, log2_rca_c, age_per_10, male]),
        ),
        (
            "m3_plus_eccdna_burden",
            [
                "intercept",
                "COVID_vs_HC",
                "log10_eligible_eccdna_count",
                "log2_RCA_per_doubling",
                "age_per_10_years",
                "male_vs_female",
            ],
            np.column_stack([ones, covid, log10_n_c, log2_rca_c, age_per_10, male]),
        ),
        (
            "m4_plus_burden_and_depth",
            [
                "intercept",
                "COVID_vs_HC",
                "log10_eligible_eccdna_count",
                "log2_mapped_alignments",
                "log2_RCA_per_doubling",
                "age_per_10_years",
                "male_vs_female",
            ],
            np.column_stack([ones, covid, log10_n_c, log2_depth_c, log2_rca_c, age_per_10, male]),
        ),
    ]


def analysis_covariate(relative: pd.DataFrame, covariates: pd.DataFrame, dirs: base.Dirs) -> pd.DataFrame:
    merged = relative.merge(covariates.drop(columns=["group"]), on="sample_id", how="left")
    if merged["rca_concentration_ng_ul"].isna().any():
        missing = sorted(merged.loc[merged["rca_concentration_ng_ul"].isna(), "sample_id"].unique())
        raise RuntimeError(f"covariates missing for samples: {missing}")

    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode"]
    for key_vals, sub in merged.groupby(keys, dropna=False):
        frame = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2_enrichment_ratio"])
        if len(frame) < 20:
            continue
        y = frame["log2_enrichment_ratio"].to_numpy(float)
        for model_name, names, x in build_models(frame):
            fit = ols_hc3(y, x)
            vifs = variance_inflation_factors(x, names)
            for j, term in enumerate(names):
                if term == "intercept":
                    continue
                rows.append(
                    {
                        **dict(zip(keys, key_vals)),
                        "model": model_name,
                        "term": term,
                        "n_samples": fit["n"],
                        "n_samples_dropped_nonfinite": int(len(sub) - len(frame)),
                        "beta": float(fit["beta"][j]),
                        "hc3_standard_error": float(fit["standard_error"][j]),
                        "t_value": float(fit["t_value"][j]),
                        "p_value": float(fit["p_value"][j]),
                        "ci_low": float(fit["ci_low"][j]),
                        "ci_high": float(fit["ci_high"][j]),
                        "vif": vifs.get(term, math.nan),
                        "model_r_squared": fit["r_squared"],
                        "degrees_of_freedom": fit["degrees_of_freedom"],
                    }
                )
    out = pd.DataFrame(rows)
    out["fdr_family"] = out["analysis_tier"] + "|" + out["mode"] + "|" + out["model"] + "|" + out["term"]
    out = bh_within_family(out, ["fdr_family"], "p_value", "bh_fdr_within_family")
    out = out.sort_values(["analysis_tier", "model", "peak_set_id", "mode", "term"])
    out.to_csv(dirs.tables / "p0_covariate_adjusted_hc3.tsv", sep="\t", index=False)
    log(f"wrote p0_covariate_adjusted_hc3.tsv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------
# analysis 4: downsampling to a common eccDNA count
# --------------------------------------------------------------------------


def per_circle_arrays(
    ecc: pd.DataFrame,
    track: base.Track,
    calculator: base.ProbabilityCalculator,
) -> Dict[str, np.ndarray]:
    """Per-eccDNA overlap indicators and matched-placement probabilities.

    Subsampling then reduces to taking means over a random index subset, because
    both the observed fraction and the exact matched expectation are means of
    per-circle quantities. This keeps the downsampled expectation matched to the
    length distribution of each individual subsample.
    """
    n = len(ecc)
    chrom_values = ecc["chrom"].to_numpy()
    starts = ecc["start"].to_numpy(dtype=np.int64)
    ends = ecc["end"].to_numpy(dtype=np.int64)
    lengths = ecc["length"].to_numpy(dtype=np.int64)
    midpoints = starts + (lengths // 2)

    hit_full = np.zeros(n, dtype=bool)
    hit_mid = np.zeros(n, dtype=bool)
    hit_start = np.zeros(n, dtype=bool)
    hit_end = np.zeros(n, dtype=bool)
    p_full = np.full(n, np.nan, dtype=float)
    p_mid = np.full(n, np.nan, dtype=float)
    p_junc = np.full(n, np.nan, dtype=float)

    for chrom in base.CANON_CHROMS:
        idx = np.where(chrom_values == chrom)[0]
        if len(idx) == 0:
            continue
        peak_starts, peak_ends = track.intervals[chrom]
        hit_full[idx] = base.interval_hits(peak_starts, peak_ends, starts[idx], ends[idx])
        hit_mid[idx] = base.point_hits(peak_starts, peak_ends, midpoints[idx])
        hit_start[idx] = base.point_hits(peak_starts, peak_ends, starts[idx])
        hit_end[idx] = base.point_hits(peak_starts, peak_ends, ends[idx] - 1)
        uniq, inverse = np.unique(lengths[idx], return_inverse=True)
        pf = np.empty(len(uniq))
        pm = np.empty(len(uniq))
        pj = np.empty(len(uniq))
        for k, length in enumerate(uniq):
            probs = calculator.probabilities(track, chrom, int(length))
            pf[k] = probs["p_full_interval"]
            pm[k] = probs["p_midpoint"]
            pj[k] = probs["p_junction_unit"]
        p_full[idx] = pf[inverse]
        p_mid[idx] = pm[inverse]
        p_junc[idx] = pj[inverse]

    return {
        "hit_full": hit_full,
        "hit_mid": hit_mid,
        "hit_start": hit_start,
        "hit_end": hit_end,
        "p_full": p_full,
        "p_mid": p_mid,
        "p_junc": p_junc,
    }


def _log2_ratio(observed: float, expected: float) -> float:
    if not np.isfinite(expected) or expected <= 0 or not np.isfinite(observed):
        return math.nan
    if observed <= 0:
        return -math.inf
    return float(math.log2(observed / expected))


def analysis_downsample(
    sample_manifest: pd.DataFrame,
    tracks: Dict[str, base.Track],
    dirs: base.Dirs,
    chrom_sizes: Mapping[str, int],
    forbidden_arrays: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
    targets: Sequence[int],
    repeats: int,
) -> pd.DataFrame:
    calculator = base.ProbabilityCalculator(allowed)
    rows: List[dict] = []
    progress_path = dirs.logs / "p0_downsample_progress.tsv"
    with open(progress_path, "wt", encoding="utf-8") as progress:
        progress.write("timestamp\tsample_id\tgroup\teligible\tseconds\n")

    samples = sample_manifest.sort_values(["group", "sample_id"]).to_dict(orient="records")
    for sample_index, sample in enumerate(samples, 1):
        started = time.time()
        sample_id = sample["sample_id"]
        log(f"[{sample_index}/{len(samples)}] downsampling {sample_id} ({sample['group']})")
        ecc, _ = base.load_sample_eccdna(Path(sample["eccdna_bed"]), chrom_sizes, forbidden_arrays)
        n_eligible = int(len(ecc))

        # One index draw per (target, repeat), shared across all peak sets so the
        # subsample is a property of the sample rather than of the track.
        draws: Dict[int, List[np.ndarray]] = {}
        for target in targets:
            if n_eligible < target:
                continue
            # crc32 rather than hash(): the built-in string hash is salted per
            # process, which would make the draws irreproducible across runs.
            rng = np.random.default_rng([P0_SEED, target, zlib.crc32(sample_id.encode("utf-8"))])
            draws[target] = [rng.choice(n_eligible, size=target, replace=False) for _ in range(repeats)]

        for track in tracks.values():
            arrays = per_circle_arrays(ecc, track, calculator)
            for target, index_sets in draws.items():
                acc = {mode: [] for mode in base.MODES}
                obs_acc = {mode: [] for mode in base.MODES}
                exp_acc = {mode: [] for mode in base.MODES}
                for idx in index_sets:
                    obs_full = float(arrays["hit_full"][idx].mean())
                    obs_mid = float(arrays["hit_mid"][idx].mean())
                    obs_junc = float(
                        (arrays["hit_start"][idx].sum() + arrays["hit_end"][idx].sum()) / (2.0 * target)
                    )
                    exp_full = float(np.nanmean(arrays["p_full"][idx]))
                    exp_mid = float(np.nanmean(arrays["p_mid"][idx]))
                    exp_junc = float(np.nanmean(arrays["p_junc"][idx]))
                    for mode, obs, exp in (
                        ("full_interval", obs_full, exp_full),
                        ("midpoint", obs_mid, exp_mid),
                        ("junction_start_end", obs_junc, exp_junc),
                    ):
                        acc[mode].append(_log2_ratio(obs, exp))
                        obs_acc[mode].append(obs)
                        exp_acc[mode].append(exp)
                for mode in base.MODES:
                    vals = np.asarray(acc[mode], dtype=float)
                    finite = vals[np.isfinite(vals)]
                    rows.append(
                        {
                            "sample_id": sample_id,
                            "group": sample["group"],
                            "peak_set_id": track.peak_set_id,
                            "dataset": track.dataset,
                            "mark": track.mark,
                            "analysis_tier": track.analysis_tier,
                            "mode": mode,
                            "downsample_target": target,
                            "downsample_repeats": repeats,
                            "eligible_eccdna_count": n_eligible,
                            "mean_observed_fraction": float(np.mean(obs_acc[mode])),
                            "mean_expected_fraction": float(np.mean(exp_acc[mode])),
                            "mean_log2_enrichment_ratio": float(finite.mean()) if len(finite) else math.nan,
                            "sd_log2_enrichment_ratio": float(finite.std(ddof=1)) if len(finite) > 1 else math.nan,
                            "n_repeats_finite": int(len(finite)),
                        }
                    )
        elapsed = time.time() - started
        with open(progress_path, "a", encoding="utf-8") as progress:
            progress.write(
                f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{sample_id}\t{sample['group']}\t{n_eligible}\t{elapsed:.1f}\n"
            )

    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(dirs.tables / "p0_downsampled_sample_enrichment.tsv", sep="\t", index=False)
    log(f"wrote p0_downsampled_sample_enrichment.tsv ({len(per_sample)} rows)")
    return per_sample


def downsample_group_stats(per_sample: pd.DataFrame, dirs: base.Dirs) -> pd.DataFrame:
    rng = np.random.default_rng(P0_SEED + 1)
    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode", "downsample_target"]
    for key_vals, sub in per_sample.groupby(keys, dropna=False):
        record = dict(zip(keys, key_vals))
        arrays = {}
        for group in GROUPS:
            arrays[group] = (
                sub.loc[sub["group"] == group, "mean_log2_enrichment_ratio"]
                .replace([np.inf, -np.inf], np.nan)
                .dropna()
                .to_numpy(float)
            )
        covid, hc = arrays["COVID"], arrays["HC"]
        if len(covid) == 0 or len(hc) == 0:
            continue
        test = stats.mannwhitneyu(covid, hc, alternative="two-sided")
        for group, vals in arrays.items():
            lo, hi = median_ci(vals, rng)
            suffix = "covid" if group == "COVID" else "hc"
            record[f"n_{suffix}"] = len(vals)
            record[f"median_{suffix}"] = float(np.median(vals))
            record[f"median_{suffix}_ci_low"] = lo
            record[f"median_{suffix}_ci_high"] = hi
            record[f"wilcoxon_p_vs_0_{suffix}"] = (
                float(stats.wilcoxon(vals, alternative="two-sided").pvalue) if len(vals) >= 6 else math.nan
            )
        record["median_difference_covid_minus_hc"] = record["median_covid"] - record["median_hc"]
        record["mannwhitneyu_p"] = float(test.pvalue)
        record["cliffs_delta_covid_vs_hc"] = base.cliffs_delta(covid, hc)
        rows.append(record)
    out = pd.DataFrame(rows)
    out["fdr_family"] = (
        out["analysis_tier"] + "|" + out["mode"] + "|" + out["downsample_target"].astype(str)
    )
    out = bh_within_family(out, ["fdr_family"], "mannwhitneyu_p", "bh_fdr_within_family")
    out = out.sort_values(["analysis_tier", "downsample_target", "peak_set_id", "mode"])
    out.to_csv(dirs.tables / "p0_downsampled_group_stats.tsv", sep="\t", index=False)
    log(f"wrote p0_downsampled_group_stats.tsv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------
# analysis 5: how the enrichment metric relates to detection burden
# --------------------------------------------------------------------------


def analysis_burden_association(relative: pd.DataFrame, dirs: base.Dirs) -> pd.DataFrame:
    """Spearman correlation of log2 enrichment with detected eccDNA count.

    Run within each cohort as well as across all samples. Within-cohort
    correlation is the key quantity: it shows whether the metric tracks
    detection burden independently of group, which is what makes the group
    contrast and the burden contrast statistically inseparable here.
    """
    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode"]
    for key_vals, sub in relative.groupby(keys, dropna=False):
        for stratum in ("COVID", "HC", "all_samples"):
            frame = sub if stratum == "all_samples" else sub[sub["group"] == stratum]
            frame = frame.replace([np.inf, -np.inf], np.nan).dropna(
                subset=["log2_enrichment_ratio", "eligible_eccdna_count"]
            )
            if len(frame) < 6:
                continue
            rho, p = stats.spearmanr(
                frame["eligible_eccdna_count"].to_numpy(float),
                frame["log2_enrichment_ratio"].to_numpy(float),
            )
            rows.append(
                {
                    **dict(zip(keys, key_vals)),
                    "stratum": stratum,
                    "n_samples": int(len(frame)),
                    "spearman_rho": float(rho),
                    "spearman_p": float(p),
                }
            )
    out = pd.DataFrame(rows)
    out["fdr_family"] = out["analysis_tier"] + "|" + out["mode"] + "|" + out["stratum"]
    out = bh_within_family(out, ["fdr_family"], "spearman_p", "bh_fdr_within_family")
    out = out.sort_values(["analysis_tier", "peak_set_id", "mode", "stratum"])
    out.to_csv(dirs.tables / "p0_burden_association.tsv", sep="\t", index=False)
    log(f"wrote p0_burden_association.tsv ({len(out)} rows)")
    return out


def analysis_downsample_agreement(
    relative: pd.DataFrame, per_sample: pd.DataFrame, dirs: base.Dirs
) -> pd.DataFrame:
    """Quantify how much downsampling changes each sample's enrichment estimate.

    The observed/expected ratio is an intensive per-sample quantity, so a random
    subsample is an unbiased estimator of the full-sample value. Documenting the
    near-perfect agreement is what licenses the (narrow) claim downsampling
    supports: the group contrast is not an artefact of the sevenfold difference
    in the number of units entering each sample's fraction. It is not a control
    for detection burden, and this table is the evidence for saying so.
    """
    merged = relative.merge(
        per_sample[
            [
                "sample_id",
                "peak_set_id",
                "mode",
                "downsample_target",
                "mean_log2_enrichment_ratio",
            ]
        ],
        on=["sample_id", "peak_set_id", "mode"],
        how="inner",
    ).replace([np.inf, -np.inf], np.nan)
    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode", "downsample_target"]
    for key_vals, sub in merged.groupby(keys, dropna=False):
        frame = sub.dropna(subset=["log2_enrichment_ratio", "mean_log2_enrichment_ratio"])
        if len(frame) < 6:
            continue
        full = frame["log2_enrichment_ratio"].to_numpy(float)
        down = frame["mean_log2_enrichment_ratio"].to_numpy(float)
        r, _ = stats.pearsonr(full, down)
        rows.append(
            {
                **dict(zip(keys, key_vals)),
                "n_samples": int(len(frame)),
                "pearson_r_full_vs_downsampled": float(r),
                "mean_absolute_difference": float(np.mean(np.abs(full - down))),
                "max_absolute_difference": float(np.max(np.abs(full - down))),
            }
        )
    out = pd.DataFrame(rows).sort_values(
        ["analysis_tier", "downsample_target", "peak_set_id", "mode"]
    )
    out.to_csv(dirs.tables / "p0_downsample_vs_full_agreement.tsv", sep="\t", index=False)
    log(f"wrote p0_downsample_vs_full_agreement.tsv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------


def make_p0_figure(dirs: base.Dirs) -> None:
    """Supplementary Figure S4, formatted to the Nature figure specification.

    183 mm double-column canvas, Arial (Liberation Sans where Arial is absent),
    5-7 pt label type, 8 pt bold lowercase panel labels, and every group separated by
    marker shape or by a step in lightness as well as by colour, so the figure survives
    greyscale and colour-blind viewing without hatching.

    Significance is annotated locally: each asterisk sits on the estimate or over the
    bracketed pair it tests, each panel is titled with the comparison and test that
    produced it, and one key at the foot of the figure defines the tiers. Every value
    and every q comes from the locked tables in `tables/`; nothing is recomputed here.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    base.setup_matplotlib()
    # setup_matplotlib only sets the body font. Without these, mathtext ($\delta$,
    # $\rho$, log$_2$) falls back to DejaVu Sans and the PDF ends up with a second,
    # non-compliant font family.
    body_font = mpl.rcParams["font.family"]
    body_font = body_font[0] if isinstance(body_font, (list, tuple)) else body_font
    mpl.rcParams.update(
        {
            "mathtext.fontset": "custom",
            "mathtext.rm": body_font,
            "mathtext.it": f"{body_font}:italic",
            "mathtext.bf": f"{body_font}:bold",
            "mathtext.sf": body_font,
            "mathtext.default": "regular",
            "hatch.linewidth": 0.4,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "legend.handlelength": 1.6,
            "legend.handletextpad": 0.5,
            "legend.labelspacing": 0.35,
            "legend.borderpad": 0.2,
        }
    )

    within = pd.read_csv(dirs.tables / "p0_within_group_enrichment.tsv", sep="\t")
    down_samples = pd.read_csv(dirs.tables / "p0_downsampled_sample_enrichment.tsv", sep="\t")
    down_stats = pd.read_csv(dirs.tables / "p0_downsampled_group_stats.tsv", sep="\t")
    full_stats = pd.read_csv(dirs.tables / "relative_enrichment_group_stats.tsv", sep="\t")
    matched = pd.read_csv(dirs.tables / "p0_burden_matched_subset.tsv", sep="\t")
    model = pd.read_csv(dirs.tables / "p0_covariate_adjusted_hc3.tsv", sep="\t")
    relative = pd.read_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t")
    burden = pd.read_csv(dirs.tables / "p0_burden_association.tsv", sep="\t")

    primary_target = int(down_stats["downsample_target"].min())

    # Palette: exact Hex values extracted from Figure 1b/c (#034E61 for HC, #D70000 for COVID-19).
    # Cohort colours are fixed once and reused in every panel.
    HC_COLOR, COVID_COLOR = "#034E61", "#D70000"          # Extracted from Figure 1b/c (HC, COVID-19)
    HC_FILL, COVID_FILL = "#A9BDC9", "#E6B0B0"            # 45% tints of the cohort inks for box faces
    # The three burden controls are one hue in three lightness steps, so the order of
    # increasing stringency is visible and the series stay separable in greyscale
    # without hatching. The nested models use the same idea in the COVID-19 hue.
    ANALYSIS = ["#9CC3C9", "#4E8E97", "#1E4F58"]
    MODEL = ["#E9A9A6", "#D70000", "#7B2A29"]
    ZERO_LINE = "#8A93A3"
    STAR_COLOR = "#1A1A1A"                                 # one ink for every significance mark
    NOTE_COLOR = "#3F3F3F"
    LABEL_SIZE, TICK_SIZE, LEGEND_SIZE, PANEL_SIZE = 6.0, 5.5, 5.0, 8.0
    NOTE_SIZE, STAR_SIZE = 5.0, 5.5

    x = np.arange(len(PRIMARY_MARK_ORDER))

    def primary(df, extra=None):
        sel = (df["analysis_tier"] == "primary") & (df["mode"] == "full_interval")
        if extra is not None:
            sel &= extra(df)
        return df[sel]

    def tidy(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=TICK_SIZE, length=2.0, pad=1.5)
        ax.yaxis.label.set_size(LABEL_SIZE)
        ax.xaxis.label.set_size(LABEL_SIZE)

    def sig_tier(q):
        """Manuscript-wide convention: * q<0.05, ** q<0.01, *** q<0.001, **** q<1e-4."""
        if q is None or pd.isna(q):
            return ""
        q = float(q)
        if q < 1e-4:
            return "****"
        if q < 1e-3:
            return "***"
        if q < 1e-2:
            return "**"
        if q < 0.05:
            return "*"
        return "NS"

    def panel_label(ax, letter, dx=-0.26):
        ax.text(
            dx, 1.075, letter, transform=ax.transAxes, fontsize=PANEL_SIZE,
            fontweight="bold", va="bottom", ha="left",
        )

    def panel_title(ax, text):
        """Name the comparison each panel's significance marks refer to."""
        ax.set_title(text, loc="left", fontsize=NOTE_SIZE, color=NOTE_COLOR, pad=2.5)

    def star(ax, x_pos, y_pos, tier, ha="center", va="bottom", rotation=0):
        if not tier:
            return
        ax.text(
            x_pos, y_pos, tier, ha=ha, va=va, rotation=rotation, color=STAR_COLOR,
            fontsize=STAR_SIZE if tier != "NS" else NOTE_SIZE,
            fontweight="bold" if tier != "NS" else "normal", zorder=6,
        )

    fig, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 130 / 25.4))

    # (a) enrichment relative to the matched expectation, within each cohort
    ax = axes[0][0]
    sub = primary(within).set_index(["mark", "group"])
    offsets = {"HC": -0.16, "COVID": 0.16}
    tops = {}
    for group, color, marker, label in (
        ("HC", HC_COLOR, "o", "HC"),
        ("COVID", COVID_COLOR, "s", "COVID-19"),
    ):
        med = np.array([sub.loc[(m, group), "median_log2_enrichment"] for m in PRIMARY_MARK_ORDER])
        lo = np.array([sub.loc[(m, group), "median_log2_ci_low"] for m in PRIMARY_MARK_ORDER])
        hi = np.array([sub.loc[(m, group), "median_log2_ci_high"] for m in PRIMARY_MARK_ORDER])
        tops[group] = hi
        ax.errorbar(
            x + offsets[group], med, yerr=[med - lo, hi - med], fmt=marker, color=color,
            markersize=3.2, capsize=1.5, elinewidth=0.7, capthick=0.7, linestyle="none",
            markeredgecolor="white", markeredgewidth=0.35, label=label, zorder=3,
        )
    ax.axhline(0.0, color=ZERO_LINE, linewidth=0.5, linestyle=(0, (3.5, 2.0)), zorder=0)
    # One mark carries one asterisk group when both cohorts fall in the same tier, and
    # a separate mark per cohort when they do not, so the label always sits over the
    # estimate it belongs to instead of in a blanket note.
    pad = 0.022
    for i, m in enumerate(PRIMARY_MARK_ORDER):
        tiers = {g: sig_tier(sub.loc[(m, g), "bh_fdr_within_family"]) for g in ("HC", "COVID")}
        if tiers["HC"] == tiers["COVID"]:
            star(ax, i, max(tops["HC"][i], tops["COVID"][i]) + pad, tiers["HC"])
        else:
            for g in ("HC", "COVID"):
                star(ax, i + offsets[g], tops[g][i] + pad, tiers[g])
    ax.set_xticks(x)
    ax.set_xticklabels(PRIMARY_MARK_ORDER, rotation=35, ha="right")
    ax.set_ylabel("log$_2$(observed / matched expected)")
    ax.set_ylim(-0.50, 0.53)
    ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left", ncol=1,
              handletextpad=0.4, borderaxespad=0.25, labelspacing=0.3)
    panel_title(ax, "vs matched expectation (Wilcoxon signed-rank)")
    tidy(ax)
    panel_label(ax, "a")

    # (b) COVID-HC effect size under each burden control
    ax = axes[0][1]
    full_sub = primary(full_stats).set_index("mark")
    down_sub = primary(down_stats, lambda d: d["downsample_target"] == primary_target).set_index("mark")
    matched_sub = primary(matched).set_index("mark")
    series = (
        ("All samples", full_sub, "bh_fdr_all_tests"),
        (f"Downsampled to {primary_target:,}", down_sub, "bh_fdr_within_family"),
        ("Burden-matched subset", matched_sub, "bh_fdr_within_family"),
    )
    height = 0.25
    for k, ((label, frame, q_col), color) in enumerate(zip(series, ANALYSIS)):
        pos = x + (1 - k) * height
        deltas = [frame.loc[m, "cliffs_delta_covid_vs_hc"] for m in PRIMARY_MARK_ORDER]
        ax.barh(pos, deltas, height=height * 0.88, color=color, edgecolor="none",
                label=label, zorder=2)
        # One ink for every asterisk: the row the mark sits on already identifies the
        # control, so the label does not have to be colour-matched to be read.
        for m, d, yy in zip(PRIMARY_MARK_ORDER, deltas, pos):
            mark_tier = sig_tier(frame.loc[m, q_col])
            if mark_tier and mark_tier != "NS":
                star(ax, d + (0.025 if d >= 0 else -0.025), yy, mark_tier,
                     ha="left" if d >= 0 else "right", va="center")
    ax.axvline(0.0, color=ZERO_LINE, linewidth=0.5, linestyle=(0, (3.5, 2.0)), zorder=0)
    ax.set_yticks(x)
    ax.set_yticklabels(PRIMARY_MARK_ORDER)
    ax.set_xlabel("Cliff's $\\delta$ (COVID-19 vs HC)")
    ax.set_xlim(-0.80, 0.92)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_ylim(-1.55, len(PRIMARY_MARK_ORDER) - 0.42)
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left", handlelength=1.2,
              handletextpad=0.45, labelspacing=0.3, borderaxespad=0.2)
    panel_title(ax, "COVID-19 vs HC (Mann–Whitney)")
    tidy(ax)
    panel_label(ax, "b")

    # (c) HC3 group coefficient across nested models
    ax = axes[0][2]
    model_order = [
        ("m1_group_only", "Group only", "o"),
        ("m2_group_rca_age_sex", "+ RCA, age, sex", "s"),
        ("m3_plus_eccdna_burden", "+ eccDNA burden", "^"),
    ]
    sub = primary(model, lambda d: d["term"] == "COVID_vs_HC").set_index(["model", "mark"])
    for k, ((name, label, marker), color) in enumerate(zip(model_order, MODEL)):
        pos = x + (1 - k) * 0.23
        betas = np.array([sub.loc[(name, m), "beta"] for m in PRIMARY_MARK_ORDER])
        lo = np.array([sub.loc[(name, m), "ci_low"] for m in PRIMARY_MARK_ORDER])
        hi = np.array([sub.loc[(name, m), "ci_high"] for m in PRIMARY_MARK_ORDER])
        qs = np.array([sub.loc[(name, m), "bh_fdr_within_family"] for m in PRIMARY_MARK_ORDER])
        ax.errorbar(
            betas, pos, xerr=[betas - lo, hi - betas], fmt="none", ecolor=color,
            elinewidth=0.7, capsize=1.3, capthick=0.7, zorder=2,
        )
        significant = qs < 0.05
        ax.plot(
            betas[significant], pos[significant], marker, color=color, markersize=3.0,
            linestyle="none", markeredgecolor="white", markeredgewidth=0.35, label=label,
            zorder=3,
        )
        ax.plot(
            betas[~significant], pos[~significant], marker, markerfacecolor="white",
            markeredgecolor=color, markersize=3.0, linestyle="none", markeredgewidth=0.7,
            zorder=3,
        )
    ax.axvline(0.0, color=ZERO_LINE, linewidth=0.5, linestyle=(0, (3.5, 2.0)), zorder=0)
    ax.set_yticks(x)
    ax.set_yticklabels(PRIMARY_MARK_ORDER)
    ax.set_xlabel("Adjusted COVID-19 − HC difference\nin log$_2$ enrichment (HC3 95% CI)")
    ax.set_xlim(-0.21, 0.185)
    ax.set_xticks([-0.1, 0.0, 0.1])
    ax.set_ylim(-1.55, len(PRIMARY_MARK_ORDER) - 0.42)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, frameon=False, fontsize=LEGEND_SIZE, loc="lower left",
              handlelength=0.9, handletextpad=0.45, labelspacing=0.3, borderaxespad=0.2)
    ax.text(0.99, 0.02, "filled, $\\it{q}$ < 0.05\nopen, NS", transform=ax.transAxes,
            fontsize=NOTE_SIZE, color=NOTE_COLOR, va="bottom", ha="right", linespacing=1.35)
    panel_title(ax, "COVID-19 − HC, HC3 robust regression")
    tidy(ax)
    panel_label(ax, "c")

    # (d, e) the group/burden collinearity for the two marks with unadjusted signal
    # These panels are about the slope, not the offset from the matched expectation,
    # so the 0 reference line is omitted; the axis label states the quantity.
    rel_p = primary(relative)
    burden_p = primary(burden).set_index(["mark", "stratum"])
    for ax, mark, letter, note_xy in (
        (axes[1][0], "H3K27me3", "d", (0.03, 0.97)),
        (axes[1][1], "H3K9me3", "e", (0.97, 0.97)),
    ):
        sub = rel_p[rel_p["mark"] == mark]
        notes = []
        for group, color, marker, label in (
            ("HC", HC_COLOR, "o", "HC"),
            ("COVID", COVID_COLOR, "^", "COVID-19"),
        ):
            ss = sub[sub["group"] == group]
            xv = np.log10(ss["eligible_eccdna_count"].to_numpy(float))
            yv = ss["log2_enrichment_ratio"].to_numpy(float)
            ax.plot(xv, yv, marker, color=color, markersize=2.6, linestyle="none",
                    markeredgecolor="white", markeredgewidth=0.3, alpha=0.95, label=label)
            if len(xv) > 2:
                slope, intercept = np.polyfit(xv, yv, 1)
                xs = np.linspace(xv.min(), xv.max(), 20)
                ax.plot(xs, slope * xs + intercept, "-", color=color, linewidth=0.9)
            # Read the correlation and its BH-adjusted q from the locked table rather
            # than recomputing it, and show both so a nominally significant P is not
            # mistaken for a significant result after correction.
            row = burden_p.loc[(mark, group)]
            # True minus sign on the coefficient, not a hyphen (the cohort label has one).
            rho = f"{row['spearman_rho']:+.2f}".replace("-", "−")
            notes.append(
                (
                    f"{label}: $\\it{{\\rho}}$ = {rho}, "
                    f"$\\it{{P}}$ = {row['spearman_p']:.3f}, "
                    f"$\\it{{q}}$ = {row['bh_fdr_within_family']:.3f}",
                    color,
                )
            )
        ax.set_xlabel("log$_{10}$ detected eccDNAs")
        ax.set_ylabel("log$_2$(observed / matched expected)")
        # Headroom for the two annotation lines, so they never sit on a data point.
        y_all = sub["log2_enrichment_ratio"].to_numpy(float)
        y_range = float(y_all.max() - y_all.min())
        ax.set_ylim(y_all.min() - 0.07 * y_range, y_all.max() + 0.34 * y_range)
        # The two lines are coloured by cohort, which removes the need for a legend.
        for line_index, (text, color) in enumerate(notes):
            ax.text(
                note_xy[0], note_xy[1] - line_index * 0.075, text, transform=ax.transAxes,
                fontsize=LEGEND_SIZE, color=color, va="top",
                ha="left" if note_xy[0] < 0.5 else "right",
                bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=0.8),
            )
        panel_title(ax, mark)
        tidy(ax)
        panel_label(ax, letter)

    # (f) downsampled per-sample enrichment at a common depth
    ax = axes[1][2]
    sub = down_samples[
        (down_samples["analysis_tier"] == "primary")
        & (down_samples["mode"] == "full_interval")
        & (down_samples["downsample_target"] == primary_target)
    ]
    box_offset = 0.185
    box_top = np.full(len(PRIMARY_MARK_ORDER), -np.inf)
    box_bottom = np.full(len(PRIMARY_MARK_ORDER), np.inf)
    for offset, group, color, fill, label in (
        (-box_offset, "HC", HC_COLOR, HC_FILL, "HC"),
        (box_offset, "COVID", COVID_COLOR, COVID_FILL, "COVID-19"),
    ):
        data = [
            sub.loc[(sub["mark"] == m) & (sub["group"] == group), "mean_log2_enrichment_ratio"]
            .replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
            for m in PRIMARY_MARK_ORDER
        ]
        bp = ax.boxplot(data, positions=x + offset, widths=0.28, patch_artist=True,
                        showfliers=False,
                        medianprops={"color": "black", "linewidth": 0.7})
        for patch in bp["boxes"]:
            patch.set_facecolor(fill)
            patch.set_edgecolor(color)
            patch.set_linewidth(0.5)
        for line in bp["whiskers"] + bp["caps"]:
            line.set_color(color)
            line.set_linewidth(0.5)
        for i, cap_pair in enumerate(zip(bp["caps"][0::2], bp["caps"][1::2])):
            lower, upper = (float(c.get_ydata()[0]) for c in cap_pair)
            box_top[i] = max(box_top[i], lower, upper)
            box_bottom[i] = min(box_bottom[i], lower, upper)
        ax.plot([], [], "s", color=fill, markeredgecolor=color, markeredgewidth=0.5,
                markersize=3.2, linestyle="none", label=label)
    ax.axhline(0.0, color=ZERO_LINE, linewidth=0.5, linestyle=(0, (3.5, 2.0)), zorder=0)
    # A bracket over the two boxes of a mark states exactly which pair each asterisk
    # tests; a row of labels along the top of the panel would not.
    f_stats = primary(down_stats, lambda d: d["downsample_target"] == primary_target).set_index("mark")
    span = float(box_top.max() - box_bottom.min())
    lift, drop = 0.045 * span, 0.022 * span
    for xi, m in zip(x, PRIMARY_MARK_ORDER):
        bracket_y = box_top[xi] + lift
        ax.plot(
            [xi - box_offset, xi - box_offset, xi + box_offset, xi + box_offset],
            [bracket_y - drop, bracket_y, bracket_y, bracket_y - drop],
            color=STAR_COLOR, linewidth=0.5, solid_capstyle="butt", zorder=5,
        )
        star(ax, xi, bracket_y + 0.1 * drop, sig_tier(f_stats.loc[m, "bh_fdr_within_family"]))
    ax.set_xticks(x)
    ax.set_xticklabels(PRIMARY_MARK_ORDER, rotation=35, ha="right")
    ax.set_xlim(-0.6, len(PRIMARY_MARK_ORDER) - 0.4)
    ax.set_ylim(box_bottom.min() - 0.08 * span, box_top.max() + 0.15 * span)
    ax.set_ylabel("log$_2$(observed / matched expected)")
    ax.legend(frameon=False, fontsize=LEGEND_SIZE, loc="lower left", ncol=1,
              handletextpad=0.4, borderaxespad=0.25, labelspacing=0.3)
    panel_title(ax, f"COVID-19 vs HC at {primary_target:,} eccDNAs (Mann–Whitney)")
    tidy(ax)
    panel_label(ax, "f")

    fig.subplots_adjust(left=0.075, right=0.99, top=0.925, bottom=0.115, wspace=0.40, hspace=0.50)
    # One key for every asterisk in the figure, so no panel depends on the caption
    # for its significance code.
    fig.text(
        0.075, 0.018,
        "Significance: *$\\it{q}$ < 0.05   **$\\it{q}$ < 0.01   ***$\\it{q}$ < 0.001   "
        "****$\\it{q}$ < 10$^{-4}$   NS or no asterisk, $\\it{q}$ $\\geq$ 0.05.   "
        "$\\it{q}$, Benjamini–Hochberg FDR within the test family.",
        fontsize=NOTE_SIZE, color=NOTE_COLOR, va="bottom", ha="left",
    )
    for ext in ("pdf", "svg", "png"):
        fig.savefig(dirs.figures / f"p0_burden_control.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    log("wrote figures/p0_burden_control.{pdf,svg,png}")

    import shutil
    revise_dir = Path("/home/gao/eccDNA/revise")
    if revise_dir.exists():
        for ext in ("pdf", "svg", "png"):
            shutil.copy2(dirs.figures / f"p0_burden_control.{ext}", revise_dir / f"Figure_S4.{ext}")
        log("copied to revise/Figure_S4.{pdf,svg,png}")
    deposit_fig = Path("/home/gao/eccDNA/analyse/chromatin/deposit/figures")
    if deposit_fig.exists():
        shutil.copy2(dirs.figures / "p0_burden_control.pdf", deposit_fig / "p0_burden_control.pdf")
        log("copied to deposit/figures/p0_burden_control.pdf")

def write_summary(dirs: base.Dirs, targets: Sequence[int], repeats: int) -> None:
    within = pd.read_csv(dirs.tables / "p0_within_group_enrichment.tsv", sep="\t")
    matched = pd.read_csv(dirs.tables / "p0_burden_matched_subset.tsv", sep="\t")
    model = pd.read_csv(dirs.tables / "p0_covariate_adjusted_hc3.tsv", sep="\t")
    down = pd.read_csv(dirs.tables / "p0_downsampled_group_stats.tsv", sep="\t")
    full_stats = pd.read_csv(dirs.tables / "relative_enrichment_group_stats.tsv", sep="\t")

    def fmt(df: pd.DataFrame, cols: Sequence[str]) -> str:
        # No rounding here: markdown_table already renders floats as %.4g, whereas
        # round(6) would collapse p values such as 3.64e-12 to 0.
        return base.markdown_table(df[list(cols)])

    primary_target = int(min(targets))
    lines: List[str] = []
    lines.append("# P0 burden-control analyses\n")
    lines.append(f"Run time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Host: {socket.gethostname()}")
    lines.append(f"Output directory: `{dirs.outdir}`")
    lines.append(f"Downsample targets: {', '.join(str(t) for t in targets)}; repeats per target: {repeats}")
    lines.append(f"Seed: {P0_SEED}\n")

    lines.append("## 1. Enrichment relative to the matched expectation, within each cohort\n")
    lines.append(
        "One-sample two-sided Wilcoxon signed-rank test of log2(observed/expected) against 0, "
        "computed separately in COVID-19 and HC. This is the formal test behind the statement "
        "that plasma eccDNA is non-randomly distributed with respect to reference histone marks.\n"
    )
    sub = within[(within["analysis_tier"] == "primary") & (within["mode"] == "full_interval")]
    lines.append(
        fmt(
            sub,
            [
                "mark",
                "group",
                "n_samples",
                "median_enrichment_ratio",
                "ratio_ci_low",
                "ratio_ci_high",
                "n_samples_above_expected",
                "wilcoxon_signed_rank_p_vs_0",
                "bh_fdr_within_family",
            ],
        )
    )

    lines.append("\n## 2. COVID-19 vs HC: all samples versus burden-matched subset\n")
    merged = full_stats[
        (full_stats["analysis_tier"] == "primary") & (full_stats["mode"] == "full_interval")
    ][["mark", "cliffs_delta_covid_vs_hc", "bh_fdr_all_tests"]].rename(
        columns={"cliffs_delta_covid_vs_hc": "delta_all_samples", "bh_fdr_all_tests": "q_all_samples"}
    ).merge(
        matched[(matched["analysis_tier"] == "primary") & (matched["mode"] == "full_interval")][
            ["mark", "n_covid", "n_hc", "cliffs_delta_covid_vs_hc", "bh_fdr_within_family"]
        ].rename(
            columns={
                "cliffs_delta_covid_vs_hc": "delta_matched_subset",
                "bh_fdr_within_family": "q_matched_subset",
            }
        ),
        on="mark",
    )
    lines.append(fmt(merged, list(merged.columns)))

    lines.append(f"\n## 3. COVID-19 vs HC after downsampling to {primary_target:,} eccDNAs per sample\n")
    sub = down[
        (down["analysis_tier"] == "primary")
        & (down["mode"] == "full_interval")
        & (down["downsample_target"] == primary_target)
    ]
    lines.append(
        fmt(
            sub,
            [
                "mark",
                "n_covid",
                "n_hc",
                "median_covid",
                "median_hc",
                "median_difference_covid_minus_hc",
                "mannwhitneyu_p",
                "cliffs_delta_covid_vs_hc",
                "bh_fdr_within_family",
            ],
        )
    )

    lines.append("\n## 4. HC3 robust regression: group effect with and without the burden term\n")
    sub = model[
        (model["analysis_tier"] == "primary")
        & (model["mode"] == "full_interval")
        & (model["term"].isin(["COVID_vs_HC", "log10_eligible_eccdna_count"]))
        & (model["model"].isin(["m1_group_only", "m2_group_rca_age_sex", "m3_plus_eccdna_burden"]))
    ]
    lines.append(
        fmt(
            sub,
            [
                "mark",
                "model",
                "term",
                "beta",
                "hc3_standard_error",
                "ci_low",
                "ci_high",
                "p_value",
                "bh_fdr_within_family",
                "vif",
            ],
        )
    )

    lines.append("\n## 5. Does downsampling change each sample's estimate?\n")
    lines.append(
        "The observed/expected ratio is an intensive per-sample quantity, so a random subsample "
        "is an unbiased estimator of the full-sample value. Near-perfect agreement below means "
        "downsampling equalizes the precision of each sample's estimate but does NOT control for "
        "detection burden, and must not be presented as if it did.\n"
    )
    agreement = pd.read_csv(dirs.tables / "p0_downsample_vs_full_agreement.tsv", sep="\t")
    sub = agreement[
        (agreement["analysis_tier"] == "primary")
        & (agreement["mode"] == "full_interval")
        & (agreement["downsample_target"] == primary_target)
    ]
    lines.append(
        fmt(
            sub,
            [
                "mark",
                "n_samples",
                "pearson_r_full_vs_downsampled",
                "mean_absolute_difference",
                "max_absolute_difference",
            ],
        )
    )

    lines.append("\n## 6. Association between the enrichment metric and detection burden\n")
    lines.append(
        "Spearman correlation of log2 enrichment with the eligible eccDNA count, within each "
        "cohort and across all samples. A within-cohort association means the metric tracks "
        "detection burden independently of group, which is what makes the group contrast and the "
        "burden contrast inseparable in this design.\n"
    )
    burden = pd.read_csv(dirs.tables / "p0_burden_association.tsv", sep="\t")
    sub = burden[(burden["analysis_tier"] == "primary") & (burden["mode"] == "full_interval")]
    lines.append(
        fmt(sub, ["mark", "stratum", "n_samples", "spearman_rho", "spearman_p", "bh_fdr_within_family"])
    )

    lines.append("\n## Output tables\n")
    for name in (
        "p0_within_group_enrichment.tsv",
        "p0_burden_matched_membership.tsv",
        "p0_burden_matched_subset.tsv",
        "p0_covariate_adjusted_hc3.tsv",
        "p0_downsampled_sample_enrichment.tsv",
        "p0_downsampled_group_stats.tsv",
        "p0_downsample_vs_full_agreement.tsv",
        "p0_burden_association.tsv",
    ):
        lines.append(f"- `tables/{name}`")
    lines.append("- `figures/p0_burden_control.{pdf,svg,png}`")
    lines.append("")

    (dirs.outdir / "README_p0_burden_control.md").write_text("\n".join(lines), encoding="utf-8")
    log("wrote README_p0_burden_control.md")


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default="/gpfs/data/gao/Covid-eccDNA")
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument(
        "--covariates",
        default="/gpfs/data/gao/Covid-eccDNA/Revise/RCA/input/rca_metadata.tsv",
        help="Sample metadata with rca_concentration_ng_ul, age_years, sex.",
    )
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument(
        "--targets",
        default=",".join(str(t) for t in DEFAULT_TARGETS),
        help="Comma-separated downsampling targets; samples below a target are skipped for it.",
    )
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument(
        "--stages",
        default="within,matched,covariate,downsample,burden,figure",
        help="Comma-separated subset of stages to run.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    stages = {s.strip() for s in args.stages.split(",") if s.strip()}
    targets = tuple(int(t) for t in args.targets.split(",") if t.strip())
    dirs = base.ensure_dirs(Path(args.outdir))

    relative = pd.read_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t")
    log(f"loaded sample_relative_enrichment.tsv ({len(relative)} rows)")

    if "within" in stages:
        analysis_within_group(relative, dirs)
    if "matched" in stages:
        analysis_matched_subset(relative, dirs)
    if "covariate" in stages:
        analysis_covariate(relative, load_covariates(Path(args.covariates)), dirs)
    if "downsample" in stages:
        chrom_sizes = base.read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
        allowed = base.build_forbidden_and_allowed(dirs, chrom_sizes)
        forbidden_df = base.read_bed3(
            dirs.reference / "hg38_forbidden.blacklist_gaps.canonical.merged.bed", chrom_sizes
        )
        forbidden_arrays = base.interval_dict_to_arrays(base.merge_intervals_from_df(forbidden_df))
        peak_manifest = base.make_peak_manifest(Path(args.project_dir), dirs)
        sample_manifest = base.make_sample_manifest(Path(args.project_dir), dirs, args.samtools)
        tracks = build_tracks_from_merged_beds(peak_manifest, dirs, chrom_sizes, allowed)
        per_sample = analysis_downsample(
            sample_manifest, tracks, dirs, chrom_sizes, forbidden_arrays, allowed, targets, args.repeats
        )
        downsample_group_stats(per_sample, dirs)
    if "burden" in stages:
        analysis_burden_association(relative, dirs)
        analysis_downsample_agreement(
            relative,
            pd.read_csv(dirs.tables / "p0_downsampled_sample_enrichment.tsv", sep="\t"),
            dirs,
        )
    if "figure" in stages:
        make_p0_figure(dirs)
        write_summary(dirs, targets, args.repeats)

    manifest_path = dirs.outdir / "run_manifest_p0.json"
    manifest_path.write_text(
        json.dumps(
            {
                "stages": sorted(stages),
                "downsample_targets": list(targets),
                "downsample_repeats": args.repeats,
                "seed": P0_SEED,
                "covariates": args.covariates,
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
