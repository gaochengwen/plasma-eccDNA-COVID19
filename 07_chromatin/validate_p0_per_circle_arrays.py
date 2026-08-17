#!/usr/bin/env python3
"""Check that the P0 per-circle path reproduces the main analysis exactly.

analysis_downsample() computes the observed fraction and the exact matched
expectation as means of per-eccDNA quantities, whereas the main analysis
aggregates over (chromosome, length) bins. Drawing a subsample of size N without
replacement from N circles returns the whole sample, so for that degenerate case
the two code paths must agree to floating-point precision. This script asserts
that on a small number of samples and prints the largest discrepancy.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402
import run_p0_burden_control as p0  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default="/gpfs/data/gao/Covid-eccDNA")
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument("--n-samples", type=int, default=3, help="Smallest N samples, to keep this quick.")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    dirs = base.ensure_dirs(Path(args.outdir))
    chrom_sizes = base.read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    allowed = base.build_forbidden_and_allowed(dirs, chrom_sizes)
    forbidden_df = base.read_bed3(
        dirs.reference / "hg38_forbidden.blacklist_gaps.canonical.merged.bed", chrom_sizes
    )
    forbidden_arrays = base.interval_dict_to_arrays(base.merge_intervals_from_df(forbidden_df))
    peak_manifest = base.make_peak_manifest(Path(args.project_dir), dirs)
    # Reuse the manifest the main analysis already wrote: rebuilding it would call
    # samtools idxstats on all 78 BAMs, which this check does not need.
    sample_manifest = pd.read_csv(dirs.inputs / "eccdna_sample_manifest.tsv", sep="\t")
    tracks = p0.build_tracks_from_merged_beds(peak_manifest, dirs, chrom_sizes, allowed)

    reference = pd.read_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t")
    qc = pd.read_csv(dirs.tables / "eccdna_sample_qc.tsv", sep="\t")
    smallest = qc.nsmallest(args.n_samples, "eligible_eccdna_count")["sample_id"].tolist()
    base.log(f"validating samples: {smallest}")

    calculator = base.ProbabilityCalculator(allowed)
    worst = 0.0
    checked = 0
    for sample_id in smallest:
        row = sample_manifest[sample_manifest["sample_id"] == sample_id].iloc[0]
        ecc, _ = base.load_sample_eccdna(Path(row["eccdna_bed"]), chrom_sizes, forbidden_arrays)
        n = len(ecc)
        for track in tracks.values():
            arrays = p0.per_circle_arrays(ecc, track, calculator)
            got = {
                "full_interval": (
                    float(arrays["hit_full"].mean()),
                    float(np.nanmean(arrays["p_full"])),
                ),
                "midpoint": (
                    float(arrays["hit_mid"].mean()),
                    float(np.nanmean(arrays["p_mid"])),
                ),
                "junction_start_end": (
                    float((arrays["hit_start"].sum() + arrays["hit_end"].sum()) / (2.0 * n)),
                    float(np.nanmean(arrays["p_junc"])),
                ),
            }
            for mode, (obs, exp) in got.items():
                ref = reference[
                    (reference["sample_id"] == sample_id)
                    & (reference["peak_set_id"] == track.peak_set_id)
                    & (reference["mode"] == mode)
                ].iloc[0]
                for label, new, old in (
                    ("observed_fraction", obs, float(ref["observed_fraction"])),
                    ("expected_fraction", exp, float(ref["mean_shuffled_overlap_fraction"])),
                ):
                    diff = abs(new - old)
                    worst = max(worst, diff)
                    checked += 1
                    if diff > args.tolerance:
                        print(
                            f"MISMATCH {sample_id} {track.peak_set_id} {mode} {label}: "
                            f"p0={new!r} main={old!r} diff={diff:.3e}"
                        )
    print(f"checked {checked} values across {len(smallest)} samples; max absolute difference {worst:.3e}")
    if worst > args.tolerance:
        print("FAIL: P0 per-circle path does not reproduce the main analysis")
        return 1
    print("PASS: P0 per-circle path reproduces the main analysis")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
