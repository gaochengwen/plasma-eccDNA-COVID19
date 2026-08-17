#!/usr/bin/env python3
"""Assert that VectorProbabilityCalculator matches base.ProbabilityCalculator.

The extensions use a numpy reimplementation of the matched-placement probability
because the original per-component Python loop is too slow for 33 peak sets. The
two must agree exactly; this checks them over many (track, chromosome, length)
combinations, including short and long lengths and every chromosome.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402
import run_chromatin_extensions as ext  # noqa: E402
import run_p0_burden_control as p0  # noqa: E402

FIELDS = (
    "placement_count",
    "p_full_interval",
    "p_midpoint",
    "p_junction_start",
    "p_junction_end",
    "p_junction_unit",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default="/gpfs/data/gao/Covid-eccDNA")
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument("--lengths", default="1,50,151,192,372,579,771,1000,2500,10000,100000")
    args = parser.parse_args()

    dirs = base.ensure_dirs(Path(args.outdir))
    chrom_sizes = base.read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    allowed = base.build_forbidden_and_allowed(dirs, chrom_sizes)
    peak_manifest = base.make_peak_manifest(Path(args.project_dir), dirs)
    tracks = p0.build_tracks_from_merged_beds(peak_manifest, dirs, chrom_sizes, allowed)

    reference = base.ProbabilityCalculator(allowed)
    fast = ext.VectorProbabilityCalculator(allowed)

    lengths = [int(v) for v in args.lengths.split(",")]
    checked = 0
    worst = 0.0
    failures = []
    for track in tracks.values():
        for chrom in base.CANON_CHROMS:
            for length in lengths:
                want = reference.probabilities(track, chrom, length)
                got = fast.probabilities(track, chrom, length)
                for field in FIELDS:
                    a, b = want[field], got[field]
                    checked += 1
                    if isinstance(a, float) and np.isnan(a):
                        if not (isinstance(b, float) and np.isnan(b)):
                            failures.append((track.peak_set_id, chrom, length, field, a, b))
                        continue
                    diff = abs(float(a) - float(b))
                    worst = max(worst, diff)
                    tolerance = 0 if field == "placement_count" else 1e-12
                    if diff > tolerance:
                        failures.append((track.peak_set_id, chrom, length, field, a, b))

    print(f"checked {checked} values across {len(tracks)} tracks x {len(base.CANON_CHROMS)} "
          f"chromosomes x {len(lengths)} lengths")
    print(f"max absolute difference {worst:.3e}")
    if failures:
        for row in failures[:20]:
            print(f"MISMATCH {row}")
        print(f"FAIL: {len(failures)} mismatches")
        return 1
    print("PASS: vectorized calculator reproduces the reference calculator")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
