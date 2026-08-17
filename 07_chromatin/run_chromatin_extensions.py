#!/usr/bin/env python3
"""Cross-cell-type, mappability and peak-centred-density extensions to the hg38 analysis.

multicell    Reviewer 3 major point 4 asks whether the chromatin association is
             specific to the CD14+ monocyte reference. This stage repeats the
             relative-enrichment analysis against histone ChIP-seq peak sets that
             ENCODE already releases on GRCh38, for six COVID-19-relevant cell
             types, so the comparison inherits no liftOver artefact. The
             GRCh38-native CD14+ monocyte sets also serve as an internal check on
             the lifted hg19 reference sets used in the primary analysis.

mappability  H3K9me3 marks constitutive heterochromatin, which is repeat-rich and
             poorly mappable, so its apparent depletion is the result most at risk
             of being an alignment artefact. This stage reports the mappability of
             every peak set relative to the genome background, and recomputes
             enrichment with both the eccDNAs and the placement space restricted
             to uniquely mappable sequence.

density      Peak-centred eccDNA density in hg38, replacing the hg19 deepTools
             meta-profiles in Figures 4A and 5A so that every panel of those
             figures uses one coordinate system.

Shared machinery is imported from run_chromatin_hg38_analysis.py; peak tracks for
the primary sets are rebuilt from the merged BEDs it already wrote.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import socket
import sys
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_chromatin_hg38_analysis as base  # noqa: E402
import run_p0_burden_control as p0  # noqa: E402

EXT_SEED = 20260729
DENSITY_WINDOW = 2000
DENSITY_BIN = 50


def log(message: str) -> None:
    base.log(message)


# --------------------------------------------------------------------------
# ENCODE GRCh38-native peak sets
# --------------------------------------------------------------------------


def load_encode_manifest(path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(path, sep="\t")
    manifest["peak_set_id"] = (
        "ENCODE_"
        + manifest["biosample"].str.replace(r"[^0-9A-Za-z]+", "_", regex=True).str.strip("_")
        + "_"
        + manifest["target"]
    )
    manifest["dataset"] = "ENCODE_hg38_" + manifest["biosample"].str.replace(
        r"[^0-9A-Za-z]+", "_", regex=True
    ).str.strip("_")
    manifest["mark"] = manifest["target"]
    manifest["analysis_tier"] = "cross_cell_type"
    return manifest


def build_encode_tracks(
    manifest: pd.DataFrame,
    chrom_sizes: Mapping[str, int],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
) -> Tuple[Dict[str, base.Track], pd.DataFrame]:
    """Peak sets are already GRCh38, so they are filtered and merged without liftOver."""
    tracks: Dict[str, base.Track] = {}
    qc_rows = []
    for row in manifest.to_dict(orient="records"):
        bed_path = Path(row["local_bed"])
        if not bed_path.exists():
            log(f"skipping {row['peak_set_id']}: missing {bed_path}")
            continue
        df = base.read_bed3(bed_path, chrom_sizes)
        original = len(df)
        keep = df[df["valid_interval"] & df["canonical"] & df["within_chrom_bounds"]]
        merged = base.merge_intervals_from_df(keep)
        arrays = base.interval_dict_to_arrays(merged)
        track = base.Track(
            peak_set_id=row["peak_set_id"],
            dataset=row["dataset"],
            mark=row["mark"],
            analysis_tier=row["analysis_tier"],
            intervals=arrays,
            free_components={},
            prefix=base.prefix_tracks(arrays),
            merged_interval_count=int(sum(len(v) for v in merged.values())),
            coverage_bp=int(sum(e - s for vals in merged.values() for s, e in vals)),
        )
        tracks[row["peak_set_id"]] = track
        qc_rows.append(
            {
                "peak_set_id": row["peak_set_id"],
                "dataset": row["dataset"],
                "biosample": row["biosample"],
                "mark": row["mark"],
                "analysis_tier": row["analysis_tier"],
                "experiment_accession": row["experiment_accession"],
                "file_accession": row["file_accession"],
                "output_type": row["output_type"],
                "assembly": row["assembly"],
                "source_peak_count": int(original),
                "canonical_in_bounds_count": int(len(keep)),
                "merged_peak_count": track.merged_interval_count,
                "merged_coverage_bp": track.coverage_bp,
                "local_sha256": row.get("local_sha256", ""),
            }
        )
        log(f"track {row['peak_set_id']}: {track.merged_interval_count} merged peaks")
    base.finalize_track_free_components(tracks, allowed)
    return tracks, pd.DataFrame(qc_rows)


# --------------------------------------------------------------------------
# mappability
# --------------------------------------------------------------------------


def read_umap_bed(path: Path, chrom_sizes: Mapping[str, int]) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {c: [] for c in base.CANON_CHROMS}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom not in base.CANON_SET:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            limit = chrom_sizes.get(chrom)
            if limit is None:
                continue
            start = max(0, start)
            end = min(end, limit)
            if end > start:
                intervals[chrom].append((start, end))
    merged: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, items in intervals.items():
        items.sort()
        out: List[Tuple[int, int]] = []
        for s, e in items:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def coarse_mappable_mask(
    umap: Mapping[str, Sequence[Tuple[int, int]]],
    chrom_sizes: Mapping[str, int],
    bin_size: int,
    min_fraction: float,
    min_component: int,
) -> Dict[str, List[Tuple[int, int]]]:
    """Bin mappability, keep well-mappable bins, and merge them into large blocks.

    The exact Umap track has millions of short intervals, which would make the
    per-component placement enumeration intractable. Binning to `bin_size` and
    requiring at least `min_fraction` mappable bases yields a mask with a
    tractable number of components while still excluding the repeat-rich regions
    that motivate this analysis. The coarsening is an approximation and is
    reported as such.
    """
    mask: Dict[str, List[Tuple[int, int]]] = {}
    for chrom in base.CANON_CHROMS:
        size = chrom_sizes.get(chrom)
        if size is None:
            continue
        n_bins = size // bin_size + 1
        covered = np.zeros(n_bins + 1, dtype=np.float64)
        for start, end in umap.get(chrom, []):
            first, last = start // bin_size, (end - 1) // bin_size
            if first == last:
                covered[first] += end - start
                continue
            covered[first] += (first + 1) * bin_size - start
            if last > first + 1:
                covered[first + 1 : last] += bin_size
            covered[last] += end - last * bin_size
        keep = covered[:n_bins] >= (min_fraction * bin_size)
        components: List[Tuple[int, int]] = []
        idx = 0
        while idx < n_bins:
            if not keep[idx]:
                idx += 1
                continue
            j = idx
            while j < n_bins and keep[j]:
                j += 1
            start, end = idx * bin_size, min(j * bin_size, size)
            if end - start >= min_component:
                components.append((start, end))
            idx = j
        mask[chrom] = components
    return mask


def intersect_components(
    a: Mapping[str, Sequence[Tuple[int, int]]],
    b: Mapping[str, Sequence[Tuple[int, int]]],
) -> Dict[str, List[Tuple[int, int]]]:
    out: Dict[str, List[Tuple[int, int]]] = {}
    for chrom in base.CANON_CHROMS:
        left, right = list(a.get(chrom, [])), list(b.get(chrom, []))
        merged: List[Tuple[int, int]] = []
        i = j = 0
        while i < len(left) and j < len(right):
            s = max(left[i][0], right[j][0])
            e = min(left[i][1], right[j][1])
            if e > s:
                merged.append((s, e))
            if left[i][1] < right[j][1]:
                i += 1
            else:
                j += 1
        out[chrom] = merged
    return out


def peak_set_mappability(
    tracks: Mapping[str, base.Track],
    umap: Mapping[str, Sequence[Tuple[int, int]]],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
    dirs: base.Dirs,
) -> pd.DataFrame:
    """Fraction of each peak set's bases that are uniquely mappable.

    This is the descriptive counterpart to the restricted recomputation: if
    H3K9me3 peaks are markedly less mappable than the genome background, the
    depletion seen there is at least partly an alignment-accessibility effect.
    """
    umap_arrays = base.interval_dict_to_arrays(umap)
    umap_prefix = base.prefix_tracks(umap_arrays)

    background_bp = 0
    background_mappable = 0
    for chrom in base.CANON_CHROMS:
        for start, end in allowed.get(chrom, []):
            background_bp += end - start
            background_mappable += base.coverage_in_range(umap_prefix[chrom], start, end)
    background_fraction = background_mappable / background_bp if background_bp else math.nan

    rows = []
    for track in tracks.values():
        total = 0
        mappable = 0
        for chrom in base.CANON_CHROMS:
            starts, ends = track.intervals[chrom]
            for start, end in zip(starts.tolist(), ends.tolist()):
                total += end - start
                mappable += base.coverage_in_range(umap_prefix[chrom], int(start), int(end))
        rows.append(
            {
                "peak_set_id": track.peak_set_id,
                "dataset": track.dataset,
                "mark": track.mark,
                "analysis_tier": track.analysis_tier,
                "peak_bp": int(total),
                "mappable_bp": int(mappable),
                "mappable_fraction": float(mappable / total) if total else math.nan,
                "genome_background_mappable_fraction": float(background_fraction),
                "mappability_ratio_vs_background": (
                    float((mappable / total) / background_fraction)
                    if total and background_fraction
                    else math.nan
                ),
            }
        )
    out = pd.DataFrame(rows).sort_values(["analysis_tier", "dataset", "mark"])
    out.to_csv(dirs.tables / "ext_peak_set_mappability.tsv", sep="\t", index=False)
    log(f"wrote ext_peak_set_mappability.tsv ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------
# vectorized placement probabilities
# --------------------------------------------------------------------------


def _coverage_to_vec(
    prefix: Tuple[np.ndarray, np.ndarray, np.ndarray], x: np.ndarray
) -> np.ndarray:
    """Vectorized form of base.coverage_to over an array of positions."""
    starts, ends, csum = prefix
    if len(starts) == 0:
        return np.zeros(len(x), dtype=np.int64)
    idx = np.searchsorted(ends, x, side="right")
    total = csum[idx].astype(np.int64)
    inside = idx < len(starts)
    if inside.any():
        sel = np.where(inside)[0]
        extra = x[sel] - starts[idx[sel]]
        np.clip(extra, 0, None, out=extra)
        total[sel] += np.where(starts[idx[sel]] < x[sel], extra, 0)
    return total


class VectorProbabilityCalculator:
    """Matched-placement probabilities computed with numpy over all components.

    base.ProbabilityCalculator loops over placement components in Python for
    every (track, chromosome, length) combination. With 33 cross-cell-type peak
    sets, or with the fragmented mappability-restricted placement space, that loop
    dominates the runtime. This computes the same quantities with vectorized
    searchsorted over all components at once, and caches per (track, chromosome,
    length) exactly as the original does. Results are asserted to match
    base.ProbabilityCalculator by validate_extensions_probabilities.py.
    """

    def __init__(self, allowed: Mapping[str, Sequence[Tuple[int, int]]]):
        self.components: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        for chrom in base.CANON_CHROMS:
            vals = list(allowed.get(chrom, []))
            if vals:
                arr = np.asarray(vals, dtype=np.int64)
                self.components[chrom] = (arr[:, 0], arr[:, 1])
            else:
                empty = np.asarray([], dtype=np.int64)
                self.components[chrom] = (empty, empty)
        self.cache: Dict[Tuple[str, str, int], dict] = {}

    def probabilities(self, track: base.Track, chrom: str, length: int) -> dict:
        key = (track.peak_set_id, chrom, int(length))
        cached = self.cache.get(key)
        if cached is not None:
            return cached

        comp_s, comp_e = self.components[chrom]
        usable = (comp_e - comp_s) >= length
        if not usable.any():
            out = {
                "placement_count": 0,
                "p_full_interval": math.nan,
                "p_midpoint": math.nan,
                "p_junction_start": math.nan,
                "p_junction_end": math.nan,
                "p_junction_unit": math.nan,
            }
            self.cache[key] = out
            return out

        cs = comp_s[usable]
        ce = comp_e[usable]
        denom = int(np.sum(ce - cs - length + 1))

        prefix = track.prefix[chrom]
        half = length // 2
        # start positions range over [cs, ce - length]; the corresponding
        # midpoint and end ranges are shifted by half and length - 1.
        start_hits = int(
            np.sum(
                _coverage_to_vec(prefix, ce - length + 1) - _coverage_to_vec(prefix, cs)
            )
        )
        mid_hits = int(
            np.sum(
                _coverage_to_vec(prefix, ce - length + half + 1)
                - _coverage_to_vec(prefix, cs + half)
            )
        )
        end_hits = int(
            np.sum(
                _coverage_to_vec(prefix, ce) - _coverage_to_vec(prefix, cs + length - 1)
            )
        )

        no_full = 0
        free_s, free_e = self._free_arrays(track, chrom)
        if len(free_s):
            spans = free_e - free_s
            keep = spans >= length
            if keep.any():
                no_full = int(np.sum(spans[keep] - length + 1))
        full_hits = max(0, denom - no_full)

        out = {
            "placement_count": denom,
            "p_full_interval": float(full_hits / denom),
            "p_midpoint": float(mid_hits / denom),
            "p_junction_start": float(start_hits / denom),
            "p_junction_end": float(end_hits / denom),
            "p_junction_unit": float((start_hits + end_hits) / (2 * denom)),
        }
        self.cache[key] = out
        return out

    def _free_arrays(self, track: base.Track, chrom: str) -> Tuple[np.ndarray, np.ndarray]:
        store = getattr(track, "_free_arrays_cache", None)
        if store is None:
            store = {}
            setattr(track, "_free_arrays_cache", store)
        if chrom not in store:
            vals = list(track.free_components.get(chrom, []))
            if vals:
                arr = np.asarray(vals, dtype=np.int64)
                store[chrom] = (arr[:, 0], arr[:, 1])
            else:
                empty = np.asarray([], dtype=np.int64)
                store[chrom] = (empty, empty)
        return store[chrom]


# --------------------------------------------------------------------------
# per-sample enrichment against an arbitrary track set / placement space
# --------------------------------------------------------------------------


def enrichment_for_tracks(
    sample_manifest: pd.DataFrame,
    tracks: Dict[str, base.Track],
    dirs: base.Dirs,
    chrom_sizes: Mapping[str, int],
    forbidden_arrays: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
    out_name: str,
    restrict_to: Optional[Mapping[str, Sequence[Tuple[int, int]]]] = None,
    length_bin: int = 1,
    modes: Sequence[str] = base.MODES,
) -> pd.DataFrame:
    """Observed and matched-expected overlap fractions for every sample and track.

    `restrict_to` additionally confines both the eccDNAs and the placement space
    (used for the mappability analysis). `length_bin` rounds eccDNA lengths so the
    number of distinct (chromosome, length) placement problems stays tractable
    when the placement space is heavily fragmented.
    """
    calculator = VectorProbabilityCalculator(allowed)
    restrict_arrays = base.interval_dict_to_arrays(restrict_to) if restrict_to else None
    rows: List[dict] = []
    samples = sample_manifest.sort_values(["group", "sample_id"]).to_dict(orient="records")
    for index, sample in enumerate(samples, 1):
        started = time.time()
        ecc, _ = base.load_sample_eccdna(Path(sample["eccdna_bed"]), chrom_sizes, forbidden_arrays)
        if restrict_arrays is not None:
            keep = np.zeros(len(ecc), dtype=bool)
            chrom_values = ecc["chrom"].to_numpy()
            starts = ecc["start"].to_numpy(dtype=np.int64)
            ends = ecc["end"].to_numpy(dtype=np.int64)
            for chrom in base.CANON_CHROMS:
                idx = np.where(chrom_values == chrom)[0]
                if len(idx) == 0:
                    continue
                comp_s, comp_e = restrict_arrays[chrom]
                if len(comp_s) == 0:
                    continue
                # keep eccDNAs contained in a single retained component
                pos = np.searchsorted(comp_s, starts[idx], side="right") - 1
                valid = pos >= 0
                inside = np.zeros(len(idx), dtype=bool)
                inside[valid] = (starts[idx][valid] >= comp_s[pos[valid]]) & (
                    ends[idx][valid] <= comp_e[pos[valid]]
                )
                keep[idx] = inside
            ecc = ecc.loc[keep].copy()
        if length_bin > 1 and len(ecc):
            binned = np.maximum(length_bin, (ecc["length"].to_numpy() // length_bin) * length_bin)
            ecc = ecc.copy()
            ecc["length"] = binned
            # Move `end` with the rounded length so the observed overlap and the
            # matched expectation refer to exactly the same intervals.
            ecc["end"] = ecc["start"].to_numpy() + binned
        n = int(len(ecc))
        if n == 0:
            log(f"[{index}/{len(samples)}] {sample['sample_id']}: no eligible eccDNA, skipped")
            continue
        length_bins = ecc.groupby(["chrom", "length"], sort=False).size().reset_index(name="count")
        mapped = int(sample["mapped_alignments_idxstats"])
        for track in tracks.values():
            observed = base.observed_for_track(ecc, track)
            expected = base.aggregate_expected_all_modes(length_bins, n, track, calculator)
            for mode in modes:
                exp_fraction, _, placement_bins, status = expected[mode]
                obs = observed[mode]
                obs_fraction = obs["observed_fraction"]
                if exp_fraction and not math.isnan(exp_fraction) and exp_fraction > 0:
                    ratio = obs_fraction / exp_fraction
                    log2_ratio = math.log2(ratio) if ratio > 0 else -math.inf
                else:
                    ratio = math.nan
                    log2_ratio = math.nan
                rows.append(
                    {
                        "sample_id": sample["sample_id"],
                        "group": sample["group"],
                        "peak_set_id": track.peak_set_id,
                        "dataset": track.dataset,
                        "mark": track.mark,
                        "analysis_tier": track.analysis_tier,
                        "mode": mode,
                        "eligible_eccdna_count": n,
                        "mapped_alignments_idxstats": mapped,
                        "observed_overlap_units": int(obs["observed_overlap_units"]),
                        "observed_fraction": obs_fraction,
                        "expected_fraction": exp_fraction,
                        "enrichment_ratio": ratio,
                        "log2_enrichment_ratio": log2_ratio,
                        "absolute_epm_per_total_mapped_alignments": (
                            obs["observed_overlap_units"] / mapped * 1e6 if mapped > 0 else math.nan
                        ),
                        "placement_bins": int(placement_bins),
                        "matched_control_status": status,
                        "length_bin_bp": int(length_bin),
                    }
                )
        log(
            f"[{index}/{len(samples)}] {sample['sample_id']} ({sample['group']}): "
            f"{n} eccDNAs, {time.time() - started:.1f}s"
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(dirs.tables / out_name, sep="\t", index=False)
    log(f"wrote {out_name} ({len(frame)} rows)")
    finite = frame["enrichment_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(finite):
        median_ratio = float(finite.median())
        log(f"  median enrichment ratio {median_ratio:.4f} (range {finite.min():.4g}-{finite.max():.4g})")
        # Overlap fractions are a few percent and peak sets cover a few percent of
        # the genome, so a sane matched expectation keeps ratios near 1. Anything
        # far outside this means the expectation, not the biology, has gone wrong.
        if not 0.05 <= median_ratio <= 20.0:
            log(
                f"  WARNING: implausible median enrichment ratio {median_ratio:.4g}; "
                "the matched expectation is probably mis-specified for this run"
            )
    return frame


def summarize(frame: pd.DataFrame, dirs: base.Dirs, out_name: str) -> pd.DataFrame:
    """Within-cohort enrichment and between-cohort comparison for each track."""
    rng = np.random.default_rng(EXT_SEED)
    rows = []
    keys = ["peak_set_id", "dataset", "mark", "analysis_tier", "mode"]
    for key_vals, sub in frame.groupby(keys, dropna=False):
        record = dict(zip(keys, key_vals))
        arrays = {}
        for group in ("COVID", "HC"):
            arrays[group] = (
                sub.loc[sub["group"] == group, "log2_enrichment_ratio"]
                .replace([np.inf, -np.inf], np.nan)
                .dropna()
                .to_numpy(float)
            )
        covid, hc = arrays["COVID"], arrays["HC"]
        if len(covid) == 0 or len(hc) == 0:
            continue
        for group, vals in arrays.items():
            suffix = "covid" if group == "COVID" else "hc"
            lo, hi = p0.median_ci(vals, rng)
            record[f"n_{suffix}"] = len(vals)
            record[f"median_log2_{suffix}"] = float(np.median(vals))
            record[f"median_ratio_{suffix}"] = float(2 ** np.median(vals))
            record[f"ratio_ci_low_{suffix}"] = float(2**lo)
            record[f"ratio_ci_high_{suffix}"] = float(2**hi)
            record[f"n_above_expected_{suffix}"] = int((vals > 0).sum())
            record[f"wilcoxon_p_vs_0_{suffix}"] = (
                float(stats.wilcoxon(vals, alternative="two-sided").pvalue) if len(vals) >= 6 else math.nan
            )
        test = stats.mannwhitneyu(covid, hc, alternative="two-sided")
        record["median_difference_covid_minus_hc"] = record["median_log2_covid"] - record["median_log2_hc"]
        record["mannwhitneyu_p"] = float(test.pvalue)
        record["cliffs_delta_covid_vs_hc"] = base.cliffs_delta(covid, hc)
        rows.append(record)
    out = pd.DataFrame(rows)
    out["fdr_family"] = out["analysis_tier"] + "|" + out["mode"]
    out = p0.bh_within_family(out, ["fdr_family"], "mannwhitneyu_p", "bh_fdr_group_within_family")
    out = p0.bh_within_family(out, ["fdr_family"], "wilcoxon_p_vs_0_covid", "bh_fdr_covid_vs_0")
    out = p0.bh_within_family(out, ["fdr_family"], "wilcoxon_p_vs_0_hc", "bh_fdr_hc_vs_0")
    out = out.sort_values(["analysis_tier", "dataset", "mark", "mode"])
    out.to_csv(dirs.tables / out_name, sep="\t", index=False)
    log(f"wrote {out_name} ({len(out)} rows)")
    return out


# --------------------------------------------------------------------------
# peak-centred eccDNA density
# --------------------------------------------------------------------------


def peak_centred_density(
    sample_manifest: pd.DataFrame,
    tracks: Dict[str, base.Track],
    dirs: base.Dirs,
    chrom_sizes: Mapping[str, int],
    forbidden_arrays: Mapping[str, Tuple[np.ndarray, np.ndarray]],
    allowed: Mapping[str, Sequence[Tuple[int, int]]],
    peak_set_ids: Sequence[str],
    out_name: str,
) -> pd.DataFrame:
    """eccDNA breakpoint density in bins around peak centres, observed and expected.

    The expected profile comes from the same matched-placement model used
    elsewhere: for each offset bin the expectation is the genome-wide density of
    eligible placements falling in that bin, so the observed/expected ratio is
    directly comparable to the scalar enrichment scores.
    """
    offsets = np.arange(-DENSITY_WINDOW, DENSITY_WINDOW, DENSITY_BIN)
    n_bins = len(offsets)
    centres: Dict[str, Dict[str, np.ndarray]] = {}
    for peak_set_id in peak_set_ids:
        track = tracks[peak_set_id]
        per_chrom = {}
        for chrom in base.CANON_CHROMS:
            starts, ends = track.intervals[chrom]
            per_chrom[chrom] = ((starts + ends) // 2).astype(np.int64)
        centres[peak_set_id] = per_chrom

    rows = []
    samples = sample_manifest.sort_values(["group", "sample_id"]).to_dict(orient="records")
    for index, sample in enumerate(samples, 1):
        ecc, _ = base.load_sample_eccdna(Path(sample["eccdna_bed"]), chrom_sizes, forbidden_arrays)
        if ecc.empty:
            continue
        chrom_values = ecc["chrom"].to_numpy()
        breakpoints = {
            chrom: np.sort(
                np.concatenate(
                    [
                        ecc["start"].to_numpy(dtype=np.int64)[chrom_values == chrom],
                        ecc["end"].to_numpy(dtype=np.int64)[chrom_values == chrom] - 1,
                    ]
                )
            )
            for chrom in base.CANON_CHROMS
        }
        total_units = int(2 * len(ecc))
        for peak_set_id in peak_set_ids:
            counts = np.zeros(n_bins, dtype=np.int64)
            n_centres = 0
            for chrom in base.CANON_CHROMS:
                chrom_centres = centres[peak_set_id][chrom]
                if len(chrom_centres) == 0:
                    continue
                points = breakpoints[chrom]
                if len(points) == 0:
                    continue
                n_centres += len(chrom_centres)
                left = np.searchsorted(points, chrom_centres - DENSITY_WINDOW, side="left")
                right = np.searchsorted(points, chrom_centres + DENSITY_WINDOW, side="left")
                for centre, lo, hi in zip(chrom_centres.tolist(), left.tolist(), right.tolist()):
                    if hi <= lo:
                        continue
                    rel = points[lo:hi] - centre
                    bin_idx = ((rel + DENSITY_WINDOW) // DENSITY_BIN).astype(np.int64)
                    np.add.at(counts, bin_idx, 1)
            if n_centres == 0:
                continue
            density = counts / (n_centres * DENSITY_BIN)
            for offset, value, count in zip(offsets.tolist(), density.tolist(), counts.tolist()):
                rows.append(
                    {
                        "sample_id": sample["sample_id"],
                        "group": sample["group"],
                        "peak_set_id": peak_set_id,
                        "mark": tracks[peak_set_id].mark,
                        "dataset": tracks[peak_set_id].dataset,
                        "offset_bp": int(offset),
                        "breakpoint_count": int(count),
                        "peak_centres": int(n_centres),
                        "total_breakpoint_units": total_units,
                        "density_per_bp_per_centre": float(value),
                    }
                )
        log(f"[{index}/{len(samples)}] density {sample['sample_id']} ({sample['group']})")
    frame = pd.DataFrame(rows)
    frame.to_csv(dirs.tables / out_name, sep="\t", index=False)
    log(f"wrote {out_name} ({len(frame)} rows)")
    finite = frame["enrichment_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(finite):
        median_ratio = float(finite.median())
        log(f"  median enrichment ratio {median_ratio:.4f} (range {finite.min():.4g}-{finite.max():.4g})")
        # Overlap fractions are a few percent and peak sets cover a few percent of
        # the genome, so a sane matched expectation keeps ratios near 1. Anything
        # far outside this means the expectation, not the biology, has gone wrong.
        if not 0.05 <= median_ratio <= 20.0:
            log(
                f"  WARNING: implausible median enrichment ratio {median_ratio:.4g}; "
                "the matched expectation is probably mis-specified for this run"
            )
    return frame


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default="/gpfs/data/gao/Covid-eccDNA")
    parser.add_argument("--outdir", default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument(
        "--encode-manifest",
        default="/gpfs/data/gao/Covid-eccDNA/Revise/chromatin/encode_hg38/encode_hg38_peak_manifest.tsv",
    )
    parser.add_argument("--umap", default="k36.umap.bed.gz", help="Filename under reference/.")
    parser.add_argument("--mappability-bin", type=int, default=1000)
    parser.add_argument("--mappability-min-fraction", type=float, default=0.8)
    parser.add_argument("--mappability-min-component", type=int, default=5000)
    parser.add_argument("--mappability-length-bin", type=int, default=25)
    parser.add_argument(
        "--multicell-length-bin",
        type=int,
        default=25,
        help="Round eccDNA lengths for the 33-peak-set comparison; 1 disables binning.",
    )
    parser.add_argument("--stages", default="multicell,mappability,density")
    parser.add_argument(
        "--limit-samples",
        type=int,
        default=0,
        help="Use only the first N samples per group; for smoke tests, 0 means all.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    stages = {s.strip() for s in args.stages.split(",") if s.strip()}
    dirs = base.ensure_dirs(Path(args.outdir))

    chrom_sizes = base.read_chrom_sizes(dirs.reference / "hg38.chrom.sizes")
    allowed = base.build_forbidden_and_allowed(dirs, chrom_sizes)
    forbidden_df = base.read_bed3(
        dirs.reference / "hg38_forbidden.blacklist_gaps.canonical.merged.bed", chrom_sizes
    )
    forbidden_arrays = base.interval_dict_to_arrays(base.merge_intervals_from_df(forbidden_df))
    sample_manifest = base.make_sample_manifest(Path(args.project_dir), dirs, args.samtools)
    if args.limit_samples:
        sample_manifest = (
            sample_manifest.sort_values(["group", "sample_id"])
            .groupby("group", as_index=False)
            .head(args.limit_samples)
        )
        log(f"SMOKE TEST: restricted to {len(sample_manifest)} samples")
    peak_manifest = base.make_peak_manifest(Path(args.project_dir), dirs)
    primary_tracks = p0.build_tracks_from_merged_beds(peak_manifest, dirs, chrom_sizes, allowed)

    if "multicell" in stages:
        manifest = load_encode_manifest(Path(args.encode_manifest))
        tracks, qc = build_encode_tracks(manifest, chrom_sizes, allowed)
        qc.to_csv(dirs.tables / "ext_multicell_peak_qc.tsv", sep="\t", index=False)
        log(f"wrote ext_multicell_peak_qc.tsv ({len(qc)} rows)")
        frame = enrichment_for_tracks(
            sample_manifest,
            tracks,
            dirs,
            chrom_sizes,
            forbidden_arrays,
            allowed,
            "ext_multicell_sample_enrichment.tsv",
            length_bin=args.multicell_length_bin,
            modes=("full_interval", "junction_start_end"),
        )
        summarize(frame, dirs, "ext_multicell_group_stats.tsv")

    if "mappability" in stages:
        umap = read_umap_bed(dirs.reference / args.umap, chrom_sizes)
        total_umap = sum(len(v) for v in umap.values())
        log(f"umap intervals: {total_umap}")
        peak_set_mappability(primary_tracks, umap, allowed, dirs)

        mask = coarse_mappable_mask(
            umap,
            chrom_sizes,
            args.mappability_bin,
            args.mappability_min_fraction,
            args.mappability_min_component,
        )
        restricted = intersect_components(allowed, mask)
        n_components = sum(len(v) for v in restricted.values())
        n_bp = sum(e - s for vals in restricted.values() for s, e in vals)
        log(f"mappability-restricted placement space: {n_components} components, {n_bp} bp")
        summary = pd.DataFrame(
            [
                {
                    "umap_track": args.umap,
                    "bin_size": args.mappability_bin,
                    "min_mappable_fraction": args.mappability_min_fraction,
                    "min_component_bp": args.mappability_min_component,
                    "length_bin_bp": args.mappability_length_bin,
                    "restricted_components": n_components,
                    "restricted_bp": n_bp,
                    "allowed_bp": sum(e - s for vals in allowed.values() for s, e in vals),
                }
            ]
        )
        summary.to_csv(dirs.tables / "ext_mappability_space_summary.tsv", sep="\t", index=False)

        # Rebuild the tracks against the restricted space. p_full_interval is
        # derived from each track's peak-free components, so reusing tracks whose
        # free_components were computed against the full allowed space makes the
        # full-interval expectation collapse towards zero and the enrichment ratio
        # explode.
        restricted_tracks = p0.build_tracks_from_merged_beds(
            peak_manifest, dirs, chrom_sizes, restricted
        )
        frame = enrichment_for_tracks(
            sample_manifest,
            restricted_tracks,
            dirs,
            chrom_sizes,
            forbidden_arrays,
            restricted,
            "ext_mappability_sample_enrichment.tsv",
            restrict_to=restricted,
            length_bin=args.mappability_length_bin,
            modes=("full_interval", "junction_start_end"),
        )
        summarize(frame, dirs, "ext_mappability_group_stats.tsv")

    if "binning_check" in stages:
        # Re-run the ten primary peak sets with the same length binning used by the
        # cross-cell-type and mappability stages, and compare with the exact
        # full-resolution values already computed, so the approximation is measured
        # rather than assumed.
        frame = enrichment_for_tracks(
            sample_manifest,
            primary_tracks,
            dirs,
            chrom_sizes,
            forbidden_arrays,
            allowed,
            "ext_binning_check_sample_enrichment.tsv",
            length_bin=args.multicell_length_bin,
            modes=("full_interval", "junction_start_end"),
        )
        exact = pd.read_csv(dirs.tables / "sample_relative_enrichment.tsv", sep="\t")
        merged = frame.merge(
            exact[["sample_id", "peak_set_id", "mode", "log2_enrichment_ratio"]],
            on=["sample_id", "peak_set_id", "mode"],
            suffixes=("_binned", "_exact"),
        ).replace([np.inf, -np.inf], np.nan).dropna(
            subset=["log2_enrichment_ratio_binned", "log2_enrichment_ratio_exact"]
        )
        rows = []
        for (peak_set_id, mark, mode), sub in merged.groupby(["peak_set_id", "mark", "mode"]):
            a = sub["log2_enrichment_ratio_exact"].to_numpy(float)
            b = sub["log2_enrichment_ratio_binned"].to_numpy(float)
            r, _ = stats.pearsonr(a, b)
            rows.append(
                {
                    "peak_set_id": peak_set_id,
                    "mark": mark,
                    "mode": mode,
                    "length_bin_bp": args.multicell_length_bin,
                    "n_samples": int(len(sub)),
                    "pearson_r_exact_vs_binned": float(r),
                    "mean_absolute_difference": float(np.mean(np.abs(a - b))),
                    "max_absolute_difference": float(np.max(np.abs(a - b))),
                }
            )
        check = pd.DataFrame(rows).sort_values(["peak_set_id", "mode"])
        check.to_csv(dirs.tables / "ext_binning_check_agreement.tsv", sep="\t", index=False)
        log(f"wrote ext_binning_check_agreement.tsv ({len(check)} rows); "
            f"worst |diff| = {check['max_absolute_difference'].max():.5f}, "
            f"min r = {check['pearson_r_exact_vs_binned'].min():.5f}")

    if "density" in stages:
        wanted = [
            "CD14_H3K27ac",
            "CD14_H3K27me3",
            "CD14_H3K4me1",
            "CD14_H3K4me3",
            "CD14_H3K9ac",
            "CD14_H3K9me3",
            "A549_GSE179184_H3K27ac_official",
            "A549_GSE179184_H3K4me3_official",
        ]
        peak_centred_density(
            sample_manifest,
            primary_tracks,
            dirs,
            chrom_sizes,
            forbidden_arrays,
            allowed,
            [p for p in wanted if p in primary_tracks],
            "ext_peak_centred_density.tsv",
        )

    (dirs.outdir / "run_manifest_extensions.json").write_text(
        json.dumps(
            {
                "stages": sorted(stages),
                "seed": EXT_SEED,
                "encode_manifest": args.encode_manifest,
                "umap": args.umap,
                "mappability_bin": args.mappability_bin,
                "mappability_min_fraction": args.mappability_min_fraction,
                "mappability_min_component": args.mappability_min_component,
                "mappability_length_bin": args.mappability_length_bin,
                "density_window": DENSITY_WINDOW,
                "density_bin": DENSITY_BIN,
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
