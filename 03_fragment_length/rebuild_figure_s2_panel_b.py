#!/usr/bin/env python3
"""Rebuild Supplementary Figure S2 so that panel b plots the stability criterion.

Panel b previously plotted `subsample_detection_frequency`, the 100 equal-count
repeats, which the Methods list as a sensitivity analysis. The prespecified
classification of a stable consensus peak is made on
`bootstrap_detection_frequency` instead (>= 90% overall and >= 80% in each
cohort; see the classification block in `run_fragment_length_peak_analysis.py`).
The two disagree for the 760-bp candidate: 90.0% under equal-count subsampling
in HC but 66.8% under the bootstrap, so the panel appeared to clear the 80%
threshold that the Results correctly report it as failing.

This script redraws the figure from the deposited fragment-length tables with
panel b showing the bootstrap frequencies. Only the plotted series and the axis
label change; no threshold rules or marks are added, because the thresholds
belong in the legend and in Supplementary Table S5b rather than on the panel.
Panels a and c are untouched. Before it writes anything it rebuilds the original panel b
and asserts that the result reproduces the submitted figure, following the
`rerun_figure_s8_locally.py` pattern, so a silent style drift cannot slip in.

Run from the repository root.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
STAGE = ROOT / "deposit" / "scripts" / "02_fragment_length"
TABLES = ROOT / "deposit" / "tables" / "02_fragment_length"
SUBMITTED = ROOT / "Submit" / "Figure_S2.pdf"
DEPOSITED = ROOT / "deposit" / "figures" / "Figure_S2.pdf"

sys.path.insert(0, str(STAGE))
import run_fragment_length_peak_analysis as fl  # noqa: E402



def read_tsv(name: str) -> list[dict[str, str]]:
    with (TABLES / name).open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_inputs():
    by_sigma: dict[float, list[tuple[int, float]]] = defaultdict(list)
    for row in read_tsv("bandwidth_sensitivity_curves.tsv"):
        by_sigma[float(row["gaussian_sigma_bp"])].append(
            (int(row["length_bp"]), float(row["overall_mean_density"])))
    sigma_curves, grid = {}, None
    for sigma, points in by_sigma.items():
        points.sort()
        lengths = np.array([p[0] for p in points], dtype=float)
        if grid is None:
            grid = lengths
        elif not np.array_equal(grid, lengths):
            raise SystemExit("bandwidth curves are not on a common length grid")
        sigma_curves[sigma] = np.array([p[1] for p in points], dtype=float)

    sigma_calls: dict[float, list[dict[str, object]]] = {s: [] for s in sigma_curves}
    for row in read_tsv("automatic_peak_calls.tsv"):
        if row["analysis"] != "sample_weighted_bandwidth" or row["dataset"] != "Overall":
            continue
        sigma_calls[float(row["gaussian_sigma_bp"])].append(
            {"position_bp": int(row["position_bp"])})
    for calls in sigma_calls.values():
        calls.sort(key=lambda c: c["position_bp"])

    candidates = read_tsv("objective_peak_identification_and_stability.tsv")
    candidates.sort(key=lambda c: int(str(c["peak_id"]).lstrip("P")))

    detected: dict[tuple[str, str], list[int]] = defaultdict(list)
    for row in read_tsv("equal_count_subsample_peak_detection.tsv"):
        detected[(row["dataset"], row["peak_id"])].append(
            int(row["detected_within_match_window"]))
    subsample = {k: {"subsample_detection_frequency": sum(v) / len(v)}
                 for k, v in detected.items()}

    bootstrap = {}
    for candidate in candidates:
        for dataset, prefix in (("Overall", "overall"), (fl.GROUP_COVID, "covid"),
                                (fl.GROUP_HC, "hc")):
            bootstrap[(dataset, str(candidate["peak_id"]))] = {
                "subsample_detection_frequency":
                    float(candidate[f"{prefix}_bootstrap_detection_frequency"])}

    sample_rows = read_tsv("peak_window_sample_proportions.tsv")
    for row in sample_rows:
        row["window_proportion"] = float(row["window_proportion"])

    return dict(grid=grid, sigma_curves=sigma_curves, sigma_calls=sigma_calls,
                candidates=candidates, subsample=subsample, bootstrap=bootstrap,
                sample_rows=sample_rows,
                comparison_rows=read_tsv("peak_window_group_comparisons.tsv"))


def draw(base: Path, data, summary, annotate: bool) -> None:
    """Call the deposited figure builder, optionally annotating panel b."""
    captured = {}
    real_save = fl.save_figure

    def capture(fig, base_path, dpi=600):
        captured["fig"] = fig
        captured["base"] = base_path
        captured["dpi"] = dpi

    fl.save_figure = capture
    try:
        fl.make_supplementary_figure(
            base, data["grid"], data["sigma_curves"], data["sigma_calls"],
            data["candidates"], summary, data["sample_rows"], data["comparison_rows"])
    finally:
        fl.save_figure = real_save

    fig = captured["fig"]
    if annotate:
        # Only the axis label changes: the bars now carry a different quantity from
        # the one the deposited builder plotted. No threshold rules, no marks — the
        # thresholds belong in the legend and in Table S5b, not on the panel.
        ax_b = next(ax for ax in fig.axes
                    if ax.get_ylabel() == "Detection frequency (%)")
        ax_b.set_ylabel("Bootstrap detection frequency (%)")
    real_save(fig, captured["base"], captured["dpi"])


def page_pixels(pdf: Path) -> np.ndarray:
    out = Path(tempfile.mkdtemp()) / "page.png"
    subprocess.run([sys.executable, "-c",
                    f"import fitz;d=fitz.open({str(pdf)!r});"
                    f"d[0].get_pixmap(dpi=110).save({str(out)!r})"], check=True)
    from PIL import Image
    return np.asarray(Image.open(out).convert("L"), dtype=np.int16)


def main() -> None:
    data = load_inputs()
    tmp = Path(tempfile.mkdtemp())

    draw(tmp / "replica", data, data["subsample"], annotate=False)
    a, b = page_pixels(SUBMITTED), page_pixels(tmp / "replica.pdf")
    if a.shape != b.shape:
        raise SystemExit(f"replica geometry differs: {a.shape} vs {b.shape}")
    delta = float(np.mean(np.abs(a - b)))
    if delta > 0.6:
        raise SystemExit("replica does not reproduce the submitted figure "
                         f"(mean |delta| = {delta:.3f} grey levels)")
    print(f"replica of the submitted panel b reproduces Submit/Figure_S2.pdf "
          f"(mean |delta| = {delta:.3f} grey levels)")

    p4 = next(c for c in data["candidates"] if str(c["peak_id"]) == "P4")
    print(f"  P4 ({p4['primary_peak_position_bp']} bp) HC support: bootstrap "
          f"{100 * float(p4['hc_bootstrap_detection_frequency']):.1f}% "
          f"(criterion, threshold 80%) vs equal-count "
          f"{100 * data['subsample'][(fl.GROUP_HC, 'P4')]['subsample_detection_frequency']:.1f}%")

    draw(tmp / "Figure_S2", data, data["bootstrap"], annotate=True)
    for target in (SUBMITTED, DEPOSITED):
        shutil.copy2(tmp / "Figure_S2.pdf", target)
        print(f"wrote {target.relative_to(ROOT)}")
    for suffix in (".svg", ".png"):
        src = (tmp / "Figure_S2").with_suffix(suffix)
        if src.exists():
            shutil.copy2(src, ROOT / "revise" / f"Figure_S2_criterion{suffix}")
    print(f"preview: {tmp / 'Figure_S2.png'}")


if __name__ == "__main__":
    main()
