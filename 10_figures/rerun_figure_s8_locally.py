#!/usr/bin/env python3
"""Regenerate the burden-control supplementary figure (Figure S8) from the updated tables.

Panel c draws the nested-model COVID-19 - HC coefficients from
p0_covariate_adjusted_hc3.tsv, which was refitted with the rounded RCA
concentrations. `make_p0_figure` recomputes nothing, it only reads the tables in
`tables/`, and all of them are deposited, so the figure is rebuilt here with the
module's own plotting code.

Note on numbering: the module's docstring and its own output name still say
"Figure S4", from two numbering rounds ago. The file this figure is submitted as
is Submit/Figure_S8.pdf -- it was Figure_S7.pdf until the clinical-variables
subsection moved to Results position 2 and shifted the series. The validation
step below confirms the identity byte-for-byte rather than trusting the name.

Validation: the figure is rebuilt from the CURRENT deposited tables and must
reproduce the submitted PDF byte-for-byte once creation timestamps are stripped.
That is the invariant worth holding from here on -- the submitted figure is
exactly what the deposited tables produce. A second render from the
pre-rounding table is reported alongside it, so the effect of the correction on
every printed annotation stays visible.

    python3 revise/scripts/fixes/rerun_figure_s8_locally.py --check
    python3 revise/scripts/fixes/rerun_figure_s8_locally.py

Run from the repository root.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MODDIR = ROOT / "deposit" / "scripts" / "06_chromatin"
T06 = ROOT / "deposit" / "tables" / "06_chromatin"
BACKUP = ROOT / "revise" / "backups" / "rca_rounding"
PREROUND = BACKUP / "deposit_tables_06_chromatin_p0_covariate_adjusted_hc3.tsv.preround"

SUBMITTED = ROOT / "Submit" / "Figure_S8.pdf"
DESTS = (ROOT / "Submit" / "Figure_S8.pdf", ROOT / "deposit" / "figures" / "Figure_S8.pdf")

NEEDED = [
    "p0_within_group_enrichment.tsv",
    "p0_downsampled_sample_enrichment.tsv",
    "p0_downsampled_group_stats.tsv",
    "relative_enrichment_group_stats.tsv",
    "p0_burden_matched_subset.tsv",
    "p0_covariate_adjusted_hc3.tsv",
    "sample_relative_enrichment.tsv",
    "p0_burden_association.tsv",
]


def load_module():
    sys.path.insert(0, str(MODDIR))
    spec = importlib.util.spec_from_file_location("p0mod", MODDIR / "run_p0_burden_control.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["p0mod"] = m
    spec.loader.exec_module(m)
    return m


def pdf_fingerprint(path: Path) -> str:
    b = path.read_bytes()
    for pat in (rb"/CreationDate\s*\([^)]*\)", rb"/ModDate\s*\([^)]*\)",
                rb"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d", rb"D:\d{14}"):
        b = re.sub(pat, b"", b)
    return hashlib.md5(b).hexdigest()


def svg_text(path: Path) -> list[str]:
    """Every glyph run matplotlib wrote, in document order."""
    tree = ET.parse(path)
    out = []
    for el in tree.iter():
        tag = el.tag.rsplit("}", 1)[-1]
        if tag in ("text", "tspan") and (el.text or "").strip():
            out.append(el.text.strip())
        if tag == "use":
            ref = el.get("{http://www.w3.org/1999/xlink}href") or el.get("href") or ""
            if ref.startswith("#Deja") or ref.startswith("#Liberation") or ref.startswith("#Arial"):
                out.append(ref)
    return out


def render(m, covariate_table: Path, workdir: Path, tag: str) -> tuple[Path, Path]:
    tmp = Path(tempfile.mkdtemp(prefix=f"p0fig_{tag}_"))
    dirs = m.base.ensure_dirs(tmp)
    for name in NEEDED:
        src = covariate_table if name == "p0_covariate_adjusted_hc3.tsv" else T06 / name
        shutil.copy2(src, dirs.tables / name)
    m.make_p0_figure(dirs)
    pdf, svg = workdir / f"{tag}.pdf", workdir / f"{tag}.svg"
    shutil.copy2(dirs.figures / "p0_burden_control.pdf", pdf)
    shutil.copy2(dirs.figures / "p0_burden_control.svg", svg)
    shutil.rmtree(tmp, ignore_errors=True)
    return pdf, svg


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    m = load_module()
    work = Path(tempfile.mkdtemp(prefix="p0fig_out_"))

    new_pdf, new_svg = render(m, T06 / "p0_covariate_adjusted_hc3.tsv", work, "rounded")
    if pdf_fingerprint(new_pdf) != pdf_fingerprint(SUBMITTED):
        print("VALIDATION FAILED - the deposited tables no longer reproduce "
              f"{SUBMITTED.relative_to(ROOT)}; nothing written.")
        sys.exit(1)
    print(f"  validation: the deposited tables reproduce {SUBMITTED.relative_to(ROOT)} "
          "byte-for-byte through the module's own plotting code")

    old_pdf, old_svg = render(m, PREROUND, work, "preround")

    a, b = svg_text(old_svg), svg_text(new_svg)
    if a == b:
        print("  annotations: identical before and after rounding "
              "(every printed value and every significance tier is unchanged)")
    else:
        print(f"  annotations: {sum(x != y for x, y in zip(a, b))} glyph runs differ "
              f"(lengths {len(a)} vs {len(b)}) - inspect before installing")
        for x, y in zip(a, b):
            if x != y:
                print(f"    {x!r} -> {y!r}")

    if pdf_fingerprint(new_pdf) == pdf_fingerprint(old_pdf):
        print("  the rendered PDF is unchanged; the coefficient shift is below "
              "the resolution of the drawing")

    if args.check:
        print(f"  renders kept in {work}")
        return

    for dest in DESTS:
        shutil.copy2(new_pdf, dest)
        print(f"  installed -> {dest.relative_to(ROOT)}")
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
