#!/usr/bin/env python3
"""Regenerate a submitted figure from the v2 (methods) tables, with a fidelity gate.

The method is borrowed from the v1 author's own rerun harness, which is the only
safe way to re-run these generators: several of them write hard-coded output
names that predate two rounds of supplementary renumbering, so trusting a
filename has already caused one figure to overwrite another.

Three steps per figure:

  1. Fidelity gate. Render with the v1 deposited tables and require the result
     to fingerprint-match the submitted PDF (creation timestamps stripped). If
     the harness cannot reproduce what was submitted, it has no business
     producing a replacement, and nothing is written.

  2. v2 render. Run the same plotting code against the v2 tables.

  3. Annotation diff. Extract every glyph run from both SVGs in document order
     and report the differences. That list is exactly the set of printed values
     the manuscript and caption must be updated to match.

Usage:
    python3 scripts/local/regen_figure_v2.py --figure S8 [--write]

Without --write nothing is installed; the diff is reported and the renders are
left in the work directory for inspection.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import re
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
V1_MODDIR = PROJECT / "deposit" / "scripts" / "06_chromatin"
V1_TABLES = PROJECT / "deposit" / "tables" / "06_chromatin"
V2_TABLES = RUN / "tables" / "chromatin"
SUBMIT = PROJECT / "Submit"

# Per-figure recipe: which module, which plotting entry point, which tables it
# reads, what the module names its own output, and which submitted file the
# result actually corresponds to after renumbering.
RECIPES = {
    "S8": {
        "module": "run_p0_burden_control.py",
        "entry": "make_p0_figure",
        "stem": "p0_burden_control",
        "submitted": "Figure_S8.pdf",
        "tables": [
            "p0_within_group_enrichment.tsv",
            "p0_downsampled_sample_enrichment.tsv",
            "p0_downsampled_group_stats.tsv",
            "relative_enrichment_group_stats.tsv",
            "p0_burden_matched_subset.tsv",
            "p0_covariate_adjusted_hc3.tsv",
            "sample_relative_enrichment.tsv",
            "p0_burden_association.tsv",
        ],
    },
}


def load_module(filename: str):
    sys.path.insert(0, str(V1_MODDIR))
    spec = importlib.util.spec_from_file_location("figmod", V1_MODDIR / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["figmod"] = mod
    spec.loader.exec_module(mod)
    return mod


def pdf_fingerprint(path: Path) -> str:
    data = path.read_bytes()
    for pattern in (
        rb"/CreationDate\s*\([^)]*\)",
        rb"/ModDate\s*\([^)]*\)",
        rb"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d",
        rb"D:\d{14}",
    ):
        data = re.sub(pattern, b"", data)
    return hashlib.md5(data).hexdigest()


def svg_text(path: Path) -> list[str]:
    """Every glyph run the backend wrote, in document order."""
    out: list[str] = []
    for el in ET.parse(path).iter():
        tag = el.tag.rsplit("}", 1)[-1]
        if tag in ("text", "tspan") and (el.text or "").strip():
            out.append(el.text.strip())
        if tag == "use":
            ref = el.get("{http://www.w3.org/1999/xlink}href") or el.get("href") or ""
            if ref.startswith(("#Deja", "#Liberation", "#Arial")):
                out.append(ref)
    return out


@contextlib.contextmanager
def sandbox(allowed: Path):
    """Block copies out of the render sandbox for the duration of the call.

    These plotting modules do more than draw. make_p0_figure, for example, ends
    by copying its output over two hard-coded absolute paths, one of which is
    `revise/Figure_S4.*` -- a name that, after two rounds of supplementary
    renumbering, no longer refers to this figure at all. Calling the module
    directly therefore silently overwrites an unrelated working figure with
    this one's content. That is the same class of incident the renumbering note
    records, and it happened here on the first run of this harness.

    Copies whose destination is inside the sandbox are allowed through; anything
    else is dropped and reported, so the module still runs unmodified.
    """
    allowed = allowed.resolve()
    originals = {name: getattr(shutil, name) for name in ("copy", "copy2", "copyfile")}
    blocked: list[str] = []

    def guard(fn):
        def wrapper(src, dst, *a, **kw):
            try:
                inside = Path(dst).resolve().is_relative_to(allowed)
            except (OSError, ValueError):
                inside = False
            if inside:
                return fn(src, dst, *a, **kw)
            blocked.append(str(dst))
            return dst
        return wrapper

    for name, fn in originals.items():
        setattr(shutil, name, guard(fn))
    try:
        yield blocked
    finally:
        for name, fn in originals.items():
            setattr(shutil, name, fn)


def render(mod, recipe: dict, source: Path, work: Path, tag: str) -> tuple[Path, Path]:
    tmp = Path(tempfile.mkdtemp(prefix=f"regen_{tag}_"))
    dirs = mod.base.ensure_dirs(tmp)
    missing = [n for n in recipe["tables"] if not (source / n).is_file()]
    if missing:
        raise SystemExit(f"{tag}: missing input tables in {source}: {missing}")
    for name in recipe["tables"]:
        shutil.copy2(source / name, dirs.tables / name)
    with sandbox(tmp) as blocked:
        getattr(mod, recipe["entry"])(dirs)
    for dst in blocked:
        print(f"  [sandbox] blocked out-of-tree write: {dst}")
    pdf = work / f"{tag}.pdf"
    svg = work / f"{tag}.svg"
    shutil.copy2(dirs.figures / f"{recipe['stem']}.pdf", pdf)
    shutil.copy2(dirs.figures / f"{recipe['stem']}.svg", svg)
    shutil.rmtree(tmp, ignore_errors=True)
    return pdf, svg


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", required=True, choices=sorted(RECIPES))
    ap.add_argument("--write", action="store_true",
                    help="install the v2 render into New_run/figures/")
    args = ap.parse_args()

    recipe = RECIPES[args.figure]
    submitted = SUBMIT / recipe["submitted"]
    work = RUN / "logs" / f"regen_{args.figure}"
    work.mkdir(parents=True, exist_ok=True)

    mod = load_module(recipe["module"])

    print(f"[1/3] fidelity gate: rendering {args.figure} from v1 deposited tables")
    v1_pdf, v1_svg = render(mod, recipe, V1_TABLES, work, "v1")
    if pdf_fingerprint(v1_pdf) != pdf_fingerprint(submitted):
        print(f"  FAIL: v1 tables do not reproduce {submitted.name}; nothing written.")
        print(f"  rendered={pdf_fingerprint(v1_pdf)} submitted={pdf_fingerprint(submitted)}")
        return 1
    print(f"  PASS: v1 tables reproduce {submitted.name} through the module's own code")

    print("[2/3] rendering from v2 (methods) tables")
    v2_pdf, v2_svg = render(mod, recipe, V2_TABLES, work, "v2")

    print("[3/3] annotation diff (every printed value that changed)")
    a, b = svg_text(v1_svg), svg_text(v2_svg)
    if a == b:
        print("  no printed annotation changed")
    else:
        n = 0
        for x, y in zip(a, b):
            if x != y:
                n += 1
                print(f"    {x!r} -> {y!r}")
        if len(a) != len(b):
            print(f"    (glyph-run count differs: v1={len(a)} v2={len(b)})")
        print(f"  {n} annotation(s) changed")

    if pdf_fingerprint(v2_pdf) == pdf_fingerprint(v1_pdf):
        print("  NOTE: the rendered PDF is unchanged; the v2 shift is below drawing resolution")

    if args.write:
        dest = RUN / "figures" / recipe["submitted"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(v2_pdf, dest)
        shutil.copy2(v2_svg, dest.with_suffix(".svg"))
        print(f"  installed: {dest.relative_to(RUN)}")
    else:
        print(f"  (dry run; renders in {work.relative_to(RUN)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
