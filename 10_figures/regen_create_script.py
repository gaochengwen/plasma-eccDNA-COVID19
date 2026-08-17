#!/usr/bin/env python3
"""Re-run a `revise/scripts/create_*.py` figure script against the v2 tables.

These scripts differ from the analysis modules: they take no arguments and
resolve their inputs and outputs from module-level Path constants derived from
__file__. They are also the layer that produced the *submitted* look -- the
analysis modules' own figure output does not match what was submitted for most
figures, because these scripts restyle it.

So they are imported, their path constants are rebound to the v2 tree and to a
sandbox output directory, and only then are they executed. The script file
itself is never edited.

Same three-step contract as regen_figure_v2.py:
  1. fidelity gate -- run against v1 inputs, require a fingerprint match with
     the submitted PDF;
  2. render against v2 inputs;
  3. diff every printed glyph run between the two.
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
V2TREE = RUN / "v2tree"
SUBMIT = PROJECT / "Submit"

# `rebind` maps a module-level constant to the v1 and v2 value it should take.
# `produces` maps the file the script writes into its output dir to the
# submitted figure it corresponds to after the two renumbering rounds.
RECIPES = {
    "S4": {
        "script": PROJECT / "revise" / "scripts" / "create_figure_s4.py",
        "rebind": {
            "TABLE_DIR": (PROJECT / "analyse" / "Circle_finder" / "tables",
                          V2TREE / "Circle_finder" / "tables"),
        },
        # OUTPUT_STEM is a full path stem, not a directory, so it is rebound
        # from the sandbox directory rather than taken from `outdir_const`.
        "outdir_const": "OUTPUT_DIR",
        "stem_const": ("OUTPUT_STEM", "Figure_S4_callset_robustness"),
        "produces": {"Figure_S4_callset_robustness.pdf": "Figure_S4.pdf"},
    },
    # RCA technical-bias figure. Its own name is Figure_S2, which after the two
    # renumbering rounds is submitted as Figure_S3.
    "S3": {
        "script": PROJECT / "revise" / "scripts" / "create_figure_s2.py",
        "rebind": {
            "INPUT_DIR": (PROJECT / "analyse" / "RCA" / "results",
                          V2TREE / "RCA" / "results"),
        },
        "outdir_const": "OUTPUT_DIR",
        "produces": {"Figure_S2.pdf": "Figure_S3.pdf"},
    },
    # Gene-assignment sensitivity. Written as Figure_S3, submitted as Figure_S6.
    # This script also emits Figure_3 and the graphical abstract; only the
    # sensitivity panel is gated and installed here.
    "S6": {
        "script": RUN / "scripts" / "hpc" / "v2_patched" / "create_eccgene_revision_figures.py",
        "rebind": {
            "ANALYSIS": (PROJECT / "analyse" / "eccGene", V2TREE / "eccGene"),
            "RESULTS": (PROJECT / "analyse" / "eccGene" / "results",
                        V2TREE / "eccGene" / "results"),
            "MATRICES": (PROJECT / "analyse" / "eccGene" / "matrices",
                         V2TREE / "eccGene" / "matrices"),
            # v1 reads its cached g:Profiler response; v2 gets a fresh cache in
            # the sandbox so the v1 response is never overwritten.
            "SOURCE_DATA": (PROJECT / "revise" / "source_data", None),
        },
        "derived": {"SOURCE_DATA": "source_data"},
        "outdir_const": "REVISE",
        # main() also rebuilds Figure 3 and the graphical abstract. Neither is
        # what the manuscript ships -- the submitted Figure 3 is an Illustrator
        # assembly -- and the Figure 3 path additionally needs a static circle
        # matrix that is not part of the v2 regeneration. Only the enrichment
        # call and the sensitivity panel are run.
        "entries": ["gprofiler_enrichment", "plot_figure_s3"],
        "produces": {"Figure_S3.pdf": "Figure_S6.pdf"},
    },
    # Age-matched sensitivity analysis. Written as Figure_S1_redesigned,
    # submitted unchanged as Figure_S1.
    "S1": {
        "script": RUN / "scripts" / "hpc" / "v2_patched" / "create_figure_s1_redesigned.py",
        "rebind": {
            # The patched copy lives in this repo, so ROOT-derived constants no
            # longer resolve to the project tree; every one of them is bound
            # explicitly. Sample metadata is static and identical in v1 and v2.
            "METADATA_PATH": (
                PROJECT / "analyse" / "EPM" / "input" / "sample_metadata.tsv",
                PROJECT / "analyse" / "EPM" / "input" / "sample_metadata.tsv"),
            "BURDEN_PATH": (
                PROJECT / "analyse" / "EPM" / "results"
                / "sample_burden_recomputed_from_qduh_bed.tsv",
                V2TREE / "EPM" / "results"
                / "sample_burden_recomputed_from_qduh_bed.tsv"),
            "LENGTH_CURVES_PATH": (
                PROJECT / "analyse" / "fragment-length" / "source_data"
                / "sample_level_primary_curves.tsv",
                V2TREE / "fragment_length" / "source_data"
                / "sample_level_primary_curves.tsv"),
            "PEAK_PATH": (
                PROJECT / "analyse" / "fragment-length" / "tables"
                / "objective_peak_identification_and_stability.tsv",
                V2TREE / "fragment_length" / "tables"
                / "objective_peak_identification_and_stability.tsv"),
            # Age matching uses age with a 1-year caliper and no burden term,
            # so the pair list is identical under v2 and is reused as-is.
            "PAIRS_PATH": (PROJECT / "revise" / "source_data" / "age_matched_pairs.json",
                           PROJECT / "revise" / "source_data" / "age_matched_pairs.json"),
            "SOURCE_DATA_DIR": (PROJECT / "revise" / "source_data", None),
            # v1 peaks are 196/365/571; v2 moves P2 by 1 bp.
            "EXPECT_STABLE_PEAKS": ([196, 365, 571], [196, 366, 571]),
        },
        "derived": {"SOURCE_DATA_DIR": "source_data"},
        "outdir_const": "SOURCE_DATA_DIR",
        "stem_const": ("OUTPUT_BASE", "Figure_S1_redesigned"),
        "produces": {"Figure_S1_redesigned.pdf": "Figure_S1.pdf"},
    },
    # Downstream robustness across call sets. This script emits Figure_S6 and
    # Figure_S7 under its own numbering; only its Figure_S7 is submitted, as
    # Figure_S7. Its Figure_S6 output matches no submitted figure.
    "S7": {
        "script": RUN / "scripts" / "hpc" / "v2_patched" / "create_circlemap_robustness_figures.py",
        "rebind": {
            "ROOT": (PROJECT / "analyse" / "Circle_finder", V2TREE / "Circle_finder"),
            "TABLES": (PROJECT / "analyse" / "Circle_finder" / "tables",
                       V2TREE / "Circle_finder" / "tables"),
            # The historical v1 figure showed four targets. The revised v2
            # figure reports the three candidates selected for wet-lab work.
            "VALIDATION_TARGETS": (
                ["CTNNA2", "CAB39", "chr22_target", "SDK1"],
                ["CTNNA2", "SDK1", "TCF7L1"]),
        },
        "outdir_const": "REVISE",
        "produces": {"Figure_S7.pdf": "Figure_S7.pdf"},
    },
}


def pdf_fingerprint(path: Path) -> str:
    data = path.read_bytes()
    for pattern in (rb"/CreationDate\s*\([^)]*\)", rb"/ModDate\s*\([^)]*\)",
                    rb"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d", rb"D:\d{14}"):
        data = re.sub(pattern, b"", data)
    return hashlib.md5(data).hexdigest()


def pdf_normalized_text(path: Path) -> str:
    """All printed text, whitespace-stripped.

    The byte fingerprint cannot match for figures that were submitted from
    matplotlib 3.11 when only 3.9 is available here: identical drawing calls
    produce different PDF bytes across versions. For those, fidelity is
    established on content instead -- every printed character, in order.
    """
    import fitz
    doc = fitz.open(path)
    text = "".join(page.get_text("text") for page in doc)
    doc.close()
    return re.sub(r"\s+", "", text.replace("\u2212", "-").replace("\xa0", " "))


def svg_text(path: Path) -> list[str]:
    out: list[str] = []
    for el in ET.parse(path).iter():
        tag = el.tag.rsplit("}", 1)[-1]
        if tag in ("text", "tspan") and (el.text or "").strip():
            out.append(el.text.strip())
    return out


@contextlib.contextmanager
def sandbox(allowed: Path):
    allowed = allowed.resolve()
    originals = {n: getattr(shutil, n) for n in ("copy", "copy2", "copyfile")}
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


def run_variant(recipe: dict, which: int, work: Path, tag: str) -> Path:
    """which: 0 = v1 inputs, 1 = v2 inputs. Returns the output directory."""
    outdir = work / tag
    outdir.mkdir(parents=True, exist_ok=True)
    path = recipe["script"]
    spec = importlib.util.spec_from_file_location(f"createmod_{tag}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"createmod_{tag}"] = mod
    spec.loader.exec_module(mod)

    for const, values in recipe["rebind"].items():
        value = values[which]
        # A None on the v2 side means "derive this from the sandbox" -- used for
        # constants the script computes from its output directory at import
        # time, such as a response cache that must not be shared between the
        # v1 and v2 runs.
        if value is None:
            value = outdir / recipe["derived"][const]
        setattr(mod, const, value)
    setattr(mod, recipe["outdir_const"], outdir)
    if "stem_const" in recipe:
        const, stem = recipe["stem_const"]
        setattr(mod, const, outdir / stem)

    with sandbox(outdir) as blocked:
        for entry in recipe.get("entries", ["main"]):
            getattr(mod, entry)()
    for dst in blocked:
        print(f"  [sandbox] blocked out-of-tree write: {dst}")
    return outdir


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", required=True, choices=sorted(RECIPES))
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    recipe = RECIPES[args.figure]
    work = RUN / "logs" / f"regen_{args.figure}"
    work.mkdir(parents=True, exist_ok=True)

    print(f"[1/3] fidelity gate: {recipe['script'].name} against v1 inputs")
    v1_dir = run_variant(recipe, 0, work, "v1")
    ok = True
    for produced, submitted_name in recipe["produces"].items():
        got, want = v1_dir / produced, SUBMIT / submitted_name
        if not got.is_file():
            print(f"  FAIL: script did not produce {produced}")
            ok = False
            continue
        if pdf_fingerprint(got) == pdf_fingerprint(want):
            print(f"  PASS (byte): {produced} reproduces {submitted_name}")
        elif pdf_normalized_text(got) == pdf_normalized_text(want):
            print(f"  PASS (content): {produced} reproduces every printed "
                  f"character of {submitted_name}; bytes differ because the "
                  f"submitted file came from a different matplotlib")
        else:
            print(f"  FAIL: {produced} does not reproduce {submitted_name}")
            ok = False
    if not ok:
        print("  nothing written.")
        return 1

    print("[2/3] rendering against v2 tables")
    v2_dir = run_variant(recipe, 1, work, "v2")

    print("[3/3] annotation diff")
    for produced, submitted_name in recipe["produces"].items():
        a_svg = (v1_dir / produced).with_suffix(".svg")
        b_svg = (v2_dir / produced).with_suffix(".svg")
        if not (a_svg.is_file() and b_svg.is_file()):
            print(f"  {produced}: no SVG pair to diff")
            continue
        a, b = svg_text(a_svg), svg_text(b_svg)
        changed = [(x, y) for x, y in zip(a, b) if x != y]
        print(f"  {submitted_name}: {len(changed)} annotation(s) changed"
              f"{'' if len(a) == len(b) else f' (glyph runs {len(a)} vs {len(b)})'}")
        for x, y in changed[:40]:
            print(f"    {x!r} -> {y!r}")

        if args.write:
            dest = RUN / "figures" / submitted_name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(v2_dir / produced, dest)
            if b_svg.is_file():
                shutil.copy2(b_svg, dest.with_suffix(".svg"))
            print(f"  installed: figures/{submitted_name}")
    if not args.write:
        print(f"  (dry run; renders in {work.relative_to(RUN)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
