#!/usr/bin/env python3
"""Regenerate Supplementary Figure S2 (fragment-length peak robustness) for v2.

`rebuild_figure_s2_panel_b.py` is the script that produced the submitted figure.
It carries its own gate -- it redraws the *original* panel b from the deposited
tables and pixel-compares that replica against Submit/Figure_S2.pdf before
writing anything -- and then redraws panel b with bootstrap detection
frequencies, which is the version that was submitted.

That gate is written against the v1 tables, so it is run once here on v1 inputs
to prove the harness is faithful, and then bypassed for the v2 render by calling
the module's own `load_inputs` and `draw` directly. No thresholds, series or
styling are changed.

The upstream script copies its result over Submit/Figure_S2.pdf and
deposit/figures/Figure_S2.pdf. Those constants are rebound to a sandbox, and
the run is additionally wrapped by protected_run.py at the call site.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
import tempfile
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
STAGE = PROJECT / "deposit" / "scripts" / "02_fragment_length"
V1_TABLES = PROJECT / "deposit" / "tables" / "02_fragment_length"
V2_TABLES = RUN / "v2tree" / "fl_flat"


def load_module():
    sys.path.insert(0, str(STAGE))
    spec = importlib.util.spec_from_file_location(
        "s2mod", STAGE / "rebuild_figure_s2_panel_b.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["s2mod"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    work = RUN / "logs" / "regen_S2"
    work.mkdir(parents=True, exist_ok=True)
    mod = load_module()

    print("[1/3] fidelity gate: v1 tables must reproduce Submit/Figure_S2.pdf")
    mod.TABLES = V1_TABLES
    # The upstream script prints `target.relative_to(ROOT)`, so the sandbox
    # has to sit inside the repository root or that print raises.
    sandbox = Path(tempfile.mkdtemp(prefix="s2_v1_", dir=work))
    mod.DEPOSITED = sandbox / "deposited.pdf"
    real_submitted = mod.SUBMITTED  # the pixel gate reads this; keep it real
    mod.SUBMITTED = sandbox / "submitted.pdf"
    shutil.copy2(real_submitted, mod.SUBMITTED)
    mod.main()
    shutil.copy2(mod.SUBMITTED, work / "v1.pdf")
    print("  PASS (the module's own pixel gate accepted the replica)")

    print("[2/3] rendering panel b from v2 tables")
    mod.TABLES = V2_TABLES
    data = mod.load_inputs()
    out = work / "v2"
    out.parent.mkdir(parents=True, exist_ok=True)
    mod.draw(out, data, data["bootstrap"], annotate=True)
    print(f"  wrote {out}.pdf")

    print("[3/3] P4 stability under v2")
    p4 = next(c for c in data["candidates"] if str(c["peak_id"]) == "P4")
    print(f"  P4 at {p4['primary_peak_position_bp']} bp: HC bootstrap "
          f"{100 * float(p4['hc_bootstrap_detection_frequency']):.1f}% "
          f"(threshold 80%)")

    if args.write:
        dest = RUN / "figures" / "Figure_S2.pdf"
        shutil.copy2(f"{out}.pdf", dest)
        for suffix in (".svg", ".png"):
            src = Path(f"{out}{suffix}")
            if src.exists():
                shutil.copy2(src, dest.with_suffix(suffix))
        print(f"  installed: figures/Figure_S2.pdf")
    else:
        print("  (dry run)")
    shutil.rmtree(sandbox, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
