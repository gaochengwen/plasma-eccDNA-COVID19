#!/usr/bin/env python3
"""Render Supplementary Figure S5 (clinical correlates) from the v2 tables.

The generator is the project's own create_figure_s9_clinical_correlates.py --
its file name still carries the pre-renumbering S9, but the figure it produces
is submitted as Figure_S5. It reads the A2-A5 tables and the merged clinical
input; it computes nothing itself, so pointing it at the v2 tables is enough.

Its COPIES constant writes straight into Submit/ and deposit/figures/, so both
that and the output stem are rebound to this repository before it is called.
"""
from __future__ import annotations
import importlib.util, sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[2]
PROJECT = RUN.parent
SCRIPT = PROJECT / "revise" / "scripts" / "create_figure_s9_clinical_correlates.py"

spec = importlib.util.spec_from_file_location("figs5", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules["figs5"] = mod
sys.path.insert(0, str(SCRIPT.parent))
spec.loader.exec_module(mod)

mod.TABLES = RUN / "v2tree" / "clinical_tables"
mod.OUT_STEM = RUN / "figures" / "Figure_S5"
mod.COPIES = ()                      # never write outside this repository

if __name__ == "__main__":
    mod.main()
    print(f"wrote {mod.OUT_STEM}.pdf")
