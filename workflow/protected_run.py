#!/usr/bin/env python3
"""Run a figure generator as a subprocess without letting it damage the tree.

Several generators in this project end by copying their output over hard-coded
absolute paths -- `revise/Figure_S4.*`, `revise/Figure_S5.*`,
`analyse/chromatin/deposit/figures/*`. Those names date from before two rounds
of supplementary renumbering, so the destination is usually a *different*
figure than the one being drawn. The project's own renumbering note records one
such incident; two more happened while building this pipeline, both recovered.

An in-process shutil guard only covers generators imported as modules. For the
argparse-style generators that have to be run as subprocesses, this wrapper
takes the belt-and-braces approach: fingerprint every file in the protected
directories first, run the command, then restore anything that changed and
report it.

    python3 scripts/local/protected_run.py -- python3 make_figure_s5.py --outdir ...
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[3]

# Directories a generator has no business writing to during a regeneration run.
PROTECTED = [
    (PROJECT / "revise", ("*.pdf", "*.svg", "*.png")),
    (PROJECT / "Submit", ("*",)),
    (PROJECT / "deposit" / "figures", ("*",)),
    (PROJECT / "deposit" / "tables", ("*.tsv",)),
    (PROJECT / "analyse" / "chromatin" / "deposit" / "figures", ("*",)),
]


def digest(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def snapshot(stash: Path) -> dict[Path, tuple[str, Path]]:
    state: dict[Path, tuple[str, Path]] = {}
    for root, patterns in PROTECTED:
        if not root.is_dir():
            continue
        for pattern in patterns:
            for path in root.glob(pattern):
                if not path.is_file():
                    continue
                keep = stash / f"{len(state):05d}_{path.name}"
                shutil.copy2(path, keep)
                state[path] = (digest(path), keep)
    return state


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        ap.error("no command given")

    stash = Path(tempfile.mkdtemp(prefix="protected_run_"))
    try:
        before = snapshot(stash)
        print(f"[protected] watching {len(before)} files", flush=True)

        result = subprocess.run(command)

        restored = []
        for path, (want, keep) in before.items():
            if not path.is_file():
                shutil.copy2(keep, path)
                restored.append(f"{path} (deleted)")
            elif digest(path) != want:
                shutil.copy2(keep, path)
                restored.append(str(path))

        if restored:
            print(f"[protected] the command modified {len(restored)} protected "
                  f"file(s); all restored:", flush=True)
            for item in restored:
                print(f"  restored {item}", flush=True)
        else:
            print("[protected] no protected file was touched", flush=True)
        return result.returncode
    finally:
        shutil.rmtree(stash, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
