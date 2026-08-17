#!/usr/bin/env python3
"""Verify a figure PDF against the Nature figure specification.

Checks the physical page size, that fonts are embedded and are an accepted
sans-serif, that text is real text rather than outlines or raster, and that no
raster image is present. Font sizes and line widths cannot be read reliably from
the PDF operator stream, so those are asserted in the plotting code and listed
here for the record.
"""

from __future__ import annotations

import argparse
import re
import sys
import zlib
from pathlib import Path

ACCEPTED_FONTS = ("Arial", "Helvetica", "LiberationSans", "Liberation Sans", "NimbusSans")
MM_PER_PT = 25.4 / 72.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf")
    parser.add_argument("--expect-width-mm", type=float, default=183.0)
    parser.add_argument("--tolerance-mm", type=float, default=1.5)
    args = parser.parse_args()

    path = Path(args.pdf)
    raw = path.read_bytes()
    ok = True

    boxes = re.findall(rb"/MediaBox\s*\[\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s*\]", raw)
    if not boxes:
        print("FAIL: no /MediaBox found")
        return 1
    x0, y0, x1, y1 = (float(v) for v in boxes[0])
    width_mm = (x1 - x0) * MM_PER_PT
    height_mm = (y1 - y0) * MM_PER_PT
    print(f"page size: {width_mm:.1f} x {height_mm:.1f} mm")
    if abs(width_mm - args.expect_width_mm) > args.tolerance_mm:
        print(f"FAIL: width {width_mm:.1f} mm differs from {args.expect_width_mm} mm")
        ok = False
    else:
        print(f"  OK width is {args.expect_width_mm} mm double column")
    if height_mm > 247.0:
        print(f"FAIL: height {height_mm:.1f} mm exceeds the 247 mm maximum")
        ok = False
    else:
        print("  OK height within 247 mm")

    base_fonts = sorted(set(m.decode() for m in re.findall(rb"/BaseFont\s*/([#\w\-+,.]+)", raw)))
    print(f"fonts referenced: {base_fonts}")
    if not base_fonts:
        print("FAIL: no fonts referenced; text may have been converted to paths")
        ok = False
    for font in base_fonts:
        if not any(name.replace(" ", "") in font.replace(" ", "") for name in ACCEPTED_FONTS):
            print(f"FAIL: {font} is not an accepted sans-serif")
            ok = False
    embedded = len(re.findall(rb"/FontFile2|/FontFile3|/FontFile\b", raw))
    print(f"embedded font programs: {embedded}")
    if embedded < len(base_fonts):
        print("FAIL: not every font is embedded")
        ok = False
    else:
        print("  OK all fonts embedded")

    subtypes = sorted(set(m.decode() for m in re.findall(rb"/Subtype\s*/(\w+)", raw)))
    print(f"object subtypes: {subtypes}")
    if "Image" in subtypes:
        print("FAIL: figure contains a raster image")
        ok = False
    else:
        print("  OK no raster image; art is vector")

    # Text-showing operators confirm live text rather than outlined glyphs.
    shown = 0
    for match in re.finditer(rb"stream\r?\n", raw):
        start = match.end()
        end = raw.find(b"endstream", start)
        if end < 0:
            continue
        try:
            content = zlib.decompress(raw[start:end])
        except zlib.error:
            continue
        shown += len(re.findall(rb"\bTj\b|\bTJ\b", content))
    print(f"text-showing operators: {shown}")
    if shown == 0:
        print("FAIL: no text operators; labels appear to be outlines")
        ok = False
    else:
        print("  OK labels are live text")

    print("\nPASS: figure meets the checked Nature requirements" if ok else "\nFAIL: see above")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
