#!/usr/bin/env python3
"""Helpers for reading and rewriting cells of the reviewer-response tables.

The three response letters are plain (untracked) three-column tables:
comment | author response | modification made in manuscript. This module
locates a row by a distinctive fragment of its comment text and replaces the
text of another cell in that row, preserving the run formatting already used in
the document so the rewritten cell looks like the ones written by hand.
"""

from __future__ import annotations

import copy
import re
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W}
# See the note in docx_tracked_edit.py: unregistered namespaces are re-serialized
# as ns1, ns2, ... which breaks the prefix names in mc:Ignorable. Run
# scripts/repair_docx_namespaces.py after the last edit of a document.
for _prefix, _uri in (
    ("w", W),
    ("mc", "http://schemas.openxmlformats.org/markup-compatibility/2006"),
    ("r", "http://schemas.openxmlformats.org/officeDocument/2006/relationships"),
    ("w14", "http://schemas.microsoft.com/office/word/2010/wordml"),
    ("wp14", "http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing"),
):
    ET.register_namespace(_prefix, _uri)


def qn(tag: str) -> str:
    prefix, local = tag.split(":", 1)
    return f"{{{W}}}{local}"


def load_document(path: Path) -> ET.ElementTree:
    with zipfile.ZipFile(path) as archive:
        data = archive.read("word/document.xml")
    return ET.ElementTree(ET.fromstring(data))


def save_document(src_path: Path, tree: ET.ElementTree, dst_path: Path) -> None:
    """Rewrite only word/document.xml, copying every other part byte-for-byte."""
    body = ET.tostring(tree.getroot(), encoding="UTF-8", xml_declaration=True)
    with zipfile.ZipFile(src_path) as src:
        items = src.infolist()
        payload = {item.filename: src.read(item.filename) for item in items}
    payload["word/document.xml"] = body
    with zipfile.ZipFile(dst_path, "w", zipfile.ZIP_DEFLATED) as dst:
        for item in items:
            dst.writestr(item, payload[item.filename])


def cell_text(cell: ET.Element) -> str:
    """Visible text of a table cell, ignoring anything marked as deleted."""
    parts: List[str] = []
    for paragraph in cell.iter(qn("w:p")):
        for node in paragraph.iter():
            if node.tag == qn("w:t"):
                parent_is_deleted = False
                parts.append(node.text or "")
            elif node.tag == qn("w:delText"):
                parent_is_deleted = True
        parts.append("\n")
    return "".join(parts).strip()


def iter_rows(tree: ET.ElementTree):
    root = tree.getroot()
    for table_index, table in enumerate(root.iter(qn("w:tbl"))):
        for row_index, row in enumerate(table.findall(qn("w:tr"))):
            cells = row.findall(qn("w:tc"))
            yield table_index, row_index, row, cells


def find_row(tree: ET.ElementTree, needle: str, column: int = 0):
    """Find the row whose `column` cell contains `needle`."""
    hits = []
    for table_index, row_index, row, cells in iter_rows(tree):
        if column < len(cells) and needle in cell_text(cells[column]):
            hits.append((table_index, row_index, row, cells))
    if not hits:
        raise KeyError(f"no row whose column {column} contains {needle!r}")
    if len(hits) > 1:
        raise KeyError(f"{len(hits)} rows contain {needle!r}; needle is not unique")
    return hits[0]


def _template_run(cell: ET.Element, document: ET.Element) -> Optional[ET.Element]:
    """A run whose formatting we can clone for new text."""
    for paragraph in cell.iter(qn("w:p")):
        for run in paragraph.findall(qn("w:r")):
            if run.find(qn("w:t")) is not None:
                return run
    # Fall back to any run with text anywhere in the document.
    for run in document.iter(qn("w:r")):
        if run.find(qn("w:t")) is not None:
            return run
    return None


def _template_paragraph(cell: ET.Element) -> Optional[ET.Element]:
    paragraphs = cell.findall(qn("w:p"))
    return paragraphs[0] if paragraphs else None


def set_cell_text(
    cell: ET.Element,
    paragraphs_text: Sequence[str],
    document: ET.Element,
) -> None:
    """Replace a cell's content with the given paragraphs, keeping its style."""
    template_p = _template_paragraph(cell)
    template_r = _template_run(cell, document)
    if template_p is None:
        raise ValueError("cell has no paragraph to use as a template")

    p_pr = template_p.find(qn("w:pPr"))
    r_pr = template_r.find(qn("w:rPr")) if template_r is not None else None

    tc_pr = cell.find(qn("w:tcPr"))
    for child in list(cell):
        if child is not tc_pr:
            cell.remove(child)

    for text in paragraphs_text:
        paragraph = ET.SubElement(cell, qn("w:p"))
        if p_pr is not None:
            paragraph.append(copy.deepcopy(p_pr))
        if text:
            run = ET.SubElement(paragraph, qn("w:r"))
            if r_pr is not None:
                run.append(copy.deepcopy(r_pr))
            t = ET.SubElement(run, qn("w:t"))
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
            t.text = text


def describe(path: Path) -> None:
    tree = load_document(path)
    print(f"### {path.name}")
    for table_index, row_index, row, cells in iter_rows(tree):
        texts = [cell_text(c).replace("\n", " ")[:70] for c in cells]
        print(f"  tbl{table_index} row{row_index} ncell={len(cells)}")
        for i, t in enumerate(texts):
            print(f"      c{i}: {t}")


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        describe(Path(arg))
