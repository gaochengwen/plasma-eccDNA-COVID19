#!/usr/bin/env python3
"""Apply tracked replacements to a Word document that already carries revisions.

The manuscript is a tracked-changes document from the previous revision round,
so some of the text being replaced is itself an unaccepted insertion. Deleting
such text correctly requires nesting a w:del inside the existing w:ins, which is
what this module does; text that is not inside a w:ins is wrapped in a plain
w:del. New text is added as a fresh w:ins attributed to the given author.

Matching is done against the "accepted" view of the document, that is, the
concatenation of w:t content with existing w:delText ignored, so a target string
is written exactly as it would read after accepting the previous round.
"""

from __future__ import annotations

import copy
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"

# ElementTree does not remember the prefixes it read, so every namespace that is
# not registered here comes back as ns1, ns2, ... Registering the ones Word uses
# keeps the rewritten document readable and, more importantly, keeps the prefix
# names that mc:Ignorable refers to. Namespaces that nothing in the body uses are
# still dropped on serialization, so run scripts/repair_docx_namespaces.py after
# the last edit of a document.
for _prefix, _uri in (
    ("w", W),
    ("mc", "http://schemas.openxmlformats.org/markup-compatibility/2006"),
    ("r", "http://schemas.openxmlformats.org/officeDocument/2006/relationships"),
    ("m", "http://schemas.openxmlformats.org/officeDocument/2006/math"),
    ("a", "http://schemas.openxmlformats.org/drawingml/2006/main"),
    ("pic", "http://schemas.openxmlformats.org/drawingml/2006/picture"),
    ("wp", "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"),
    ("w14", "http://schemas.microsoft.com/office/word/2010/wordml"),
    ("wp14", "http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing"),
    ("a14", "http://schemas.microsoft.com/office/drawing/2010/main"),
    ("o", "urn:schemas-microsoft-com:office:office"),
    ("v", "urn:schemas-microsoft-com:vml"),
    ("w10", "urn:schemas-microsoft-com:office:word"),
):
    ET.register_namespace(_prefix, _uri)


def qn(tag: str) -> str:
    return f"{{{W}}}{tag.split(':', 1)[1]}"


@dataclass
class RunSlot:
    run: ET.Element
    parent: ET.Element
    text_node: ET.Element
    start: int
    end: int


def load_document(path: Path) -> ET.ElementTree:
    with zipfile.ZipFile(path) as archive:
        return ET.ElementTree(ET.fromstring(archive.read("word/document.xml")))


def save_document(src_path: Path, tree: ET.ElementTree, dst_path: Path) -> None:
    body = ET.tostring(tree.getroot(), encoding="UTF-8", xml_declaration=True)
    with zipfile.ZipFile(src_path) as src:
        items = src.infolist()
        payload = {item.filename: src.read(item.filename) for item in items}
    payload["word/document.xml"] = body
    with zipfile.ZipFile(dst_path, "w", zipfile.ZIP_DEFLATED) as dst:
        for item in items:
            dst.writestr(item, payload[item.filename])


def _index_parents(root: ET.Element) -> dict:
    return {child: parent for parent in root.iter() for child in parent}


def build_index(root: ET.Element) -> Tuple[str, List[RunSlot]]:
    """Concatenated accepted text plus the run slots that produced it."""
    parents = _index_parents(root)
    slots: List[RunSlot] = []
    cursor = 0
    pieces: List[str] = []
    for run in root.iter(qn("w:r")):
        text_node = run.find(qn("w:t"))
        if text_node is None:
            continue
        text = text_node.text or ""
        if not text:
            continue
        parent = parents.get(run)
        if parent is None:
            continue
        slots.append(RunSlot(run, parent, text_node, cursor, cursor + len(text)))
        pieces.append(text)
        cursor += len(text)
    return "".join(pieces), slots


def _clone_run_with_text(template: ET.Element, text: str) -> ET.Element:
    run = ET.Element(qn("w:r"))
    r_pr = template.find(qn("w:rPr"))
    if r_pr is not None:
        run.append(copy.deepcopy(r_pr))
    t = ET.SubElement(run, qn("w:t"))
    t.set(XML_SPACE, "preserve")
    t.text = text
    return run


def _split_run(slot: RunSlot, local_offsets: Sequence[int]) -> List[ET.Element]:
    """Split one run into consecutive runs at the given local offsets."""
    text = slot.text_node.text or ""
    bounds = [0] + sorted(set(o for o in local_offsets if 0 < o < len(text))) + [len(text)]
    if len(bounds) == 2:
        return [slot.run]
    pieces = [text[bounds[i] : bounds[i + 1]] for i in range(len(bounds) - 1)]
    position = list(slot.parent).index(slot.run)
    new_runs = [_clone_run_with_text(slot.run, piece) for piece in pieces]
    slot.parent.remove(slot.run)
    for offset, run in enumerate(new_runs):
        slot.parent.insert(position + offset, run)
    return new_runs


def _mark_deleted(run: ET.Element, parent: ET.Element, author: str, date: str, rev_id: int) -> int:
    """Replace a run in place with a tracked deletion of the same run."""
    text_node = run.find(qn("w:t"))
    if text_node is None:
        return rev_id
    del_text = ET.Element(qn("w:delText"))
    del_text.set(XML_SPACE, "preserve")
    del_text.text = text_node.text
    run.remove(text_node)
    run.append(del_text)

    position = list(parent).index(run)
    parent.remove(run)
    del_el = ET.Element(qn("w:del"))
    del_el.set(qn("w:id"), str(rev_id))
    del_el.set(qn("w:author"), author)
    del_el.set(qn("w:date"), date)
    del_el.append(run)
    parent.insert(position, del_el)
    return rev_id + 1


def replace_tracked(
    tree: ET.ElementTree,
    old_text: str,
    new_paragraph_texts: Sequence[str],
    author: str,
    date: str,
    rev_id_start: int = 90000,
) -> int:
    """Delete `old_text` and insert `new_paragraph_texts` as tracked changes.

    Returns the next free revision id. Raises if `old_text` is absent or occurs
    more than once in the accepted view of the document.
    """
    root = tree.getroot()
    accepted, slots = build_index(root)
    occurrences = accepted.count(old_text)
    if occurrences == 0:
        raise KeyError(f"target text not found: {old_text[:90]!r}")
    if occurrences > 1:
        raise KeyError(f"target text occurs {occurrences} times, not unique: {old_text[:90]!r}")

    start = accepted.index(old_text)
    end = start + len(old_text)

    # Split the runs straddling the boundaries so the target maps to whole runs.
    for slot in slots:
        cuts = []
        if slot.start < start < slot.end:
            cuts.append(start - slot.start)
        if slot.start < end < slot.end:
            cuts.append(end - slot.start)
        if cuts:
            _split_run(slot, cuts)

    accepted, slots = build_index(root)
    start = accepted.index(old_text)
    end = start + len(old_text)
    targets = [s for s in slots if s.start >= start and s.end <= end and s.start < end]
    if not targets:
        raise RuntimeError("no runs mapped to the target span after splitting")

    rev_id = rev_id_start
    template_run = targets[0].run
    last = targets[-1]
    anchor_parent = last.parent
    anchor_index = list(anchor_parent).index(last.run)

    ins = ET.Element(qn("w:ins"))
    ins.set(qn("w:id"), str(rev_id))
    ins.set(qn("w:author"), author)
    ins.set(qn("w:date"), date)
    rev_id += 1
    for offset, piece in enumerate(new_paragraph_texts):
        if offset > 0:
            br = ET.SubElement(ins, qn("w:r"))
            r_pr = template_run.find(qn("w:rPr"))
            if r_pr is not None:
                br.append(copy.deepcopy(r_pr))
            ET.SubElement(br, qn("w:br"))
        ins.append(_clone_run_with_text(template_run, piece))
    anchor_parent.insert(anchor_index + 1, ins)

    for slot in targets:
        rev_id = _mark_deleted(slot.run, slot.parent, author, date, rev_id)

    return rev_id


def accepted_text(tree: ET.ElementTree) -> str:
    return build_index(tree.getroot())[0]
