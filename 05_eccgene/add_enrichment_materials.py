#!/usr/bin/env python3
"""Add the gene-set over-representation materials to the deposit bundle.

The bundle promised "every processed table needed to re-derive the reported
statistics", but neither over-representation analysis was represented in it: the
Results report a g:Profiler confirmation against a 26,642-gene custom background
and Figure 3c plots a Metascape run, and no file backed either claim.

This script deposits the g:Profiler materials, which do exist in full, and adds
a README that states plainly what is and is not reproducible for each analysis.
It does not reconstruct the Metascape input: 157 genes tie at the q value on the
3,000-gene boundary of Supplementary Table S14, so the stated selection rule does
not identify a unique set and any reconstruction would be a guess.

Run from the repository root.
"""

from __future__ import annotations

import hashlib
import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DEP = ROOT / "deposit"
SRC = ROOT / "revise" / "source_data"
DST = DEP / "tables" / "04_eccgene"

COPIES = {
    "eccGene_gProfiler_query_genes.txt": "gprofiler_query_genes.txt",
    "eccGene_gProfiler_background_genes.txt": "gprofiler_background_genes.txt",
    "eccGene_gProfiler_request.json": "gprofiler_request.json",
    "eccGene_gProfiler_response.json": "gprofiler_response.json",
    "eccGene_gProfiler_all_significant_terms.tsv": "gprofiler_all_significant_terms.tsv",
    "Figure_3C_gProfiler_representative_terms.tsv": "gprofiler_representative_terms.tsv",
}

REVISION_NOTE = (
    "Reassembled after depositing the g:Profiler over-representation materials in "
    "tables/04_eccgene (query genes, custom background, request, response and the "
    "55 significant terms), adding tables/04_eccgene/README_enrichment.md, and "
    "correcting two stale supplementary references in README.md. See "
    "scripts/04_eccgene/add_enrichment_materials.py."
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_inputs() -> None:
    """Refuse to deposit anything that does not match what the paper claims."""
    query = [l.strip() for l in (SRC / "eccGene_gProfiler_query_genes.txt")
             .read_text(encoding="utf-8").splitlines() if l.strip()]
    background = [l.strip() for l in (SRC / "eccGene_gProfiler_background_genes.txt")
                  .read_text(encoding="utf-8").splitlines() if l.strip()]
    terms = (SRC / "eccGene_gProfiler_all_significant_terms.tsv") \
        .read_text(encoding="utf-8").strip().splitlines()
    request = json.loads((SRC / "eccGene_gProfiler_request.json").read_text(encoding="utf-8"))

    assert len(background) == 26642, len(background)
    assert len(terms) - 1 == 55, len(terms) - 1
    assert request["domain_scope"] == "custom", request["domain_scope"]
    assert set(request["sources"]) == {"GO:BP", "REAC", "KEGG"}, request["sources"]
    assert len(query) == 3094, len(query)
    print(f"inputs verified: query {len(query)}, background {len(background)}, "
          f"{len(terms) - 1} significant terms, custom domain scope")


def main() -> None:
    check_inputs()

    for src_name, dst_name in COPIES.items():
        shutil.copy2(SRC / src_name, DST / dst_name)
    print(f"copied {len(COPIES)} g:Profiler files into {DST.relative_to(ROOT)}")

    script_dst = DEP / "scripts" / "04_eccgene" / Path(__file__).name
    shutil.copy2(Path(__file__), script_dst)
    print(f"copied this script into {script_dst.relative_to(ROOT)}")

    files = sorted(p for p in DEP.rglob("*")
                   if p.is_file()
                   and p.name not in ("checksums.sha256", "bundle_manifest.json"))
    lines, records, total = [], [], 0
    for p in files:
        rel = p.relative_to(DEP).as_posix()
        digest = sha256(p)
        size = p.stat().st_size
        total += size
        lines.append(f"{digest}  {rel}")
        records.append({"path": rel, "sha256": digest, "bytes": size})

    (DEP / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    mf = DEP / "bundle_manifest.json"
    old = json.loads(mf.read_text(encoding="utf-8")) if mf.exists() else {}
    stage_counts: dict[str, int] = {}
    for r in records:
        stage = r["path"].split("/")[0]
        stage_counts[stage] = stage_counts.get(stage, 0) + 1
    manifest = {
        "assembled_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "n_files": len(records),
        "total_bytes": total,
        "stage_counts": stage_counts,
        "files": records,
    }
    if old.get("assembled_utc"):
        manifest["previous_assembly_utc"] = old["assembled_utc"]
        manifest["revision_note"] = REVISION_NOTE
    mf.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"rebuilt checksums.sha256 and bundle_manifest.json over {len(records)} files "
          f"({total / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
