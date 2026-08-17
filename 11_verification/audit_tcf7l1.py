#!/usr/bin/env python3
"""Audit TCF7L1-circle in every original QDUH Circle-Map BED, with CTNNA2/SDK1 controls."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

OUT = Path("/tmp/claude-1000/-home-gao-eccDNA/7676d647-c0d3-42e5-9100-781aadabec79/scratchpad")

REMOTE = r'''
import csv, json, socket
from datetime import datetime, timezone
from pathlib import Path

METADATA = Path("/gpfs/data/gao/Covid-eccDNA/Revise/eccGene/sample_metadata.tsv")
TOL = 10

TARGETS = [
    ("TCF7L1circle", "chr2", 85213855, 85214235),
    ("CTNNA2circle", "chr2", 79885063, 79885444),
    ("SDK1circle",   "chr7", 3619742,  3620122),
]
CHROMS = {t[1] for t in TARGETS}

REFS = {
 "blacklist_gap": "/gpfs/data/gao/Covid-eccDNA/Revise/chromatin/reference/hg38_forbidden.blacklist_gaps.canonical.merged.bed",
 "umap_k100":     "/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/reference/normalized/umap_k100.canonical.merged.bed",
 "segdup":        "/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/reference/normalized/segmental_duplication.canonical.merged.bed",
 "repeatmasker":  "/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/reference/normalized/repeatmasker.canonical.merged.bed",
}
WINDOW = 100


def passes(split, score, mean_cov, sd, si, ei, unc):
    return (split >= 3 and score > 200 and mean_cov > sd
            and si >= 0.3 and ei >= 0.3 and unc < 0.1)


with METADATA.open(newline="") as fh:
    samples = list(csv.DictReader(fh, delimiter="\t"))
if len(samples) != 78:
    raise ValueError("expected 78 samples, got %d" % len(samples))

rows = []
for s in samples:
    bed = Path(s["bed_path"])
    if not bed.is_file():
        raise FileNotFoundError(bed)
    hits = {t[0]: {"exact": [], "tol": []} for t in TARGETS}
    total_rows = 0
    with bed.open() as fh:
        for line in fh:
            total_rows += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or f[0] not in CHROMS:
                continue
            try:
                st, en = int(f[1]), int(f[2])
            except ValueError:
                continue
            for name, chrom, ts, te in TARGETS:
                if f[0] != chrom:
                    continue
                if abs(st - ts) <= TOL and abs(en - te) <= TOL:
                    rec = {"start": st, "end": en, "fields": f}
                    hits[name]["tol"].append(rec)
                    if st == ts and en == te:
                        hits[name]["exact"].append(rec)

    for name, chrom, ts, te in TARGETS:
        h = hits[name]
        def best(recs):
            if not recs:
                return None
            # highest score wins if several calls fall inside the tolerance
            return max(recs, key=lambda r: float(r["fields"][5]))
        b_exact, b_tol = best(h["exact"]), best(h["tol"])
        row = {
            "circle": name,
            "coordinate_bed0": "%s:%d-%d" % (chrom, ts, te),
            "length_bp": te - ts,
            "sample_id": s["sample_id"],
            "group": s["group"],
            "exact_detected": bool(h["exact"]),
            "exact_row_count": len(h["exact"]),
            "tol10_detected": bool(h["tol"]),
            "tol10_row_count": len(h["tol"]),
            "tol10_coordinates": sorted({"%d-%d" % (r["start"], r["end"]) for r in h["tol"]}),
            "bed_rows_total": total_rows,
        }
        for tag, b in (("exact", b_exact), ("tol10", b_tol)):
            if b is None:
                row.update({tag + "_" + k: None for k in
                            ("discordant", "split", "score", "mean_cov", "cov_sd",
                             "start_inc", "end_inc", "uncovered")})
                row[tag + "_filter_pass"] = False
                continue
            v = [float(x) for x in b["fields"][3:11]]
            d, sp, sc, mc, sd_, si, ei, un = v
            row.update({tag + "_discordant": d, tag + "_split": sp, tag + "_score": sc,
                        tag + "_mean_cov": mc, tag + "_cov_sd": sd_,
                        tag + "_start_inc": si, tag + "_end_inc": ei,
                        tag + "_uncovered": un,
                        tag + "_filter_pass": passes(sp, sc, mc, sd_, si, ei, un)})
        rows.append(row)

# ---- breakpoint mask audit (same rules as build_masked_callsets.py) ----
windows = []
for name, chrom, ts, te in TARGETS:
    for side, pos in (("start", ts), ("end", te)):
        windows.append({"circle": name, "chrom": chrom, "side": side,
                        "pos": pos, "ws": pos - WINDOW, "we": pos + WINDOW,
                        "len": 2 * WINDOW})
for w in windows:
    for ref in REFS:
        w[ref + "_bp"] = 0
    w["repeat_annotations"] = []

for ref, path in REFS.items():
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(path)
    with p.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in CHROMS:
                continue
            try:
                a, b = int(f[1]), int(f[2])
            except ValueError:
                continue
            for w in windows:
                if f[0] != w["chrom"]:
                    continue
                ov = min(b, w["we"]) - max(a, w["ws"])
                if ov > 0:
                    w[ref + "_bp"] += ov
                    if ref == "repeatmasker" and len(f) > 3:
                        w["repeat_annotations"].append("%s:%d-%d" % (f[3], a, b))

payload = {
    "remote_host": socket.getfqdn(),
    "run_time_utc": datetime.now(timezone.utc).isoformat(),
    "metadata_path": str(METADATA),
    "tolerance_bp": TOL,
    "targets": [{"circle": n, "chrom": c, "start": a, "end": b} for n, c, a, b in TARGETS],
    "sample_count": len(samples),
    "rows": rows,
    "mask_windows": windows,
}
print(json.dumps(payload, separators=(",", ":")))
'''


def main() -> None:
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", "QDUH", "python3", "-"],
                       input=REMOTE, text=True, capture_output=True)
    if r.returncode != 0:
        raise SystemExit("remote failed:\n" + r.stderr[-4000:])
    payload = json.loads(r.stdout)
    (OUT / "tcf7l1_audit_payload.json").write_text(json.dumps(payload, indent=1))

    rows = payload["rows"]
    print("host=%s  samples=%d  tolerance=%dbp\n" % (
        payload["remote_host"], payload["sample_count"], payload["tolerance_bp"]))

    hdr = ("%-14s %-9s %7s %7s %7s %7s" %
           ("circle", "group", "exact", "ex+pass", "tol10", "t10+pass"))
    print(hdr)
    print("-" * len(hdr))
    for circle in ("TCF7L1circle", "CTNNA2circle", "SDK1circle"):
        for grp in ("COVID-19", "HC"):
            sub = [x for x in rows if x["circle"] == circle and x["group"] == grp]
            n = len(sub)
            e = sum(1 for x in sub if x["exact_detected"])
            ep = sum(1 for x in sub if x["exact_filter_pass"])
            t = sum(1 for x in sub if x["tol10_detected"])
            tp = sum(1 for x in sub if x["tol10_filter_pass"])
            print("%-14s %-9s %7s %7s %7s %7s" % (
                circle, grp, "%d/%d" % (e, n), ep, "%d/%d" % (t, n), tp))

    print("\n== breakpoint mask audit (±100 bp windows, 200 bp each) ==")
    print("%-14s %-6s %9s %9s %8s %9s  %s" % (
        "circle", "side", "umap_frac", "repeat_f", "segdup", "blacklist", "repeat annotations"))
    for w in payload["mask_windows"]:
        print("%-14s %-6s %9.3f %9.3f %8d %9d  %s" % (
            w["circle"], w["side"], w["umap_k100_bp"] / w["len"],
            w["repeatmasker_bp"] / w["len"], w["segdup_bp"], w["blacklist_gap_bp"],
            "; ".join(w["repeat_annotations"]) or "-"))

    print("\n== TCF7L1circle positive samples ==")
    pos = [x for x in rows if x["circle"] == "TCF7L1circle" and x["tol10_detected"]]
    if not pos:
        print("  none")
    for x in sorted(pos, key=lambda r: (r["group"] != "COVID-19", r["sample_id"])):
        print("  %-8s %-9s coords=%s split=%s score=%s meancov=%s sd=%s "
              "si=%s ei=%s unc=%s exact=%s pass=%s" % (
                  x["sample_id"], x["group"], ",".join(x["tol10_coordinates"]),
                  x["tol10_split"], x["tol10_score"], x["tol10_mean_cov"], x["tol10_cov_sd"],
                  x["tol10_start_inc"], x["tol10_end_inc"], x["tol10_uncovered"],
                  x["exact_detected"], x["tol10_filter_pass"]))


if __name__ == "__main__":
    main()
