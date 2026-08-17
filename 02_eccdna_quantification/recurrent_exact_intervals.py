#!/usr/bin/env python3
"""Recurrent COVID-specific exact eccDNA intervals (Figure 3d input).

REIMPLEMENTED per pre-declaration amendment A3: the original code for this
artefact is not in the repository or on the cluster, so the definition is
taken verbatim from the Figure 3d legend. There is no v1 code to gate
against; results must be labelled as a re-implementation.

Kept as a standalone file rather than a heredoc inside the PBS script --
the scheduler mangled the heredoc's indentation at run time even though the
submitted file was correct on disk.
"""

import csv, sys
from collections import Counter
from pathlib import Path

adapter, outdir = Path(sys.argv[1]), Path(sys.argv[2])
CANON = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def samples(group_dir):
    return sorted(p.name.replace("_circle_site.bed", "")
                  for p in group_dir.glob("*_circle_site.bed"))


covid = samples(adapter / "covid")
hc = samples(adapter / "normal")
print(f"COVID={len(covid)} HC={len(hc)}", flush=True)
assert len(covid) == 39 and len(hc) == 39


def intervals(bed):
    """Distinct exact (chrom, start, end) keys present in one sample."""
    seen = set()
    with bed.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in CANON:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            if end <= start:
                continue
            seen.add(f"{f[0]}:{start}_{end}")
    return seen


# Pass 1 -- recurrence across COVID-19 samples only.
per_sample = {}
counts = Counter()
for i, s in enumerate(covid, 1):
    keys = intervals(adapter / "covid" / f"{s}_circle_site.bed")
    per_sample[s] = keys
    counts.update(keys)
    print(f"[{i}/39] {s}: {len(keys):,} distinct intervals", flush=True)

for min_n in (10, 6):
    candidates = {k for k, c in counts.items() if c >= min_n}
    print(f"min{min_n}: {len(candidates):,} intervals in >= {min_n} COVID samples",
          flush=True)

    # Pass 2 -- drop anything seen in any HC sample.
    for j, s in enumerate(hc, 1):
        if not candidates:
            break
        candidates -= intervals(adapter / "normal" / f"{s}_circle_site.bed")
    print(f"min{min_n}: {len(candidates):,} remain after HC exclusion", flush=True)

    def sort_key(key):
        chrom, rest = key.split(":")
        start = int(rest.split("_")[0])
        order = (int(chrom[3:]) if chrom[3:].isdigit()
                 else {"X": 23, "Y": 24}[chrom[3:]])
        return (-counts[key], order, start)

    rows = sorted(candidates, key=sort_key)
    out = outdir / f"covid_specific_exact.matrix.min{min_n}.tsv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([""] + covid)
        for key in rows:
            w.writerow([key] + [1 if key in per_sample[s] else 0 for s in covid])
    print(f"wrote {out.name}: {len(rows)} rows x {len(covid)} COVID samples",
          flush=True)
