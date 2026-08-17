#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
OUTPUT="$ROOT/tables/software_versions.tsv"
PYTHON=/public/opt/miniforge3_2024/bin/python
SAMBLASTER=/gpfs/softwares/STRetch/tools/bwa.kit/samblaster

module load bwa-0.7.17 >/dev/null 2>&1 || true
module load samtools-1.17 >/dev/null 2>&1 || true
module load bedtools2-2.28.0 >/dev/null 2>&1 || true

printf 'component\tversion_or_commit\tpath_or_note\n' > "$OUTPUT"
printf 'Circle_finder\t%s\t%s\n' \
  '3eb333db2ea6277dde36cbf640be9afeb710c717' \
  "$ROOT/tools/Circle_finder" >> "$OUTPUT"
printf 'Circle_finder_pipeline_sha256\t%s\t%s\n' \
  "$(sha256sum "$ROOT/tools/Circle_finder/circle_finder-pipeline-bwa-mem-samblaster.sh" | awk '{print $1}')" \
  "$ROOT/tools/Circle_finder/circle_finder-pipeline-bwa-mem-samblaster.sh" >> "$OUTPUT"
printf 'BWA\t%s\t%s\n' \
  "$(bwa 2>&1 | awk -F': ' '/Version:/ {print $2; exit}' || true)" \
  "$(command -v bwa)" >> "$OUTPUT"
printf 'samtools\t%s\t%s\n' \
  "$(samtools --version | awk 'NR==1 {print $2}')" \
  "$(command -v samtools)" >> "$OUTPUT"
printf 'bedtools\t%s\t%s\n' \
  "$(bedtools --version | awk '{print $2}')" \
  "$(command -v bedtools)" >> "$OUTPUT"
printf 'samblaster\t%s\t%s\n' \
  "$("$SAMBLASTER" --version 2>&1 | head -n 1 | tr '\t' ' ' || true)" \
  "$SAMBLASTER" >> "$OUTPUT"
"$PYTHON" - <<'PY' >> "$OUTPUT"
import platform
import matplotlib
import numpy
import scipy
print(f"Python\t{platform.python_version()}\t{platform.python_implementation()}")
print(f"NumPy\t{numpy.__version__}\tPython package")
print(f"SciPy\t{scipy.__version__}\tPython package")
print(f"Matplotlib\t{matplotlib.__version__}\tPython package")
PY
