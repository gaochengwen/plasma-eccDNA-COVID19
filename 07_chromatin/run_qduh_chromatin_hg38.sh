#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=${PROJECT_DIR:-/gpfs/data/gao/Covid-eccDNA}
OUTDIR=${OUTDIR:-/gpfs/data/gao/Covid-eccDNA/Revise/chromatin}
PYTHON_BIN=${PYTHON_BIN:-/public/opt/miniforge3_2024/bin/python}
SAMTOOLS_BIN=${SAMTOOLS_BIN:-samtools}
LIFTOVER_BIN=${LIFTOVER_BIN:-/gpfs/softwares/ricopili202511/dependencies/liftover/liftOver}

mkdir -p "$OUTDIR"/{scripts,logs}

unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY

module load Miniforge3_2024 >/dev/null 2>&1 || true
module load samtools-1.9 >/dev/null 2>&1 || true
module load bedtools2-2.28.0 >/dev/null 2>&1 || true

cd "$OUTDIR"

"$PYTHON_BIN" scripts/run_chromatin_hg38_analysis.py \
  --stage prepare \
  --project-dir "$PROJECT_DIR" \
  --outdir "$OUTDIR" \
  --liftover "$LIFTOVER_BIN" \
  --samtools "$SAMTOOLS_BIN" \
  2>&1 | tee "$OUTDIR/logs/prepare_reference.$(date +%Y%m%d_%H%M%S).log"

PBS="$OUTDIR/scripts/chromatin_hg38_analysis.pbs"
cat > "$PBS" <<'PBS'
#!/usr/bin/env bash
#PBS -N chromatin_hg38
#PBS -q batch
#PBS -l nodes=1:ppn=4
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -o /gpfs/data/gao/Covid-eccDNA/Revise/chromatin/logs/chromatin_hg38.o
#PBS -e /gpfs/data/gao/Covid-eccDNA/Revise/chromatin/logs/chromatin_hg38.e

set -euo pipefail

PROJECT_DIR=/gpfs/data/gao/Covid-eccDNA
OUTDIR=/gpfs/data/gao/Covid-eccDNA/Revise/chromatin
PYTHON_BIN=/public/opt/miniforge3_2024/bin/python
LIFTOVER_BIN=/gpfs/softwares/ricopili202511/dependencies/liftover/liftOver

unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY

module load Miniforge3_2024 >/dev/null 2>&1 || true
module load samtools-1.9 >/dev/null 2>&1 || true
module load bedtools2-2.28.0 >/dev/null 2>&1 || true

cd "$OUTDIR"
"$PYTHON_BIN" scripts/run_chromatin_hg38_analysis.py \
  --stage analyze \
  --project-dir "$PROJECT_DIR" \
  --outdir "$OUTDIR" \
  --liftover "$LIFTOVER_BIN" \
  --samtools samtools
PBS

qsub "$PBS"
