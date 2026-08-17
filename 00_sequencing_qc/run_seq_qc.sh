#!/bin/bash
# Per-sample sequencing QC for the Covid-eccDNA revision.
# Emits flagstat (read pairs, mapping rate) and samtools markdup stats (duplicate rate).
set -euo pipefail

SAMPLE="$1"
GROUP="$2"
ROOT=/gpfs/data/gao/Covid-eccDNA
OUT="$ROOT/Revise/seq_qc"
if [ "$GROUP" = "COVID-19" ]; then BAMDIR="$ROOT/covid_bam"; else BAMDIR="$ROOT/normal_bam"; fi
BAM="$BAMDIR/sorted_${SAMPLE}_circle.bam"

module load samtools-1.17

mkdir -p "$OUT/flagstat" "$OUT/markdup" "$OUT/tmp"
TMP="$OUT/tmp/$SAMPLE"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

samtools flagstat -@ 8 "$BAM" > "$OUT/flagstat/${SAMPLE}.flagstat.txt"

# Duplicates were never marked in these BAMs, so mark them here.
# Stream collate -> fixmate -m -> sort -> markdup; discard the marked BAM, keep only stats.
samtools collate -@ 8 -O -u "$BAM" "$TMP/collate" \
  | samtools fixmate -@ 4 -m -u - - \
  | samtools sort -@ 8 -u -m 2G -T "$TMP/sort" - \
  | samtools markdup -@ 4 -f "$OUT/markdup/${SAMPLE}.markdup.txt" - /dev/null

echo "DONE ${SAMPLE}"
