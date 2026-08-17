#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 SAMPLE_ID" >&2
    exit 2
fi

SAMPLE_ID=$1
ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
PROJECT=/gpfs/data/gao/Covid-eccDNA
MANIFEST=/gpfs/data/gao/Covid-eccDNA/Revise/chromatin/inputs/eccdna_sample_manifest.tsv
REFERENCE="$ROOT/reference/hg38.canonical.fa"
UPSTREAM="$ROOT/tools/Circle_finder/circle_finder-pipeline-bwa-mem-samblaster.sh"
SAMBLASTER=/gpfs/softwares/STRetch/tools/bwa.kit/samblaster
THREADS=${PBS_NP:-12}
JOB_TAG=${PBS_JOBID:-manual_$(date -u +%Y%m%dT%H%M%SZ)}
RUN_DIR="$ROOT/work/$SAMPLE_ID/run_$JOB_TAG"
RESULT_DIR="$ROOT/results/circlefinder_raw/$SAMPLE_ID"
LOG_DIR="$ROOT/logs/samples"

mkdir -p "$RUN_DIR" "$RESULT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/${SAMPLE_ID}.${JOB_TAG}.log"
exec > >(tee "$LOG") 2>&1

module load bwa-0.7.17 >/dev/null 2>&1
module load samtools-1.17 >/dev/null 2>&1
module load bedtools2-2.28.0 >/dev/null 2>&1
module load parallel-20170206 >/dev/null 2>&1 || true

# The samblaster installation directory also contains an obsolete samtools
# binary. Expose only samblaster through an isolated shim so that samtools 1.17
# and the module-provided BWA cannot be shadowed.
SHIM_DIR="$RUN_DIR/command_shims"
mkdir -p "$SHIM_DIR"
ln -s "$SAMBLASTER" "$SHIM_DIR/samblaster"
export PATH="$SHIM_DIR:$PATH"
hash -r
export LC_ALL=C
export TMPDIR="$RUN_DIR/tmp"
mkdir -p "$TMPDIR"

echo "sample_id=$SAMPLE_ID"
echo "job_tag=$JOB_TAG"
echo "started_utc=$(date -u +%FT%TZ)"
echo "host=$(hostname)"
echo "threads=$THREADS"
echo "bwa_path=$(command -v bwa)"
echo "samtools_path=$(command -v samtools)"
echo "bedtools_path=$(command -v bedtools)"
echo "samblaster_path=$(command -v samblaster)"

BWA_VERSION=$(bwa 2>&1 | awk -F': ' '/Version:/ {print $2; exit}' || true)
SAMTOOLS_VERSION=$(samtools --version | awk 'NR==1 {print $2}')
BEDTOOLS_VERSION=$(bedtools --version | awk '{print $2}')
SAMBLASTER_VERSION=$(samblaster --version 2>&1 | head -n 1 || true)
echo "bwa_version=$BWA_VERSION"
echo "samtools_version=$SAMTOOLS_VERSION"
echo "bedtools_version=$BEDTOOLS_VERSION"
echo "samblaster_version=$SAMBLASTER_VERSION"
[[ "$BWA_VERSION" == 0.7.17* ]]
[[ "$SAMTOOLS_VERSION" == 1.17* ]]
[[ "$BEDTOOLS_VERSION" == v2.28.0* ]]

test -s "$MANIFEST"
test -s "$REFERENCE"
test -s "$REFERENCE.bwt"
test -s "$UPSTREAM"
test -x "$SAMBLASTER"

ROW=$(awk -F '\t' -v sample="$SAMPLE_ID" '
    NR > 1 && $1 == sample {print; found=1; exit}
    END {if (!found) exit 1}
' "$MANIFEST")
GROUP=$(awk -F '\t' '{print $2}' <<< "$ROW")
BED=$(awk -F '\t' '{print $3}' <<< "$ROW")
BAM=$(awk -F '\t' '{print $4}' <<< "$ROW")

echo "group=$GROUP"
echo "bam=$BAM"
echo "circlemap_bed=$BED"
test -s "$BAM"
test -s "$BAM.bai"
test -s "$BED"

samtools quickcheck -v "$BAM"
samtools idxstats "$BAM" |
    awk 'BEGIN{m=0;u=0} {m+=$3;u+=$4}
         END{printf "mapped=%d\nunmapped=%d\ntotal_alignments=%d\n",m,u,m+u}' |
    tee "$RUN_DIR/input_alignment_counts.txt"

R1="$RUN_DIR/${SAMPLE_ID}.R1.fq.gz"
R2="$RUN_DIR/${SAMPLE_ID}.R2.fq.gz"
export BAM R1 R2

cd "$RUN_DIR"

set +e
/usr/bin/time -v -o "$RUN_DIR/fastq_recovery.time.txt" \
  bash -o pipefail -c '
    samtools collate -@ 4 -u -O -T "$TMPDIR/collate" "$BAM" |
        samtools fastq -@ 4 -N \
            -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null -
' > "$RUN_DIR/fastq_recovery.stdout.log" \
  2> "$RUN_DIR/fastq_recovery.stderr.log"
PRODUCER_STATUS=$?
set -e

echo "fastq_recovery_exit_status=$PRODUCER_STATUS"
if [[ "$PRODUCER_STATUS" -ne 0 ]]; then
    exit 1
fi
test -s "$R1"
test -s "$R2"
gzip -t "$R1" "$R2"
echo "r1_fastq_gzip_bytes=$(stat -c %s "$R1")"
echo "r2_fastq_gzip_bytes=$(stat -c %s "$R2")"

set +e
/usr/bin/time -v -o "$RUN_DIR/circle_finder.time.txt" \
    bash "$UPSTREAM" "$THREADS" "$REFERENCE" "$R1" "$R2" \
        10 "$SAMPLE_ID" hg38
CALLER_STATUS=$?
set -e

echo "circle_finder_exit_status=$CALLER_STATUS"
if [[ "$CALLER_STATUS" -ne 0 || "$PRODUCER_STATUS" -ne 0 ]]; then
    exit 1
fi

RAW="$RUN_DIR/${SAMPLE_ID}-hg38.microDNA-JT.txt"
test -e "$RAW"

awk 'BEGIN{FS=OFS="\t"}
     NF == 4 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ &&
         $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ &&
         $4 ~ /^[0-9]+$/ && $2 < $3 {
         key=$1 FS $2 FS $3
         support[key]+=$4
     }
     END {
         for (key in support) print key, support[key]
     }' "$RAW" |
    sort -k1,1V -k2,2n -k3,3n > "$RESULT_DIR/${SAMPLE_ID}.circle_finder.deduplicated.bed"

gzip -c "$RAW" > "$RESULT_DIR/${SAMPLE_ID}.circle_finder.raw.tsv.gz"
cp "$RUN_DIR/input_alignment_counts.txt" "$RESULT_DIR/"
cp "$RUN_DIR/fastq_recovery.stderr.log" "$RESULT_DIR/"
cp "$RUN_DIR/fastq_recovery.time.txt" "$RESULT_DIR/"
cp "$RUN_DIR/circle_finder.time.txt" "$RESULT_DIR/"

{
    printf 'sample_id\tgroup\traw_rows\tdeduplicated_calls\tbam_size_bytes\tbed_size_bytes\tr1_fastq_gzip_bytes\tr2_fastq_gzip_bytes\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE_ID" "$GROUP" \
        "$(wc -l < "$RAW")" \
        "$(wc -l < "$RESULT_DIR/${SAMPLE_ID}.circle_finder.deduplicated.bed")" \
        "$(stat -c %s "$BAM")" "$(stat -c %s "$BED")" \
        "$(stat -c %s "$R1")" "$(stat -c %s "$R2")"
} > "$RESULT_DIR/${SAMPLE_ID}.run_summary.tsv"

sha256sum \
    "$RESULT_DIR/${SAMPLE_ID}.circle_finder.deduplicated.bed" \
    "$RESULT_DIR/${SAMPLE_ID}.circle_finder.raw.tsv.gz" \
    "$RESULT_DIR/${SAMPLE_ID}.run_summary.tsv" \
    > "$RESULT_DIR/${SAMPLE_ID}.sha256"

echo "completed_utc=$(date -u +%FT%TZ)"
touch "$RESULT_DIR/.complete"

# The upstream workflow creates very large sample-specific SAM/BAM/TXT
# intermediates. Results, checksums, timing, and diagnostic logs have already
# been copied above. Clean only this successful job's generated work directory;
# failed runs retain their work directory for diagnosis.
cd "$ROOT"
find "$RUN_DIR" \( -type f -o -type l \) -delete
find "$RUN_DIR" -depth -type d -empty -delete
