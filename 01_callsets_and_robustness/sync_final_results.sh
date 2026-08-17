#!/usr/bin/env bash
set -euo pipefail

QDUH_ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
LOCAL_ROOT=/home/gao/eccDNA/analyse/Circle_finder
MEI_ROOT=/home/gao/eccDNA/analyse/Circle_finder

FILE_LIST=$(mktemp)
FAILED_LIST=$(mktemp)
CHECK_LOG=$(mktemp)
trap 'rm -f "$FILE_LIST" "$FAILED_LIST" "$CHECK_LOG"' EXIT

mkdir -p "$LOCAL_ROOT"

ssh QDUH "cd '$QDUH_ROOT' && tar -cf - final_output_manifest.tsv" \
  | tar -C "$LOCAL_ROOT" -xf -
{
    echo "final_output_manifest.tsv"
    awk 'NR > 1 {print $1}' "$LOCAL_ROOT/final_output_manifest.tsv"
} | awk '!seen[$0]++' > "$FILE_LIST"
if [[ ! -s "$FILE_LIST" ]]; then
    echo "No files listed in QDUH final_output_manifest.tsv" >&2
    exit 1
fi

for iter in 1 2 3 4 5; do
    if (cd "$LOCAL_ROOT" && awk 'NR > 1 {print $3 "  " $1}' final_output_manifest.tsv | sha256sum -c - >"$CHECK_LOG" 2>&1); then
        break
    fi
    awk '
      /^sha256sum: .*: No such file/ {line=$0; sub(/^sha256sum: /, "", line); sub(/: No such file.*/, "", line); print line; next}
      /^sha256sum: .*: cannot open/ {line=$0; sub(/^sha256sum: /, "", line); sub(/: cannot open.*/, "", line); print line; next}
      /: FAILED open or read$/ {line=$0; sub(/: FAILED open or read$/, "", line); print line; next}
      /: FAILED$/ {line=$0; sub(/: FAILED$/, "", line); print line; next}
    ' "$CHECK_LOG" | sort -u > "$FAILED_LIST"
    if [[ ! -s "$FAILED_LIST" ]]; then
        cat "$CHECK_LOG" >&2
        exit 1
    fi
    ssh QDUH "cd '$QDUH_ROOT' && tar -cf - -T -" \
      < "$FAILED_LIST" \
      | tar -C "$LOCAL_ROOT" -xf -
    if [[ "$iter" == 5 ]]; then
        (cd "$LOCAL_ROOT" && awk 'NR > 1 {print $3 "  " $1}' final_output_manifest.tsv | sha256sum -c - >"$CHECK_LOG" 2>&1) || {
            cat "$CHECK_LOG" >&2
            exit 1
        }
    fi
done

ssh mei "mkdir -p '$MEI_ROOT'"
tar -C "$LOCAL_ROOT" -cf - -T "$FILE_LIST" \
  | ssh mei "tar -C '$MEI_ROOT' -xf -"

(cd "$LOCAL_ROOT" && awk 'NR > 1 {print $3 "  " $1}' final_output_manifest.tsv | sha256sum -c - >/dev/null)
ssh mei "cd '$MEI_ROOT' && awk 'NR > 1 {print \$3 \"  \" \$1}' final_output_manifest.tsv | sha256sum -c - >/dev/null"

manifest_files=$(awk 'END {print NR - 1}' "$LOCAL_ROOT/final_output_manifest.tsv")
manifest_bytes=$(awk 'NR > 1 {sum += $2} END {print sum}' "$LOCAL_ROOT/final_output_manifest.tsv")

{
  echo "PASS"
  echo "qduh_root=$QDUH_ROOT"
  echo "local_root=$LOCAL_ROOT"
  echo "mei_root=$MEI_ROOT"
  echo "manifest_files=$manifest_files"
  echo "manifest_bytes=$manifest_bytes"
  echo "verified_utc=$(date -u +%FT%TZ)"
  echo "verification=final_output_manifest_sha256_all_pass"
} > "$LOCAL_ROOT/synchronization_validation.txt"
rsync -a "$LOCAL_ROOT/synchronization_validation.txt" \
  "mei:${MEI_ROOT}/synchronization_validation.txt"
