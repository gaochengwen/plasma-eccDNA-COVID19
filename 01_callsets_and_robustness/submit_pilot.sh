#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
PILOT="$ROOT/config/pilot_samples.tsv"
JOB_MANIFEST="$ROOT/manifests/pilot_job_ids.tsv"
REFERENCE_JOB_ID=${1:-}

test -s "$PILOT"
mkdir -p "$ROOT/manifests"

printf 'sample_id\tjob_id\treference_dependency\tsubmitted_utc\n' > "$JOB_MANIFEST"

tail -n +2 "$PILOT" | while IFS=$'\t' read -r order sample group reason; do
    if [[ -n "$REFERENCE_JOB_ID" ]]; then
        job_id=$(qsub -N "CF_${sample}" \
            -v "SAMPLE_ID=$sample" \
            -W "depend=afterok:$REFERENCE_JOB_ID" \
            "$ROOT/scripts/run_one_sample.pbs")
    else
        job_id=$(qsub -N "CF_${sample}" \
            -v "SAMPLE_ID=$sample" \
            "$ROOT/scripts/run_one_sample.pbs")
    fi
    printf '%s\t%s\t%s\t%s\n' \
        "$sample" "$job_id" "${REFERENCE_JOB_ID:-none}" \
        "$(date -u +%FT%TZ)" | tee -a "$JOB_MANIFEST"
done

