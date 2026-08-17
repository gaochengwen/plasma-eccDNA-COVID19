#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
MANIFEST=/gpfs/data/gao/Covid-eccDNA/Revise/chromatin/inputs/eccdna_sample_manifest.tsv
PILOT="$ROOT/config/pilot_samples.tsv"
DEPENDENCY=${1:-}
OUTPUT="$ROOT/manifests/full_remaining_job_ids.tsv"

declare -A pilot
while IFS=$'\t' read -r order sample group reason; do
    [[ "$order" == "pilot_order" ]] && continue
    pilot["$sample"]=1
done < "$PILOT"

printf 'sample_id\tjob_id\tdependency\tsubmitted_utc\ttarget_node\n' > "$OUTPUT"
while IFS=$'\t' read -r sample group bed bam bam_exists mapped; do
    [[ "$sample" == "sample_id" ]] && continue
    [[ -n "${pilot[$sample]:-}" ]] && continue
    if [[ -f "$ROOT/results/circlefinder_raw/$sample/.complete" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\n' \
          "$sample" "already_complete" "none" "$(date -u +%FT%TZ)" "none" | tee -a "$OUTPUT"
        continue
    fi
    if [[ -n "$DEPENDENCY" ]]; then
        job=$(qsub \
          -N "CF_${sample}" \
          -v "SAMPLE_ID=$sample" \
          -l "nodes=1:ppn=20" \
          -W "depend=afterok:$DEPENDENCY" \
          "$ROOT/scripts/run_one_sample.pbs")
    else
        job=$(qsub \
          -N "CF_${sample}" \
          -v "SAMPLE_ID=$sample" \
          -l "nodes=1:ppn=20" \
          "$ROOT/scripts/run_one_sample.pbs")
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$sample" "$job" "${DEPENDENCY:-none}" "$(date -u +%FT%TZ)" "scheduler" | tee -a "$OUTPUT"
done < "$MANIFEST"
