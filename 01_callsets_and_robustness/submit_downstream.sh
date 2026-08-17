#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
DEPENDENCY=${1:-}
MANIFEST="$ROOT/manifests/downstream_job_ids.tsv"
CALLSETS=(circlemap_methods circlemap_high_support_masked consensus_t10)

printf 'callset\tjob_id\tdependency\tsubmitted_utc\n' > "$MANIFEST"
for callset in "${CALLSETS[@]}"; do
    if [[ -n "$DEPENDENCY" ]]; then
        job=$(qsub \
          -N "DS_${callset:0:11}" \
          -v "CALLSET=$callset" \
          -W "depend=afterok:$DEPENDENCY" \
          "$ROOT/scripts/run_downstream_callset.pbs")
    else
        job=$(qsub \
          -N "DS_${callset:0:11}" \
          -v "CALLSET=$callset" \
          "$ROOT/scripts/run_downstream_callset.pbs")
    fi
    printf '%s\t%s\t%s\t%s\n' \
      "$callset" "$job" "${DEPENDENCY:-none}" "$(date -u +%FT%TZ)" | tee -a "$MANIFEST"
done
