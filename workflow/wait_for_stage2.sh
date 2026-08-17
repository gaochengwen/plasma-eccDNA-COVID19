#!/usr/bin/env bash
# Poll QDUH until every Stage 2 job has produced its completion marker.
# Exits 0 when all markers are present, 1 on timeout, 2 if a job disappeared
# from the queue without leaving a marker (i.e. it failed).
#
# Each job writes its own marker only after verifying its output artefact, so a
# marker means "the table is on disk and has the expected shape", not merely
# "the process exited 0".

set -uo pipefail

L=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/downstream/circlemap_methods/logs
ART=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/downstream/circlemap_artifact_masked

MARKERS=(
    "$L/extensions.multicell.done"
    "$L/extensions.mappability.done"
    "$L/extensions.density.done"
    "$L/extensions.binning_check.done"
    "$L/p0_burden_control.done"
    "$L/fragment_length.done"
    "$L/rca_bias.done"
    "$L/rank_candidate_loci.done"
    "$L/target_mask_audit.done"
    "$L/display_candidate_annotation.done"
    "$ART/downstream_validation.txt"
)

INTERVAL=${INTERVAL:-300}
MAX_HOURS=${MAX_HOURS:-10}
deadline=$(( $(date +%s) + MAX_HOURS * 3600 ))

while true; do
    report=$(ssh -o BatchMode=yes QDUH "
        for m in ${MARKERS[*]}; do
            if [ -f \"\$m\" ]; then echo \"DONE \$(basename \$(dirname \$m))/\$(basename \$m)\";
            else echo \"WAIT \$(basename \$m)\"; fi
        done
        echo \"QUEUE \$(qstat -u gaochengwen 2>/dev/null | grep -cE ' [QRH] ')\"
    " 2>/dev/null)

    n_wait=$(printf '%s\n' "$report" | grep -c '^WAIT' || true)
    n_queue=$(printf '%s\n' "$report" | awk '/^QUEUE/{print $2}')

    printf '[%s] pending=%s running_or_queued=%s\n' \
        "$(date -u +%FT%TZ)" "$n_wait" "${n_queue:-?}"

    if [[ "$n_wait" -eq 0 ]]; then
        echo "ALL STAGE 2 MARKERS PRESENT"
        printf '%s\n' "$report"
        exit 0
    fi

    # Nothing left in the queue but markers still missing => something failed.
    if [[ "${n_queue:-1}" -eq 0 ]]; then
        echo "QUEUE EMPTY BUT MARKERS MISSING -- job failure"
        printf '%s\n' "$report" | grep '^WAIT'
        exit 2
    fi

    if [[ $(date +%s) -ge $deadline ]]; then
        echo "TIMEOUT after ${MAX_HOURS}h"
        printf '%s\n' "$report" | grep '^WAIT'
        exit 1
    fi

    sleep "$INTERVAL"
done
