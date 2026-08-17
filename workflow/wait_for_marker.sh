#!/usr/bin/env bash
# Poll QDUH for a single completion marker. Exits 0 when it appears, 2 if the
# queue drains without it (job failed), 1 on timeout.
set -uo pipefail
MARKER=${1:?usage: wait_for_marker.sh <absolute marker path>}
INTERVAL=${INTERVAL:-120}
MAX_HOURS=${MAX_HOURS:-6}
deadline=$(( $(date +%s) + MAX_HOURS * 3600 ))
while true; do
    state=$(ssh -o BatchMode=yes QDUH "test -f '$MARKER' && echo DONE || echo WAIT; \
        qstat -u gaochengwen 2>/dev/null | grep -cE ' [QRH] '" 2>/dev/null)
    have=$(printf '%s\n' "$state" | head -1)
    nq=$(printf '%s\n' "$state" | tail -1)
    printf '[%s] %s queue=%s\n' "$(date -u +%FT%TZ)" "$have" "$nq"
    [[ "$have" == DONE ]] && exit 0
    [[ "${nq:-1}" -eq 0 ]] && { echo "queue drained without marker"; exit 2; }
    [[ $(date +%s) -ge $deadline ]] && { echo TIMEOUT; exit 1; }
    sleep "$INTERVAL"
done
