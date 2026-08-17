#!/usr/bin/env bash
# S2 orchestrator -- submit the whole methods-call-set recompute.
#
# Run this on mu01. It only submits; nothing long-running is tied to the SSH
# session. Compute nodes have no outbound network, so every input is verified
# to exist on /gpfs before any job is queued.
#
#   ssh QDUH 'cd /gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder/v2_scripts && ./submit_v2_stage2.sh'
#
# Branches A, B and C are independent. The only cross-branch dependency is
# A4 (Figure 5 locus tables) -> C3 (locus ranking) -> C2 (artifact-masked
# downstream), expressed with PBS afterok.

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/v2_env.sh"

MANIFEST=$V2ROOT/logs/stage2_job_ids.tsv
mkdir -p "$V2ROOT/logs"

echo "== preflight =="
test -d "$ADAPTER/covid" && test -d "$ADAPTER/normal"
n_bed=$(find "$ADAPTER/covid" "$ADAPTER/normal" -type l | wc -l)
n_bam=$(find "$ADAPTER/covid_bam" "$ADAPTER/normal_bam" -type l | wc -l)
n_broken=$(find "$ADAPTER" -xtype l | wc -l)
echo "adapter bed symlinks : $n_bed (expect 78)"
echo "adapter bam symlinks : $n_bam (expect >= 78)"
echo "broken symlinks      : $n_broken (expect 0)"
test "$n_bed" -eq 78
test "$n_broken" -eq 0
test -s "$V2ROOT/chromatin/tables/sample_relative_enrichment.tsv"
test -s "$V2ROOT/eccGene/results/junction_abundance_wilcoxon.tsv"
v2_prepare_dirs
echo "preflight=PASS"
echo

printf 'branch\tjob_name\tjob_id\tdepends_on\tsubmitted_utc\n' > "$MANIFEST"
record() {
    printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-none}" "$(date -u +%FT%TZ)" \
        | tee -a "$MANIFEST"
}

echo "== branch A: chromatin =="
# One stage per job: qsub -v truncates comma-separated values.
for stage in multicell mappability density; do
    jid=$(qsub -v "EXT_STAGES=$stage" -N "v2ext_$stage" "$HERE/v2_chromatin_extensions.pbs")
    record A2 "ext_$stage" "$jid"
done
jid_p0=$(qsub "$HERE/v2_p0_burden_control.pbs")
record A3 p0_burden "$jid_p0"

echo
echo "== branch B: fragment length and RCA =="
jid_len=$(qsub "$HERE/v2_fragment_length.pbs")
record B1 fragment_length "$jid_len"
jid_rca=$(qsub "$HERE/v2_rca_bias.pbs")
record B2 rca_bias "$jid_rca"

echo
echo "== branch C: artifact-masked downstream, then locus ranking =="
jid_c2=$(qsub -N v2_ds_artmask -v CALLSET=circlemap_artifact_masked \
    "$CF_ROOT/scripts/run_downstream_callset.pbs")
record C2 downstream_artifact_masked "$jid_c2"

jid_c3=$(qsub -W "depend=afterok:$jid_c2" "$HERE/v2_rank_candidate_loci.pbs")
record C3 rank_candidate_loci "$jid_c3" "$jid_c2"

echo
echo "Submitted. Job id manifest: $MANIFEST"
echo "Figure 5 locus tables (A4) are held until Gate 1 is decided:"
echo "  inspect $V2ROOT/BCL3/tables/selected_loci.tsv once C3 finishes,"
echo "  compare against $V1_BCL3/tables/selected_loci.tsv, then submit A4."
