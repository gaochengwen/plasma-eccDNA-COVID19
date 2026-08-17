# 01 — call sets and caller robustness

Builds the seven-call-set strictness ladder and the robustness synthesis reported in
**Tables S6–S8, S15, S16** and **Figures S4, S7**. See `docs/call_sets.md` for the
definition of each call set.

## Circle_finder as an independent second caller

| File | Role |
|---|---|
| `prepare_reference.pbs` | Build the reference and indices Circle_finder needs. |
| `prepare_strict_resources.py` / `.pbs` | ENCODE blacklist and assembly-gap resources for the artifact masks. |
| `run_circle_finder_sample.sh`, `run_one_sample.pbs` | Per-sample Circle_finder run. |
| `submit_pilot.sh`, `submit_remaining_samples.sh` | Pilot first, then the remaining samples. |
| `validate_pilot.pbs`, `validate_all.pbs`, `validate_circlefinder_runs.py` | Completion and integrity checks before anything downstream is allowed to start. |

## Call-set construction

| File | Role |
|---|---|
| `build_circlemap_sensitivity.py` / `.pbs` | The methods-consistent set and the two stricter split-read/score sets. |
| `build_masked_callsets.py` / `.pbs` | Breakpoint-level artifact masking (full-circle and ±100 bp window). |
| `build_consensus_callsets.py`, `build_pilot_consensus.pbs`, `build_all_consensus.pbs` | Circle-Map / Circle_finder two-caller consensus at 10 bp tolerance. |

## Downstream and synthesis

| File | Role |
|---|---|
| `run_downstream_callset.pbs`, `submit_downstream.sh` | Run the full downstream stack against one named call set. This is how stages 02 and 05 are executed for each arm of the ladder. |
| `run_core_robustness_analysis.py` + `v2_core_robustness.pbs` | Statistical synthesis across call sets. **Note:** `--root` is both input and output directory, so the job stages a private copy rather than pointing at the live tree. |
| `summarize_downstream_robustness.py` + `v2_summarize_downstream.pbs` | Cross-call-set eccGene and chromatin concordance (feeds Figure S7). |
| `audit_validation_target_masks.py` + `v2_target_mask_audit.pbs` | Audit the five validation targets against every prespecified mask; three go forward to the wet lab. |
| `generate_final_report.py` / `.pbs` | Human-readable robustness report plus a checksum manifest. |
| `build_circlemap_supplementary_tables.py` | Format Tables S6–S8 from the above. |
| `collect_software_versions.sh` | Record tool versions at run time. |
| `sync_final_results.sh` | Copy finalized results out of the scratch tree. |
