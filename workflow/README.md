# workflow — shared environment and submission

| File | Role |
|---|---|
| `v2_env.sh` | The single source of truth for cluster paths, the primary call set, random seeds and resampling counts. Every PBS job sources it. Edit paths here, not in individual jobs. |
| `submit_v2_stage2.sh` | Stage-2 orchestrator. Gates on every input before queueing anything (78 call-set symlinks, no broken links, the two prerequisite tables), then submits branches A, B and C and writes a job-id manifest. |
| `wait_for_marker.sh` | Block until one named `*.done` marker appears. |
| `wait_for_stage2.sh` | Block until all stage-2 markers appear, then report the table counts each job recorded. |

Run `submit_v2_stage2.sh` from the cluster login node. Compute nodes have no outbound
network, so all reference files must already be on shared storage.
