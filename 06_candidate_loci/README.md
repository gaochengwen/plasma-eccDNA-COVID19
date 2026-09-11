# 06 — candidate locus eligibility and ranking

The prespecified screen that chooses which loci are displayed and which go to the wet
lab. The ranking is computed first and the display loci follow from it — not the other
way round. Feeds **Table S16**.

| File | Role |
|---|---|
| `rank_candidate_loci.py` + `v2_rank_candidate_loci.pbs` | The eligibility audit and ranking. Primary statistics are read from the locked junction-based eccGene tables; this script does not recompute them. |
| `run_candidate_locus_ranking.pbs`, `run_ranking_only.pbs` | Job wrappers. |
| `locus_candidate_screen.pbs` | The upstream screen. |
| `annotate_chromatin_display_candidates.py` + `v2_display_candidate_annotation.pbs` | Apply a separate, explicitly recorded feasibility filter for display. Produces the ranking and selection-criteria tables. |
| `build_ranking_workbook.py` | Assemble the auditable ranking workbook. |
| `prepare_artifact_links.py` | Link each candidate to its artifact-mask evidence. |

Under v2 the display loci changed from {HLA-E, MCAM} to {TFE3, MCAM}; the validation
targets are CTNNA2, SDK1 and TCF7L1.

Dependency: `07_chromatin/figure5_*.pbs` → this stage →
`01_callsets_and_robustness/run_downstream_callset.pbs` for the artifact-masked arm
(PBS `afterok`).
