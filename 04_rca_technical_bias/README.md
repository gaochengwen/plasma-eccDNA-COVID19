# 04 — rolling-circle amplification bias

Whether post-RCA yield explains the cohort difference in eccDNA burden. Feeds
**Figure S3**.

| File | Role |
|---|---|
| `run_rca_bias_analysis.py` | Compare RCA concentration between cohorts (Wilcoxon rank-sum with Cliff's δ), correlate it with detected eccDNA count, Circle-Map support events and EPM (Spearman, BH-corrected), and fit HC3 covariate models with group × RCA interaction contrasts. |
| `v2_rca_bias.pbs` | Job wrapper. |
| `rerun_module02_locally.py` | Local re-run for figure regeneration. |

5,000 bootstrap resamples under a fixed seed.
