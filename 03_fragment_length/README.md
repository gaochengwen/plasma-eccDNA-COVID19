# 03 — fragment-length structure

Discovery-based peak calling on the eccDNA length distribution. Feeds **Figure 1d**,
**Figure S2** and **Table S5**.

| File | Role |
|---|---|
| `run_fragment_length_peak_analysis.py` | The analysis engine. Deduplicates per sample by chromosome/start/end, restricts to unique eccDNAs of 50–1,000 bp, bins at 1 bp, smooths with a Gaussian kernel (σ = 15 bp, reflect mode, 4σ truncation), renormalizes to unit area and averages sample curves with equal weight. Neither the number nor the location of peaks is assumed. |
| `run_fragment_length_peak_analysis.pbs`, `v2_fragment_length.pbs` | Job wrappers. |
| `rebuild_figure_s2_panel_b.py` | Redraw the peak-robustness panel; imports the engine as a library. |

Seed `20260729`, 1,000 bootstrap resamples, 100 subsampling repeats. The equal-count
subsampling target is **2,500** unique eccDNAs per sample under v2, down from 4,000: the
smallest sample retains 2,598 calls after filtering, so 4,000 is infeasible. Recorded as
pre-declaration amendment A1.
