# 02 — eccDNA burden and genomic distribution

Normalized abundance and where eccDNAs fall in the genome. Feeds **Figures 1b, 1c, 1e,
2b, 2d, 3d** and **Tables S10, S11, S12**.

| File | Role |
|---|---|
| `run_epm_reanalysis.py` | Per-sample eccDNA counts and EPM (eccDNA per million mapped reads); repeat-class ratios. The EPM denominator is total mapped alignments from `samtools idxstats`. |
| `calc_epm_in_peaks_fast.py` | EPM restricted to a peak set, using the same global denominator as the study-wide EPM. |
| `genomic_distribution.py` + `v2_genomic_distribution.pbs` | Per-sample chromosome distribution and gene-element observed/expected. **Re-implemented** from the Methods under pre-declaration amendments A3/A4; the original code differed on start-coordinate vs interval-intersection assignment and on the genome-occupied denominator. |
| `recurrent_exact_intervals.py` + `v2_recurrent_exact_intervals.pbs` | Recurrent COVID-19-specific exact intervals (Figure 3d input). **Re-implemented** under amendment A3, following the Figure 3d legend verbatim. |
| `export_supplementary_tables.py` | Emit the tables above in supplementary format. |

Run via `01_callsets_and_robustness/run_downstream_callset.pbs` to execute against a
named call set.
