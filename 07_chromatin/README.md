# 07 — chromatin association and burden controls

All reference-chromatin analyses, plus the chromatin figure generators. Feeds
**Figure 4**, **Figure 5**, **Figures S8, S9** and **Tables S17, S18, S19**.

The figure scripts live here rather than in `10_figures/` because they import
`run_chromatin_hg38_analysis` as a library via `sys.path.insert(<own directory>)` and
must stay adjacent to it.

## Analysis

| File | Role |
|---|---|
| `fetch_encode_hg38_peaks.py` | Retrieve the ENCODE hg38 peak sets. Run where there is outbound network — compute nodes have none. |
| `run_chromatin_hg38_analysis.py` + `chromatin_hg38_analysis.pbs`, `run_qduh_chromatin_hg38.sh` | The core analysis in a single hg38 coordinate system: eccDNA BEDs stay in hg38, legacy hg19 ChIP-seq peaks are lifted over. Uses the exact chromosome- and length-matched placement expectation. |
| `run_chromatin_extensions.py` + `chromatin_extensions.pbs`, `v2_chromatin_extensions.pbs` | Cross-cell-type, mappability and peak-centred-density extensions. Three stages, submitted one per job (`qsub -v` truncates comma-separated values). |
| `run_p0_burden_control.py` + `p0_burden_control.pbs`, `v2_p0_burden_control.pbs` | Burden controls: within-group enrichment, the age-matched subset, HC3 covariate adjustment, depth downsampling, and the downsample-vs-full agreement check. |
| `validate_extensions_probabilities.py`, `validate_p0_per_circle_arrays.py` | Independent recomputation of the probability arrays. |
| `rerun_chromatin_covariate_locally.py` | Local re-run of the covariate-adjusted table. |
| `build_chromatin_supplementary_tables.py` | Format Tables S17–S19. |
| `figure5_statistical_de.pbs`, `figure5_hlae_mcam.pbs`, `figure5_krt1_gp5.pbs` | Locus-level jobs feeding Figure 5. |

## Figures

| File | Produces |
|---|---|
| `make_figure4new.py` | Figure 4 |
| `make_main_figures.py`, `make_locus_panels.py` | Shared plotting library for the locus panels |
| `make_figure5_reference_layout.py`, `make_figure5_revised.py` | Figure 5 layout |
| `make_figure5_locked_abc_reference_de.py`, `make_figure5_locked_abc_statistical_de.py` | Figure 5 panels a–c (locked, approved) |
| `make_figure_s5.py` | Figure S9 (**note**: the script's own output name is `figure_s5_crosscell_mappability`) |
| `make_figure_chromatin_merged.py` | Merged chromatin overview |
| `run_p0_burden_control.py` | Figure S8 (**note**: writes `Figure_S4`) |

## Superseded

An earlier relative-enrichment score was used for the COVID-19 vs HC histone-mark
comparisons. It was replaced during revision by the exact chromosome- and length-matched
placement expectation. Anything computed with the old score does not support any
statement in the manuscript.
