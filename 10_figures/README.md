# 10 — figure generation

Non-chromatin figures. The chromatin figure generators live in `07_chromatin/` because
they import that module as a library.

**Read `docs/figure_source_map.tsv` before running anything here.** Supplementary
numbering changed twice during revision, so a script's own output filename is often not
the figure it produces — `create_figure_s2.py` makes Figure S3, `create_eccgene_revision_figures.py`
makes Figure S6, and so on.

## Harnesses — use these, not the generators directly

| File | Role |
|---|---|
| `regen_figure_v2.py` | Re-render a submitted figure from the v2 tables behind a fidelity gate: render with the v1 tables first and require that the submitted file is reproduced, only then render with v2. A framework that cannot reproduce the submitted figure is not allowed to produce a replacement. |
| `regen_create_script.py` | The same contract for the `create_*.py` generators, which take no arguments and hard-code their paths. Holds the recipe table: which module constants to rebind, and which output file corresponds to which submitted figure. |
| `regen_figure_s2.py` | Figure S2 specifically; carries its own pixel-level gate. |
| `../workflow/protected_run.py` | For subprocess-style generators: fingerprints 89 protected files before the run and restores anything overwritten. |

Gates come in three levels, recorded per figure: **byte-level** (PDF content md5 after
stripping timestamps), **pixel-level** (mean grey difference < 0.6 at 110 dpi), and
**content-level** (all extracted text identical after normalization). Content-level is
used where the submitted PDF was produced by matplotlib 3.11.0, which is present in
neither environment, making byte-level reproduction unattainable.

## Generators

| File | Produces (submitted numbering) |
|---|---|
| `rebuild_illustrator_panels_v2.py` | Replacement data panels for Figures 1b, 1c, 1e, 2b, 2d, 3a, 3b |
| `build_illustrator_panel_data.py` | Their source data, v1 beside v2 |
| `rebuild_figure_3d_panel_v2.py`, `build_figure_3d_panel_data_v2.py` | Figure 3d panel and its data |
| `create_figure_s1_redesigned.py` | Figure S1 |
| `create_figure_s2.py` | Figure **S3** |
| `create_figure_s4.py` | Figure S4 |
| `create_eccgene_revision_figures.py` | Figure **S6** |
| `create_circlemap_robustness_figures.py` | Figure S7 |
| `rerun_figure_s8_locally.py` | Figure S8 |
| `create_figure_s9_clinical_correlates.py`, `build_figure_s5_v2.py` | Figure **S5** |
| `create_figure3_ab_replica.py`, `create_figure_3d_frequency_sorted.py` | Figure 3 working panels |
| `create_graphical_abstract.py` | Graphical abstract |
| `check_figure_nature_spec.py` | Check dimensions, fonts and resolution against the journal specification |

Figures 1, 2 and 3 were assembled by hand in Adobe Illustrator 29.0. Their page layout
has no scripted generator; only the data panels listed above are scripted.
