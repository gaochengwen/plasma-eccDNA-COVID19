# 05 — eccGene assignment and endpoints

An eccGene is annotated whenever a 1-bp Circle-Map junction coordinate falls within a
merged UCSC hg38 refGene gene body. It is a **positional label**: it does not imply that
an intact or functional gene is carried by the circle. Feeds **Figures 3a, 3b, S6** and
**Table S13**.

| File | Role |
|---|---|
| `run_eccgene_reanalysis.py` | Abundance and detection-frequency endpoints under three assignment definitions (junction, interval, midpoint). Of 28,278 merged gene bodies on canonical chromosomes, 26,087 carry at least one junction and are tested. |
| `add_enrichment_materials.py` | Assemble the g:Profiler over-representation request: query genes, the custom background, the exact request and the unmodified response. |
| `build_eccgene_results_package.py`, `validate_eccgene_results_package.py` | Package the results and check the package against its sources. |

**Figure 3c** plots a Metascape run made through the Metascape web service. Its
term-level export is not scriptable here; the g:Profiler confirmation reported alongside
it is fully re-derivable from `add_enrichment_materials.py`.

Run via `01_callsets_and_robustness/run_downstream_callset.pbs` for each call set.
