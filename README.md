# Plasma eccDNA in COVID-19 — analysis code

Analysis code for the manuscript

> **Plasma eccDNA profiles in COVID-19: increased abundance and associations with
> genomic and reference-chromatin features**
> Manuscript ID **JARE-D-26-04883**

This repository contains **code only**: the analysis, job-submission, quality-control
and figure-generation scripts that produce every number, table and figure reported in
the manuscript and its supplement. It contains no sequencing data, no per-sample eccDNA
call sets and no result tables — see [Data availability](#data-availability).

---

## Contents

| Directory | Files | Contents |
|---|---:|---|
| [`workflow/`](workflow/) | 5 | Shared environment, the stage-2 submission orchestrator, and the write-protection harness used when re-running figure generators |
| [`00_sequencing_qc/`](00_sequencing_qc/) | 4 | Per-sample `flagstat`, duplicate marking and sequencing-covariate models (Table S2) |
| [`01_callsets_and_robustness/`](01_callsets_and_robustness/) | 30 | Circle_finder second caller, the seven-call-set strictness ladder, consensus building, artifact masking and the robustness synthesis |
| [`02_eccdna_quantification/`](02_eccdna_quantification/) | 7 | EPM burden, chromosome distribution, gene-element observed/expected, and recurrent COVID-19-specific exact intervals |
| [`03_fragment_length/`](03_fragment_length/) | 4 | Objective fragment-length peak calling, bootstrap stability and subsampling robustness |
| [`04_rca_technical_bias/`](04_rca_technical_bias/) | 3 | Post-RCA yield, its correlates and the HC3 covariate models |
| [`05_eccgene/`](05_eccgene/) | 4 | eccGene assignment under three definitions, abundance/detection endpoints, over-representation materials |
| [`06_candidate_loci/`](06_candidate_loci/) | 10 | Prespecified eligibility screen and ranking used to choose display and validation loci |
| [`07_chromatin/`](07_chromatin/) | 26 | hg38 liftOver, chromosome- and length-matched placement expectation, cross-cell-type and mappability extensions, burden controls, and the chromatin figure generators (they import the analysis module as a library, so they must stay beside it) |
| [`08_age_matching/`](08_age_matching/) | 2 | The age-matched subset and its sensitivity analysis |
| [`09_clinical_correlates/`](09_clinical_correlates/) | 3 | Clinical annotation, correlates, strata and outcome analyses |
| [`10_figures/`](10_figures/) | 19 | Non-chromatin figure generators and the fidelity-gated regeneration harnesses, plus the replacement data panels for the three hand-assembled main figures |
| [`11_verification/`](11_verification/) | 10 | Cross-checks that every reported number agrees between tables, figures, manuscript and response letters |
| [`12_supplementary_tables/`](12_supplementary_tables/) | 5 | Assembly of the 21-sheet supplementary workbook from the locked tables |
| [`docs/`](docs/) | — | Source maps, call-set definitions, reference-file manifest, provenance and the pre-analysis declaration |
| [`env/`](env/) | — | Interpreter and package versions for both execution environments |

Every file records where it came from in the working analysis tree, with a SHA-256,
in [`docs/provenance.tsv`](docs/provenance.tsv).

---

## The two call sets

The manuscript Methods declare that all downstream analyses use high-confidence calls
satisfying five filters simultaneously:

1. ≥ 3 split reads supporting the junction
2. Circle-Map score > 200
3. mean per-base coverage across the circle exceeding its standard deviation
4. ≥ 0.3 coverage increase at both the start and the end
5. < 10% uncovered bases

During revision it was found that the originally submitted downstream analyses had run
on the **unfiltered** archived Circle-Map output. All downstream analyses were therefore
re-run on the call set that actually satisfies the declared filters
(`circlemap_methods`, 10,099,219 of 16,857,160 archived calls).

**The code in this repository is the v2 state — the state that backs the manuscript.**
Analysis engines reused verbatim from the first round are included unchanged; where an
engine needed a patch for v2 the patched copy is the one shipped. The pre-analysis
declaration written before any v2 job was submitted is archived at
[`docs/PRE_DECLARATION.md`](docs/PRE_DECLARATION.md) (original Chinese); an English
summary of the two call sets and the resulting differences is in
[`docs/call_sets.md`](docs/call_sets.md).

---

## Execution order

The upstream Circle-seq calling (FastQC → fastp → BWA-MEM → samtools → Circle-Map
Realign) was run before this repository's scope and is described in the Methods; the
call sets it produced are the input to stage 00 below.

```
workflow/v2_env.sh                 # shared paths, seeds and resampling counts
        │
        ├── 00_sequencing_qc       # independent; produces Table S2
        │
        └── workflow/submit_v2_stage2.sh
                │                  # preflight-gates every input, then submits:
                ├── branch A  07_chromatin  (extensions ×3 stages, burden controls)
                ├── branch B  03_fragment_length, 04_rca_technical_bias
                └── branch C  01_callsets_and_robustness (artifact-masked downstream)
                                  └── 06_candidate_loci (ranking; PBS afterok)
                                            └── 07_chromatin figure5_*.pbs

02_eccdna_quantification, 05_eccgene   # run on each call set via
                                       # 01_.../run_downstream_callset.pbs
09_clinical_correlates                 # local, from the locked burden tables
12_supplementary_tables                # assemble the 21-sheet workbook
08_age_matching                        # AFTER stage 12: reads Table S2a
10_figures                             # local, from the locked tables
11_verification                        # last; must pass before anything ships
```

Branches A, B and C are independent. The only cross-branch dependency is
`figure5_*` → `rank_candidate_loci` → artifact-masked downstream, expressed with PBS
`afterok`. Every job writes its own log and a `*.done` marker containing the number of
tables it produced.

Long jobs are submitted with `qsub` only. No analysis is run inside an interactive SSH
session, `screen`, `tmux` or `nohup`. Compute nodes have no outbound network, so every
reference file must exist on shared storage before submission; each job gates on its
inputs and fails fast rather than writing half a result set.

---

## Environment

Two environments were used and both are recorded in
[`env/software_versions.tsv`](env/software_versions.tsv):

- **HPC (PBS cluster)** — Python 3.10.13, NumPy 2.2.6, SciPy 1.13.1, Matplotlib 3.9.0,
  on Linux 3.10.0 / glibc 2.17. All compute-heavy analyses.
- **Local workstation** — Python 3.9.25, NumPy 2.0.2, SciPy 1.13.1, Matplotlib 3.9.4,
  pandas 2.3.3, statsmodels 0.14.6, openpyxl 3.1.5. Figure regeneration, auditing and
  consistency verification.

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r env/requirements.txt
# or
conda env create -f env/environment.yml
```

### Third-party software (not redistributed)

Install these separately; the versions used are pinned in
`env/software_versions.tsv`.

| Tool | Version | Source |
|---|---|---|
| FastQC | 0.11.3 | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| fastp | 0.23.x | https://github.com/OpenGene/fastp |
| BWA | 0.7.17-r1188 | https://github.com/lh3/bwa |
| samtools | 1.9 and 1.17 | https://github.com/samtools/samtools |
| bedtools | 2.28.0 | https://github.com/arq5x/bedtools2 |
| samblaster | 0.1.20 | https://github.com/GregoryFaust/samblaster |
| Circle-Map | Realign | https://github.com/iprada/Circle-Map |
| Circle_finder | commit `3eb333db2ea6277dde36cbf640be9afeb710c717` | https://github.com/pk7zuva/Circle_finder |
| liftOver | UCSC | https://hgdownload.soe.ucsc.edu/admin/exe/ |
| GraphPad Prism | 10.4.1 | commercial |
| R | 4.1.0 | https://www.r-project.org/ |

Reference files are **not** redistributed either. Their download URLs, byte sizes and
SHA-256 checksums are in [`docs/reference_manifest.tsv`](docs/reference_manifest.tsv);
`07_chromatin/fetch_encode_hg38_peaks.py` retrieves the ENCODE peak sets.

---

## Finding the code behind a figure or table

Figure and supplementary numbering changed twice during revision, so **a script's own
output filename is not a reliable guide to which submitted figure it makes**. The
verified mapping is in:

- [`docs/figure_source_map.tsv`](docs/figure_source_map.tsv) — every main and
  supplementary figure, its generator, and how the mapping was confirmed
- [`docs/table_source_map.tsv`](docs/table_source_map.tsv) — every supplementary table
  and the script that builds it

---

## Known limits of scripted reproduction

These are stated so that no reader mistakes a gap for an error.

- **Figures 1, 2 and 3** were assembled by hand in Adobe Illustrator 29.0 and have no
  scripted generator for the page layout. The seven data panels (1b, 1c, 1e, 2b, 2d,
  3a, 3b) *are* scripted, in `10_figures/rebuild_illustrator_panels_v2.py`, which
  redraws them to the submitted geometry. Figure 1d comes from the fragment-length
  module directly.
- **Figure S5** cannot be rebuilt end to end. Its generator only reads the
  pre-computed clinical tables A2–A5; the original code that produced A1–A6 exists
  neither in this repository nor on the cluster. `09_clinical_correlates/clinical_correlates_v2.py`
  is a documented re-implementation from the Methods, written under pre-declaration
  amendments A3/A4, and is the authoritative generator going forward.
- **Figure 3c** plots a Metascape run made through the Metascape web service; its
  term-level export is not scriptable here. The g:Profiler confirmation reported
  alongside it is fully re-derivable from `05_eccgene/add_enrichment_materials.py`.
- Several generators end by copying output to hard-coded absolute paths that later
  renumbering made wrong. `workflow/protected_run.py` and the `sandbox()` helpers in
  `10_figures/regen_figure_v2.py` and `regen_create_script.py` intercept those writes.
  **Re-run figure generators through those harnesses, not directly.**
- Some scripts hard-code v1 results as fatal assertions. This is deliberate — it is how
  the v2 differences were found. The patched copies make each guard configurable and
  record the observed v2 value; see the guard table in `docs/call_sets.md`.

---

## Data availability

Original raw FASTQ files were not retained. Patient-derived aligned BAM files and the
paired-read records reconstructed from them are held under controlled institutional
access. Per-sample eccDNA BED call sets are likewise patient-derived and are not
distributed here.

Processed results sufficient to re-derive every reported statistic are published as
Supplementary Tables S1–S21 with the manuscript.

Sample identifiers appearing in the code (`NS*`, `SP*`, `XS*`, `ZXS*`) are
de-identified study codes; they carry no personal information and match the codes used
in Supplementary Table S1.

Absolute paths of the form `/gpfs/data/gao/Covid-eccDNA/...` are the cluster locations
the analyses read from. They are retained verbatim rather than rewritten, so that the
deposited code is byte-identical to the code that produced the results.

---

## Verifying this repository

```bash
sha256sum -c checksums.sha256
```

---

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the manuscript when using this code.

## License

[MIT](LICENSE). Third-party tools listed above are distributed under their own
licenses and are not included here.
