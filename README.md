# Plasma eccDNA in COVID-19 — analysis code

Analysis code for the manuscript

> **Plasma eccDNA profiling in COVID-19 reveals higher abundance and
> immune-associated genomic patterns**
> Manuscript ID **JARE-D-26-04883**

This repository contains **code only** — the analysis, job-submission and
quality-control scripts behind every number reported in the manuscript and its
supplement. The manuscript figure generators and the supplementary-workbook
assembly scripts are not deposited. There is no sequencing data, no per-sample
eccDNA call set and no result table here; see [Data availability](#data-availability).

---

## Contents

| Directory | Files | Contents |
|---|---:|---|
| [`workflow/`](workflow/) | 4 | Shared environment, the stage-2 submission orchestrator and the marker-waiting helpers |
| [`00_sequencing_qc/`](00_sequencing_qc/) | 4 | Per-sample `flagstat`, duplicate marking and sequencing-covariate models |
| [`01_callsets_and_robustness/`](01_callsets_and_robustness/) | 29 | Circle_finder second caller, the seven-call-set strictness ladder, consensus building, artifact masking and the robustness synthesis |
| [`02_eccdna_quantification/`](02_eccdna_quantification/) | 7 | EPM burden, chromosome distribution, gene-element observed/expected, and recurrent COVID-19-specific exact intervals |
| [`03_fragment_length/`](03_fragment_length/) | 3 | Objective fragment-length peak calling, bootstrap stability and subsampling robustness |
| [`04_rca_technical_bias/`](04_rca_technical_bias/) | 3 | Post-RCA yield, its correlates and the HC3 covariate models |
| [`05_eccgene/`](05_eccgene/) | 4 | eccGene assignment under three definitions, abundance and detection endpoints, over-representation materials |
| [`06_candidate_loci/`](06_candidate_loci/) | 9 | Prespecified eligibility screen and ranking used to choose the validation loci |
| [`07_chromatin/`](07_chromatin/) | 13 | hg38 liftOver, the chromosome- and length-matched placement expectation, cross-cell-type and mappability extensions, and burden controls |
| [`08_age_matching/`](08_age_matching/) | 2 | The age-matched subset and its caliper sensitivity analysis |
| [`09_clinical_correlates/`](09_clinical_correlates/) | 3 | Clinical annotation, correlates, strata and outcome analyses |
| [`docs/`](docs/) | — | Table source map, call-set definitions, reference-file manifest, provenance and the pre-analysis declaration |
| [`env/`](env/) | — | Interpreter and package versions for both execution environments |

[`docs/table_source_map.tsv`](docs/table_source_map.tsv) maps every supplementary
table to the analysis that produces it. [`docs/provenance.tsv`](docs/provenance.tsv)
records where each file came from in the working analysis tree, with a SHA-256.

---

## The two call sets

The Methods declare that all downstream analyses use high-confidence calls satisfying
five filters simultaneously: ≥ 3 split reads; Circle-Map score > 200; mean per-base
coverage exceeding its standard deviation; ≥ 0.3 coverage increase at both breakpoints;
< 10% uncovered bases.

During revision it was found that the originally submitted downstream analyses had run
on the **unfiltered** archived Circle-Map output. Every downstream analysis was
re-run on the call set that actually satisfies those filters (`circlemap_methods`,
10,099,219 of 16,857,160 archived calls). **The code here is that state — the one that
backs the manuscript.**

Some scripts hard-code first-round results as fatal assertions. That is deliberate: it
is how the differences were found. In the shipped copies each guard is configurable and
records the observed value; the guard table is in
[`docs/call_sets.md`](docs/call_sets.md), and the declaration written before any job was
submitted is archived at [`docs/PRE_DECLARATION.md`](docs/PRE_DECLARATION.md).

---

## Execution order

Upstream Circle-seq calling (FastQC → fastp → BWA-MEM → samtools → Circle-Map Realign)
ran before this repository's scope and is described in the Methods; its call sets are
the input to stage 00.

```
workflow/v2_env.sh                 # shared paths, seeds and resampling counts
        │
        ├── 00_sequencing_qc       # independent
        │
        └── workflow/submit_v2_stage2.sh
                │                  # preflight-gates every input, then submits:
                ├── branch A  07_chromatin  (extensions ×3 stages, burden controls)
                ├── branch B  03_fragment_length, 04_rca_technical_bias
                └── branch C  01_callsets_and_robustness (artifact-masked downstream)
                                  └── 06_candidate_loci (ranking; PBS afterok)

02_eccdna_quantification, 05_eccgene   # run on each call set via
                                       # 01_.../run_downstream_callset.pbs
09_clinical_correlates                 # local, from the locked burden tables
08_age_matching                        # reads the per-sample burden table
```

Branches A, B and C are independent; the only cross-branch dependency is
`rank_candidate_loci` → artifact-masked downstream, expressed with PBS `afterok`. Each
job writes its own log and a `*.done` marker.

Long jobs are submitted with `qsub` only. Compute nodes have no outbound network, so
every reference file must be on shared storage before submission; each job gates on its
inputs and fails fast rather than writing half a result set.

---

## Environment

Both environments are recorded in
[`env/software_versions.tsv`](env/software_versions.tsv):

- **HPC (PBS cluster)** — Python 3.10.13, NumPy 2.2.6, SciPy 1.13.1, Matplotlib 3.9.0.
  All compute-heavy analyses.
- **Local workstation** — Python 3.9.25, NumPy 2.0.2, SciPy 1.13.1, pandas 2.3.3,
  statsmodels 0.14.6, openpyxl 3.1.5.

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r env/requirements.txt   # or: conda env create -f env/environment.yml
```

### Third-party software (not redistributed)

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
| R | 4.1.0 | https://www.r-project.org/ |
| GraphPad Prism | 10.4.1 | commercial |

Reference files are not redistributed either. Their download URLs, byte sizes and
SHA-256 checksums are in
[`docs/reference_manifest.tsv`](docs/reference_manifest.tsv);
`07_chromatin/fetch_encode_hg38_peaks.py` retrieves the ENCODE peak sets.

---

## Data availability

Patient-derived sequence data — the aligned BAM files and the paired FASTQ records
derived from those alignments — are deposited in the Genome Sequence Archive under
BioProject accession **PRJCA073440**. Per-sample eccDNA BED call sets are not
distributed here. Processed results sufficient to re-derive every reported statistic
are published as Supplementary Tables S1–S20 with the manuscript.

Sample identifiers in the code (`NS*`, `SP*`, `XS*`, `ZXS*`) are de-identified study
codes carrying no personal information; they match Supplementary Table S1.

Absolute paths of the form `/gpfs/data/gao/Covid-eccDNA/...` are the cluster locations
the analyses read from. They are kept verbatim so that the deposited code is
byte-identical to the code that produced the results.

---

## Verifying this repository

```bash
sha256sum -c checksums.sha256
```

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the manuscript when using this code.

## License

[MIT](LICENSE). Third-party tools are distributed under their own licenses and are not
included here.
