# Call sets, and what changed between the two analysis rounds

English summary of [`PRE_DECLARATION.md`](PRE_DECLARATION.md), which is the archived
original (Chinese) and was written and frozen **before any recompute job was submitted**.

## The declared filters

The manuscript Methods state that all downstream analyses use high-confidence calls
meeting all five of the following:

| # | Filter |
|---|---|
| i | ≥ 3 split reads supporting the junction |
| ii | Circle-Map score > 200 |
| iii | mean per-base coverage across the circle exceeds its standard deviation |
| iv | ≥ 0.3 coverage increase at both the start and the end |
| v | < 10% uncovered bases |

## v1 and v2

Applying those five filters line by line to the 16,857,160 archived Circle-Map records
leaves 10,099,219 — that is, the originally submitted downstream analyses had run on the
**unfiltered** archive, not on the set the Methods describe. All downstream analyses were
re-run on the filtered set.

| | v1 (original submission) | v2 (current, and what this repository contains) |
|---|---:|---:|
| primary call set | archived Circle-Map BED, filters not applied | `circlemap_methods` |
| total calls | 16,857,160 | 10,099,219 |
| COVID-19 | 15,006,670 | 8,830,916 |
| HC | 1,850,490 | 1,268,303 |
| median per-sample retention | 1.000 (by definition) | 0.652 |

Retention is 59.9% against the archive denominator, or 60.1% against the 16,807,187
records that parse cleanly.

The change also removed two internal inconsistencies: the remaining six call sets in the
robustness ladder were all rooted in the methods set while the primary set alone sat
outside it, and the support counts for the three wet-lab validation targets had always
been computed on the methods set.

## Seeds and resampling, carried over unchanged

Held fixed so that any v1/v2 difference is attributable to the call set, not to
resampling (`workflow/v2_env.sh`):

```
SEED_FRAGLEN=20260729     BOOTSTRAP_FRAGLEN=1000    BOOTSTRAP_RCA=5000
SUBSAMPLE_REPEATS=100     P0_REPEATS=100            P0_TARGETS=4000,20000
```

One forced change: the fragment-length equal-count subsampling target dropped from 4,000
to 2,500 unique eccDNAs per sample. Under the methods call set the smallest sample
(ZXS181) retains 2,598 calls, so 4,000 is arithmetically infeasible and the analysis
aborts. 2,500 is the largest round value every sample can supply. This is recorded as
pre-declaration amendment A1 — a feasibility constraint, not a post-hoc choice.

## The seven call sets

| Call set | Definition |
|---|---|
| `circlemap_current` | archived Circle-Map output, no filters (the v1 primary; retained as the permissive sensitivity arm) |
| `circlemap_methods` | the five declared filters applied (**the v2 primary**) |
| strict split/score ×2 | split reads ≥ 10 with score > 1,000; split reads ≥ 20 with score > 2,000 |
| `circlemap_artifact_masked` | additionally removes circles with any full-circle or ±100 bp breakpoint-window overlap with the ENCODE blacklist or assembly gaps |
| `circlemap_high_support_masked` | high-support calls, artifact-masked |
| `consensus_t10` | Circle-Map / Circle_finder two-caller consensus at 10 bp breakpoint tolerance |

## Re-implemented analyses

Three analyses had no surviving original code — neither in the working tree nor on the
cluster — and were re-implemented from the Methods text under pre-declaration amendments
A3/A4. Each re-implementation carries a docstring stating exactly what it assumed and
where it may differ from the original:

| Analysis | Script | Amendment |
|---|---|---|
| chromosome distribution and gene-element O/E (Figure 2b, 2d) | `02_eccdna_quantification/genomic_distribution.py` | A3/A4 |
| recurrent COVID-19-specific exact intervals (Figure 3d) | `02_eccdna_quantification/recurrent_exact_intervals.py` | A3 |
| clinical correlates (Figure S5, Tables S9/A1–A6) | `09_clinical_correlates/clinical_correlates_v2.py` | A3/A4 |

Amendment A4 records two definitional points where the original code differed from the
Methods: whether an eccDNA is assigned by its **start coordinate** or by **interval
intersection**, and how the genome-occupied denominator is formed. The re-implementations
follow the Methods.

## v1 guard values

Several scripts hard-code v1 results as fatal assertions. This is deliberate: it is how
the v2 differences were detected. In the patched copies each guard became configurable,
and the observed v2 value is recorded:

| Script | Guard | v1 | v2 |
|---|---|---|---|
| `10_figures/create_eccgene_revision_figures.py` | enrichment set sizes | 3,094 / 26,642 | 3,189 / 26,087 |
| `10_figures/create_figure_s1_redesigned.py` | stable peak positions | 196 / 365 / 571 bp | 196 / 366 / 571 bp |
| `06_candidate_loci/annotate_chromatin_display_candidates.py` | display loci | HLA-E, MCAM | TFE3, MCAM |
| `10_figures/create_circlemap_robustness_figures.py` | validation targets | CTNNA2, CAB39, chr22_target, SDK1 | CTNNA2, SDK1, TCF7L1 |

The last one is the most easily missed: a hard-coded target list had silently dropped
TCF7L1. Promoting it to a module-level `VALIDATION_TARGETS` constant made the figure, its
caption and the wet-lab plan use the same three loci.
