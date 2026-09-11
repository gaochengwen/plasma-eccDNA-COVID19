# 09 — clinical correlates

**Figure S5** and **Tables S9 (A1–A6), S20**.

| File | Role |
|---|---|
| `clinical_correlates_v2.py` | The authoritative generator for sub-tables A1–A6 and the clinical block of Table S9. **Re-implemented** under pre-declaration amendments A3/A4: the original code exists neither in the working tree nor on the cluster, so the definitions are taken from the Methods and from the column structure of the v1 artefacts. Confidence intervals are 5,000-resample percentile bootstrap under a fixed seed, not a Fisher *z* transform — the v1 intervals do not reproduce under Fisher *z*, and the Methods specify percentile bootstrap for the parallel RCA correlations. |
| `derive_v2_burden_table.py` | Derive the per-sample v2 burden table from the locked robustness metrics. Not a new computation. |
| `build_table_s20.py` | Table S20, the non-COVID pneumonia comparator samples. |

Only the eccDNA burden columns differ between v1 and v2. Every clinical variable is
carried through unchanged, and `mapped_reads` is asserted identical row by row because
the BAM files are shared.

Note that outside lymphocyte count, the clinical covariates are available for the
COVID-19 group only, so cohort comparisons cannot be adjusted for them.
