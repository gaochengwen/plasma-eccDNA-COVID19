# 12 — supplementary workbook assembly

Writes the v2 results back into the 21-sheet supplementary workbook, preserving titles,
notes, headers, column widths and number formats. Rows are matched by header name, in
the sheet's own column order; columns no longer in the v2 data model are removed
explicitly, and unmapped cells are left blank rather than carrying stale values.

| File | Sheets |
|---|---|
| `build_table_s3_v2.py` | S3 — the master statistics table: eight named FDR families, each BH-corrected within itself, two-sided Wilcoxon rank-sum with Cliff's δ |
| `build_supplementary_v2.py` | S3, S10, S12, S13 |
| `build_s16_s20_v2.py` | S15, S16, S20 |
| `rebuild_table_s14_v2_20260813.py` | S14 — the only sheet whose row count changes (3,402 → 3,487 significant genes) |
| `repair_v1_supplementary_tables_20260813.py` | S2a, S5b–e, S17b–g, S18a/c/e/f, S19b/d — blocks a cell-by-cell audit found still on v1 |

**Order matters.** `repair_v1_supplementary_tables_20260813.py` must run before
`08_age_matching/rebuild_age_matched_pairs_v2_20260813.py`, which reads Table S2a.

Sheets deliberately left on their original values, because the call-set migration does
not reach them: **S11** (repeat classes — computed from BAM reads, identical under both
call sets) and **S17a / S19a / S19c** (peak-set provenance and mappability, independent
of the eccDNA call set).
