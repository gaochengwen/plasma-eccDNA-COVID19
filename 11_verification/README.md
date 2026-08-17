# 11 — verification

Run last. These scripts exist because the failure mode that matters is a **partial**
migration: a number corrected in the manuscript but not in the tables, or in one
response letter but not the other two.

| File | Role |
|---|---|
| `verify_v2_consistency.py` | Cross-check every v2 number across the manuscript, the three response letters and the locked tables. Must report PASS before anything ships. |
| `v2_manuscript_numbers.py` | Compute every Results number the manuscript prints, from the v2 locked tables. |
| `extract_key_numbers.py` | Extract every manuscript-facing number beside its v1 value. |
| `verify_workbook_against_sources.py` | Check every data block of the supplementary workbook against its locked TSV. |
| `audit_tcf7l1.py` | Audit the TCF7L1 circle in every original Circle-Map BED, with CTNNA2 and SDK1 as controls. |
| `validate_circlemap_revision.py`, `validate_chromatin_consistency.py`, `validate_rca_revision.py` | Per-section checks: each assertion recomputes the value from the workbook or reads it out of the figure's SVG rather than trusting the prose. |
| `docx_tracked_edit.py`, `docx_table_tools.py` | Helper modules the three validators import. Kept adjacent because they are resolved by `sys.path.insert(<own directory>)`. |

The three `validate_*` scripts read the manuscript `.docx` and the supplementary
workbook, which are not part of this repository. They are included because they document
exactly which numbers were checked and how.
