# 00 — sequencing quality control

Per-sample alignment statistics and the sequencing-covariate models behind
**Table S2**. Independent of the call-set choice: it runs on BAM files, which are
identical under v1 and v2.

| File | Role |
|---|---|
| `run_seq_qc.sh` | Driver; submits the array job. |
| `seq_qc_array.pbs` | PBS array: `samtools flagstat` and duplicate marking per sample. |
| `collect_seq_qc.py` | Collate per-sample output into the QC table. |
| `run_seq_qc_covariate_analysis.py` | Test whether cohort differences in burden track sequencing depth, duplication or mapping rate. |

Order: `run_seq_qc.sh` → `collect_seq_qc.py` → `run_seq_qc_covariate_analysis.py`.
