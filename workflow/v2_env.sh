#!/usr/bin/env bash
# Shared environment for the v2 rerun (primary call set = circlemap_methods).
#
# v1 (archived)  primary call set = the un-filtered archived Circle-Map BEDs
#                under /gpfs/data/gao/Covid-eccDNA/{covid,normal}.
# v2 (this run)  primary call set = circlemap_methods, i.e. the calls that
#                satisfy all five filters declared in the manuscript Methods:
#                  (i)   >= 3 split reads
#                  (ii)  Circle-Map score > 200
#                  (iii) mean per-base coverage > its standard deviation
#                  (iv)  >= 0.3 coverage increase at both start and end
#                  (v)   < 10% uncovered bases
#
# Nothing here rebuilds a call set. circlemap_methods and the strictness ladder
# derived from it already exist and are verified by SHA-256 in stage S0.

set -euo pipefail

export CF_ROOT=/gpfs/data/gao/Covid-eccDNA/Revise/Circle_finder
export V1_CHROMATIN=/gpfs/data/gao/Covid-eccDNA/Revise/chromatin
export V1_BCL3=/gpfs/data/gao/Covid-eccDNA/Revise/BCL3
export V1_RCA=/gpfs/data/gao/Covid-eccDNA/Revise/RCA
export V1_EPM=/gpfs/data/gao/Covid-eccDNA/Revise/EPM

# Symlink adapter that presents the methods call set with the directory layout
# the analysis scripts expect ({covid,normal}/<sample>_circle_site.bed).
export ADAPTER=$CF_ROOT/work/downstream_projects/circlemap_methods

# Everything this rerun produces lands under V2ROOT. v1 outputs are never
# overwritten.
export V2ROOT=$CF_ROOT/downstream/circlemap_methods
export V2_ARTIFACT=$CF_ROOT/downstream/circlemap_artifact_masked

export PYTHON_BIN=/public/opt/miniforge3_2024/bin/python
export SAMTOOLS_19=/opt/samtools-1.9/samtools

# Seeds are carried over from v1 unchanged so that any difference between the
# two runs is attributable to the call set and not to resampling.
export SEED_FRAGLEN=20260729
export BOOTSTRAP_FRAGLEN=1000
export BOOTSTRAP_RCA=5000

# Equal-count subsampling target for the fragment-length analysis.
# v1 used 4,000 unique eccDNAs per sample. Under the methods call set the
# smallest sample (ZXS181) retains 2,598, so 4,000 is arithmetically
# infeasible and the analysis aborts. 2,500 is the largest round value that
# every sample can supply. This is a forced feasibility change, not a
# post-hoc choice; see PRE_DECLARATION.md amendment A1.
export SUBSAMPLE_COUNT=2500
export SUBSAMPLE_REPEATS=100
export P0_REPEATS=100
export P0_TARGETS=4000,20000

v2_prepare_dirs() {
    mkdir -p "$V2ROOT"/{fragment_length,RCA,BCL3,logs}
    # run_chromatin_extensions.py resolves --umap relative to <outdir>/reference.
    # The core downstream job copied only four reference files, so link the rest.
    ln -sfn "$V1_CHROMATIN/reference/k36.umap.bed.gz" \
            "$V2ROOT/chromatin/reference/k36.umap.bed.gz"
}
