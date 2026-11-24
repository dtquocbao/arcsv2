#!/bin/bash
#SBATCH --job-name=arcsv2-cov
#SBATCH --account=132702112504
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-40
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# Backfill coverage for the 40 ML-0 samples from a CRAM URL. One core per
# sample. Reads only the windows ml-coverage asks for (samtools region queries
# through the Grace proxy). Does not rewrite merged.tab.
#
# Do not submit until the pilot (step 3) has been reviewed. Step 4 of the
# ML2 task is the run that reports wall time.
#
#   cd $SCRATCH/arcsv2/repo && /usr/bin/sbatch scripts/hpc/backfill_coverage.sbatch

set -euo pipefail
REPO="${ARCSV2_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${REPO}/scripts/hpc/run_label.sh" ]] || REPO="${SCRATCH}/arcsv2/repo"
SIF="${ARCSV2_SIF:-${SCRATCH}/arcsv2/sif/arcsv2_1.6.1.sif}"
SAMPLES="${ARCSV2_SAMPLES:-${REPO}/workflows/config/samples.ml.all40.tsv}"
MERGED_DIR="${ARCSV2_MERGED_DIR:-${REPO}/export/ml0/merged}"
REF="${ARCSV2_REFERENCE:-${SCRATCH}/arcsv2/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa}"
GAPS="${ARCSV2_GAPS:-${SCRATCH}/arcsv2/ref/gaps.bed}"
OUTDIR="${ARCSV2_COVERAGE_DIR:-${SCRATCH}/arcsv2/ml/coverage_backfill}"
mkdir -p "${OUTDIR}"

export http_proxy="${http_proxy:-http://10.73.132.63:8080}"
export https_proxy="${https_proxy:-http://10.73.132.63:8080}"
export PYTHONPATH="${REPO}/python${PYTHONPATH:+:${PYTHONPATH}}"

singularity exec -B /scratch "${SIF}" \
  python "${REPO}/scripts/hpc/backfill_coverage.py" \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --samples "${SAMPLES}" \
  --merged-dir "${MERGED_DIR}" \
  --reference "${REF}" \
  --gaps "${GAPS}" \
  --out-dir "${OUTDIR}"
