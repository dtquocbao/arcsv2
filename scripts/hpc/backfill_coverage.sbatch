#!/bin/bash
#SBATCH --job-name=arcsv2-cov
#SBATCH --account=132702112504
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-40%4
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# Backfill coverage for the 40 ML-0 samples from a CRAM URL. One core per
# sample, at most four at a time. Reads only the windows ml-coverage asks for
# (samtools region queries through the Grace proxy). Does not rewrite merged.tab.
# A sample whose coverage.tsv already exists is skipped, so a resubmit of the
# full array retries only the failures (34 of 40 died mid-stream with
# "Container header CRC32 failure" / "retry_reconnect ... samtools bedcov").
#
#   cd $SCRATCH/arcsv2/repo && /usr/bin/sbatch scripts/hpc/backfill_coverage.sbatch
#
# Array id = data-row number in samples.ml.all40.tsv (header is line 1, so
# task 1 is the first sample). To submit specific indices:
#
#   /usr/bin/sbatch --array=1-4,6-9%4 scripts/hpc/backfill_coverage.sbatch
#
# Indices that do not yet have a coverage file:
#
#   awk -F'\t' 'NR>1 { print (NR-1) "\t" $1 }' workflows/config/samples.ml.all40.tsv |
#     while read -r i sample; do
#       [[ -s "$SCRATCH/arcsv2/ml/coverage_backfill/${sample}.coverage.tsv" ]] || printf '%s\n' "$i"
#     done

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

# samples.ml.all40.tsv has a header on line 1 and no comment line.
sample="$(awk -F'\t' -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n+1 { print $1; exit }' "${SAMPLES}")"
out="${OUTDIR}/${sample}.coverage.tsv"
if [[ -s "${out}" ]]; then
  echo "skip ${sample}; ${out} exists"
  exit 0
fi

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
