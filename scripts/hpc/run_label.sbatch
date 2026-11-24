#!/bin/bash
#SBATCH --job-name=arcsv2-label
#SBATCH --account=132702112504
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Controller for workflows/label.smk. One core on the long partition so the
# process survives a disconnect. Snakemake submits the call/label jobs from
# this compute node. run_label.sh puts /usr/bin/sbatch on PATH (the
# /sw/local/bin/sbatch wrapper exits 0 when a job is rejected).
#
# Do not submit a batch until the pilot has been reviewed.
#
#   cd $SCRATCH/arcsv2/repo
#   /usr/bin/sbatch scripts/hpc/run_label.sbatch workflows/config/label.ml2.batch01.yaml
#
# Pass extra snakemake args after the config, or set them in ARCSV2_SNAKE_ARGS.

set -euo pipefail
REPO="${ARCSV2_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${REPO}/scripts/hpc/run_label.sh" ]] || REPO="${SCRATCH}/arcsv2/repo"
CONFIG="${1:-${ARCSV2_CONFIG:-}}"
if [[ -z "${CONFIG}" ]]; then
  echo "usage: sbatch scripts/hpc/run_label.sbatch <config.yaml>" >&2
  exit 2
fi
shift || true
cd "${REPO}"
bash scripts/hpc/run_label.sh "${CONFIG}" --rerun-incomplete "$@"
