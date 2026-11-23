#!/bin/bash
#SBATCH --job-name=arcsv2-ml1
#SBATCH --account=132702112504
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Submit with /usr/bin/sbatch, not /sw/local/bin/sbatch.
# The 1.6.1 image has the old trainer, so PYTHONPATH points at this repo.
# Do not submit until the laptop report has been reviewed.
#
#   cd $SCRATCH/arcsv2/repo && /usr/bin/sbatch scripts/hpc/ml1_eval.sbatch

set -euo pipefail
REPO="${ARCSV2_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${REPO}/scripts/ml1_eval.py" ]] || REPO="${SCRATCH}/arcsv2/repo"
SIF="${ARCSV2_SIF:-${SCRATCH}/arcsv2/sif/arcsv2_1.6.1.sif}"
LABELS="${ARCSV2_LABELS:-${SCRATCH}/arcsv2/ml/dataset_v0.tsv.gz}"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="${SCRATCH}/arcsv2/ml/eval/ml1_${STAMP}"
mkdir -p "${OUT}"
export PYTHONPATH="${REPO}/python${PYTHONPATH:+:${PYTHONPATH}}"
/usr/bin/time -v singularity exec -B /scratch "${SIF}" \
  python "${REPO}/scripts/ml1_eval.py" \
  --labels "${LABELS}" \
  --out "${OUT}"
echo "report ${OUT}"
