#!/bin/bash
#SBATCH --job-name=arcsv2-ml2
#SBATCH --account=132702112504
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Sample holdout and chromosome holdout on dataset_v1. The frozen test split
# is dropped before folds. Four models: logistic and XGBoost, each with f1
# and f1+context. Logistic fills a missing coverage ratio with the training-fold
# median and adds a missing indicator. XGBoost keeps NaN.
#
# ML-0 chrom eval peaked near 2.5 GiB. This table is about 4x that, so the
# request is 32G. Medium allows one day.
#
# Submit with /usr/bin/sbatch, not /sw/local/bin/sbatch.
# Do not submit until dataset_v1 has been reviewed.
#
#   cd $SCRATCH/arcsv2/repo && /usr/bin/sbatch scripts/hpc/ml2_eval.sbatch

set -euo pipefail
REPO="${ARCSV2_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${REPO}/scripts/ml1_eval.py" ]] || REPO="${SCRATCH}/arcsv2/repo"
SIF="${ARCSV2_SIF:-${SCRATCH}/arcsv2/sif/arcsv2_1.6.1.sif}"
LABELS="${ARCSV2_LABELS:-${REPO}/export/ml2/dataset_v1.tsv.gz}"
export PYTHONPATH="${REPO}/python${PYTHONPATH:+:${PYTHONPATH}}"

for group in sample chrom; do
  OUT="${REPO}/export/ml2/eval/${group}"
  mkdir -p "${OUT}"
  echo "ml-eval outer-group=${group} labels=${LABELS} out=${OUT}"
  /usr/bin/time -v singularity exec -B /scratch "${SIF}" \
    env PYTHONPATH="${REPO}/python${PYTHONPATH:+:${PYTHONPATH}}" \
    python "${REPO}/scripts/ml1_eval.py" \
    --labels "${LABELS}" \
    --out "${OUT}" \
    --outer-group "${group}" \
    --exclude-split test \
    --include-rf none \
    --no-ablation
  echo "report ${OUT}"
done
