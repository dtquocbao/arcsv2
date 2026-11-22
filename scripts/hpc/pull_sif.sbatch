#!/bin/bash
#SBATCH --job-name=arcsv2-pull
#SBATCH --account=132702112504
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# From Grace login (repo already on $SCRATCH or $HOME):
#   mkdir -p $SCRATCH/arcsv2/logs
#   cd $SCRATCH/arcsv2/repo && sbatch scripts/hpc/pull_sif.sbatch
#
# Or copy-paste the body into: srun ... --pty bash

set -euo pipefail
# Resolve script dir whether submitted from repo or via absolute path.
# SLURM runs a spooled copy of this file, so BASH_SOURCE is not the repo path.
REPO="${ARCSV2_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "${REPO}/scripts/hpc/pull_sif.sh" ]] || REPO="${SCRATCH}/arcsv2/repo"
bash "${REPO}/scripts/hpc/pull_sif.sh"
