#!/usr/bin/env bash
# Controller env for the Snakemake SLURM plugin. Install on $SCRATCH only.
# $HOME has a 10k-file quota; do not create this env there.
#
# Run inside a compute allocation (the solve is heavy for a login node):
#   srun -A 132702112504 -p short -c 2 --mem=8G -t 01:00:00 --pty bash
#   bash scripts/hpc/setup_controller.sh
#
# Pins (also docs/CLUSTER_TAMU.md):
#   micromamba 2.0.5
#   python 3.11          (snakemake-minimal requires >=3.11,<3.13)
#   snakemake-minimal 8.30.0 (bioconda pyhdfd78af_0)
#   snakemake-executor-plugin-slurm 2.8.0 (bioconda pyhdfd78af_0)
#   snakemake-executor-plugin-slurm-jobstep 0.6.1
set -euo pipefail

: "${SCRATCH:?SCRATCH is unset; run this on Grace}"

PREFIX="${SCRATCH}/arcsv2"
BIN="${PREFIX}/bin"
ENV_PREFIX="${PREFIX}/envs/smk"
MAMBA_VERSION="${MAMBA_VERSION:-2.0.5}"

export MAMBA_ROOT_PREFIX="${PREFIX}/mamba"
mkdir -p "${BIN}" "${MAMBA_ROOT_PREFIX}" "${PREFIX}/envs"

if command -v module >/dev/null 2>&1; then
  module load WebProxy
fi

if [[ ! -x "${BIN}/micromamba" ]]; then
  tmp="$(mktemp -d)"
  curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/${MAMBA_VERSION}" \
    | tar -xj -C "${tmp}" bin/micromamba
  install -m 0755 "${tmp}/bin/micromamba" "${BIN}/micromamba"
  rm -rf "${tmp}"
fi

export PATH="${BIN}:${PATH}"

if [[ -x "${ENV_PREFIX}/bin/snakemake" ]]; then
  echo "Controller env already present: ${ENV_PREFIX}"
  "${ENV_PREFIX}/bin/snakemake" --version
  exit 0
fi

micromamba create -y -p "${ENV_PREFIX}" \
  -c conda-forge -c bioconda \
  python=3.11 \
  "snakemake-minimal=8.30.0=pyhdfd78af_0" \
  "snakemake-executor-plugin-slurm=2.8.0=pyhdfd78af_0" \
  "snakemake-executor-plugin-slurm-jobstep=0.6.1"

echo "Controller env: ${ENV_PREFIX}"
"${ENV_PREFIX}/bin/snakemake" --version
