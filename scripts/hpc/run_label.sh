#!/usr/bin/env bash
# Launch the label workflow on Grace with the SLURM executor plugin.
# Usage: bash scripts/hpc/run_label.sh <configfile> [extra snakemake args]
#
# The controller is one process. Run it in tmux on a login node, or as a
# 1-core job on the long partition so it stays up for the whole DAG.
# TODO(Bao): confirm HPRC policy on long-running login processes.
#
#   tmux new -s arcsv2
#   bash scripts/hpc/run_label.sh workflows/config/label.yaml
#
# Compute nodes have no internet. http_proxy / https_proxy are exported here
# (SLURM's default --export=ALL forwards them) and the profile's apptainer-args
# pass them into the container.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <configfile> [extra snakemake args]" >&2
  exit 2
fi
CONFIGFILE="$1"
shift

: "${SCRATCH:?SCRATCH is unset; run this on Grace}"

PREFIX="${SCRATCH}/arcsv2"
ENV_PREFIX="${PREFIX}/envs/smk"
if [[ ! -x "${ENV_PREFIX}/bin/snakemake" ]]; then
  echo "Controller env missing (${ENV_PREFIX}). Run scripts/hpc/setup_controller.sh first." >&2
  exit 1
fi

export MAMBA_ROOT_PREFIX="${PREFIX}/mamba"
if [[ -x "${PREFIX}/bin/micromamba" ]]; then
  # shellcheck disable=SC1090
  eval "$("${PREFIX}/bin/micromamba" shell hook -s bash)"
  micromamba activate "${ENV_PREFIX}"
else
  export PATH="${ENV_PREFIX}/bin:${PATH}"
fi

export ARCSV2_SIF="${ARCSV2_SIF:-${PREFIX}/sif/arcsv2_1.6.1.sif}"
# Grace compute nodes have no direct internet. Set these before snakemake
# starts: the profile expands $http_proxy into apptainer-args, and SLURM's
# default --export=ALL forwards the variables into each job.
export http_proxy=http://10.73.132.63:8080
export https_proxy=http://10.73.132.63:8080
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-${PREFIX}/singularity-cache}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SINGULARITY_CACHEDIR}}"
mkdir -p "${SINGULARITY_CACHEDIR}"

# Snakemake 8's method name is apptainer. Grace only ships singularity.
export APPTAINER_CMD="${APPTAINER_CMD:-singularity}"
export SINGULARITY_CMD="${SINGULARITY_CMD:-singularity}"

# Grace puts a wrapper at /sw/local/bin/sbatch ahead of the real binary. It
# exits 0 when a submission is rejected and mangles the long quoted --wrap
# command the SLURM plugin sends, so Snakemake waits forever on a job that was
# never queued. Call the real sbatch. Site job_submit checks still run on the
# SLURM server.
if [[ -x /usr/bin/sbatch ]]; then
  SBATCH_DIR="${PREFIX}/bin/slurm-direct"
  mkdir -p "${SBATCH_DIR}"
  ln -sfn /usr/bin/sbatch "${SBATCH_DIR}/sbatch"
  export PATH="${SBATCH_DIR}:${PATH}"
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

LOGDIR="${PREFIX}/logs"
mkdir -p "${LOGDIR}"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
MANIFEST="${LOGDIR}/run_manifest.${STAMP}.txt"
{
  printf 'time_utc=%s\n' "${STAMP}"
  printf 'config=%s\n' "${CONFIGFILE}"
  printf 'git_describe=%s\n' "$(git describe --always --dirty 2>/dev/null || echo unknown)"
  printf 'commit=%s\n' "$(git rev-parse HEAD 2>/dev/null || echo unknown)"
  printf 'sif=%s\n' "${ARCSV2_SIF}"
  if [[ -f "${ARCSV2_SIF}" ]]; then
    printf 'sif_sha256=%s\n' "$(sha256sum "${ARCSV2_SIF}" | awk '{print $1}')"
  else
    printf 'sif_sha256=missing\n'
  fi
  printf 'snakemake=%s\n' "$(snakemake --version)"
  printf 'sbatch=%s\n' "$(command -v sbatch)"
} > "${MANIFEST}"
echo "run manifest: ${MANIFEST}"

snakemake -s workflows/label.smk \
  --profile workflows/profiles/slurm-grace \
  --configfile "${CONFIGFILE}" \
  "$@"
