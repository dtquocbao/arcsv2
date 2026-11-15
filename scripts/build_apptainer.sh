#!/usr/bin/env bash
# Build arcsv2.sif for SLURM. Prefer GHCR pull on clusters that cannot build
# (e.g. Texas A&M Grace — use scripts/hpc/pull_sif.sh instead).
#
# Runtime: singularity if present, else apptainer.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-${ROOT}/arcsv2.sif}"

if command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
elif command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
else
  echo "Neither singularity nor apptainer found." >&2
  exit 1
fi

if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
  echo "Building Docker image arcsv2:dev …"
  docker build -t arcsv2:dev "${ROOT}"
  echo "Converting to SIF via ${RUNTIME}: ${OUT}"
  "$RUNTIME" build "${OUT}" "docker-daemon://arcsv2:dev"
  exit 0
fi

echo "Docker not available; building from apptainer.def with repo bind …"
"$RUNTIME" build --bind "${ROOT}:/opt/arcsv2-src" "${OUT}" "${ROOT}/apptainer.def"
