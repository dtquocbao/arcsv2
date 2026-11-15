#!/usr/bin/env bash
# H1: pull frozen ARC-SV2 SIF on a Texas A&M Grace *compute* node.
# Do NOT run on a login node (singularity pull is heavy + needs WebProxy).
#
# Interactive (from Grace login):
#   srun --account=132702112504 --partition=short --nodes=1 \
#     --ntasks=1 --cpus-per-task=4 --mem=30G --time=01:00:00 \
#     --pty bash scripts/hpc/pull_sif.sh
#
# Or: sbatch scripts/hpc/pull_sif.sbatch
#
# Env overrides:
#   ARCSV2_GHCR   default ghcr.io/dtquocbao/arcsv2:1.6.1
#   ARCSV2_SIFDIR default $SCRATCH/arcsv2/sif
#   ARCSV2_SIF    default $ARCSV2_SIFDIR/arcsv2_1.6.1.sif
set -euo pipefail

ARCSV2_GHCR="${ARCSV2_GHCR:-ghcr.io/dtquocbao/arcsv2:1.6.1}"
ARCSV2_SIFDIR="${ARCSV2_SIFDIR:-${SCRATCH:?SCRATCH unset}/arcsv2/sif}"
ARCSV2_SIF="${ARCSV2_SIF:-${ARCSV2_SIFDIR}/arcsv2_1.6.1.sif}"

# Prefer singularity (Grace); fall back to apptainer only if needed.
if command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
elif command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
else
  echo "Neither singularity nor apptainer found on PATH." >&2
  exit 1
fi

echo "host=$(hostname) runtime=$RUNTIME"
echo "image=$ARCSV2_GHCR"
echo "out=$ARCSV2_SIF"

# Grace compute nodes: no outbound net until WebProxy is loaded.
if command -v module >/dev/null 2>&1; then
  module load WebProxy || true
fi

# Keep Singularity cache off $HOME (10 GB / 10k-file quota).
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-${TMPDIR:-/tmp}/.singularity}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SINGULARITY_CACHEDIR}"
mkdir -p "$SINGULARITY_CACHEDIR" "$ARCSV2_SIFDIR"
cd "$ARCSV2_SIFDIR"

if [[ -f "$ARCSV2_SIF" ]]; then
  echo "SIF already present: $ARCSV2_SIF"
  ls -lh "$ARCSV2_SIF"
else
  echo "Pulling (this can take a while)…"
  "$RUNTIME" pull "$ARCSV2_SIF" "docker://${ARCSV2_GHCR}"
fi

echo "Smoke:"
"$RUNTIME" exec "$ARCSV2_SIF" arcsv2 --version
"$RUNTIME" exec "$ARCSV2_SIF" python -c "import arcsv2_native; print('native ok')"
echo "H1 pull OK → $ARCSV2_SIF"
