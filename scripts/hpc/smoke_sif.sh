#!/usr/bin/env bash
# Smoke-test an already-pulled ARC-SV2 SIF (safe on login node; light).
set -euo pipefail
ARCSV2_SIF="${ARCSV2_SIF:-${SCRATCH:?SCRATCH unset}/arcsv2/sif/arcsv2_1.6.1.sif}"

if command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
elif command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
else
  echo "Neither singularity nor apptainer found." >&2
  exit 1
fi

if [[ ! -f "$ARCSV2_SIF" ]]; then
  echo "Missing SIF: $ARCSV2_SIF" >&2
  echo "Pull first with scripts/hpc/pull_sif.sh inside an srun/sbatch job." >&2
  exit 1
fi

"$RUNTIME" exec "$ARCSV2_SIF" arcsv2 --version
"$RUNTIME" exec "$ARCSV2_SIF" python -c "import arcsv2_native; print('native ok')"
echo "smoke OK ($RUNTIME) $ARCSV2_SIF"
