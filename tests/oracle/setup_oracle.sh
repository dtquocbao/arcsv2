#!/usr/bin/env bash
# Install original ARC-SV (commit dd6f795) with the NumPy 2 fix np.Inf -> np.inf.
# Usage: ARCSV_ORIG=/path/to/arcsv ARCSV_VENV=/path/to/venv bash tests/oracle/setup_oracle.sh
set -euo pipefail

ORIG="${ARCSV_ORIG:-}"
if [[ -z "${ORIG}" ]]; then
  echo "Set ARCSV_ORIG to a checkout of SUwonglab/arcsv at dd6f795" >&2
  exit 1
fi
VENV="${ARCSV_VENV:-$(pwd)/.oracle-venv}"
python3 -m venv "${VENV}"
# shellcheck disable=SC1091
source "${VENV}/bin/activate"
python -m pip install --upgrade pip
python -m pip install "numpy<2" "scipy" "scikit-learn<1.6" "pysam" "pyinter"
# One-line NumPy fix required by the roadmap.
grep -RIl "np.Inf" "${ORIG}" | while read -r f; do
  sed -i.bak "s/np\\.Inf/np.inf/g" "$f"
done
python -m pip install --no-deps "${ORIG}"
echo "oracle ready in ${VENV}"
