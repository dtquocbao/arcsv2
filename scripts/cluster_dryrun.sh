#!/usr/bin/env bash
# U5: dry-run label workflow with SLURM profile (no jobs submitted).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"
mkdir -p tests/fixtures/label-smoke haps logs/slurm
printf '>chrSim\nACGTACGT\n' > tests/fixtures/label-smoke/ref.fa
cp tests/fixtures/label-smoke/ref.fa haps/sim.hap1.fa
cp tests/fixtures/label-smoke/ref.fa haps/sim.hap2.fa
: > tests/fixtures/label-smoke/gaps.bed
: > tests/fixtures/label-smoke/unreliable.bed
: > tests/fixtures/label-smoke/sim.cram

export ARCSV2_SIF="${ARCSV2_SIF:-}"
snakemake -s workflows/label.smk -n \
  --profile workflows/profiles/slurm \
  --configfile workflows/config/label.smoke.yaml

echo "OK: SLURM profile dry-run"
