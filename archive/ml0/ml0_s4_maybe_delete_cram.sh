#!/usr/bin/env bash
set -euo pipefail
cd /work
CHROMS=(chr{1..22} chrX)
for s in $(awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' workflows/config/samples.ml.tsv); do
  ok=1
  for c in "${CHROMS[@]}"; do
    if [[ ! -f "/data/work/${s}/${c}.bam" || ! -f "/data/work/${s}/${c}.bam.bai" ]]; then
      ok=0
      break
    fi
  done
  if (( ok )); then
    cram="/data/samples/${s}/reads/${s}.cram"
    if [[ -f "$cram" ]]; then
      echo "maybe_delete_cram ${s}"
      rm -f "$cram" "${cram}.crai" || true
    fi
  fi
done
find /data/work -type f \( -name '*.bam.tmp.*' -o -name '*.tmp' -o -name '*.partial' \) -print -delete 2>/dev/null || true
df -h /data
