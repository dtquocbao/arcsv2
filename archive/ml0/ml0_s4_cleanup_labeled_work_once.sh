#!/usr/bin/env bash
# One-off: free disk by removing BAM work dirs (and stray CRAMs) for labeled samples only.
set -euo pipefail
LABEL_DIR=/data/labels_ml0_s3
FREE_BEFORE=$(df -BG /data | awk 'NR==2{gsub(/G/,"",$4); print $4}')
echo "free_gb_before=${FREE_BEFORE}"
for f in "${LABEL_DIR}"/*.labeled.tsv; do
  [[ -f "$f" ]] || continue
  s=$(basename "$f" .labeled.tsv)
  echo "=== sample=${s} ==="
  if [[ -d "/data/work/${s}" ]]; then
    du -sh "/data/work/${s}"
    rm -rf "/data/work/${s}"
    echo "deleted work/${s}"
  else
    echo "no work dir"
  fi
  cram="/data/samples/${s}/reads/${s}.cram"
  if [[ -f "$cram" ]]; then
    du -sh "$cram"
    rm -f "$cram" "${cram}.crai" || true
    echo "deleted cram"
  else
    echo "no cram"
  fi
done
df -h /data
FREE_AFTER=$(df -BG /data | awk 'NR==2{gsub(/G/,"",$4); print $4}')
echo "free_gb_after=${FREE_AFTER}"
echo "freed_approx_gb=$((FREE_AFTER - FREE_BEFORE))"
