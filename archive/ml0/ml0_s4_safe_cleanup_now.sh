#!/usr/bin/env bash
set -euo pipefail
cd /work
CHROMS=(chr{1..22} chrX)
FREE_BEFORE=$(df -BG /data | awk 'NR==2{gsub(/G/,"",$4); print $4}')
echo "=== CLEANUP START free_gb_before=${FREE_BEFORE} ==="

for f in /data/labels_ml0_s3/*.labeled.tsv; do
  [[ -f "$f" ]] || continue
  s=$(basename "$f" .labeled.tsv)
  echo "=== labeled sample=${s} ==="
  if [[ -d "/data/work/${s}" ]]; then
    du -sh "/data/work/${s}"
    rm -rf "/data/work/${s}"
    echo "deleted work/${s}"
  fi
  cram="/data/samples/${s}/reads/${s}.cram"
  if [[ -f "$cram" ]]; then
    du -sh "$cram" "${cram}.crai" 2>/dev/null || du -sh "$cram"
    rm -f "$cram" "${cram}.crai" || true
    echo "deleted cram for ${s}"
  fi
done

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
      du -sh "$cram" 2>/dev/null || true
      rm -f "$cram" "${cram}.crai" || true
    fi
  fi
done

while IFS= read -r -d '' f; do
  echo "delete 0-byte tab $f"
  rm -f "$f"
done < <(find /data/results_ml0_s3 -name 'arcsv_out.tab' -size 0 -print0 2>/dev/null)

find /data/work -type f \( -name '*.bam.tmp.*' -o -name '*.tmp' -o -name '*.partial' \) -print -delete 2>/dev/null || true

FREE_AFTER=$(df -BG /data | awk 'NR==2{gsub(/G/,"",$4); print $4}')
echo "=== CLEANUP END free_gb_after=${FREE_AFTER} freed_gb=$((FREE_AFTER - FREE_BEFORE)) ==="
df -h /data
n_lab=$(find /data/labels_ml0_s3 -maxdepth 1 -name '*.labeled.tsv' 2>/dev/null | wc -l)
n_call=$(find /data/results_ml0_s3 -name arcsv_out.tab 2>/dev/null | wc -l)
n_cram=$(find /data/samples -path '*/reads/*.cram' 2>/dev/null | wc -l)
echo "S4_PROGRESS labeled=${n_lab}/40 calls=${n_call} crams=${n_cram} free_gb=${FREE_AFTER}"
