#!/usr/bin/env bash
# Capped CRAM prefetch for ML-0 S4: keep about 2 unlabeled CRAMs ahead.
# Skips samples that already have a complete CRAM (coordinates with prior rest37).
# Haplotypes are fetched by ml0_s4_run.sh (DAG requirement); this script is CRAM-only.
set -euo pipefail
export PYTHONPATH=
cd /work
mkdir -p /work/bench/logs/ml0
LOG=/work/bench/logs/ml0/s4_prefetch.log
exec > >(tee -a "$LOG") 2>&1

TSV=workflows/config/samples.ml.tsv
LABEL_DIR=/data/labels_ml0_s3
MAX_AHEAD=2
MIN_FREE_GB=50

echo "START_S4_PREFETCH $(date -u +%Y-%m-%dT%H:%M:%SZ)"

sample_list() {
  awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' "$TSV"
}

CHROMS=(chr{1..22} chrX)

cram_ready() {
  local s="$1"
  [[ -f "/data/samples/${s}/reads/${s}.cram" ]]
}

# True when all 23 chrom BAM+BAI exist — CRAM is no longer needed for fetch_region.
bams_complete() {
  local s="$1" c
  for c in "${CHROMS[@]}"; do
    if [[ ! -f "/data/work/${s}/${c}.bam" || ! -f "/data/work/${s}/${c}.bam.bai" ]]; then
      return 1
    fi
  done
  return 0
}

labeled() {
  local s="$1"
  [[ -f "${LABEL_DIR}/${s}.labeled.tsv" ]]
}

bam_count() {
  local s="$1" c n=0
  for c in "${CHROMS[@]}"; do
    if [[ -f "/data/work/${s}/${c}.bam" && -f "/data/work/${s}/${c}.bam.bai" ]]; then
      n=$((n + 1))
    fi
  done
  echo "$n"
}

free_gb() {
  df -BG /data | awk 'NR==2{gsub(/G/,"",$4); print $4}'
}

ahead_count() {
  local n=0 s
  for s in $(sample_list); do
    # Count only in-flight CRAMs (still need chrom extraction), not post-BAM leftovers.
    if cram_ready "$s" && ! labeled "$s" && ! bams_complete "$s"; then
      n=$((n + 1))
    fi
  done
  echo "$n"
}

next_missing() {
  # Prefer samples with fewest BAMs so S4 gets CRAM for fresh samples, not only tail chroms.
  local s best="" best_n=999 n
  for s in $(sample_list); do
    if ! cram_ready "$s" && ! labeled "$s" && ! bams_complete "$s"; then
      n=$(bam_count "$s")
      if (( n < best_n )); then
        best_n=$n
        best=$s
      fi
    fi
  done
  if [[ -n "$best" ]]; then
    echo "$best"
    return 0
  fi
  return 1
}

# Wait until S4 hap phase has at least the pilot assemblies (DAG can build for some).
# Prefer waiting briefly so we do not race the hap fetch for bandwidth.
sleep 30

while true; do
  if [[ -f "${LABEL_DIR}/dataset.tsv" ]]; then
    echo "dataset present; prefetch exiting"
    break
  fi
  ahead=$(ahead_count)
  free=$(free_gb)
  echo "prefetch_status ahead=${ahead} free_gb=${free} $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if (( ahead >= MAX_AHEAD )); then
    sleep 120
    continue
  fi
  if (( free < MIN_FREE_GB )); then
    echo "disk low (free_gb=${free} < ${MIN_FREE_GB}); waiting"
    sleep 180
    continue
  fi
  if ! s=$(next_missing); then
    echo "no missing CRAMs; waiting for S4 to finish"
    sleep 300
    continue
  fi
  echo "prefetch_one_cram sample=${s}"
  tmp_tsv=$(mktemp /tmp/ml0_s4_one.XXXXXX.tsv)
  # Emit real header (first non-comment line) plus the matching sample row.
  awk -F'\t' -v s="$s" '/^#/ {next} $1=="sample" {print; next} $1==s{print}' "$TSV" > "$tmp_tsv"
  # --fetch-cram also refreshes haps if missing; safe and idempotent.
  python scripts/fetch_data.py --samples "$tmp_tsv" --dest /data --fetch-cram --yes
  rm -f "$tmp_tsv"
  echo "prefetch_done sample=${s}"
done

echo "END_S4_PREFETCH $(date -u +%Y-%m-%dT%H:%M:%SZ)"
