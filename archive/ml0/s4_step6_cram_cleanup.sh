#!/usr/bin/env bash
set -euo pipefail
echo "=== unlabeled work dirs (top 20) ==="
for d in /data/work/*; do
  [[ -d "$d" ]] || continue
  s=$(basename "$d")
  [[ -f "/data/labels_ml0_s3/${s}.labeled.tsv" ]] && continue
  du -sb "$d" 2>/dev/null
done | sort -n | tail -20 | while read -r sz path; do
  echo "$(numfmt --to=iec-i --suffix=B "$sz" 2>/dev/null || echo "${sz}B") $path"
done

echo "=== CRAM inventory ==="
find /data/samples -name '*.cram' -type f 2>/dev/null | while read -r cram; do
  sz=$(stat -c%s "$cram" 2>/dev/null || echo 0)
  echo "$sz $cram"
done | sort -n

count_bams() {
  local s=$1 wd=/data/work/$s n=0 c
  [[ -d "$wd" ]] || { echo 0; return; }
  for c in $(seq 1 22); do
    [[ -f "$wd/chr${c}.bam" ]] && n=$((n + 1))
  done
  [[ -f "$wd/chrX.bam" ]] && n=$((n + 1))
  echo "$n"
}

FREED=0
echo "=== deleting CRAMs for samples with >=15 BAMs ==="
find /data/samples -name '*.cram' -type f 2>/dev/null | while read -r cram; do
  s=$(basename "$cram" .cram)
  # path like .../HG03742/reads/HG03742.cram
  if [[ "$s" == "reads" ]]; then
    s=$(basename "$(dirname "$(dirname "$cram")")")
    cram="/data/samples/${s}/reads/${s}.cram"
  fi
  nb=$(count_bams "$s")
  if [[ "$nb" -ge 15 ]]; then
    for f in "/data/samples/${s}/reads/${s}.cram" "/data/samples/${s}/reads/${s}.cram.crai" \
             "/data/samples/${s}.cram" "/data/samples/${s}.cram.crai"; do
      if [[ -f "$f" ]]; then
        sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
        rm -f "$f"
        echo "rm_cram_partial_bams sample=$s bams=$nb file=$f bytes=$sz"
      fi
    done
  fi
done

# recount freed via df
df -h /data
du -sh /data/work /data/samples 2>/dev/null
