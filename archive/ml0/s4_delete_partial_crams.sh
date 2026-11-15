#!/usr/bin/env bash
set -euo pipefail
AVAIL_BEFORE=$(df -B1 /data | awk 'NR==2{print $4}')
FREED=0

count_bams() {
  local s=$1 wd=/data/work/$s n=0 c
  [[ -d "$wd" ]] || { echo 0; return; }
  for c in $(seq 1 22); do
    [[ -f "$wd/chr${c}.bam" ]] && n=$((n + 1))
  done
  [[ -f "$wd/chrX.bam" ]] && n=$((n + 1))
  echo "$n"
}

for s in NA21102 HG03742 HG03784 NA20850; do
  nb=$(count_bams "$s")
  if [[ "$nb" -ge 15 ]]; then
    for f in "/data/samples/${s}/reads/${s}.cram" "/data/samples/${s}/reads/${s}.cram.crai"; do
      if [[ -f "$f" ]]; then
        sz=$(stat -c%s "$f")
        rm -f "$f"
        FREED=$((FREED + sz))
        echo "rm_cram sample=$s bams=$nb file=$f bytes=$sz"
      fi
    done
  else
    echo "keep_cram sample=$s bams=$nb"
  fi
done

AVAIL_AFTER=$(df -B1 /data | awk 'NR==2{print $4}')
echo "step6_freed_bytes=$FREED"
echo "avail_after_bytes=$AVAIL_AFTER"
df -h /data
