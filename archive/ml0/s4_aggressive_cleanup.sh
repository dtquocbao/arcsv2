#!/usr/bin/env bash
set -euo pipefail
AVAIL_BEFORE=$(df -B1 /data | awk 'NR==2{print $4}')
FREED=0

add_freed_dir() {
  local d=$1
  [[ -d "$d" ]] || return 0
  local sz
  sz=$(du -sb "$d" | awk '{print $1}')
  rm -rf "$d"
  FREED=$((FREED + sz))
  echo "rm_dir $d bytes=$sz"
}

add_freed_file() {
  local f=$1
  [[ -f "$f" ]] || return 0
  local sz
  sz=$(stat -c%s "$f")
  rm -f "$f"
  FREED=$((FREED + sz))
  echo "rm_file $f bytes=$sz"
}

# 4a again: labeled work
for tsv in /data/labels_ml0_s3/*.labeled.tsv; do
  [[ -f "$tsv" ]] || continue
  s=$(basename "$tsv" .labeled.tsv)
  add_freed_dir "/data/work/$s"
done

count_bams() {
  local s=$1 wd=/data/work/$s n=0 c
  [[ -d "$wd" ]] || { echo 0; return; }
  for c in $(seq 1 22); do
    [[ -f "$wd/chr${c}.bam" ]] && n=$((n + 1))
  done
  [[ -f "$wd/chrX.bam" ]] && n=$((n + 1))
  echo "$n"
}

# 6: CRAMs for >=15 bams (unlabeled partial progress)
while IFS= read -r cram; do
  [[ -n "$cram" ]] || continue
  s=$(basename "$cram" .cram)
  nb=$(count_bams "$s")
  if [[ "$nb" -ge 15 ]]; then
    add_freed_file "/data/samples/${s}/reads/${s}.cram"
    add_freed_file "/data/samples/${s}/reads/${s}.cram.crai"
  fi
done < <(find /data/samples -name '*.cram' -type f 2>/dev/null)

AVAIL_AFTER=$(df -B1 /data | awk 'NR==2{print $4}')
echo "total_freed_bytes=$FREED"
echo "avail_before_bytes=$AVAIL_BEFORE"
echo "avail_after_bytes=$AVAIL_AFTER"
df -h /data
du -sh /data/work /data/samples 2>/dev/null
