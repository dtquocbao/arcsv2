#!/usr/bin/env bash
set -euo pipefail
AVAIL_BEFORE=$(df -B1 /data | awk 'NR==2{print $4}')
echo "avail_before_bytes=$AVAIL_BEFORE"
FREED=0
for tsv in /data/labels_ml0_s3/*.labeled.tsv; do
  [[ -f "$tsv" ]] || continue
  s=$(basename "$tsv" .labeled.tsv)
  wd=/data/work/$s
  if [[ -d "$wd" ]]; then
    sz=$(du -sb "$wd" 2>/dev/null | awk '{print $1}')
    rm -rf "$wd"
    FREED=$((FREED + sz))
    echo "rm_labeled_work $s $sz"
  fi
done

maybe_delete_cram() {
  local s=$1
  local wd=/data/work/$s
  local n=0 c
  for c in $(seq 1 22); do
    [[ -f "$wd/chr${c}.bam" && -f "$wd/chr${c}.bam.bai" ]] && n=$((n + 1))
  done
  [[ -f "$wd/chrX.bam" && -f "$wd/chrX.bam.bai" ]] && n=$((n + 1))
  if [[ "$n" -eq 23 ]]; then
    for f in "/data/samples/${s}.cram" "/data/samples/${s}.cram.crai"; do
      if [[ -f "$f" ]]; then
        sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
        rm -f "$f"
        FREED=$((FREED + sz))
        echo "rm_cram_complete $f $sz"
      fi
    done
  fi
}

for tsv in /data/labels_ml0_s3/*.labeled.tsv; do
  [[ -f "$tsv" ]] || continue
  s=$(basename "$tsv" .labeled.tsv)
  maybe_delete_cram "$s"
done

TMP_CNT=0
while IFS= read -r -d '' f; do
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  rm -f "$f"
  FREED=$((FREED + sz))
  TMP_CNT=$((TMP_CNT + 1))
done < <(find /data/work -name '*.bam.tmp*' -print0 2>/dev/null)
echo "rm_bam_tmp count=$TMP_CNT"

ZCNT=0
while IFS= read -r -d '' f; do
  rm -f "$f"
  ZCNT=$((ZCNT + 1))
done < <(find /data/results_ml0_s3 -name arcsv_out.tab -empty -print0 2>/dev/null)
echo "rm_zero_arcsv count=$ZCNT"

AVAIL_AFTER=$(df -B1 /data | awk 'NR==2{print $4}')
echo "freed_bytes=$FREED"
echo "avail_after_bytes=$AVAIL_AFTER"
df -h /data
du -sh /data/work /data/samples 2>/dev/null
