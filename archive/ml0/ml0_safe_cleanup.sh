#!/bin/sh
set -eu
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
FREED=0
DELETED=""
for s in $(awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' /work/workflows/config/samples.ml.tsv); do
  ok=1
  for c in $CHROMS; do
    if [ ! -f "/data/work/${s}/${c}.bam" ] || [ ! -f "/data/work/${s}/${c}.bam.bai" ]; then
      ok=0
      break
    fi
  done
  if [ "$ok" = 1 ]; then
    cram="/data/samples/${s}/reads/${s}.cram"
    crai="/data/samples/${s}/reads/${s}.crai"
    if [ -f "$cram" ]; then
      sz=$(stat -c %s "$cram" 2>/dev/null || echo 0)
      crai_sz=0
      [ -f "$crai" ] && crai_sz=$(stat -c %s "$crai")
      rm -f "$cram" "$crai"
      total=$((sz + crai_sz))
      FREED=$((FREED + total))
      DELETED="${DELETED}${s}:$((total/1024/1024))MB "
    fi
  fi
done
echo "CRAM_DELETED: ${DELETED:-none}"
echo "FREED_BYTES: $FREED"
awk -v b="$FREED" 'BEGIN{printf "FREED_GB: %.2f\n", b/1024/1024/1024}'

TAB_FREED=0
TAB_LIST=""
for f in $(find /data/results_ml0_s3 -name 'arcsv_out.tab' -size -1k 2>/dev/null); do
  sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
  rm -f "$f"
  TAB_FREED=$((TAB_FREED + sz))
  TAB_LIST="${TAB_LIST}${f} "
done
echo "TABS_DELETED: ${TAB_LIST:-none}"
echo "TAB_FREED_BYTES: $TAB_FREED"

df -h /data
