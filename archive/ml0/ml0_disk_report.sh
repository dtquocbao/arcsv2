#!/bin/sh
set -eu
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
echo "=== COMPLETE_BAMS_STILL_CRAM ==="
for s in $(awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' /work/workflows/config/samples.ml.tsv); do
  nb=0
  for c in $CHROMS; do
    if [ -f "/data/work/${s}/${c}.bam" ] && [ -f "/data/work/${s}/${c}.bam.bai" ]; then
      nb=$((nb + 1))
    fi
  done
  cram="/data/samples/${s}/reads/${s}.cram"
  if [ "$nb" = 23 ] && [ -f "$cram" ]; then
    sz=$(du -sh "$cram" | cut -f1)
    echo "$s cram=$sz"
  fi
done
echo "=== PARTIAL_BAMS_WITH_CRAM ==="
for s in $(awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' /work/workflows/config/samples.ml.tsv); do
  nb=0
  for c in $CHROMS; do
    if [ -f "/data/work/${s}/${c}.bam" ] && [ -f "/data/work/${s}/${c}.bam.bai" ]; then
      nb=$((nb + 1))
    fi
  done
  cram="/data/samples/${s}/reads/${s}.cram"
  if [ "$nb" -gt 0 ] && [ "$nb" -lt 23 ] && [ -f "$cram" ]; then
    ws=$(du -sh "/data/work/$s" 2>/dev/null | cut -f1)
    cs=$(du -sh "$cram" 2>/dev/null | cut -f1)
    echo "$s bams=${nb}/23 work=$ws cram=$cs"
  fi
done
echo "=== CRAM inventory ==="
find /data/samples -name '*.cram' -exec du -sh {} \; 2>/dev/null | sort -hr
