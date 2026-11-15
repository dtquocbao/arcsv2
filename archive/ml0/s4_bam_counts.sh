#!/usr/bin/env bash
for s in NA21102 HG03742 HG03784 NA20850 HG00738 HG01928 HG02602 HG03688; do
  n=0
  wd=/data/work/$s
  for c in $(seq 1 22); do
    [[ -f "$wd/chr${c}.bam" ]] && n=$((n + 1))
  done
  [[ -f "$wd/chrX.bam" ]] && n=$((n + 1))
  labeled=0
  [[ -f "/data/labels_ml0_s3/${s}.labeled.tsv" ]] && labeled=1
  echo "$s bams=$n labeled=$labeled"
done
