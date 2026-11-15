#!/usr/bin/env bash
set -euo pipefail
export PATH=/opt/conda/bin:$PATH
export PYTHONPATH=/work/python${PYTHONPATH:+:$PYTHONPATH}
mkdir -p /work/scratch/gap_analysis/loci
centers_true=(2883362 5962348 32311430 57219054 61995575)
centers_false=(28549495 61390333 64235350)
modes=(default precedence no_minblock no_absorb precedence_no_minblock_no_absorb)

for c in "${centers_true[@]}" "${centers_false[@]}"; do
  for mode in "${modes[@]}"; do
    out=/work/scratch/gap_analysis/loci/v2_${c}_${mode}.txt
    echo "=== V2 $c mode=$mode ==="
    python /work/scripts/diagnostics/trace_breakpoint_stages.py \
      --bam /data/work/HG01884/chr20.bam \
      --reference /data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
      --gaps /data/ref/gaps.bed \
      --chrom chr20 --center "$c" --pad 20000 --mode "$mode" \
      --out "$out" >/dev/null
    # summary line (empty match must not abort under pipefail)
    grep -nE "^CALL |^CALLS |final breakpoints|after merge|after absorb|after precedence|blocks n=|locus blocks" "$out" | head -40 || true
  done
done
echo V2_DONE
