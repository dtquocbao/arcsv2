#!/usr/bin/env bash
set -euo pipefail
OA=/data/envs/arcsv-orig/bin
export PATH="$OA:$PATH"
export PYTHONPATH=/data/orig-src${PYTHONPATH:+:$PYTHONPATH}
python -c "import arcsv; print('arcsv', arcsv.__file__)"
mkdir -p /work/scratch/gap_analysis/loci

# Original -v 2 on ±20kb around each center
centers=(2883362 5962348 32311430 57219054 61995575 28549495 61390333 64235350)
for c in "${centers[@]}"; do
  start=$((c - 20000))
  end=$((c + 20000))
  if [[ $start -lt 0 ]]; then start=0; fi
  out=/work/scratch/gap_analysis/loci/orig_${c}
  mkdir -p "$out"
  echo "=== ORIG $c chr20:${start}-${end} ==="
  "$OA/python" /data/orig-src/bin/arcsv call \
    -i /data/work/HG01884/chr20.bam \
    -r "chr20:${start}-${end}" \
    -R /data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -G /data/ref/gaps.bed \
    -o "$out" --overwrite -v 2 \
    >"$out/verbose.log" 2>&1 || echo "ORIG FAIL $c exit $?"
  # keep merge sections
  grep -nE "all_bp|just merged|merged and filtered|Beginning breakpoint|passing|softclip|cluster" "$out/verbose.log" | head -80 >"$out/merge_excerpt.txt" || true
  if [[ -f "$out/arcsv_out.tab" ]]; then
    wc -l "$out/arcsv_out.tab" | tee "$out/ncalls.txt"
  fi
done
echo ORIG_DONE
