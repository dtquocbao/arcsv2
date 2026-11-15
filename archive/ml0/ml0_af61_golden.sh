#!/usr/bin/env bash
# Fresh full-chr20 re-call with host-mounted AF.6f caller, then check_golden --fast.
set -euo pipefail
LOG=/work/bench/logs/ml0/af61_chr20.log
OUT=/data/results_ml0_af61
REF=/data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
GAPS=/data/ref/gaps.bed
mkdir -p /work/bench/logs/ml0 "$OUT/HG01884" "$OUT/NA19700"
exec > >(tee -a "$LOG") 2>&1

echo "=== AF.6f golden re-call $(date -Is) ==="
python - <<'PY'
from pathlib import Path
p = Path("/work/python/arcsv2/engine/caller.py")
text = p.read_text()
needle = 'f"{af:.6f}"'
assert needle in text, "mounted caller missing AF.6f"
print("mounted_af_fmt OK", needle)
import arcsv2.engine.caller as c
print("loaded", c.__file__)
assert needle in Path(c.__file__).read_text()
import arcsv2_native
print("native", arcsv2_native.__file__)
PY

for sample in HG01884 NA19700; do
  echo "--- call $sample chr20 ---"
  t0=$(date +%s)
  arcsv2 call \
    -i "/data/work/${sample}/chr20.bam" \
    -R "$REF" \
    -G "$GAPS" \
    --region chr20 \
    -o "$OUT/${sample}" \
    --threads 8 \
    --engine rust
  t1=$(date +%s)
  # --fast reads {sample}/merged.tab; raw chr20 table is what prior S1 used
  cp "$OUT/${sample}/arcsv_out.tab" "$OUT/${sample}/merged.tab"
  rows=$(tail -n +2 "$OUT/${sample}/merged.tab" | wc -l)
  echo "TIME wall_s=$((t1-t0)) sample=$sample rows=$rows"
  python - <<PY
from pathlib import Path
lines = Path("$OUT/${sample}/merged.tab").read_text().splitlines()
hdr = lines[0].split("\t")
ai = hdr.index("af")
afs = sorted({r.split("\t")[ai] for r in lines[1:]})[:8]
print("sample=$sample af_samples", afs)
ok = any(("." in a) and (len(a.split(".", 1)[1]) == 6) for a in afs)
assert ok, afs
PY
done

echo "=== check_golden --fast ==="
python /work/scripts/check_golden.py --fast --results "$OUT" | tee /work/bench/logs/ml0/golden_fast_af61.txt
echo "=== done $(date -Is) ==="
