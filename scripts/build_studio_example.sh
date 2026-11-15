#!/usr/bin/env bash
# Build the bundled Studio example BAM + mini reference (≈100 kb, ~2–5 min call).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/python/arcsv2/dashboard/demo/example_bam"
mkdir -p "$OUT"
python "${ROOT}/scripts/simulate_haplotypes.py" \
  --out "$OUT" --length 100000 --events 8 --coverage 20 --seed 42
# Rename for Studio discovery
cp -f "$OUT/reads.bam" "$OUT/example.bam"
cp -f "$OUT/reads.bam.bai" "$OUT/example.bam.bai" 2>/dev/null || true
# empty gaps
: > "$OUT/gaps.bed"
echo "Wrote $OUT/example.bam and ref.fa"
