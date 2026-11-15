set -euo pipefail
export PYTHONPATH=
cd /work
mkdir -p /work/bench/logs/ml0
echo "START_PILOT3_DOWNLOAD $(date -u +%Y-%m-%dT%H:%M:%SZ)"
START=$(date +%s)
python scripts/fetch_data.py --samples workflows/config/samples.ml.pilot3.tsv --dest /data --fetch-cram
END=$(date +%s)
ELAPSED=$((END-START))
echo "END_PILOT3_DOWNLOAD $(date -u +%Y-%m-%dT%H:%M:%SZ) elapsed_s=$ELAPSED"
# sizes
for s in HG01884 HG01150 HG00126; do
  echo "=== $s ==="
  du -sh /data/samples/$s 2>/dev/null || true
  ls -lh /data/samples/$s/reads/ 2>/dev/null || true
  ls -lh /data/haps/${s}.hap1.fa.gz /data/haps/${s}.hap2.fa.gz 2>/dev/null || true
done
# rough bandwidth from CRAM+hap bytes / elapsed
python - <<PY
import os, time
samples=["HG01884","HG01150","HG00126"]
total=0
for s in samples:
  base=f"/data/samples/{s}"
  for root,dirs,files in os.walk(base):
    for f in files:
      total += os.path.getsize(os.path.join(root,f))
elapsed=float(os.environ.get("ELAPSED","0") or 0)
print(f"pilot3_bytes={total} elapsed_s_env_note=see_shell")
PY
echo "ELAPSED_S=$ELAPSED"
if [ "$ELAPSED" -gt 0 ]; then
  BYTES=$(du -sb /data/samples/HG01884 /data/samples/HG01150 /data/samples/HG00126 2>/dev/null | awk '{s+=$1} END {print s}')
  python -c "b=int('$BYTES'); e=int('$ELAPSED'); print(f'pilot3_total_bytes={b} elapsed_s={e} MBps={b/e/1e6:.2f}')"
fi