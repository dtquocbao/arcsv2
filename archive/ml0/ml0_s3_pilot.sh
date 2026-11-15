set -euo pipefail
export PYTHONPATH=
cd /work
mkdir -p /work/bench/logs/ml0 /data/results_ml0_s3 /data/labels_ml0_s3 /work/scratch/ml0_s3

echo "START_S3_PILOT $(date -u +%Y-%m-%dT%H:%M:%SZ)"
START=$(date +%s)

# Background 2 Hz RSS sampler for the docker PID tree (self)
SAMPLER_PID=""
python - <<'PY' &
import os, time, json
from pathlib import Path
out = Path("/work/bench/logs/ml0/s3_rss_samples.tsv")
out.write_text("ts\tpid\tname\tvmrss_kb\n", encoding="utf-8")
start = time.time()
self = os.getpid()
# sample all processes; filter arcsv2/snakemake/samtools/minimap2
while True:
    rows = []
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\x00", b" ").decode("utf-8", "replace")[:200]
            if not any(k in cmd for k in ("snakemake", "arcsv2", "samtools", "minimap2", "ml-label", "ml-dataset")):
                continue
            rss = 0
            for line in Path(f"/proc/{pid}/status").read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("VmRSS:"):
                    rss = int(line.split()[1])
                    break
            name = cmd.split()[0].rsplit("/", 1)[-1] if cmd.strip() else pid
            rows.append((int(time.time()-start), pid, name, rss))
        except Exception:
            pass
    if rows:
        with out.open("a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
    time.sleep(2.0)
PY
SAMPLER_PID=$!

cleanup() {
  if [ -n "$SAMPLER_PID" ]; then kill "$SAMPLER_PID" 2>/dev/null || true; fi
}
trap cleanup EXIT

snakemake -s workflows/label.smk \
  --configfile workflows/config/label.ml0_s3.yaml \
  --cores 22 --resources mem_mb=100000 \
  --rerun-incomplete --rerun-triggers mtime \
  --report /work/bench/logs/ml0/s3_snakemake_report.html --report-after-run \
  -p 2>&1 | tee /work/bench/logs/ml0/s3_pilot_snakemake.log

END=$(date +%s)
echo "END_S3_PILOT $(date -u +%Y-%m-%dT%H:%M:%SZ) elapsed_s=$((END-START))" | tee -a /work/bench/logs/ml0/s3_pilot_snakemake.log

echo "=== outputs ==="
ls -lh /data/labels_ml0_s3/ || true
du -sh /data/results_ml0_s3 /data/labels_ml0_s3 2>/dev/null || true
for s in HG01884 HG01150 HG00126; do
  echo "--- $s disk ---"
  du -sh /data/work/$s /data/samples/$s /data/results_ml0_s3/$s 2>/dev/null || true
done
if [ -f /data/labels_ml0_s3/dataset.summary.json ]; then
  python -c "import json; print(json.dumps(json.load(open('/data/labels_ml0_s3/dataset.summary.json')), indent=2)[:4000])"
fi
# peak RSS summary
python - <<'PY'
from pathlib import Path
from collections import defaultdict
p = Path("/work/bench/logs/ml0/s3_rss_samples.tsv")
if not p.is_file():
    print("no rss samples")
    raise SystemExit
peak = defaultdict(int)
for i, line in enumerate(p.read_text(encoding="utf-8").splitlines()):
    if i == 0 or not line.strip():
        continue
    ts, pid, name, rss = line.split("\t")
    key = name
    peak[key] = max(peak[key], int(rss))
print("peak_rss_kb_by_cmd:")
for k,v in sorted(peak.items(), key=lambda kv: -kv[1])[:30]:
    print(f"  {k}: {v} KB ({v/1024/1024:.2f} GB)")
print(f"global_peak_kb={max(peak.values()) if peak else 0}")
PY