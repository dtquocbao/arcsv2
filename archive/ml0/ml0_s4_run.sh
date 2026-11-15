#!/usr/bin/env bash
# ML-0 S4: 40-sample whole-genome run (same settings as S3). Resumes pilot outputs.
# Phase 0 downloads all haplotype assemblies (needed for DAG build), then snakemake
# loops with --keep-going while capped CRAM prefetch runs in a sibling container.
set -euo pipefail
export PYTHONPATH=
cd /work
mkdir -p /work/bench/logs/ml0 /data/results_ml0_s3 /data/labels_ml0_s3 /work/scratch/ml0_s4

echo "START_S4 $(date -u +%Y-%m-%dT%H:%M:%SZ)"
START=$(date +%s)

SAMPLER_PID=""
DISK_CLEAN_PID=""
python - <<'PY' &
import os, time
from pathlib import Path
out = Path("/work/bench/logs/ml0/s4_rss_samples.tsv")
out.write_text("ts\tpid\tname\tvmrss_kb\n", encoding="utf-8")
start = time.time()
while True:
    rows = []
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\x00", b" ").decode("utf-8", "replace")[:200]
            if not any(k in cmd for k in ("snakemake", "arcsv2", "samtools", "minimap2", "ml-label", "ml-dataset", "fetch_data")):
                continue
            rss = 0
            for line in Path(f"/proc/{pid}/status").read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("VmRSS:"):
                    rss = int(line.split()[1])
                    break
            name = cmd.split()[0].rsplit("/", 1)[-1] if cmd.strip() else pid
            rows.append((int(time.time() - start), pid, name, rss))
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
  if [ -n "$DISK_CLEAN_PID" ]; then kill "$DISK_CLEAN_PID" 2>/dev/null || true; fi
}
trap cleanup EXIT

CHROMS=(chr{1..22} chrX)

free_gb() {
  df -BG /data | awk 'NR==2{gsub(/G/,"",$4); print $4}'
}

progress() {
  local n_lab n_call n_cram
  n_lab=$(find /data/labels_ml0_s3 -maxdepth 1 -name '*.labeled.tsv' 2>/dev/null | wc -l || true)
  n_call=$(find /data/results_ml0_s3 -name arcsv_out.tab 2>/dev/null | wc -l || true)
  n_cram=$(find /data/samples -path '*/reads/*.cram' 2>/dev/null | wc -l || true)
  echo "S4_PROGRESS labeled=${n_lab}/40 calls=${n_call} crams=${n_cram} free_gb=$(free_gb)"
}

bam_count() {
  local s="$1" c n=0
  for c in "${CHROMS[@]}"; do
    if [[ -f "/data/work/${s}/${c}.bam" && -f "/data/work/${s}/${c}.bam.bai" ]]; then
      n=$((n + 1))
    fi
  done
  echo "$n"
}

maybe_delete_cram() {
  # Drop CRAM when all 23 BAMs exist, or earlier under low disk (>=18 BAMs)
  # so BAM+CRAM+sort-temps cannot fill the ~1 TB volume again.
  local s cram n fg thresh=23
  fg=$(free_gb)
  if (( fg < 60 )); then
    thresh=18
  fi
  for s in $(awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' workflows/config/samples.ml.tsv); do
    cram="/data/samples/${s}/reads/${s}.cram"
    [[ -f "$cram" ]] || continue
    n=$(bam_count "$s")
    if (( n >= thresh )); then
      echo "maybe_delete_cram ${s} bams=${n}/23 free_gb=${fg} thresh=${thresh}"
      rm -f "$cram" "${cram}.crai" || true
    fi
  done
}

cleanup_sort_temps() {
  # Only remove STALE sort leftovers (mtime > 3h). Never delete live
  # samtools *.bam.tmp.* mid-merge — that caused "tmp.0000.bam: No such file".
  local n
  n=$(find /data/work -name '*.bam.tmp*' -mmin +180 2>/dev/null | wc -l || true)
  if (( n > 0 )); then
    echo "cleanup_sort_temps stale_count=${n}"
    find /data/work -name '*.bam.tmp*' -mmin +180 -delete 2>/dev/null || true
  fi
}

if [[ -f /data/manifest.lock.json ]]; then
  echo "S4_HAPS_SKIP manifest present $(date -u +%Y-%m-%dT%H:%M:%SZ)"
else
  echo "S4_HAPS_FETCH $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  # Assemblies only (no CRAM) so label.smk can build the full 40-sample DAG.
  python -u scripts/fetch_data.py --samples workflows/config/samples.ml.tsv --dest /data \
    2>&1 | tee -a /work/bench/logs/ml0/s4_haps_fetch.log
  echo "S4_HAPS_DONE $(date -u +%Y-%m-%dT%H:%M:%SZ)"
fi

maybe_delete_cram || true
progress

cleanup_labeled_crams() {
  local s cram
  for s in $(awk -F'\t' '/^#/ {next} $1=="sample" {next} NF{print $1}' workflows/config/samples.ml.tsv); do
    if [[ -f "/data/labels_ml0_s3/${s}.labeled.tsv" ]]; then
      if [[ -d "/data/work/${s}" ]]; then
        echo "cleanup_work ${s}"
        rm -rf "/data/work/${s}" || true
      fi
      cram="/data/samples/${s}/reads/${s}.cram"
      if [[ -f "$cram" ]]; then
        echo "cleanup_cram ${s}"
        rm -f "$cram" "${cram}.crai" || true
      fi
    fi
  done
}

# Reclaim labeled BAMs / finished CRAMs every 10 min during long snakemake attempts
# (cleanup only at attempt end was too late and caused ENOSPC).
(
  while true; do
    sleep 600
    echo "S4_PERIODIC_CLEAN $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    cleanup_labeled_crams || true
    maybe_delete_cram || true
    cleanup_sort_temps || true
    progress || true
  done
) &
DISK_CLEAN_PID=$!

attempt=0
while true; do
  attempt=$((attempt + 1))
  echo "S4_SNAKEMAKE_ATTEMPT ${attempt} $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  set +e
  snakemake -s workflows/label.smk \
    --configfile workflows/config/label.ml0_s4.yaml \
    --cores 22 --resources mem_mb=100000 \
    --rerun-incomplete --rerun-triggers mtime \
    --keep-going \
    --report /work/bench/logs/ml0/s4_snakemake_report.html --report-after-run \
    -p 2>&1 | tee -a /work/bench/logs/ml0/s4_snakemake.log
  rc=${PIPESTATUS[0]}
  set -e
  echo "S4_SNAKEMAKE_RC rc=${rc} attempt=${attempt}"
  cleanup_labeled_crams || true
  maybe_delete_cram || true

  if [[ -f /data/labels_ml0_s3/dataset.tsv ]]; then
    echo "S4_DATASET_READY"
    break
  fi
  if [[ "$rc" -eq 0 ]]; then
    break
  fi
  progress
  fg=$(free_gb)
  if (( fg < 5 )); then
    echo "waiting for disk/prefetch (free_gb=${fg}) retry in 90s"
    sleep 90
  else
    echo "waiting for CRAM prefetch / retry in 90s"
    sleep 90
  fi
done

END=$(date +%s)
echo "END_S4 $(date -u +%Y-%m-%dT%H:%M:%SZ) elapsed_s=$((END-START))" | tee -a /work/bench/logs/ml0/s4_snakemake.log
ls -lh /data/labels_ml0_s3/ || true
du -sh /data/results_ml0_s3 /data/labels_ml0_s3 /data/samples /data/work 2>/dev/null || true
if [[ -f /data/labels_ml0_s3/dataset.summary.json ]]; then
  python -c "import json; print(json.dumps(json.load(open('/data/labels_ml0_s3/dataset.summary.json')), indent=2)[:4000])"
fi
