set -euo pipefail
export PYTHONPATH=
cd /work
echo "START_FILTER_BEDS $(date -u +%Y-%m-%dT%H:%M:%SZ)"
START=$(date +%s)
python scripts/fetch_filter_beds.py --dest /data/ref
END=$(date +%s)
echo "END_FILTER_BEDS $(date -u +%Y-%m-%dT%H:%M:%SZ) elapsed_s=$((END-START))"
ls -lh /data/ref/simple_repeats.bed /data/ref/segdups.bed /data/ref/blacklist.bed