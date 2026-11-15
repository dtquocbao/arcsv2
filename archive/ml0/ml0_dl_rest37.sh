set -euo pipefail
export PYTHONPATH=
cd /work
mkdir -p /work/bench/logs/ml0
LOG=/work/bench/logs/ml0/s3_dl_rest37.log
exec > >(tee -a "$LOG") 2>&1
echo "START_REST37_DOWNLOAD $(date -u +%Y-%m-%dT%H:%M:%SZ)"
START=$(date +%s)
python scripts/fetch_data.py --samples workflows/config/samples.ml.rest37.tsv --dest /data --fetch-cram
END=$(date +%s)
echo "END_REST37_DOWNLOAD $(date -u +%Y-%m-%dT%H:%M:%SZ) elapsed_s=$((END-START))"
du -sh /data/samples
ls /data/samples | wc -l