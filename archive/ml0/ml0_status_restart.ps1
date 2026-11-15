$ErrorActionPreference = "Continue"
Write-Host "=== docker ps ==="
docker ps -a --filter name=arcsv2-ml0
Write-Host "=== progress counts ==="
docker run --rm -v arcsv2-data:/data --entrypoint bash arcsv2:ml0 -lc "echo calls:; ls /data/results_ml0_s3/*/chr*/arcsv_out.tab 2>/dev/null | wc -l; echo bams:; ls /data/work/*/chr*.bam 2>/dev/null | wc -l; echo samples:; ls /data/samples | wc -l; echo incomplete:; find /data/results_ml0_s3 -name '*.tab.incomplete' 2>/dev/null | wc -l; find /data/work -name '*.bam.incomplete' 2>/dev/null | wc -l"
