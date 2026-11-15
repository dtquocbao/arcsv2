@echo off
cd /d C:\Users\dtquocbao\arcsv2-work
echo === containers ===
docker ps -a --filter name=arcsv2-ml0
echo === sample count ===
docker run --rm -v arcsv2-data:/data --entrypoint bash arcsv2:ml0 -lc "ls /data/samples | wc -l; du -sh /data/samples 2>/dev/null; ls /data/results_ml0_s3/*/chr*/arcsv_out.tab 2>/dev/null | wc -l"
echo === rest37 log tail ===
powershell -NoProfile -Command "Get-Content -Tail 30 bench\logs\ml0\s3_dl_rest37_host.log -ErrorAction SilentlyContinue"
echo === snakemake log tail ===
powershell -NoProfile -Command "Get-Content -Tail 15 bench\logs\ml0\s3_pilot_snakemake.log -ErrorAction SilentlyContinue"
