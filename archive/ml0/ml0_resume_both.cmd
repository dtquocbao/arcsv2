@echo off
REM Resume rest37 prefetch + S3 pilot after host SIGKILL/SIGTERM (~02:45Z).
cd /d C:\Users\dtquocbao\arcsv2-work
mkdir bench\logs\ml0 2>nul

echo === status ===
docker ps -a --filter name=arcsv2-ml0
docker run --rm -v arcsv2-data:/data --entrypoint bash arcsv2:ml0 -lc "echo samples:; ls /data/samples | wc -l; echo call_tabs:; ls /data/results_ml0_s3/*/chr*/arcsv_out.tab 2>/dev/null | wc -l; du -sh /data/samples 2>/dev/null"

docker rm -f arcsv2-ml0-dl-rest37 arcsv2-ml0-s3-pilot 2>nul

echo === start rest37 prefetch ===
start "ml0-rest37" /B cmd /c "docker run --rm -e PYTHONPATH= -v arcsv2-data:/data -v C:\Users\dtquocbao\arcsv2-work:/work -w /work --name arcsv2-ml0-dl-rest37 --entrypoint bash arcsv2:ml0 /work/scratch/ml0_dl_rest37.sh >> bench\logs\ml0\s3_dl_rest37_host.log 2>&1"

timeout /t 5 /nobreak >nul

echo === start S3 resume ===
docker run --rm -e PYTHONPATH= -v arcsv2-data:/data -v C:\Users\dtquocbao\arcsv2-work:/work -w /work --name arcsv2-ml0-s3-pilot --entrypoint bash arcsv2:ml0 /work/scratch/ml0_s3_resume.sh >> bench\logs\ml0\s3_pilot_host.log 2>&1
echo s3 resume exit %ERRORLEVEL%
