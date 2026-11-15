@echo off
REM Resume ML-0 S3 pilot (SIGTERM recovery). From arcsv2-work.
cd /d C:\Users\dtquocbao\arcsv2-work
echo SIGTERM_NOTE %date% %time%>> bench\logs\ml0\s3_pilot_snakemake.log
echo RESUME_HOST %date% %time%>> bench\logs\ml0\s3_pilot_snakemake.log
docker rm -f arcsv2-ml0-s3-pilot 2>nul
docker run --rm -e PYTHONPATH= -v arcsv2-data:/data -v C:\Users\dtquocbao\arcsv2-work:/work -w /work --name arcsv2-ml0-s3-pilot --entrypoint bash arcsv2:ml0 /work/scratch/ml0_s3_resume.sh >> bench\logs\ml0\s3_pilot_host.log 2>&1
echo s3 resume exit %ERRORLEVEL%
