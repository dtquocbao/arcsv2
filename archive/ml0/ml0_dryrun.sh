export PYTHONPATH=
snakemake -s workflows/label.smk --configfile workflows/config/label.ml0_s3.yaml -n 2>&1 | tail -40