# ARC-SV2 on VISION

Status: not supported for the caller or for labeling.

VISION is the wrong machine for this workflow:

- The container runtime is Pyxis/Enroot. Grace uses Singularity. The Snakemake profile does not submit Pyxis jobs.
- Compute nodes are GPU-only DGX H200s. ARC-SV2 does not use a GPU.
- The default is 14 CPUs per GPU.
- Accounting is in GPU-hours. Tier One is 970 GPU-h over 6 months.
- The maximum walltime is 48 h.
- CPU-only jobs are not documented.

Use VISION later for a GPU model (ML-2 and after). Do not put the chromosome caller or the haplotype labeler here.

## Ad-hoc Pyxis recipe (one sample, reference only)

Do not send this as the supported path. Enroot ignores `ENTRYPOINT`, so the command has to name `arcsv2`. `$HOME` is not mounted unless `--container-mount-home` is set. `$PWD` cannot be used in `#SBATCH` lines.

```bash
#SBATCH --time=12:00:00 --qos=standard --gres=gpu:1 --cpus-per-task=14
srun --container-image=ghcr.io#dtquocbao/arcsv2:1.6.1 \
     --container-mounts=/scratch/user/<vision-id>/arcsv2:/scratch/user/<vision-id>/arcsv2 \
     arcsv2 run --input ... --out ... --threads 14
```

## Draft note to VISION support

Do not send until the questions below are filled in.

1. Are CPU-only jobs (`--gres=gpu:0`) allowed?
2. What is the policy on controlled-access human genomic data (dbGaP / institutional DUA)?
3. What tier is our project, and what is its GPU-hour balance?
