# ARC-SV2 on Texas A&M clusters

The copy-paste runbook (pull, golden chr20, one genome, storage, sync) is
added with the rest of the cluster docs. This file records the controller
pins that `scripts/hpc/setup_controller.sh` installs.

## Controller environment

Create it on `$SCRATCH`, never in `$HOME` (10k-file quota). Run the installer
from a short compute job:

```bash
srun -A 132702112504 -p short -c 2 --mem=8G -t 01:00:00 --pty bash
bash /path/to/arcsv2/scripts/hpc/setup_controller.sh
```

The Snakemake controller stays on the host. `sbatch` is not in the SIF, so
the controller cannot run inside the container. Rules still run in the image
(`container:` plus `--software-deployment-method apptainer`).

| Piece | Pin |
|:--|:--|
| micromamba binary | 2.0.5 in `$SCRATCH/arcsv2/bin` (image base tag; package cache under `$SCRATCH/arcsv2/mamba`) |
| python | 3.11 (snakemake-minimal 8.30.0 requires `>=3.11,<3.13`; the patch is whatever conda-forge solves) |
| snakemake-minimal | 8.30.0 bioconda build `pyhdfd78af_0` |
| snakemake-executor-plugin-slurm | 2.8.0 bioconda build `pyhdfd78af_0` |
| snakemake-executor-plugin-slurm-jobstep | 0.6.1 (required by plugin 2.8.0: `>=0.6.0,<1`) |

Checked against Snakemake 8.30.0 and plugin 2.8.0 source:

- Profile keys `executor`, `jobs`, `restart-times`, `latency-wait`, `printshellcmds`, `software-deployment-method`, `apptainer-args`, `default-resources`, and `set-resources` match the 8.30 CLI. A mapping under `default-resources` becomes `name=value`. A two-level `set-resources` mapping becomes `rule:resource=value`.
- `runtime` is minutes. Snakemake stores it as an integer number of minutes, and plugin 2.8.0 passes that integer to `sbatch -t`. `mem_mb` is passed to `sbatch --mem` (megabytes). Resource names are `slurm_account`, `slurm_partition`, `runtime`, and `mem_mb`.
- `--software-deployment-method` accepts `conda`, `apptainer`, and `env-modules`. `--singularity-args` is only an alias of `--apptainer-args`.

`--cleanenv` stays in the Grace and FASTER profiles. The image sets `ENV PATH="/opt/conda/bin:$PATH"` in the Dockerfile, and Singularity applies image environment variables after `--cleanenv`, so `arcsv2` remains on `PATH`. The proxy is passed back with `--env` because `--cleanenv` would drop `http_proxy`. This was not executed here (no `singularity` binary on the laptop). If a Grace job reports `arcsv2: command not found`, remove `--cleanenv` from `apptainer-args`.

Launch with `scripts/hpc/run_label.sh`. The controller is one process: tmux on a login node, or a 1-core `long` job. TODO(Bao): confirm HPRC policy on long-running login processes.
