# Running ARC-SV2 on a SLURM cluster (U5 / HPC0)

Use this when Docker is not available on compute nodes. Prefer a **pre-built
image from GHCR** (required on Texas A&M Grace: no image builds there).

**Runtime name:** use **`singularity`** (SingularityCE). Keep `apptainer` only
as a fallback alias on sites that provide it. On Grace the binary is
`/usr/bin/singularity` (no `apptainer`).

Frozen caller image: **`arcsv2:ml0`** at git tag **`caller-v1.6.1`**.

| Tag on GHCR | Meaning |
|:--|:--|
| `ghcr.io/dtquocbao/arcsv2:1.6.1` | Frozen caller-v1.6.1 (use this) |
| `ghcr.io/dtquocbao/arcsv2:ml0` | Same digest, convenience tag |
| `ghcr.io/dtquocbao/arcsv2:latest-dev` | Same freeze for now; may move later |

Publish from a laptop (once): `bash scripts/publish_ghcr.sh` or
`.\scripts\publish_ghcr.ps1` (needs Docker + GHCR token). Make the package
**Public** so cluster pulls need no GitHub token.

## Texas A&M Grace (verified 2026-10-16)

- Account: `132702112504`. Login nodes: **≤8 cores** total — never pull or run
  heavy work on login nodes.
- **Pull only inside a compute job.** Compute nodes need internet via
  `module load WebProxy` (proxy `10.73.132.63:8080`).
- Cache: `export SINGULARITY_CACHEDIR=$TMPDIR/.singularity` (keep `$HOME` clean).
- Singularity binds `/scratch`, `$HOME`, `/tmp`, `$PWD` automatically.

### H1 — pull + smoke (copy-paste)

On a Grace **login** node, create dirs then submit a short job:

```bash
mkdir -p $SCRATCH/arcsv2/{sif,ref,data,haps,results,labels,logs}
```

**Interactive `srun` (recommended for first pull):**

```bash
srun --account=132702112504 --partition=short --nodes=1 \
  --ntasks=1 --cpus-per-task=4 --mem=30G --time=01:00:00 --pty bash -lc '
set -euo pipefail
module load WebProxy
export SINGULARITY_CACHEDIR=$TMPDIR/.singularity
mkdir -p "$SINGULARITY_CACHEDIR" $SCRATCH/arcsv2/sif
cd $SCRATCH/arcsv2/sif
singularity pull arcsv2_1.6.1.sif docker://ghcr.io/dtquocbao/arcsv2:1.6.1
singularity exec arcsv2_1.6.1.sif arcsv2 --version
singularity exec arcsv2_1.6.1.sif python -c "import arcsv2_native; print(\"native ok\")"
'
```

**Or batch:**

```bash
# from a checkout that has scripts/hpc/
sbatch -D $SCRATCH/arcsv2/logs /path/to/arcsv2/scripts/hpc/pull_sif.sbatch
```

Light re-check later (login node OK):

```bash
export ARCSV2_SIF=$SCRATCH/arcsv2/sif/arcsv2_1.6.1.sif
bash /path/to/arcsv2/scripts/hpc/smoke_sif.sh
```

If the package is private, create `$HOME/.singularity/docker-config.json` (or
export `SINGULARITY_DOCKER_USERNAME` / `SINGULARITY_DOCKER_PASSWORD`) before
pulling. Prefer making `arcsv2` public on GHCR.

## 1. Image on other clusters (optional local build)

Grace cannot build images. On a laptop or build node with Docker:

```bash
cd /path/to/arcsv2
bash scripts/build_apptainer.sh ~/arcsv2.sif   # uses singularity if present
scp ~/arcsv2.sif login.cluster:/path/
```

Or skip building and `singularity pull` the GHCR tag as above.

Smoke-test:

```bash
singularity exec arcsv2_1.6.1.sif arcsv2 --version
```

## 2. Data layout on the cluster

Match the ML-0 layout (recommended). On Grace use `$SCRATCH/arcsv2/...`:

| Path | Content |
|:--|:--|
| `…/ref/` | GRCh38 FASTA, `gaps.bed`, filter BEDs |
| `…/data/samples/<ID>/reads/<ID>.cram` | CRAM + `.crai` |
| `…/haps/<ID>.hap{1,2}.fa.gz` | HPRC haplotypes |
| `…/haps/<ID>.hap{1,2}.mmi` | minimap2 indexes (build just-in-time; delete after label) |

Copy or generate `workflows/config/label.yaml` with **absolute host paths** under
`/scratch/user/dtquocbao/arcsv2/...` for `data_root`, `reference`, `gaps`,
`cram`, and the filter BEDs. Snakemake checks inputs on the host and the shell
runs inside the container, so those paths must be identical in both places.
Grace Singularity auto-binds `/scratch`. Do not remap the tree to `/data`.

## 3. Dry-run (no SLURM jobs)

From the repo checkout:

```bash
bash scripts/cluster_dryrun.sh
# or:
snakemake -s workflows/label.smk -n \
  --profile workflows/profiles/slurm \
  --configfile workflows/config/label.cluster_smoke.yaml
```

You should see jobs with `--mem`, `--cpus-per-task`, and `--time` derived from
`workflows/label.smk`.

## 4. Submit whole-genome labeling

```bash
export ARCSV2_SIF=$SCRATCH/arcsv2/sif/arcsv2_1.6.1.sif
# Snakemake 8 may say "apptainer"; point it at Grace's singularity binary:
export APPTAINER_CMD="${APPTAINER_CMD:-singularity}"
export SINGULARITY_CMD="${SINGULARITY_CMD:-singularity}"
cd /path/to/arcsv2

snakemake -s workflows/label.smk \
  --profile workflows/profiles/slurm \
  --configfile workflows/config/label.yaml \
  --jobs 200 \
  --software-deployment-method apptainer
```

Snakemake 8 accepts `conda`, `apptainer`, and `env-modules` (not `singularity`).
Point the apptainer method at Grace's binary with `APPTAINER_CMD=singularity`.
The image comes from `container:` in `workflows/label.smk` (`ARCSV2_SIF`, or
`container:` in the config). Leave `container:` empty and export `ARCSV2_SIF`
to the `.sif` path. Do not pass `--bind $SCRATCH/arcsv2:/data`. Cluster configs
use the same absolute `/scratch/user/dtquocbao/arcsv2/...` paths on the host and
inside the container.

Edit `workflows/profiles/slurm/config.yaml` defaults:

- `slurm_partition` — e.g. `short` / `medium` / `long` on Grace
- `slurm_extra` — e.g. `--account=132702112504`

Retries: profile sets `restart-times: 2` and `attempts: 3`.

Logs: `logs/slurm/{rule}-<jobid>.out|err`

## 5. Alternative: `arcsv2 batch` without Snakemake

For call-only cohorts (no hap labeling):

```bash
singularity exec "$ARCSV2_SIF" \
  arcsv2 batch --manifest samples.tsv \
  --out /scratch/user/dtquocbao/arcsv2/results/batch1 --chroms chr20
```

Singularity already binds `/scratch`. Point `--out` at that absolute path; do
not bind the tree onto `/data`.

Status table: `…/results/batch1/batch_status.tsv` (see U4).

## 6. Checklist for new lab members

1. Module load / PATH: `snakemake`, `sbatch`, **`singularity`**.
2. Pull `arcsv2_1.6.1.sif` from GHCR **on a compute node** with `WebProxy` (Grace).
3. Fill `label.yaml` (`reference`, `gaps`, `cram` dict, `samples`, `chroms`).
4. Run `scripts/cluster_dryrun.sh` or `snakemake -n` with your config.
5. Submit with `--jobs` capped to your fair-share; keep mem/cores ≤ ~7.5 GB/core on Grace to avoid extra SU charges.
6. After labels finish, `ml-dataset` output is under `labels/` (see config).

## 7. Resource formulas

| Rule | Memory | Time |
|:--|:--|:--|
| `call` | `max(4000, ceil(chrom_Mb × 30.68 × 1.2))` MiB | ~12 min/Mb, floor 3 h |
| `label` | `label_mem_mb` from config (default 32000) | 8 h |
| `index_haps` | 16000 MiB | 6 h |

Tune `mem_mb_per_mb`, `label_mem_mb`, and `threads` in `label.yaml` if jobs OOM or timeout.

## 8. Related scripts

| Script | Role |
|:--|:--|
| `scripts/publish_ghcr.sh` / `.ps1` | Push `arcsv2:ml0` → GHCR `:1.6.1` / `:latest-dev` / `:ml0` |
| `scripts/hpc/pull_sif.sh` | Compute-node pull + smoke |
| `scripts/hpc/pull_sif.sbatch` | Same as `#SBATCH` job |
| `scripts/hpc/smoke_sif.sh` | Light version/native check |
| `docs/CURSOR_TASK_HPC0.md` | Full Grace/FASTER task plan |
