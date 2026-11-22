# ARC-SV2 on Texas A&M clusters

The copy-paste runbook (pull, golden chr20, one genome, storage, sync) is
added with the rest of the cluster docs. This file records the controller
pins that `scripts/hpc/setup_controller.sh` installs.

## Controller environment

Create it on `$SCRATCH`, never in `$HOME` (10k-file quota). Run the installer
from a short compute job:

```bash
srun -A 132702112504 -p short -N 1 -n 1 -c 2 --mem=8G -t 01:00:00 --pty bash
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

## Network

`scripts/hpc/run_label.sh` sets `http_proxy` and `https_proxy` to
`http://10.73.132.63:8080` before Snakemake starts. SLURM forwards them
(`--export=ALL`), and the profile's `apptainer-args` pass them into the
container so `fetch_region` can stream remote CRAMs.

To avoid streaming, prefetch CRAM + `.crai` in one `short` job (Aspera, or
`samtools view` after `module load WebProxy`) under
`$SCRATCH/arcsv2/data/samples/<ID>/reads/` and point `cram:` at those files.
The commands are in `docs/CLUSTER.md`. After `merge`, delete that sample's
BAMs and CRAM (H5). Delete each `.mmi` after `label`.

## Grace gotchas (found during H3, 2026-11-22)

- `srun`/`sbatch` need `-N 1 -n 1` (job_submit rejects jobs without `--nodes` or `--ntasks`) and `--mem`.
- `/sw/local/bin/sbatch` is a site wrapper ahead of `/usr/bin/sbatch`. It exits 0 when a submission is rejected and breaks the plugin's quoted `--wrap`, so Snakemake hangs with "No active jobs". `run_label.sh` puts `/usr/bin/sbatch` first on PATH via `$SCRATCH/arcsv2/bin/slurm-direct`; the manifest records `sbatch=`. Server-side job_submit checks still apply. `sbatch --wrap` typed by hand also goes through the wrapper, so write job files instead.
- The image has micromamba, not `conda`. Snakemake 8 runs `conda info` inside the container for any rule with `conda:`, so `label.smk` sets `conda:` to None when a SIF is in use.
- Interactive `singularity exec` does not bind `/scratch` by default on Grace; pass `-B /scratch` (the Snakemake jobs were fine).
- After an interrupted run: `--unlock`, then `--rerun-incomplete`.
- Long work goes in `sbatch` jobs, not interactive `srun`: portal/SSH disconnects kill `srun --pty`. tmux sessions live on one login node (`ssh grace5.hprc.tamu.edu` to return).
- The OnDemand web shell inserts blank lines on paste; use SSH from a terminal for heredocs.

H3 result on Grace (2026-11-22): fast and full golden PASS (HG01884 280 rows, NA19700 271, dataset 376, precision 0.5904), 14 steps in ~25 min.

## Hand-off: H1, then H3, then H4

Run these on Grace. Clone the pushed `master` (do not check out only the tag):

```bash
git clone https://github.com/dtquocbao/arcsv2.git $SCRATCH/arcsv2/repo
cd $SCRATCH/arcsv2/repo
```

Tag `caller-v1.6.1` is commit `203ee14`, which already contains
`scripts/hpc/run_label.sh`. This hand-off and the Grace profile are on
`master` after that clone. `run_label.sh` records the commit it actually
runs (`git rev-parse HEAD` and `git describe --always --dirty`) in
`$SCRATCH/arcsv2/logs/run_manifest.<timestamp>.txt`. Do not remap paths to
`/data`. Singularity on Grace binds `/scratch`.

`scripts/hpc/pull_sif.sh` pulls `ghcr.io/dtquocbao/arcsv2:1.6.1`. An anonymous
pull returns 403 while the package is private. Export `GHCR_TOKEN` in the
job (the script passes it to Singularity and does not write it down). If the
pull still fails, stop and paste the error with the token redacted. Do not
build an image on Grace.

### H1. Pull the SIF on a compute node

```bash
cd $SCRATCH/arcsv2/repo
mkdir -p $SCRATCH/arcsv2/sif
# Export GHCR_TOKEN in this shell first when the package is private.
# pull_sif.sh passes it to Singularity and does not print or save it.
srun -A 132702112504 -p short -N 1 -n 1 -c 4 --mem=30G -t 01:00:00 --pty \
  bash scripts/hpc/pull_sif.sh
bash scripts/hpc/smoke_sif.sh
```

Then install the controller, still on a compute node. Login nodes are capped
at 8 cores.

```bash
srun -A 132702112504 -p short -N 1 -n 1 -c 2 --mem=8G -t 01:00:00 --pty bash
bash $SCRATCH/arcsv2/repo/scripts/hpc/setup_controller.sh
```

Place the caller-v1 reference and BEDs at:

```text
/scratch/user/dtquocbao/arcsv2/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
/scratch/user/dtquocbao/arcsv2/ref/gaps.bed
/scratch/user/dtquocbao/arcsv2/ref/unreliable.bed
```

and the haplotypes at
`$SCRATCH/arcsv2/data/haps/<ID>.hap{1,2}.fa` (or `.fa.gz`).

### H3. chr20 golden (HG01884 and NA19700)

`threads: 1` matches `workflows/config/label.caller-v1.yaml`.
`label_mem_mb` and `label_threads` come from
`workflows/config/label.grace.example.yaml` and only change the SLURM
request. Write this on Grace:

```bash
cat > $SCRATCH/arcsv2/repo/workflows/config/label.grace.h3.yaml << 'EOF'
container: ""
data_root: /scratch/user/dtquocbao/arcsv2/data
results_dir: results_caller_v1
labels_dir: labels_caller_v1
reference: /scratch/user/dtquocbao/arcsv2/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
gaps: /scratch/user/dtquocbao/arcsv2/ref/gaps.bed
unreliable: /scratch/user/dtquocbao/arcsv2/ref/unreliable.bed
simple_repeats: ""
segdups: ""
blacklist: ""
chroms:
  - chr20
samples:
  - HG01884
  - NA19700
cram:
  HG01884: https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988940/HG01884.final.cram
  NA19700: https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3239913/NA19700.final.cram
sequential: true
threads: 1
label_mem_mb: 30000
label_threads: 4
EOF

cd $SCRATCH/arcsv2/repo
tmux new -s arcsv2
bash scripts/hpc/run_label.sh workflows/config/label.grace.h3.yaml
```

`run_label.sh` exports the Grace proxy and submits with
`workflows/profiles/slurm-grace`. Remote CRAMs stream through that proxy.
To prefetch instead, use the `samtools view` block in `docs/CLUSTER.md` and
point `cram:` at the absolute `.cram` paths.

After the DAG finishes, from the repo:

```bash
export ARCSV2_SIF=$SCRATCH/arcsv2/sif/arcsv2_1.6.1.sif
singularity exec -B /scratch "$ARCSV2_SIF" python scripts/check_golden.py --fast \
  --results /scratch/user/dtquocbao/arcsv2/data/results_caller_v1
singularity exec -B /scratch "$ARCSV2_SIF" python scripts/check_golden.py \
  --results /scratch/user/dtquocbao/arcsv2/data/results_caller_v1 \
  --labels /scratch/user/dtquocbao/arcsv2/data/labels_caller_v1
```

Both must print PASS against `export/v2fix11` before H4.

### H4. One genome, HG01884, chr1–22 and chrX

```bash
cat > $SCRATCH/arcsv2/repo/workflows/config/label.grace.h4.yaml << 'EOF'
container: ""
data_root: /scratch/user/dtquocbao/arcsv2/data
results_dir: results_h4
labels_dir: labels_h4
reference: /scratch/user/dtquocbao/arcsv2/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
gaps: /scratch/user/dtquocbao/arcsv2/ref/gaps.bed
unreliable: /scratch/user/dtquocbao/arcsv2/ref/unreliable.bed
simple_repeats: ""
segdups: ""
blacklist: ""
chroms:
  - chr1
  - chr2
  - chr3
  - chr4
  - chr5
  - chr6
  - chr7
  - chr8
  - chr9
  - chr10
  - chr11
  - chr12
  - chr13
  - chr14
  - chr15
  - chr16
  - chr17
  - chr18
  - chr19
  - chr20
  - chr21
  - chr22
  - chrX
samples:
  - HG01884
cram:
  HG01884: https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/data/ERR3988940/HG01884.final.cram
sequential: true
threads: 4
label_mem_mb: 30000
label_threads: 4
EOF

cd $SCRATCH/arcsv2/repo
# Call only. merged.tab is the production product; index_haps and label stay out.
bash scripts/hpc/run_label.sh workflows/config/label.grace.h4.yaml \
  /scratch/user/dtquocbao/arcsv2/data/results_h4/HG01884/merged.tab
# Labeling is a second run. Call jobs are already done, so this adds index_haps and label.
bash scripts/hpc/run_label.sh workflows/config/label.grace.h4.yaml
```

Save the Snakemake log from the call-only run, then:

```bash
bash scripts/hpc/su_report.sh path/to/snakemake.log
```

That prints per-rule SU from `sacct`. "SU per genome" is fetch + call + merge.

### Paste back

From a login node, after H3 and again after H4:

```bash
seff <jobid>          # one seff per call, label, and index job you care about
myproject -j 132702112504
du -sh $SCRATCH/arcsv2
```

Also paste whether `check_golden.py` printed PASS, any failed rule name, and
the `sbatch` job ids from the Snakemake log.
