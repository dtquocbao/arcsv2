# Cursor task: HPC1 — fix the cluster path, then deploy on Grace (VISION: docs only)

Context: review on 2026-11-15 (Claude). The HPC0 plan (`archive/docs/CURSOR_TASK_HPC0.md`) is still the
right plan, but the cluster code in this repo will not run correctly on its first submission. Fix it
here, test locally, then I (Bao) run H1–H4 on Grace by hand (Duo login).

**VISION is not a target for the caller.** It supports only Pyxis/Enroot containers (no Singularity or
Apptainer), it has only 8×H200 DGX nodes, and usage is counted in GPU-hours. ARC-SV2 uses no GPU. Do
not port the Snakemake workflow to VISION. Step F6 only documents it.

## Rules
- Do not change caller logic, numeric defaults, or golden outputs. This task covers workflow, profile,
  docs, scripts and tests only.
- `scripts/check_golden.py --fast` must still PASS locally when you finish.
- Keep CI green: CI runs Snakemake *inside* the Docker image without a deployment method, so any
  `container:` directive must be harmless there. It is: Snakemake ignores `container:` unless
  `--sdm apptainer` is passed.
- Use `singularity` in scripts and docs. Grace has no `apptainer` binary.
- Stop and report after each step (F1–F7), with the diff summary and test output.

---

## F1. Make rules actually run inside the container
Problem: `workflows/label.smk` and `workflows/Snakefile` have no `container:` directive. In
`docs/CLUSTER.md`:
- `--container-image` is Snakemake's Kubernetes/remote-executor flag;
- `--software-deployment-method singularity` is not a valid value in Snakemake 8 (the choices are
  `conda`, `apptainer` and `env-modules`);
- `--singularity-args` is only an alias.

As written, every job runs on the bare host and fails with `arcsv2: command not found`.

Do:
1. Near the top of both `workflows/label.smk` and `workflows/Snakefile`, after `configfile:`, add:
   ```python
   _SIF = str(config.get("container") or os.environ.get("ARCSV2_SIF") or "").strip()
   if _SIF:
       container: _SIF
   ```
   (Add `import os` to `Snakefile` if it is missing.) This is a global container that applies to
   every rule.
2. Keep the `conda:` directives. They are only used with `--sdm conda`.
3. Add `container: ""` (empty) to `workflows/config/label.yaml.example` / `label.yaml` with a comment.
4. Paths: Snakemake checks input existence **on the host**, and the shell runs **in the container**, so
   paths must be identical in both places. On Grace, Singularity auto-binds `/scratch`. Cluster configs
   must therefore use absolute host paths under `/scratch/user/dtquocbao/arcsv2/...` (`data_root`,
   `reference`, `gaps`, `cram`, BEDs). **Do not remap to `/data`** in cluster mode. Remove the
   `--bind $SCRATCH/arcsv2:/data` advice from the docs.

Test: `snakemake -n --sdm apptainer` on the smoke config prints `singularity exec` lines when
`ARCSV2_SIF` is set. This works without the binary only if Snakemake skips the version check in a
dry-run. If it doesn't, assert on `snakemake --list-rules` plus a unit test that parses the
`container:` line.

## F2. Replace the `cluster-generic` executor with the SLURM executor plugin
Problems:
- `executor: cluster-generic` needs `snakemake-executor-plugin-cluster-generic`, which is not
  installed anywhere.
- `sbatch` is not in the image, so the Snakemake controller cannot run inside the container.

Do:
1. Create `workflows/profiles/slurm-grace/config.yaml`:
   ```yaml
   executor: slurm
   jobs: 100                 # stay under fair-share; raise only after H4
   restart-times: 2
   latency-wait: 120
   printshellcmds: true
   software-deployment-method: apptainer
   apptainer-args: "--cleanenv --env http_proxy=$http_proxy,https_proxy=$https_proxy"  # see F4
   default-resources:
     slurm_account: "132702112504"
     slurm_partition: "medium"   # 1 day; every rule fits after F3
     runtime: 60                 # minutes (int) for the slurm plugin
     mem_mb: 4000
   set-resources:
     dataset:
       slurm_partition: "short"
     merge:
       slurm_partition: "short"
   ```
   Verify each key against the installed plugin version (`snakemake --help` after installing; `runtime`
   is in minutes). If `--cleanenv` breaks `PATH` inside the image, drop it and document why.
2. Create `workflows/profiles/slurm-faster/config.yaml`: the same file with `slurm_account: "142702113847"`.
   Look up FASTER's partition names and leave a `TODO(Bao): confirm partitions` comment instead of
   guessing.
3. Keep `workflows/profiles/slurm/` as a generic template, but switch it to `executor: slurm`, remove
   `slurm_partition=normal` (no such partition on Grace) and the stale singularity comments, and delete
   `slurm-status.py`. The plugin handles job status itself.
4. New `scripts/hpc/setup_controller.sh`. It creates the controller env on `$SCRATCH`, never in `$HOME`
   (10k-file quota):
   ```bash
   # run inside: srun -A 132702112504 -p short -c 2 --mem=8G -t 01:00:00 --pty bash
   module load WebProxy
   export MAMBA_ROOT_PREFIX=$SCRATCH/arcsv2/mamba
   # install micromamba binary into $SCRATCH/arcsv2/bin if missing
   micromamba create -y -p $SCRATCH/arcsv2/envs/smk -c conda-forge -c bioconda \
     python=3.11 snakemake-minimal=8.30 snakemake-executor-plugin-slurm
   ```
   Pin exact versions, and record them in `docs/CLUSTER_TAMU.md`.
5. New `scripts/hpc/run_label.sh <configfile> [extra snakemake args]`. It:
   - activates the env;
   - exports `ARCSV2_SIF`, the proxy vars (F4) and `SINGULARITY_CACHEDIR`;
   - runs `snakemake -s workflows/label.smk --profile workflows/profiles/slurm-grace --configfile "$1"`.

   Document that the controller runs in `tmux` on the login node (it is light: 1 process) **or** as a
   1-core `long` job. Leave a `TODO(Bao): confirm HPRC policy on long-running login processes`.

## F3. Fix walltime and memory requests
Problem: `call_runtime_slurm` is `max(180, Mb*12)` minutes, so chr1 requests about 50 h. The measured
chr1 run takes 31 min at 8 threads (`bench/RESULTS.md`, PERF3). That forces `long`, hurts backfill, and
any request over 48 h would be rejected outright on VISION.

Do, in both `label.smk` and `Snakefile`:
1. Call runtime: `max(call_runtime_floor_min, ceil(Mb * call_min_per_mb))`, with config defaults
   `call_min_per_mb: 0.5` and `call_runtime_floor_min: 60`. That gives chr1 ≈ 124 min.
2. Scale runtime and memory by attempt so retries get more:
   `runtime=lambda wc, attempt: base * attempt` and the same for `mem_mb`. Use Snakemake's `attempt`
   argument. Do not change the formulas' base values.
3. `Snakefile` `call_chrom` has no `runtime`. Add the same one.
4. Make `index_haps` (360 min, 16 GB) and `label` (480 min) configurable (`index_runtime_min`,
   `label_runtime_min`). Keep the current defaults.
5. Grace SU efficiency: Grace charges `max(cores, mem/7.5 GB)`. In a new
   `workflows/config/label.grace.example.yaml`, set:
   - `label_mem_mb: 30000`, `label_threads: 4`
   - `threads: 4` for call
   - comments explaining the 7.5 GB/core rule

   Don't change the laptop defaults.

Tests (`tests/test_cluster_profile.py`, rewrite):
- runtime(chr1) is between 60 and 1440 min; runtime(chr20) ≥ 60;
- attempt 2 doubles runtime and memory;
- the `slurm-grace` profile has `executor: slurm`, the account `132702112504`, and no `normal` partition;
- a dry-run with `--profile workflows/profiles/slurm-grace` on `label.cluster_smoke.yaml` succeeds
  (skip if snakemake or the plugin is missing).

## F4. Network on Grace compute nodes
Compute nodes have no internet. `fetch_region` streams CRAMs (often remote URLs), and that fails
without a proxy.
- In `run_label.sh`, before launching:
  `export http_proxy=http://10.73.132.63:8080 https_proxy=http://10.73.132.63:8080`. SLURM passes the
  environment to jobs (`--export=ALL` is the default), and F2's `apptainer-args` passes it into the
  container.
- Add a `prefetch` alternative to the docs: download CRAM+CRAI to `$SCRATCH/arcsv2/data/` in one
  `short` job (Aspera or `samtools view` with WebProxy) and point `cram:` at the local files. Keep the
  "delete per-sample BAM/CRAM after merge" rule (H5).
- Add a pre-label check to `index_haps`/`label`: run `arcsv2`'s existing `mmi_validate` (or a size
  check plus a smoke `minimap2` query) and fail loudly. This is the HG00140 lesson. If there's already
  a CLI hook for it, wire it in; if not, add `arcsv2 mmi-validate <mmi>` as a thin CLI over
  `python/arcsv2/mmi_validate.py`, and call it at the end of `index_haps`.

## F5. Reproducibility: tag, remote, image digest
The repo was re-initialized: one `init` commit and no tags. Remote `arcsv2` has a broken URL (the text
`git remote add origin https://...` was pasted in as the URL).
1. Run `git remote remove arcsv2` and `git remote rename arcsv2Repo origin`. Confirm with `git remote -v`.
2. Do **not** recreate `caller-v1.6.1` blindly. Run `check_golden.py --fast` and the full check on the
   current tree:
   - if both PASS, tag the commit `caller-v1.6.1` with the annotation "re-tagged after repo re-init;
     golden PASS";
   - otherwise tag it `hpc1-wip` and report.
3. Get the digest of `ghcr.io/dtquocbao/arcsv2:1.6.1` (`docker buildx imagetools inspect`).
   - Record it in `docs/CLUSTER_TAMU.md`.
   - Make `scripts/hpc/pull_sif.sh` pull by digest: `docker://ghcr.io/dtquocbao/arcsv2@sha256:...`,
     with `ARCSV2_GHCR` override kept.
   - Confirm the GHCR package is **Public**.
   - If the image is missing, stop and tell me. Do not publish without asking.
4. `run_label.sh` writes `git describe --always --dirty`, the SIF sha256 and the snakemake version to
   `$SCRATCH/arcsv2/logs/run_manifest.<timestamp>.txt`.

## F6. Docs
1. New `docs/CLUSTER_TAMU.md` (copy-paste guide), in this order:
   - setup_controller;
   - pull SIF (H1);
   - golden chr20 (H3) with `check_golden.py`;
   - one whole genome (H4) with `sacct`/`seff` and the `myproject` SU readout;
   - storage rules (H5);
   - rsync/Globus sync (H6).
2. Rewrite `docs/CLUSTER.md` §4 to use the new profile and `run_label.sh`. Delete `--container-image`,
   `--singularity-args`, `APPTAINER_CMD`/`SINGULARITY_CMD` and the `/data` bind advice.
3. New `docs/CLUSTER_VISION.md`, short:
   - Status: not supported for the caller or labeling. Give the reasons: Pyxis/Enroot only; GPU-only
     DGX H200 nodes; 14 CPUs per GPU by default; GPU-hour accounting (Tier One = 970 GPU-h / 6 mo);
     48 h max walltime; CPU-only jobs not documented.
   - When to use it: a future GPU model (ML-2+).
   - An ad-hoc Pyxis recipe for one sample, for reference only:
     ```bash
     #SBATCH --time=12:00:00 --qos=standard --gres=gpu:1 --cpus-per-task=14
     srun --container-image=ghcr.io#dtquocbao/arcsv2:1.6.1 \
          --container-mounts=/scratch/user/<vision-id>/arcsv2:/scratch/user/<vision-id>/arcsv2 \
          arcsv2 run --input ... --out ... --threads 14
     ```
     Note: Enroot ignores `ENTRYPOINT`, so pass the full command. `$HOME` is not mounted unless
     `--container-mount-home` is set. `$PWD` can't be used in `#SBATCH` lines.
   - A draft email to VISION support (do not send) asking:
     1. Are CPU-only jobs (`--gres=gpu:0`) allowed?
     2. What is the policy on controlled-access human genomic data (dbGaP / institutional DUA)?
     3. What tier is our project, and what is its GPU-hour balance?
4. README:
   - Fix the stale line saying the Rust core is "not wired in": `engine=rust` is the default via
     `arcsv2_native`.
   - Point "Quick start (c) Cluster" at `docs/CLUSTER_TAMU.md`.

## F7. Local verification, then hand-off
1. `pytest -m "not io and not label"` passes, plus the container test job.
2. `snakemake -n --profile workflows/profiles/slurm-grace --configfile workflows/config/label.cluster_smoke.yaml`
   shows `slurm_account=132702112504`, partition `medium`/`short`, sane runtime and memory, and
   `container` set when `ARCSV2_SIF` is exported.
3. `check_golden.py --fast` PASS.
4. Write a hand-off block at the end of `docs/CLUSTER_TAMU.md`: the exact commands I run on Grace for
   H1 → H3 → H4, and what to paste back to you (`seff`, `myproject -j 132702112504`, `du -sh $SCRATCH/arcsv2`).

## Deliverables
- `workflows/label.smk`, `workflows/Snakefile`: container directive, runtime/memory fixes
- `workflows/profiles/slurm-grace/`, `workflows/profiles/slurm-faster/`, updated `workflows/profiles/slurm/`
- `workflows/config/label.grace.example.yaml`
- `scripts/hpc/setup_controller.sh`, `scripts/hpc/run_label.sh`, digest-pinned `pull_sif.sh`
- `tests/test_cluster_profile.py` (rewritten)
- `docs/CLUSTER_TAMU.md`, `docs/CLUSTER_VISION.md`, updated `docs/CLUSTER.md` and README
- git: fixed remote, tag per F5

Then stop.
