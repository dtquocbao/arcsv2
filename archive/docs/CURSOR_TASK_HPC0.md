# Cursor task: HPC0 — set up ARC-SV2 on Texas A&M HPRC (Grace + FASTER)

## Allocation (approved 2026-10)
| System | Account | SUs | Valid |
|:--|:--|--:|:--|
| **Grace** | `132702112504` | 20,000 | 2026-08-01 → 2027-08-31 |
| **FASTER** | `142702113847` | 20,000 | 2026-08-01 → 2027-08-31 |

**Charging:** 1 SU = 1 *effective core*-hour.
- On a non-exclusive job, effective cores = max(cores requested, memory requested as a share of the node's cores). On a 48-core / ~360 GB node that is ~7.5 GB per core.
- **Exclusive jobs are charged the whole node.**
- **GPU on FASTER costs 64–128 SU/hour.**

Ground rules:
- Run no exclusive jobs and no GPU jobs unless I approve.
- Check spend with `myproject` and `myproject -j <acct>`, and use `maxconfig` to preview a job's charge.

**Storage**
| Area | Quota | Notes |
|:--|:--|:--|
| `$HOME` | 10 GB / 10k files | no big files and no conda envs here |
| `$SCRATCH` | 1 TB / 250k files | not backed up; shared with FASTER |

**Cite HPRC in publications:** https://hprc.tamu.edu/research/citations.html

**Who runs what:** login requires Duo 2FA, so I (Bao) run the commands, or Cursor uses Remote-SSH with my approval. Write every step as copy-paste commands plus a script under `scripts/hpc/`.


## Verified on Grace (2026-10-16, by Bao)
- Login works (grace5). Default account set to `132702112504`. **The PI on the allocation is Zhou, Bo.** 20,000 SU, 0 used.
- Partitions: short 2 h · medium 1 day · long 7 days · xlong 21 days · bigmem 2 days · gpu / gpu-a40 4 days. Compute nodes `c001–c800` have ~375 GB RAM.
- **Login nodes: at most 8 cores across all sessions, or the account is suspended.** Run nothing heavy on login nodes. That includes `singularity pull`, which HPRC says must run on compute nodes.
- **The container runtime is `singularity` (SingularityCE 4.5.0 at `/usr/bin/singularity`), not `apptainer`.** Use `singularity` in every script and doc; keep `apptainer` only as a fallback alias.
- **Images can't be built on Grace.** Only `singularity pull docker://...` works, so the GHCR image is required (H1).
- **Compute nodes have no internet by default.** The test job's `curl` to ghcr.io returned nothing. For pulls and downloads, run `module load WebProxy` inside the job. Inside containers, export `SINGULARITYENV_http_proxy` / `SINGULARITYENV_https_proxy` (Grace proxy `10.73.132.63:8080`).
- Set `export SINGULARITY_CACHEDIR=$TMPDIR/.singularity` before pulling, so the cache doesn't fill `$HOME`.
- Singularity binds `/scratch`, `$HOME`, `/tmp` and `$PWD` automatically.
- HPRC documents Aspera for 1000 Genomes downloads (kb/Software/Bioinformatics/Download_data). Consider it for CRAM fetches.

H1 pull job (run as an `sbatch` or `srun` on `short`, 4 cores / 30 GB / 1 h):
```
module load WebProxy
export SINGULARITY_CACHEDIR=$TMPDIR/.singularity
cd $SCRATCH/arcsv2/sif
singularity pull arcsv2_1.6.1.sif docker://ghcr.io/dtquocbao/arcsv2:1.6.1
singularity exec arcsv2_1.6.1.sif arcsv2 --version
singularity exec arcsv2_1.6.1.sif python -c "import arcsv2_native; print('native ok')"
```

## Rules
- Code must be identical to the laptop: same git tag, same container.
- The golden check must PASS on Grace before any real run.
- Never keep per-sample BAM/CRAM after its sample finishes (1 TB quota).
- Snakemake and its envs live in the container or in `$SCRATCH`, never in `$HOME`, because of the file-count quota.
- Stop and report after each step.

---

## H1. Container on Grace
1. Publish the frozen image `arcsv2:ml0` (caller-v1.6.1) to GHCR as
   `ghcr.io/dtquocbao/arcsv2:1.6.1` (also tag `latest-dev` and `ml0`).
   Laptop: `bash scripts/publish_ghcr.sh` or `.\scripts\publish_ghcr.ps1`.
   Make the package **Public** so Grace can pull without a token.
2. On Grace: create dirs on the **login** node, then **pull only inside a
   compute job** (`short`, 4 cores / 30 GB / 1 h) with `module load WebProxy`
   and `SINGULARITY_CACHEDIR=$TMPDIR/.singularity`. See the “Verified on Grace”
   block above and `docs/CLUSTER.md`. Helper: `scripts/hpc/pull_sif.sh` /
   `pull_sif.sbatch`.
3. Smoke test (after pull; light checks OK on login via `scripts/hpc/smoke_sif.sh`):
   ```
   singularity exec arcsv2_1.6.1.sif arcsv2 --version
   singularity exec arcsv2_1.6.1.sif python -c "import arcsv2_native; print('native ok')"
   ```
4. Clone the repo at tag `caller-v1.6.1` into `$SCRATCH/arcsv2/repo` (or `$HOME`, since code is small).

## H2. Snakemake + SLURM profile for Grace
- Run Snakemake from the login node, or from a tiny long-running job if login-node policy requires it (check HPRC policy). Every rule runs inside the `.sif`.
- Update `workflows/profiles/slurm/config.yaml` with a `grace` variant:
  - `--account=132702112504`
  - partition by walltime: `short` ≤2 h, `medium` ≤1 day, `long` ≤7 days
  - per-rule `cpus-per-task` and `mem`; keep **mem/cores ≤ 7.5 GB** to avoid extra SU charges. Labeling at 11 GB → 2 cores is fine.
  - jobs and latency settings
- Add `--account=142702113847` as a `faster` profile variant for later.
- `snakemake -n` dry-run must show the correct account, partition, mem and time.

## H3. Reference + golden check on Grace
1. Fetch reference, BEDs and the two golden samples' chr20 BAMs (HG01884, NA19700) into `$SCRATCH/arcsv2`, using `fetch_data.py` with checksum verification.
2. Run the chr20 golden config through the Grace profile.
3. Run `check_golden.py` (fast + full) vs `export/v2fix11`. **It must PASS.** Copy the `export/v2fix11` golden fixtures to Grace with `scp`/`rsync`.
4. Record per-job wall time, peak RSS (`seff <jobid>`) and **SUs used**.

## H4. One whole genome: measure real cost
- HG01884, chr1–22 + X, threads 4 vs 8 per chromosome job, using the Grace profile.
- Report:
  - total SUs per genome (call + label)
  - wall time with chromosomes in parallel
  - `$SCRATCH` peak usage
  - any failures
- Compare with the laptop run (caller-v1.6.1 output must be identical for the same sample).
- From this, estimate genomes per 20k SUs, and the SUs needed for (a) 100 more HPRC samples and (b) the 10k Mt. Sinai cohort, as input for a larger allocation request.

## H5. Storage plan
1 TB doesn't hold 40+ samples of assemblies (~15 GB) plus `.mmi` indexes (~14 GB):
- build `.mmi` per sample just before labeling and delete it after
- keep only results and labels (~12 MB/sample)
- stream CRAM → per-chromosome BAM → delete

Also draft (do not send) an email to HPRC help asking for a `$SCRATCH` increase to ~5 TB, justified by the H4 numbers. I will review and send it.

## H6. Data sync
- Copy `export/ml0/` (the ML dataset, ~200 MB) to `$SCRATCH/arcsv2/ml/` for ML-2 training on Grace CPU nodes. Tree models on CPU are cheap; no GPU.
- Document the laptop ↔ Grace sync commands (rsync over ssh; Globus if HPRC offers it).

## Deliverables
- `scripts/hpc/` (setup, pull, dryrun, golden, wg_one)
- `workflows/profiles/slurm-grace/`
- `docs/CLUSTER_TAMU.md` (copy-paste guide)
- `bench/RESULTS.md` section "HPC0 Grace": golden PASS, SUs per genome, wall times
- the draft HPRC storage email

Then stop.
