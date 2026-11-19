# ARC-SV2

Clean-room reimplementation of ARC-SV (Zhou, Arthur, Guo et al., Cell 2024) and its machine-learning filter. The method spec is `../ARC-SV2 - From-Scratch Implementation Roadmap.md`.

Paper mode is the default. Every numeric default in `python/arcsv2/config.py` comes from roadmap §4.12. Extensions in §8 stay off unless a flag turns them on. C5 optional smoove-style noise filters (`--filter-high-nm`, `--filter-both-mates-clipped`) also default off so golden scores stay identical.

## Layout

- `python/arcsv2/` — CLI, configuration, caller, features, labels, training, VCF, and the results dashboard
- `core/arcsv-core/` — Rust kernels used by the default `engine=rust` path through `arcsv2_native`
- `tests/` — pytest
- `scripts/simulate.py` — deterministic paired-end simulator
- `scripts/simulate_haplotypes.py` — two haplotypes, truth, and reads for label tests
- `workflows/` — Snakemake cohort skeleton, including `label.smk`
- `tests/oracle/` — install and run the original ARC-SV as the parity oracle

The Python engine is the reference implementation of §4. `arcsv2 call` defaults to `engine=rust` via `arcsv2_native`. The Rust crate implements the same kernels and is checked by `cargo test`. SAM input, features, and training run as pure Python on Windows. BAM/CRAM (`pysam`) and haplotype labeling (`mappy`) need the Linux container in `docs/CROSS_PLATFORM.md`.

## Install

```bash
pip install -e ".[dev]"
cargo test --manifest-path core/arcsv-core/Cargo.toml
pytest -m "not io and not label"
arcsv2 --version
```

Install the dashboard in its own environment, not Anaconda `base`. Pip can pull NumPy 2.x into `base` while conda's `pyarrow` is still built for NumPy 1.x, and Streamlit then fails to import.

```powershell
conda create -n arcsv2 python=3.11 -y
conda activate arcsv2
pip install -e ".[dashboard]"
```

On Windows, macOS, or Linux, `docker compose -f compose.yaml build cli` (or `.\scripts\arcsv2.ps1 build`) then `scripts/arcsv2.ps1` / `scripts/arcsv2.sh` runs the full pipeline. Compose project name is `arcsv2`; data lives in the `arcsv2-data` volume at `/data`. Two haplotype indexes need well over 16 GB of RAM.

### Quick start (a) ARC-SV Studio app (no git needed)

1. Install [Docker Desktop](https://docs.docker.com/desktop/) and start it.
2. Download the **Studio release zip** (launchers at the top level; no source checkout).
   Maintainers build it with:
   `.\scripts\package_studio_release.ps1 -IncludeImage` (Windows) or
   `INCLUDE_IMAGE=1 bash scripts/package_studio_release.sh` (macOS/Linux).
   That also writes a Docker image `.tar` to send alongside the zip when GHCR is not available.
3. Unzip. If you received an image `.tar`, load it once:
   **Windows** `Load Studio Image.bat` · **macOS** `Load Studio Image.command` · **Ubuntu** `./Load Studio Image.sh`
4. Start: **Windows** `Start ARC-SV.bat` · **macOS** `Start ARC-SV.command` · **Ubuntu** `./Start ARC-SV.sh`
   (On Ubuntu after unzip from Windows: `chmod +x "Start ARC-SV.sh" "Load Studio Image.sh" scripts/*.sh`.)
5. Accept the default data folder (**Documents/ARC-SV**) or pick another in the folder dialog (no typing).
6. Wait for the browser to open → **Try example** (runs a real ~2–5 min call) → **Results**.

The launcher uses a local `arcsv2:studio` image if present, otherwise loads the sibling `.tar`, otherwise **pulls** `ghcr.io/arcsv2/arcsv2:studio` (override with `ARCSV2_STUDIO_IMAGE`). On Apple Silicon the launcher runs the `linux/amd64` image with Docker emulation.  
Developers with a checkout: `Start ARC-SV.bat --build` or `./scripts/start-arcsv.sh --build`.

Command line: **`arcsv2 run`**. Cluster: **`docs/CLUSTER.md`**. Fixes: **`docs/TROUBLESHOOTING.md`**.

### Quick start (b) Command line

```bash
# inside the arcsv2 container / env
arcsv2 run --input sample.cram --reference GRCh38 --out results/run1 --chroms chr20 --threads auto
arcsv2 batch --manifest samples.tsv --out results/batch1 --chroms chr20
```

### Quick start (c) Cluster (SLURM + Singularity)

See **`docs/CLUSTER_TAMU.md`**. On Texas A&M Grace: pull
`ghcr.io/dtquocbao/arcsv2:1.6.1` inside a compute job with `module load WebProxy`
(no image builds on Grace).

### Chr20 benchmark (original vs ARC-SV2)

| Caller | Validated / labeled | Precision |
|:--|--:|--:|
| Original ARC-SV | 243 / 449 | 54.1% |
| ARC-SV2 (`caller-v1.6.1` / `v2fix11`) | 222 / 376 | 59.0% |

Overall Fisher *p* ≈ 0.16; after Dr. Bo's XGBoost: 89.8% vs 82.1% (*p* ≈ 0.033).
Details: `docs/ORIGINAL_COMPARISON.md`, report
[artifact](https://claude.ai/artifact/9tCLAfkKAsnH4ShswQpuZW), scripts under
`scripts/analysis/`. Large tables: `docs/DATA.md`.

### Citation

Zhou et al., Cell 2024 (ARC-SV). Please cite the original ARC-SV paper when using
this software. Cluster work used Texas A&M High Performance Research Computing
(HPRC); cite HPRC as required: https://hprc.tamu.edu/research/citations.html

## Quick start (M2 on a laptop)

Heavy tools run in Docker. From the repo root:

**Windows (PowerShell)** — `scripts/arcsv2.ps1`

```powershell
./scripts/arcsv2.ps1 build
./scripts/arcsv2.ps1 fetch --dry-run --n-samples 2 --chroms chr20
./scripts/arcsv2.ps1 fetch --n-samples 2 --chroms chr20 --yes
./scripts/arcsv2.ps1 --entrypoint snakemake `
  -s workflows/label.smk `
  --configfile workflows/config/label.auto.yaml -j 2
./scripts/arcsv2.ps1 dashboard   # defaults to /data (Snakemake outputs)
```

**macOS / Linux (bash)** — `scripts/arcsv2.sh` (make it executable once: `chmod +x scripts/arcsv2.sh`)

```bash
./scripts/arcsv2.sh build
./scripts/arcsv2.sh fetch --dry-run --n-samples 2 --chroms chr20
./scripts/arcsv2.sh fetch --n-samples 2 --chroms chr20 --yes
./scripts/arcsv2.sh --entrypoint snakemake \
  -s workflows/label.smk \
  --configfile workflows/config/label.auto.yaml -j 2
./scripts/arcsv2.sh dashboard   # defaults to /data (Snakemake outputs)
```

On Apple Silicon, if a dependency lacks arm64: `docker compose build --platform linux/amd64 cli`.

CRAMs are streamed by chromosome unless you pass `--fetch-cram`. Generated
`workflows/config/samples.tsv` and `label.auto.yaml` are gitignored. Details:
`docs/CROSS_PLATFORM.md`. Record wall time in `bench/RESULTS.md` after
Snakemake. The dashboard **Guide** tab has the same steps with a
Windows / macOS·Linux switch.

## Chr20 benchmark (caller-v1 golden)

One-command re-run against the frozen baseline (`export/v2fix11/`, tag
`caller-v1`). Reuses BAMs/haps already under `/data`:

```bash
# inside arcsv2:dev (or via scripts/arcsv2.ps1 --entrypoint snakemake ...)
snakemake -s workflows/label.smk \
  --configfile workflows/config/label.caller-v1.yaml -j 1
python scripts/check_golden.py \
  --results /data/results_caller_v1 \
  --labels /data/labels_caller_v1
```

Fast locus check (minutes): `python scripts/check_golden.py --fast`.
Default config `workflows/config/label.yaml` uses the same caller settings
(writes to `results/` / `labels/`). Method notes: `docs/METHOD.md`,
deviations: `docs/PARITY.md`, history: `docs/HISTORY.md`.
SLURM + Apptainer: `docs/CLUSTER.md`.

## Call

```bash
# Whole CRAM/BAM (primary chroms + chrX; caches GRCh38 under ~/.cache/arcsv2):
arcsv2 run --input sample.cram --reference GRCh38 --out results/run1 --chroms chr20 --threads auto
# Many samples (isolated failures, batch_status.tsv, per-sample qc/ and FAILED.txt):
arcsv2 batch --manifest workflows/config/batch_run.tsv --out results/batch1 --chroms chr20

arcsv2 call -i sample.sam -R ref.fa -G gaps.bed -o out/ --region chr20 --seed 204178949
# Paper BEDs (once): python scripts/fetch_filter_beds.py --dest /data/ref
arcsv2 merge-filter out/*/arcsv_out.tab \
  --simple-repeats /data/ref/simple_repeats.bed \
  --segdups /data/ref/segdups.bed \
  --blacklist /data/ref/blacklist.bed \
  --gaps /data/ref/gaps.bed \
  -o merged.tab
arcsv2 ml-label --calls calls.tab -R ref.fa --hap1 hap1.fa --hap2 hap2.fa --sample SAMPLE -o labeled.tsv
# On a laptop, add --sequential so only one haplotype index is in RAM.
arcsv2 ml-annotate --calls merged.tab --segdups segdups.bed -o annot.tsv
arcsv2 ml-train --labels labels.tsv -o model_v1/
arcsv2 ml-predict --model model_v1 --calls merged.tab -o final.tab
```

`--mode paper` sets the §4.12 defaults, including a path cap of 1,000. `-R` supplies chromosome length and `-G` marks assembly gaps. `call` writes `graph_log.txt` beside the table and prints stage progress to stderr (BAM load, discover, build_graph, genotype). The table appends `blocks` (locus intervals, breakpoint gaps left uncovered), `blocks_point` (the same blocks closed at the soft-clip base, or at the interval midpoint when there is no clip), and `bp_point_source` (`clip` or `mid` per internal junction).

## Labels

`ml-label` builds the alt and reference sequences from `blocks_point` and keeps 150 bp of each outer block next to the event. A single soft clip inside an uncovered interval refines the point estimate. For `mid` junctions, labeling tries positions inside that interval and records `label_bp_shift` when a shift validates the call. It aligns both to the two haplotype assemblies with minimap2. A call is 1 when the alt passes on at least one haplotype and the alt and reference do not both pass on the same haplotype, so a heterozygote whose other haplotype is reference is 1. Insertions (`_` or `=`) are NA. The dashboard block drawing still uses `blocks`. Details are in `docs/LABELS_PLAN.md` and `docs/METHOD.md`.

Labeling needs the Linux container (`docs/CROSS_PLATFORM.md`). Inside it, tools are on `PATH`, including for `bash -lc` and a replaced entrypoint, so omit `--use-conda`. The wrappers set `PYTHONPATH=/work/python`.

CI runs a 200 kb, 5-event simulation through `workflows/label.smk` and `scripts/check_label_pilot.py`. Every called non-insertion event must be label 1, and at least 3 of the 4 non-insertion events must be called. A miss is printed. On seed 7 the deletion, inversion, and inv.del are label 1; the insertion is NA; the heterozygous duplication is reported as missed.

A local dress rehearsal, 40 events on 2 Mb, is `make pilot-large` (or `python scripts/pilot_large.py`) inside the container. It writes `labels/dataset.tsv`. Point the dashboard at the repo to fill the ML tab: `arcsv2 dashboard --results .`. That run is not an M2 timing record.

## What is checked in tests

Library KDE and robustness π, discordant-pair cliques and FDR, split-read types, soft-clip order invariance, uniform breakpoint intervals, affected length, the complex-type table, a constructed deletion call, `blocks_point`, 150 bp label flanks, the pilot check (label 1 on called events, misses reported), feature rows, label criteria, merge-filter allele frequency, VCF records, a deterministic simulator checksum, and dashboard truth-matching (correct, wrong genotype, wrong type, missed, false call). The container job also runs the mappy and pysam tests and the 200 kb label workflow.

Oracle parity against `example/expected_output.tab` needs the original install and the example BAM. See `docs/PARITY.md`.

## Data (M2 pilot)

Same steps as [Quick start](#quick-start-m2-on-a-laptop). Public catalogs →
config → download into Compose volume `/data`. Keep `samples.tsv.example` as
the committed template.

## ML-0 S4 long run (pause / resume / unstick)

S4 uses two detached containers (`arcsv2:ml0`, scripts under `scratch/`):
`arcsv2-ml0-s4` (Snakemake) and `arcsv2-ml0-s4-prefetch` (capped CRAM download).
Data stays in the `arcsv2-data` volume. On Windows, run these from a repo path
**without an apostrophe** (e.g. the `C:\Users\dtquocbao\arcsv2-work` junction).

**Never** `docker volume prune` or `docker compose down --volumes` — that
destroys `arcsv2-data`.

**Pause** (safe; finished BAMs/calls/labels remain on `/data`):

```powershell
docker stop arcsv2-ml0-s4 arcsv2-ml0-s4-prefetch
```

**Resume** (same containers still exist):

```powershell
docker start arcsv2-ml0-s4-prefetch arcsv2-ml0-s4
```

**Clear stuck** (`LockException` on `/work`, or wrapper retrying every 90s
without progress). Stop containers, remove stale locks, then recreate:

```powershell
docker rm -f arcsv2-ml0-s4 arcsv2-ml0-s4-prefetch
docker run --rm --entrypoint bash `
  -v arcsv2-data:/data -v ${PWD}:/work -w /work arcsv2:ml0 `
  -lc "rm -rf .snakemake/locks; mkdir -p .snakemake/locks; snakemake -s workflows/label.smk --configfile workflows/config/label.ml0_s4.yaml --unlock || true"
docker run -d --name arcsv2-ml0-s4-prefetch --restart unless-stopped --entrypoint bash `
  -v arcsv2-data:/data -v ${PWD}:/work -w /work -e PYTHONPATH=/work/python arcsv2:ml0 `
  -lc "bash archive/ml0/ml0_s4_prefetch.sh >> bench/logs/ml0/s4_prefetch_host.log 2>> bench/logs/ml0/s4_prefetch_host.err.log"
docker run -d --name arcsv2-ml0-s4 --restart unless-stopped --entrypoint bash `
  -v arcsv2-data:/data -v ${PWD}:/work -w /work -e PYTHONPATH=/work/python arcsv2:ml0 `
  -lc "bash archive/ml0/ml0_s4_run.sh >> bench/logs/ml0/s4_host.log 2>> bench/logs/ml0/s4_host.err.log"
```

S4 finished (40/40 labeled). One-shot helpers live under `archive/ml0/`.
Dataset: `docs/ML_DATASET.md`, `/data/ml/dataset_v0.parquet`.

Status: `docker ps --filter name=arcsv2-ml0-s4` and
`Get-Content .\bench\logs\ml0\s4_host.log -Tail 20`. Disk:
`docker run --rm -v arcsv2-data:/data alpine df -h /data`.

**Progress snapshot** (BAMs/calls/labels, SV type & genotype tallies, chrom
coverage, live call stages from the host log):

```powershell
docker run --rm --entrypoint bash `
  -v arcsv2-data:/data -v ${PWD}:/work -w /work -e PYTHONPATH=/work/python `
  arcsv2:ml0 -lc "python scripts/ml0_progress.py --matrix"
# optional: python scripts/ml0_progress.py --by-sample HG01884
```

## Local scratch you can delete

These directories are gitignored and safe to remove from a local checkout
(~175 MB when fully populated). They are not needed to rebuild or re-run
against `/data`:

- `pilot/`, `pilot-large/`
- `work/`, `results/`, `labels/`, `haps/`, `target/`
- `.snakemake/`

Do not delete `/data` contents in the Compose volume.

## Dashboard

`arcsv2 dashboard` opens a local web app at <http://localhost:8501>. It only reads result files, so it runs on Windows, macOS, and Linux without `pysam` or `mappy`. The Compose wrappers default to `--results /data` (where fetch/Snakemake write). A bare `arcsv2 dashboard` with no `--results` shows the bundled 2 Mb demo. A run is any directory that contains `arcsv_out.tab` or `merged.tab`. Optional `truth.tsv` and `run.json` beside the table score accuracy and fill the sidebar. `dataset.tsv`, `dataset.summary.json`, and `manifest.json` anywhere under that folder fill the ML tab.

Tabs: **Guide** (step-by-step fetch → call → label → dashboard), Overview, How it works, Call explorer, ML & labels, and Project status. File layout and truth matching are in `archive/docs/DASHBOARD.md`.

```powershell
./scripts/arcsv2.ps1 dashboard
# or: docker compose --profile dashboard up dashboard
```
