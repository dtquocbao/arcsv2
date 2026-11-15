# Running ARC-SV2 on Windows, macOS and Linux

Status: container runtime is in the tree (2026-09-29). `conda-lock.yml` is not generated, the base image is pinned by tag (`mambaorg/micromamba:2.0.5`), and the M2 sample pilot has not been run.

## 1. The problem

| Component | Linux | macOS | Windows (native) |
|:--|:--|:--|:--|
| Core caller on SAM, features, training (pure Python + numpy/sklearn) | yes | yes | yes |
| `pysam` (BAM/CRAM input) | yes | yes | **no** (no Windows builds) |
| `mappy` / `minimap2` (labels) | yes | yes | **no** |
| `samtools` | yes | yes | **no** (not native) |
| Snakemake shell rules (bash, pipes) | yes | yes | fragile |

A native Windows install can never run the full pipeline. There are three possible routes:

- a Linux container everywhere;
- WSL2 on Windows plus conda on Linux/macOS;
- two separate native setups.

**Decision: one Linux container image is the supported runtime on every OS.** Windows uses it through Docker Desktop with the WSL2 backend, and macOS and Linux use the same image. Cursor opens the repo inside the same image as a Dev Container. The pure-Python core stays installable natively everywhere (`pip install -e .`) for quick edits and unit tests.

Why this route:
- It's the only one that is identical on all three OSes.
- Tool versions are pinned in one place.
- It matches the HPC use later: the same image converts to Apptainer/Singularity.

## What is in the tree

- `environment.yml` on conda-forge and bioconda: Python 3.11, the scientific and ML packages, `pysam`, `mappy`, `minimap2`, `samtools`, Snakemake, and the dashboard libraries. It does not `pip install -e .` during the conda solve. The image installs the package with pip after the solve. `conda-lock.yml` is still absent; the lock command is a comment at the top of `environment.yml`.
- `Dockerfile` from `mambaorg/micromamba:2.0.5` (tag, not digest). Includes
  conda-forge `rust` + `maturin` and bakes `arcsv2_native` into the image
  (`OPENSSL_NO_VENDOR=1`, conda libssl; build needs `pkg-config`, `binutils`,
  `clangdev`/`libclang`, `patchelf`). `CallConfig.engine` defaults to `rust`
  after R1/R2 golden PASS; missing native raises (use `--engine python` or
  rebuild `arcsv2:dev`). Full BAM/call stays in Docker; host-only Python is for
  edits/tests.
- `compose.yaml` with fixed project name `arcsv2` (cli / dashboard / dev services) and named volume `arcsv2-data`. Wrappers call Compose; Docker Desktop groups containers under **arcsv2**. Optional host bind: copy `compose.override.yaml.example` → `compose.override.yaml`.
- **ML-0 freeze:** `compose.ml0.yaml` overrides the image to `arcsv2:ml0` (built from `caller-v1.6.1`). All ML-0 S3/S4 runs must use `arcsv2:ml0` via `docker compose -f compose.yaml -f compose.ml0.yaml ...`. Default `compose.yaml` still uses `arcsv2:dev` for day-to-day work.
- `.devcontainer/devcontainer.json` uses the Compose `dev` service (`/work` + `/data`).
- `scripts/arcsv2.ps1` and `scripts/arcsv2.sh` call Compose (`build`, `dashboard`, `shell`, `fetch`, `--entrypoint`, or default `cli`).
- `.gitattributes` forces LF for text and CRLF for PowerShell. `.github/workflows/ci.yml` runs the core tests on Ubuntu, macOS, and Windows, and the full pytest suite plus the 200 kb `label.smk` pilot inside the image.
- Inside the container, omit `--use-conda`. Outside it, `snakemake --use-conda` builds `environment.yml` per rule.

## 2. Original spec

The sketches below are what the container was built from. Where they differ from the tree, "What is in the tree" above is the current behavior.

### 2.1 `environment.yml` (single source of versions)

```yaml
name: arcsv2
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - numpy
  - scikit-learn
  - xgboost
  - lightgbm
  - pysam>=0.22
  - mappy>=2.28
  - minimap2>=2.28
  - samtools>=1.20
  - snakemake-minimal>=8
  - pip
  - pip:
      - -e .[dev]
```

Lock it with `conda-lock -p linux-64 -p linux-aarch64 -p osx-arm64 -p osx-64` and commit `conda-lock.yml`. The container builds from the lock file, so every machine gets the same versions.

### 2.2 `Dockerfile` (replace the current one)

- **Base image:** a conda base (e.g. `condaforge/miniforge3` or `mambaorg/micromamba`), pinned by digest. Create the env from `conda-lock.yml`.
- **Rust:** drop `rustup` and the `cargo test` from the runtime image. The Rust crate isn't called from Python yet, and building it makes the image slow and large. Keep a separate `rust` stage or CI job for `cargo test`.
- **Run as a non-root user** with `WORKDIR /work`.
- **Entrypoint:** `ENTRYPOINT ["arcsv2"]`, with `docker run ... --entrypoint bash` for shell work.
- **Apple Silicon:** build `linux/arm64` if every package has an aarch64 build. Otherwise build with `--platform linux/amd64`; it runs under emulation, slower but correct.

### 2.3 `.devcontainer/devcontainer.json` (Cursor / VS Code)

```json
{
  "name": "arcsv2",
  "build": { "dockerfile": "../Dockerfile", "context": ".." },
  "workspaceFolder": "/work",
  "workspaceMount": "source=${localWorkspaceFolder},target=/work,type=bind",
  "mounts": ["source=arcsv2-data,target=/data,type=volume"],
  "overrideCommand": true,
  "postCreateCommand": "pip install -e .[dev] && pytest -q"
}
```

Cursor supports Dev Containers: "Reopen in Container" gives the same shell, Python and tools on every OS.

### 2.4 Wrapper scripts (one command per OS)

- **`scripts/arcsv2.ps1` (Windows PowerShell):**
  - Mounts the repo at `/work` and `$env:ARCSV2_DATA` (default: a named volume `arcsv2-data`) at `/data`.
  - Passes the arguments through to `docker run --rm arcsv2:<version> ...`.
  - Always quotes the paths. This project lives under `D:\...\Bao's windows\...`, so the path has spaces and an apostrophe.
- **`scripts/arcsv2.sh`:** the same for macOS and Linux.
- **Inside the container, paths are always `/work/...` and `/data/...`.** Configs (`label.yaml`) use those container paths, never `D:\...`.

### 2.5 Data location (performance)

- Bind-mounting a Windows drive (`D:\...`) into a WSL2 container is slow for random I/O. CRAM fetches, `.mmi` loading and FASTA indexing all do random I/O.
- Keep large inputs and outputs in the Docker **named volume** `arcsv2-data` or on the WSL2 ext4 filesystem (`\\wsl$\...`): CRAMs or streamed regions, assemblies, minimap2 indexes, results.
- Only the code lives in the bind-mounted repo.
- Keep data out of synced folders (OneDrive, Obsidian vault sync). Large binary files and indexes there cause sync churn and file-lock errors.
- **Memory:** two haplotype minimap2 indexes plus the reference need well over 16 GB. On Windows set `memory=` in `%UserProfile%\.wslconfig`, e.g. 24-32 GB if the laptop allows, or label on an HPC node with the same image through Apptainer.

### 2.6 `.gitattributes`

```
* text=auto eol=lf
*.ps1 text eol=crlf
*.bat text eol=crlf
*.bam binary
*.cram binary
*.mmi binary
```

This stops Windows checkouts from turning `.sh`, `.smk` and test fixtures to CRLF, which breaks shell rules and byte-identical output tests.

### 2.7 Code hygiene (already mostly done; keep it that way)

- Use `pathlib` for every path. Never build paths with `/` string joins or split on `/`.
- Open text with `encoding="utf-8", newline="\n"`. The writers already do this.
- Call no shell from Python (`subprocess` with argument lists only, no `shell=True`).
- Keep multiprocessing spawn-safe: top-level worker functions, picklable arguments, the CLI entry under `if __name__ == "__main__"`. Checked in round 7: `threads=3` under `spawn` matches `threads=1`.
- Optional imports (`pysam`, `mappy`) must fail with a clear message that points to the container. Tests `importorskip` them; that's already done.
- In Snakemake, quote every path (`{input:q}`, `{output:q}`), and run rules inside the container or with `--use-conda`.

### 2.8 CI matrix (`.github/workflows/ci.yml`)

| Job | OS | What |
|:--|:--|:--|
| core | ubuntu, macos, **windows** | `pip install -e ".[dev,ml,dashboard]"`; `pytest -m "not io and not label and not slow"`; ruff |
| full | ubuntu | build the image; `pytest` inside the container; `label.smk` dry run; 200 kb simulated label workflow |
| rust | ubuntu | `cargo test`, `cargo clippy -D warnings` |

Pytest markers `io` (needs pysam) and `label` (needs mappy) keep those tests off the Windows core job. The container job runs them.

## 3. Daily use on Windows

```powershell
# one-time
./scripts/arcsv2.ps1 build
# or: docker compose -f compose.yaml build cli

# run a command (Compose project name arcsv2)
./scripts/arcsv2.ps1 call -i /data/HG002/chr20.bam -R /data/ref/GRCh38.fa -G /data/ref/gaps.bed -o /data/out/HG002/chr20 --region chr20

# auto-fetch catalogs → samples.tsv + label.auto.yaml (streams CRAMs by region)
# Both config files are gitignored; catalogs and assemblies live in volume arcsv2-data.
./scripts/arcsv2.ps1 fetch --dry-run --n-samples 2 --chroms chr20
./scripts/arcsv2.ps1 fetch --n-samples 2 --chroms chr20 --yes

# After fetch prints "next: snakemake …":
./scripts/arcsv2.ps1 --entrypoint snakemake -s workflows/label.smk --configfile workflows/config/label.auto.yaml -j 2
# Record wall time / peak RSS in bench/RESULTS.md
```

`label.yaml` / `label.auto.yaml` set `data_root: /data` so `haps/`, `work/`, `results/`, and `labels/` land on the Docker volume (not the NTFS bind mount). Pilot configs keep `data_root: "."`.

Named volumes are created as root. The wrappers run `init-data` (chown to `mambauser`) before `fetch` and other commands. If you hit `PermissionError: /data`, run:

```powershell
docker compose --profile init run --rm init-data
```

To bind a host folder instead of the named volume, copy `compose.override.yaml.example` to `compose.override.yaml` and set e.g. `D:/arcsv2-data:/data` on each service.
The wrappers set `PYTHONPATH=/work/python`, so the mounted checkout is what `arcsv2` imports. `--use-conda` is for a run outside the container, where each rule builds `environment.yml`.

The 200 kb label check, from the repo root inside the container:

```bash
python scripts/simulate_haplotypes.py --out pilot --length 200000 --events 5 --coverage 20 --seed 7
mkdir -p haps
cp pilot/hap1.fa haps/sim.hap1.fa
cp pilot/hap2.fa haps/sim.hap2.fa
: > pilot/gaps.bed
: > pilot/unreliable.bed
snakemake -s workflows/label.smk -j 2 --configfile workflows/config/label.pilot.yaml
python scripts/check_label_pilot.py
```

`make pilot-large` is the 2 Mb local run. It is not CI and not the M2 timing.

In Cursor: open the repo and choose "Reopen in Container". The terminal, tests and debugger all run in Linux.

## 4. Still open

1. Generate `conda-lock.yml` and pin the base image by digest.
2. Run the M2 pilot inside the container (HG002 + 5 HPRC samples, chr20–22) and record wall time and peak RAM in `bench/RESULTS.md`.
