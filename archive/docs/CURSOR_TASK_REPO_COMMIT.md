# Cursor task: REPO-COMMIT — clean the repo and commit all code and results

Goal: one clean commit series that someone else, or Dr. Bo, can clone and reproduce from:
- code, docs and configs in git
- small benchmark evidence in git
- big data published outside git, with checksums

No caller behavior changes: `check_golden.py --fast` must PASS at the end.
**Stop after step 1 and after step 6 for my review. Never push or rewrite history without my OK.**

## Current state (2026-11-07)
| Path | Size | Plan |
|:--|--:|:--|
| `dist/` | 16 GB | Studio image tar. Never in git; publish the image to GHCR instead (already gitignored) |
| `target/` | 549 MB | Rust build output; gitignored |
| `export/ml0/` | 200 MB | 40-genome dataset; release asset (step 4) |
| `bench/` | 131 MB | Two logs are 64 MB (`bench/logs/ml0/s4_host.log`, `s4_snakemake.log`); compress or trim (step 3) |
| `pilot-large/`, `pilot/`, `haps/`, `work/`, `results/`, `labels/`, `.snakemake/` | ~210 MB | Local scratch; gitignored; leave on disk |
| `export/v2fix … v2fix10`, `export/results`, `export/labels` | ~15 MB | Superseded intermediate results (step 4) |
| `scratch/` | 46 files | Gitignored; empty it except the 2 helpers |
| Root launchers (`Start ARC-SV.*`, `Load Studio Image.*`) | small | Keep; check that they're referenced in the README |

## 1. Inspect git first (report, change nothing)
- `git status --short | head -100`, `git count-objects -vH`, current branch, remotes and remote URL.
  Is the GitHub repo **private or public**? It should stay private until Dr. Bo agrees to release.
- Find big blobs already in history:
  `git rev-list --objects --all | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' | awk '$1=="blob" && $3>20000000' | sort -k3 -n`
  Anything over 100 MB blocks a GitHub push. Report it. Rewriting history (git filter-repo) needs my OK.
- Secrets and personal-path scan of tracked and to-be-tracked files:
  - GitHub/GHCR tokens, `.env`, passwords
  - `D:\WORKSPACE`, `C:\Users\dtquocbao`
  - HPRC account numbers are fine in docs; tokens are not
- **Stop and report.**

## 2. Code and docs layout
- **Move to `archive/docs/`:** all `docs/CURSOR_TASK_*.md` (this one too, once done), `docs/R1_STATUS.md`, `docs/R2_STATUS.md`, `docs/ARCSV_CORE_AUDIT.md`, `docs/QUESTIONS_FOR_DR_BO_SHORT.md`.
- **Keep in `docs/`:** METHOD, PARITY, HISTORY, ORIGINAL_COMPARISON, OUTPUT, CLUSTER (+ `CLUSTER_TAMU.md` if created), CROSS_PLATFORM, LABELS_PLAN, ML_DATASET, ML_SAMPLES, TROUBLESHOOTING, USABILITY_TEST.
- **Update `docs/ORIGINAL_COMPARISON.md`** with the new benchmark (numbers below), and link the report doc: https://claude.ai/artifact/9tCLAfkKAsnH4ShswQpuZW
- **New analysis scripts** were added at `scripts/analysis/compare_orig_vs_new_chr20.py` and `scripts/analysis/eval_drbo_xgb.py`. They were written in a scratch environment.
  - Make them run from the repo: take input paths as arguments (`export/orig/labels`, `export/v2fix11/labels`) and remove the hard-coded `/mnt/user-data/...` and `/tmp/...` paths.
  - Add a short `scripts/analysis/README.md`.
  - They must reproduce these numbers:
    - before the filter: orig 243/449 (54.1%) vs new 222/376 (59.0%), p = 0.16; DUP p = 0.011; 221/243 orig true calls matched
    - after Dr. Bo's XGBoost: orig 193/215 (89.8%) vs new 165/201 (82.1%), p = 0.033
- `scratch/`: delete everything except `chrom_lens.py` and `measure_rss.py`. Move those two into `scripts/hpc/` or `scripts/bench/` if useful.
- README: quick starts (app, CLI, cluster), the benchmark table, and the citation block (Zhou et al., Cell 2024; HPRC acknowledgement for cluster work).

## 3. bench/ and logs
- Keep `bench/RESULTS.md` and small logs (< 1 MB).
- For logs > 5 MB: keep the first and last 500 lines as `*.head_tail.log`, and put the full gzipped log in the release asset (step 4). Don't commit the 64 MB originals.

## 4. Results: what goes where
**In git (small, the evidence behind the report)**
- `export/orig/` (calls, labels, time files)
- `export/v2fix11/` (the golden reference; `check_golden` needs it)
- `export/ml0/dataset_v0.summary.json`, `export/ml0/merged/` (if < 20 MB total; else release)
- `tests/golden/`

**Release asset, not in git:** a GitHub Release named `benchmark-2026-11`, or Zenodo, with a `SHA256SUMS` file and a `docs/DATA.md` that lists each asset, its size, its checksum and the command to download it.
- `export/ml0/dataset_v0.parquet` + `dataset_v0.tsv.gz`
- `export/ml0/dataset.tsv`, regenerated after the HG00140 fix, or dropped
- full gzipped logs
- `export_history.tar.gz` holding `export/v2fix … v2fix10`, `export/results`, `export/labels` (superseded history)

**Not anywhere public:** Dr. Bo's `model_replication` package (pickled models, HG002 demo).
- It is unpublished lab material. Keep it at a gitignored path, `third_party/drbo_model/`, and document in `docs/ML_DATASET.md` how to obtain it from Dr. Bo.
- `eval_drbo_xgb.py` takes its location as an argument.

Update `.gitignore` for: `third_party/drbo_model/`, `export/ml0/*.parquet`, `export/ml0/*.tsv*`, `bench/logs/**/*.log` > 5 MB (or explicit names), `dist/`, `*.tar`, `*.sif`.

## 5. Verify
- Run `pytest -q` (default suite) and `check_golden.py --fast` vs `export/v2fix11`: both must PASS.
- Run `scripts/analysis/compare_orig_vs_new_chr20.py`: it must reproduce the numbers in step 2.
- Fresh-clone test: clone to a temp dir, then `docker compose build` (or pull the GHCR image) → `arcsv2 --version` → `check_golden --fast` on the committed fixtures.
- Repo size after cleanup: `git count-objects -vH`. Target < 100 MB.

## 6. Commit plan (do not push yet)
Logical commits, each with the attribution lines from our convention:
1. `chore: archive task docs and status notes`
2. `chore: trim logs, ignore large artifacts`
3. `feat(analysis): orig-vs-new benchmark and XGBoost evaluation scripts`
4. `docs: benchmark results, data manifest (DATA.md), README quick starts`
5. `data: commit benchmark evidence (orig, v2fix11, ml0 summary)`

Tag `benchmark-2026-11` on the last commit.

**Stop and report:**
- the `git log --oneline` of the new commits
- repo size
- the list of release assets with checksums
- the output of the verify step

I'll approve the push and the release upload.
