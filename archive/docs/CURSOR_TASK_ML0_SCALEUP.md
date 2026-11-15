# Cursor task: ML-0 — scale up the labeled dataset (whole genome, many samples)

Caller: tag `caller-v1.6.1` (engine=rust, frozen image `arcsv2:ml0` with native baked in;
AF call-table format is six decimals). Day-to-day `compose.yaml` still defaults to
`arcsv2:dev`; ML-0 S3/S4 must add `-f compose.ml0.yaml` so every service uses
`arcsv2:ml0`.
Goal: produce the training data for the ML filter. This task does not design or train any model.
Stop and report after each step.

## Machine (ThinkPad P16, Windows + Docker Desktop/WSL2)
- i9-12950HX: 16 cores / 24 threads
- 128 GB RAM
- 7.28 TB disk, 5.7 TB free
- RTX A5500 16 GB (not needed for ML-0)

## S0. Host setup (report before running anything big)
1. **WSL2 limits.** WSL defaults to 50% of RAM. Create or update
   `%UserProfile%\.wslconfig` with:
   ```
   [wsl2]
   memory=112GB
   processors=24
   swap=16GB
   ```
   Then run `wsl --shutdown`, restart Docker Desktop, and confirm that
   `docker run --rm arcsv2:dev nproc` and `free -g` show 24 CPUs and about 110 GB.
2. **Disk location.** Report which Windows drive holds Docker's WSL disk image
   (`docker-desktop-data` / `ext4.vhdx`) and how much is free there. This step
   needs about **1.2 TB**:
   - CRAM-derived BAMs: about 18 GB × 40
   - assemblies: about 6 GB × 40
   - minimap2 indexes
   - results

   If that drive lacks space, move the disk image with Docker Desktop
   Settings → Resources → Disk image location. **Do not** bind-mount a Windows
   folder as `/data` (9p I/O is slow).
   **The `arcsv2-data` volume must survive any move.** Back up or verify it first,
   and ask me before moving it.
3. **Laptop sustained load.** Plugged in, best-performance power plan, sleep
   disabled while jobs run.

## S1. Image sanity check
From a **fresh** `docker compose run` container (not the long-lived `arcsv2-dev`):
- `python -c "import arcsv2_native"`
- golden `--fast` on a re-call of the golden windows

This also closes the open caveat from PERF3: the golden check after the image
rebuild was run on existing outputs, not on a fresh call.

## S2. Sample selection
- Use `fetch_data.py --auto` (HPRC Release 2 ∩ NYGC 1KGP high-coverage):
  `--dry-run` first, and list the candidates with population / superpopulation / sex.
- Pick **40 samples**:
  - balanced across superpopulations (AFR, AMR, EAS, EUR, SAS) as far as the
    intersection allows
  - both sexes
  - include HG01884 and NA19700
  - exclude HG002 (not in 1KGP)
- Write `workflows/config/samples.ml.tsv` with pinned URLs and checksums, plus
  `docs/ML_SAMPLES.md` (table + selection rule).
- **Ask me to confirm the list before downloading.**

## S3. Pilot: 3 samples, whole genome (chr1–22 + chrX)
- Pick 3 samples from different superpopulations.
- Run the full Snakemake flow per (sample, chrom): fetch → call → merge →
  merge-filter (repeats/segdups/blacklist/gaps BEDs) → `ml-label` → dataset.
- **Scheduler settings.**
  - Use `--cores 22 --resources mem_mb=100000`.
  - Run calls at `threads=4` each, so several chromosome jobs run at once. Most
    stages are single-threaded and genotype is the only one that scales, so many
    4-thread jobs beat a few 8-thread ones. Benchmark 4 vs 8 on one sample's chr1–3.
  - Keep `call_chrom` `mem_mb` from the length-based formula.
- Record in `bench/RESULTS.md`:
  - per sample: wall time, CPU-hours, and peak RSS per job type (call, label/minimap2, index_haps)
  - disk per sample
  - download time and bandwidth
- Report per sample:
  - calls per type
  - labeled / validated / NA counts (no_sequence, query_too_long, …)
  - precision per type, compared with the chr20 numbers
- Then **extrapolate the 40-sample time and disk, and stop for my OK**.

## S4. Full run: 40 samples
- Use the same settings as the S3 pilot. Download ahead of calling (a prefetch
  queue of about 2 samples), but cap disk use. Delete per-sample intermediates
  you don't need (the CRAM after the per-chrom BAMs are made, temp files) only
  after that sample's `dataset` step succeeds.
- Snakemake must be resumable: `--rerun-incomplete`, plus a log per sample.
- Report progress by sample count.
- **Golden first:** run `check_golden --fast` once at the start of S4.

## S5. Dataset assembly (no modeling)
Write `/data/ml/dataset_v0.parquet` (also a `.tsv.gz`) with one row per merged call.

**Identity**
- sample, superpopulation, population, sex, chrom, start, end, call id

**Label**
- `label` (1 / 0 / NA) and the NA reason
- labeler details (identity, max gap, MAPQ, hap hit)

**Call features**
Every column the caller already writes, kept raw. Do not engineer features yet.
- SV type / class, rearrangement string, GT, `score_vs_ref`, AF
- sr / pe supports per junction, number of blocks / junctions
- affected length, `length_is_lower_bound`
- `bp_uncertainty`, `bp_point_source` (clip/split/drp)

**Context**
- overlap fraction with simple repeats, segdups, blacklist, gaps
- reference GC in the SV ± 500 bp flanks
- local coverage ratio (in vs flanks) if cheap

**Also**
- `docs/ML_DATASET.md`: the data dictionary (column, source, type)
- summary tables of label balance per SV type × superpopulation × chromosome

## S6. Optional, low priority: original ARC-SV at scale
If CPU is idle at any point, run the original ARC-SV (`src/arcsv`) on **5** of the
40 samples, whole genome, using the same filter and labeler. This gives the
original-vs-ours comparison at scale instead of only chr20. It's single-threaded
and needs about 2 GB; run it at low priority.

## Not in scope (next tasks)
- Insertion labeling (flank-based; insertions are NA today) → ML-1
- Feature engineering (Table S1-like) and models → ML-2, after Dr. Bo replies
  about Table S1 and the train/test split

Keep `scratch/` clean, put logs in `bench/logs/ml0/`, and keep the golden rule
for any caller code change (there should be none in this task).
