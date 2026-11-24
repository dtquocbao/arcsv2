# CURSOR_TASK_ML2_SCALE: more labeled genomes on Grace, with the dataset gaps fixed

## Why
ML-1 (docs/ML1_RESULTS.md) showed the models are limited by data, not by model choice:
- inversion has 241 labeled calls in 40 genomes; complex 9,625;
- on unseen chromosomes (locus holdout) complex XGBoost drops from AP 0.959 to 0.914 and inversion
  from ~0.98 to 0.86: the models need more distinct loci and more samples;
- `coverage_ratio` is NA on every ML-0 row (BAMs were deleted before the dataset build) and
  `frac_blacklist` is 0 everywhere (the 10x SV blacklist URL returns HTTP 403).

Grace is the platform (docs/CLUSTER_TAMU.md, "Grace gotchas"). Measured cost: ~46 SU per genome
for calling with two retries; ~39 expected with `call_min_per_mb: 2.0` and `mem_mb_floor: 28000`;
labeling 6.5 SU; haplotype indexing not yet measured. Allocation: 20,000 SU on Grace, <150 used.

HPRC Release 2 (`assemblies_release2_v1.0.index.csv`, pinned commit in `scripts/fetch_data.py`)
has 232 samples with both haplotypes; ML-0 used 40, so 192 are unused. How many have a 1KGP
high-coverage CRAM is decided by `fetch_data.py`'s intersection on Grace, where the catalogs are
already cached under `$SCRATCH/arcsv2/data/`.

## Rules
- Nothing here may change ML-0 results or the default `ml-train`. The 40-sample `dataset_v0` stays
  as is; the new dataset is `dataset_v1`.
- The workflow must still run on one machine (laptop, Docker) with the same Snakefile and configs;
  Grace specifics stay in the profile, job scripts and `run_label.sh`. Single-machine whole-genome
  speed is a separate future track; do not trade it away here.
- Every Grace job: `-N 1 -n 1`, explicit `--mem`, submitted with `/usr/bin/sbatch` (job files, not
  `--wrap` by hand). Repo scripts that run on Grace run inside the SIF (login-node python3 is 3.6).
- No new dependencies beyond what the image already has (samtools, pysam if present — check;
  otherwise samtools only).

## Changes

### 1. Sample selection (`scripts/select_ml2_samples.py`)
- Candidates = HPRC R2 samples with both haplotypes AND a 1KGP high-coverage CRAM (reuse
  `fetch_data.parse_hprc_csv`, `parse_1kgp_index`, `intersect_candidates`), minus the 40 ML-0
  samples (`workflows/config/samples.ml.all40.tsv`).
- Drop relatives: if a candidate shares `family_id` with an ML-0 sample or another chosen sample,
  keep at most one per family. Report how many were dropped.
- Split the remaining candidates, with a fixed seed, balanced by superpopulation:
  - **`test` (frozen):** 20 samples, never used for training, model selection or thresholds.
    Written to `workflows/config/samples.ml2.test.tsv` and listed in the dataset `split` column.
  - **`train_add`:** the rest, written in batches of 10 to `samples.ml2.batchNN.tsv`.
- Print the counts per superpopulation for each split and the total. Do not exceed 120 new samples
  in this task (20 test + up to 100 train_add); stop early if fewer qualify.
- Output a label config per batch (`workflows/config/label.ml2.batchNN.yaml`) with absolute Grace
  paths, `call_min_per_mb: 2.0`, `mem_mb_floor: 28000`, `threads: 4`, `label_threads: 4`,
  `label_mem_mb: 30000`, chroms chr1–22 + chrX, CRAMs streamed (no prefetch).

### 2. Coverage feature while the BAM exists (new rule `coverage`, command `arcsv2 ml-coverage`)
- Per sample and chromosome, after `call` and before the BAM is deleted: for each call in
  `<chrom>/arcsv_out.tab`, compute
  - `cov_inside` = mean depth over [minbp, maxbp] (for spans > 50 kb, sample 20 evenly spaced
    1 kb windows),
  - `cov_flank` = mean depth over 1 kb immediately left of minbp and 1 kb right of maxbp,
  - `cov_chrom_median` = median depth of 200 random 1 kb windows on that chromosome (fixed seed,
    excluding gaps.bed),
  - `coverage_ratio` = cov_inside / cov_flank (NA if cov_flank == 0), plus
    `coverage_ratio_chrom` = cov_inside / cov_chrom_median.
  - Use `samtools bedcov` (or pysam) with MAPQ >= 20; same thresholds everywhere.
- Output `<sample>/<chrom>/coverage.tsv` keyed by the call id; `merge` (or the dataset build) joins
  it. The merged call table gains the columns; `merged.tab` format for existing consumers is
  unchanged (write the joined table as a separate file if needed).
- **Backfill for the 40 ML-0 samples:** the same command must accept a CRAM URL + reference instead
  of a BAM and only read the needed windows (region queries through the Grace WebProxy). Add
  `scripts/hpc/backfill_coverage.sbatch` that runs it for all 40 samples from
  `export/ml0/merged/*.merged.tab` (1 core per sample is fine; a job array of 40). Report wall time
  and SU.

### 3. Storage: delete what is no longer needed
- Mark per-chromosome BAMs/BAIs `temp()` so Snakemake removes them once `call` and `coverage` are
  done. Per-sample `.mmi` files: remove after `label` (temp() or an explicit cleanup rule).
- Process batches of 10 samples; document the measured peak scratch use per batch. Hard stop:
  if `$SCRATCH` free space < 150 GB, do not start the next batch.

### 4. Blacklist replacement
- Add the ENCODE blacklist v2 for hg38 (Boyle-Lab/Blacklist on GitHub, pinned commit + sha256) to
  `scripts/fetch_filter_beds.py` as `blacklist_encode_v2`. Keep the 10x URL as the first choice;
  fall back to ENCODE only when the 10x download fails, and record which one was used in the
  dataset summary and docs/PARITY.md. It is a substitute (genomic-assay blacklist, not the 10x SV
  blacklist); say so.
- Recompute `frac_blacklist` for ML-0 calls too (needs only coordinates, not BAMs).

### 5. Controller as a job (multi-day run)
- `scripts/hpc/run_label.sbatch`: 1 core, 4 GB, `long` partition, up to 7 days, runs
  `run_label.sh <config> --rerun-incomplete` so the controller survives disconnects (no tmux).
  Snakemake submits from the compute node; verify `/usr/bin/sbatch` works from inside a job.

### 6. Dataset v1 (`scripts/build_ml2_dataset_v1.py`)
- Rows from the 40 ML-0 samples (with backfilled coverage and recomputed blacklist) plus every new
  sample; same columns as `dataset_v0` plus `coverage_ratio_chrom` and `split`
  (`ml0` | `train_add` | `test`). Write `export/ml2/dataset_v1.parquet` + `.tsv.gz` + summary JSON
  (label counts per class and split, superpopulation counts, which blacklist was used).
- `ml-eval` must be able to exclude `split == test` rows (`--exclude-split test`), so model
  selection never sees the frozen test set.

### 7. Measurements to report (per batch, from `su_report.sh` and `sacct`)
- SU per genome for calling (fetch + call + merge) with the new sizing: does it land near 39?
- Any call retries (TIMEOUT / OUT_OF_MEMORY) and on which chromosomes.
- index_haps SU per sample, label SU per sample, coverage SU per sample.
- Wall time per batch, peak scratch.

## Order of work
1. Selection script; run it on Grace (in the SIF) and paste the split counts. Do not fetch anything
   yet.
2. Coverage rule + command, blacklist fallback, temp() cleanup, controller job script. Tests:
   coverage on the label-smoke fixture BAM (known depths), temp() markers present, selection
   script determinism and family exclusion on a synthetic catalog, default `ml-train` unchanged
   (`scripts/compare_manifest.py` PASS vs the 0d7d3bf reference).
3. **Pilot:** 3 samples from batch 01, end to end on Grace (fetch, call, coverage, merge,
   index_haps, label). Report the measurements above. Stop for review.
4. After Bao approves the pilot: backfill coverage for the 40 ML-0 samples, then batches.
5. Dataset v1 + a rerun of `ml-eval` (sample holdout and `--outer-group chrom`), excluding the test
   split, with the four ML-1 models. Report in `docs/ML2_RESULTS.md` next to the ML-1 numbers.

## Budget guard
Stop and report if cumulative SU for this task passes 8,000 before the dataset build.

## Out of scope
No image rebuild/push, no caller or labeler algorithm changes, no commit without review. Insertion
labeling stays out (no inserted sequence).
