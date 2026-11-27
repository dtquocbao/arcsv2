# Cursor task: ML2 step 6: dataset_v1, test-split exclusion, coverage retry

Follow-up to `docs/CURSOR_TASK_ML2_SCALE.md`. Commit 030066d covered steps 1–5. The labeling run on
Grace is done: 120 new genomes (100 train_add + 20 test) plus the 3-genome pilot, with no failures.
This task builds what's needed to train and evaluate on them.

**Do not commit.** Leave all changes uncommitted for Bao to review.

## What exists on Grace (all under `$SCRATCH/arcsv2`)

| Path | Contents |
|---|---|
| `data/labels_ml2_batch01` … `batch10` | 10 samples each: `<sample>.labeled.tsv`, `dataset.tsv`, `dataset.summary.json` |
| `data/labels_ml2_test` | 20 samples (frozen test split) |
| `data/labels_ml2_pilot` | 3 samples (HG03195, HG01192, NA18960) |
| `data/results_ml2_<run>/<sample>/merged.tab` and `chr*/coverage.tsv` | per-chromosome coverage: `id chrom cov_inside cov_flank cov_chrom_median coverage_ratio coverage_ratio_chrom` |
| `ml/ml0/merged/` | ML-0 merged tables (40 samples) |
| `ml/coverage_backfill/<sample>.coverage.tsv` | ML-0 backfilled coverage, **only 6 of 40 so far** (see part C) |
| repo `workflows/config/samples.ml2.{batch01..10,test}.tsv` | sample lists (first line is a `#` comment, then a header) |

Locate the existing ML-0 dataset (`dataset_v0`) and the code that built it, and reuse it. Don't
write a second feature pipeline.

## A. `scripts/build_ml2_dataset_v1.py`

- Rows: all ML-0 rows + every row from the 12 ML2 label folders. Same columns as `dataset_v0`, plus
  `coverage_ratio`, `coverage_ratio_chrom`, `split` (`ml0` | `train_add` | `test`) and
  `coverage_source` (`computed` | `backfill` | `missing`).
- Join coverage on `(sample, id)`. First check whether the ML2 `dataset.tsv` already carries the
  coverage columns. If not, join from `results_ml2_*/<sample>/chr*/coverage.tsv`.
- Split is decided by **sample-list membership**, not by folder: in `samples.ml2.test.tsv` →
  `test`; ML-0 → `ml0`; anything else → `train_add`. Pilot samples may also appear in a batch or test
  list: deduplicate by sample, and fail loudly if the same sample has conflicting labels.
- Missing coverage stays NA. Never impute it in the builder.
- Blacklist: ML2 used `blacklist: ""` (same as ML-0, where `frac_blacklist` was 0 everywhere). Record
  that in the summary; don't recompute it.
- Outputs: `export/ml2/dataset_v1.tsv.gz`, plus `.parquet` if pyarrow is in the SIF (skip quietly
  if it isn't), plus `dataset_v1.summary.json` with:
  - label counts (0/1/NA) per class × split;
  - superpopulation counts per split;
  - `coverage_source` counts per split;
  - blacklist used;
  - input paths.
- Must run inside the SIF (`singularity exec -B /scratch $SIF env PYTHONPATH=python python ...`):
  the login node's python is 3.6. Add `--data-root` (default `$SCRATCH/arcsv2`) and `--out`.

## B. `ml-eval --exclude-split test`

- New option `--exclude-split <name>` (repeatable). Rows with that `split` value are dropped before
  any fold is built. If the dataset has no `split` column and the option is given, fail with a clear
  error.
- Log, and write to the eval output, the number of samples and rows excluded.
- Works with both `--outer-group sample` and `--outer-group chrom`, and with the four ML-1 models
  (logistic, logistic+ctx, xgb, xgb+ctx). XGBoost handles NA coverage natively. For logistic, use
  the existing NA handling, and say in the docstring what it does.
- `scripts/hpc/ml2_eval.sbatch`: medium partition, 1 core, enough memory for the dataset, which is
  about 4× ML-0. Runs both holdouts with `--exclude-split test` into `export/ml2/eval/{sample,chrom}`.

## C. Coverage retry (`python/arcsv2/coverage.py`, `scripts/hpc/backfill_coverage.*`)

The ML-0 backfill streams CRAMs from S3 through Grace's proxy. 34 of 40 tasks failed partway with
errors like `Container header CRC32 failure` and `retry_reconnect ... samtools bedcov: error reading
from input file`.

- In `bedcov()`, split regions into chunks (default 500), run `samtools bedcov` per chunk, and retry a
  failed chunk up to 5 times with backoff (5 s, 15 s, 45 s, …). Raise only after the last attempt.
  Results must be identical to the single-call version (add a test).
- Cache finished chunks under the output dir so a rerun of the same sample resumes, then clean up
  once `coverage.tsv` is written.
- `backfill_coverage.sbatch`: throttle the array to `%4`, and skip a sample whose output already
  exists. Document how to resubmit only failed indices, e.g. `--array=1-4,6-9,...%4`.

## D. Small fixes

- `su_report.sh`: the "SU per genome" line currently sums call + fetch over the whole run. Divide by
  the number of samples in the log, and label it clearly.
- Its SU must use Grace billing: max(cores, mem_GB / 7.5) × hours. The 120-genome run used about 8,000
  SU per `myproject`, but sacct's `billing=` field suggested about 2,330.

## Tests

- Builder: tiny fixture with 2 ML-0 + 2 train_add + 1 test sample. Check split assignment, the
  coverage join, `coverage_source`, deduplication and the conflict error.
- `ml-eval --exclude-split test` removes exactly those samples from every fold.
- `bedcov` chunk + retry: a fake samtools that fails twice then succeeds gives the same totals.
- `pytest -q` passes.

## Report back

Changed files, test output, and the exact commands Bao runs on Grace for:
1. building dataset_v1;
2. the eval sbatch;
3. resubmitting the 34 failed backfills.
