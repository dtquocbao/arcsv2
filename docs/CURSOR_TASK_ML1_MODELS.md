# CURSOR_TASK_ML1_MODELS: stronger per-class models (XGBoost + context features)

## Why
The logistic baseline on the 40-genome ML-0 set (Grace, 2026-11-22; leave-one-chromosome-out,
thresholds at out-of-fold recall >= 0.90) is already strong for simple SVs and weak for the rest:

| Class | Labeled rows (1 / 0) | Mean AP | Precision at threshold |
|:--|--:|--:|--:|
| deletion | 103,426 (89,009 / 14,417) | 0.975 | 0.967 |
| duplication | 16,717 (13,117 / 3,600) | 0.970 | 0.962 |
| complex | 9,625 (4,276 / 5,349) | 0.852 | 0.696 |
| inversion | 241 (145 / 96) | 0.626 | 0.362 |

Dr. Bo's published XGBoost models scored on the same data (docs and project note of 2026-11-05):
complex 84.0% precision / 79% of true calls kept at 0.5; 92.7% / 65% at prob >= 0.9; the simple
`sr > 0` rule gave complex 82.1% / 81% and inversion 94.7% / 74%.

Goal: better complex and inversion filtering without losing deletion/duplication, measured honestly.
Target from the project: >= 95% precision per class, at the highest recall that allows.

The trainer is now fast (featurize-once, commit 0d7d3bf: 2.5 min on Grace for all classes), so a
proper cross-validated comparison is affordable.

## Data
- Input: `export/ml0/dataset_v0.parquet` (or `.tsv.gz`, same rows): 210,141 rows, 40 samples.
  Training rows: `label` in {0, 1} and `exclude` != 1. Insertions stay out (no sequence to label).
- Class: `class_of_row` from `train.py` (same as today).
- Grouping columns: `sample`, `superpopulation`. All 40 samples are from different families
  (`family_id` has no family with more than one sample), so grouping by sample is enough.
- Context columns already in `dataset_v0` and not used today: `frac_simple_repeat`, `frac_segdup`,
  `frac_gap`, `gc_flank_500`, `coverage_ratio`, `length_is_lower_bound`. Check `frac_blacklist`:
  the blacklist BED was not fetched (HTTP 403, `docs/PARITY.md`), so if it is constant, leave it out
  and say so.

## Rules
- **Do not change the current default.** `arcsv2 ml-train --labels X` with no new flags must produce
  a manifest identical to today's: `scripts/compare_manifest.py` against
  `$SCRATCH/arcsv2/ml/models/ml0_perf/manifest.json` (or a fresh laptop run of 0d7d3bf) must PASS.
- New features get a new schema name (`arcsv2-f2`). `ml-predict` must refuse a manifest whose schema
  does not match the input's, with a clear error, and must keep scoring old `arcsv2-f1` manifests.
- No pickle. XGBoost models via `save_model` JSON (embedded in the manifest as a string, or a side
  file listed in the manifest with its sha256). RandomForest, if used, via skops. Record library
  versions in the manifest.
- Deterministic: fixed seeds, `tree_method="hist"`, fixed `n_jobs`; two runs on the same machine give
  identical manifests (check with `compare_manifest.py`).
- No leakage: every fit (PCA, scaling, imputation, model, hyperparameters, threshold) happens inside
  the training folds only. Labels' helper columns (`reason`, `identity`, `mapq`, `max_indel`,
  `label_bp_shift`, `best_hap`, `gt_concordant`, `ref_validated`, `strand`) are labeler outputs and
  must never be features. Assert this in code.

## Changes
1. **Features.** Add a `context` feature group (the columns above) behind
   `--features f1` (default, today's) | `f1+context`. Missing values: XGBoost gets NaN natively;
   logistic/RF get median imputation fit on the training fold plus a `_missing` indicator.
2. **Models.** `--model logistic` (default, today's) | `xgboost` | `rf` | `auto` (per class, the best
   of the candidates by the selection rule below). XGBoost: `objective="binary:logistic"`,
   `scale_pos_weight` = neg/pos of the training fold, a small grid chosen by inner CV:
   `max_depth` {3, 6, 10}, `n_estimators` {200, 500}, `learning_rate` {0.05, 0.1},
   `min_child_weight` {1, 5}, `subsample` 0.8, `colsample_bytree` 0.8. Include the paper presets
   (`PAPER_PRESETS`) as extra candidates but note that `max_depth` 143 / 84 look like copy errors
   for XGBoost; report how they do, do not default to them.
3. **Evaluation protocol** (new `scripts/ml1_eval.py`, also callable as `arcsv2 ml-eval`):
   - Outer: `StratifiedGroupKFold(n_splits=5)` grouped by `sample`, stratified by
     `superpopulation`, fixed seed. Every number in the report comes from outer test folds.
   - Inner (inside each outer training set): `GroupKFold(4)` by sample for hyperparameters and for
     the threshold.
   - Per class report, pooled over outer test folds: AP; ROC AUC; precision and recall at the
     inner-chosen threshold for recall >= 0.90; recall at the inner-chosen threshold for
     precision >= 0.95 (or "not reachable"); n.
   - 95% confidence intervals by bootstrap over samples (1,000 resamples, fixed seed) for AP and the
     two operating points. Required for inversion, which has 241 rows.
   - Per-superpopulation precision at the recall-0.90 point (as the 2026-11-05 eval did).
   - Baselines in the same table, same folds: the `sr > 0` rule (support from `sr_support`), current
     logistic f1, logistic f1+context, xgboost f1, xgboost f1+context, rf (inversion only, if it helps).
   - Feature importance (XGBoost gain, top 15) for the chosen complex and inversion models.
   - **PCA ablation:** for logistic and xgboost (f1+context), also score the same folds with the
     four PCA features (`rearr_pc1/2`, `ref_pc1/2`) removed. Report the AP and operating-point
     differences with their bootstrap CIs, per class. If removing PCA costs nothing within the CI,
     say so: it would simplify the feature schema and drop the numerically delicate SVD step.
     (In the 2026-11-22 logistic baseline the PCA coefficients rank 13th to 46th of 50 by size.)
4. **Selection rule** for `--model auto`, per class: highest outer-fold recall at precision >= 0.95;
   if no candidate reaches 0.95, highest AP. Ties within the bootstrap CI go to the simpler model
   (logistic < rf < xgboost). Report the decision and why.
5. **Inversion honesty.** With 241 rows, state the CI widths. If no model beats `sr > 0` outside its
   CI, say so plainly and recommend the rule plus more labeled inversions, not a fragile model.
6. **Grace job script** `scripts/hpc/ml1_eval.sbatch`: 4 cores, 15 GB, `medium`, 6 h limit,
   `/usr/bin/sbatch`, runs in the SIF with `PYTHONPATH=$SCRATCH/arcsv2/repo/python` (the 1.6.1 image
   has the old trainer), writes to `$SCRATCH/arcsv2/ml/eval/ml1_<timestamp>/`, wraps the run in
   `/usr/bin/time -v`.

## Tests
- Default path unchanged (fixture manifest from `tests/fixtures/train_perf` still matches).
- Leakage guard: training with a labeler column in the feature list raises.
- Schema guard: `ml-predict` on an f2 manifest with f1 input raises; old f1 manifests still score.
- XGBoost manifest round-trip: train on the small fixture, save, load, predict -> same probabilities.
- Determinism: two `ml1_eval` runs on the fixture give identical reports.

## Acceptance (report all of these)
1. `compare_manifest.py` PASS for the default path vs 0d7d3bf.
2. `docs/ML1_RESULTS.md`: the per-class table with CIs and all baselines, the per-superpopulation
   table, top features for complex and inversion, the `auto` decision per class, and the wall time
   and peak RSS of the Grace run.
3. Run `ml1_eval` on Grace via the sbatch script and paste the job id, `sacct` line and report path.
4. One paragraph on what would most improve complex and inversion next (more samples, features,
   labeling), backed by the numbers.

## Out of scope
- Do not rebuild or push the image. Do not change the caller or labeler. Do not commit until Bao
  reviews.
- Overlap with Dr. Bo's HPRC training samples is still unknown (his sample list was not in the
  package). Do not present the comparison with his published models as an independent test.
