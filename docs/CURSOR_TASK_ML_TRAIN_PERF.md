# CURSOR_TASK_ML_TRAIN_PERF: make `ml-train` fast without changing its output

## Why
On Grace, `arcsv2 ml-train --labels export/ml0/dataset.tsv` (130,009 labeled rows, 40 samples)
did not finish in 30 min (job 20653417, TIMEOUT). A rerun with 6 h on 1 core is job 20653724;
it finished in **38 min 02 s wall (2,237 s user), peak RSS 2.37 GB**, exit 0 (the 30 min limit was
just too short). Its `manifest.json` (209,422 bytes) is the reference output; a copy is at
`export/ml0/perf_ref/grace_20653724_manifest.json`.

The work is pure Python and repeated. For each class, `fit_labeled_calls` builds the design matrix
1 (final fit) + 4 (OOF folds) + ~24 (chromosome hold-out) times, plus a PCA fit each time. Every
build calls `featurize_row` on every row, and every PCA fit rebuilds `rearrangement_matrix` /
`reference_indicator` per row. Only the four PCA columns depend on the fold.

## Rules
- **The manifest must not change.** Same keys, same thresholds, same `grouping`, same
  `chromosome_holdout` entries. Floats may differ only in the last bits: abs diff ≤ 1e-9 and
  rel diff ≤ 1e-9 per number. Thresholds must be identical (they are chosen from sorted scores;
  check that no score crosses a threshold).
- No change to `featurize_row`, `feature_names`, `FEATURE_SCHEMA`, the PCA definition
  (first two right singular vectors of the centered matrix, same sign), scaling, or the
  `LogisticRegression(max_iter=500)` settings. `ml-predict` must keep reading old manifests.
- No new dependencies. numpy / scikit-learn only. No pickle.
- Keep it single-process. (Parallel folds can be a later task.)

## Changes
1. **Featurize once.** In `fit_labeled_calls`, per class, compute once:
   - the base feature matrix `X_base` (all `feature_names()` columns except the four PCA columns),
     by calling `featurize_row(row, None)` once per row;
   - `R` = stacked `rearrangement_matrix(...).ravel()` and `F` = stacked `reference_indicator(...)`.
   Then every fit (final, each OOF fold, each held-out chromosome) uses index arrays into these.
2. **PCA from cached arrays.** Add `fit_pca_arrays(R, F) -> PcaModel` in `features.py` that does
   exactly what `fit_pca` does on the given rows. Make `fit_pca(rows)` call it, so there is one
   implementation. Fix `fit_pca` computing each SVD twice (`_center_components(...)` is called
   separately for `[0]` and `[1]`).
3. **PCA columns from cached arrays.** For a fitted `pca`, the four columns are
   `(R[idx] - rearr_mean) @ rearr_components.T` and `(F[idx] - ref_mean) @ ref_components.T`.
   Put them in the same column positions `_design` uses today (`feature_names()` order).
   Check against `featurize_row(row, pca)` on a sample of rows (≤ 1e-9).
4. **Vectorize scoring.** Replace the per-row `_logistic_prob` loops in `_oof_probabilities` and
   `_chromosome_holdout` with a numpy version: `z = intercept + ((X - mean) / scale_safe) @ coef`,
   where `scale_safe` maps 0 to 1 exactly as `_logistic_prob` does, and the same ±60 clamp.
   Keep `_logistic_prob` itself for `ml-predict`, unchanged.
5. **Chromosome hold-out indices.** Build the chromosome → indices map once instead of scanning
   all rows per held-out chromosome.

## Tests
- New `tests/test_train_perf.py`:
  - `fit_pca_arrays` equals `fit_pca` on the same rows (≤ 1e-12).
  - Cached PCA columns equal `featurize_row(row, pca)` values.
  - Vectorized scores equal `_logistic_prob` per row (≤ 1e-12), including the clamp.
  - `fit_labeled_calls` on a small fixture (a few hundred rows, ≥ 3 samples, ≥ 3 chromosomes,
    all four classes) gives a manifest equal to the one from the current code. Generate that
    reference with the current code first and commit it as a fixture.
- Existing train/predict tests still pass.

## Acceptance (report all of these)
1. Laptop, subset: run old and new `ml-train` on the rows of 3 samples from
   `export/ml0/dataset.tsv`; a compare script (`scripts/compare_manifest.py old new`) prints
   PASS under the rule above, and the max abs/rel float diff.
2. Laptop, full 40 samples: wall time and peak RSS for the new code
   (`/usr/bin/time -v`, or the PowerShell equivalent).
3. Grace: run the new code in the SIF on the full table (job script:
   `$SCRATCH/arcsv2/jobs/ml_smoke.sh`, 1 core, 7.5 GB) and compare its manifest with the
   reference from job 20653724 using `scripts/compare_manifest.py`. Expected: PASS, and wall time
   in minutes rather than hours.
   The image has the old code; run the new code without rebuilding by putting the repo's
   `python/` first on the path: `singularity exec -B /scratch --env PYTHONPATH=$SCRATCH/arcsv2/repo/python ...`
   (check `python -c 'import arcsv2, sys; print(arcsv2.__file__)'` inside the SIF points at the repo).
4. Note in `docs/ML_DATASET.md` the old and new wall times.

Do not rebuild or push the image as part of this task; Bao decides when to cut a new image.
