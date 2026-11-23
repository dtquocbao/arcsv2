# ML-1 results

Laptop evaluation of `export/ml0/dataset_v0.tsv.gz` (40 samples, training rows with `label` in {0, 1} and `exclude` != 1). Outer `StratifiedGroupKFold` (5) by sample, stratified by superpopulation. Hyperparameters and both thresholds use inner `GroupKFold` (4) on the training samples only. Each outer test row is scored by the model refit on that fold and judged by that fold's threshold. Intervals are 95% bootstrap-over-samples percentiles (1,000 resamples, seed 204178949).

`frac_blacklist` is 0.0 on every ML-0 row, so it is not a feature. `coverage_ratio` is NA on every row; logistic and random forest see a missing indicator, and XGBoost sees NaN. It cannot separate calls.

The paper presets were extra candidates. `max_depth` 143 (complex) and 84 (deletion) never won an outer fold. Inversion's random-forest preset (`max_depth` 215) won 2 of 5 folds of `rf f1+context`; the other folds kept the unrestricted grid forest. Those depths are not defaults.

The `auto` choice uses the outer-fold numbers in the tables below. That is the selection set, not a second held-out test. Overlap with Dr. Bo's HPRC training samples is unknown, so these numbers are not a comparison with his published models.

Grace was not run. The job script is `scripts/hpc/ml1_eval.sbatch` and was not submitted. Raw output: `export/ml0/ml1_eval/full/`.

## Wall time

Laptop, one process, `n_jobs` 1, `tree_method=hist`.

| | |
|:--|--:|
| Wall time | 11,563 s (3 h 12 min 43 s) |
| Peak RSS | 2,668,859,392 bytes (2.49 GiB) |

Libraries: scikit-learn 1.9.0, xgboost 3.2.0, skops 0.14.0.

A 5-sample AFR smoke (HG01884, HG02055, HG02572, HG02922, HG03098; 18,505 training rows) finished first, in 2,185 s with peak RSS 554 MiB. The tables below are the 40-sample run.

## Per class

Recall @ precision ≥ 0.95 is "not reachable" when at least one outer fold found no threshold with precision ≥ 0.95 on its inner out-of-fold scores.

### deletion

n = 103,426 (89,009 positive / 14,417 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 103426 | 0.949 | [0.947, 0.951] | 0.833 | 0.953 | 0.955 | 0.955 |
| logistic f1 | 103426 | 0.978 | [0.977, 0.979] | 0.909 | 0.969 | 0.900 | 0.958 |
| logistic f1+context | 103426 | 0.978 | [0.976, 0.979] | 0.909 | 0.969 | 0.900 | 0.958 |
| xgboost f1 | 103426 | 0.995 | [0.994, 0.996] | 0.974 | 0.992 | 0.899 | 0.985 |
| xgboost f1+context | 103426 | 0.996 | [0.995, 0.997] | 0.980 | 0.994 | 0.901 | 0.991 |

Auto choice: `xgboost f1+context` — highest recall at precision ≥ 0.95 (0.991). Every outer fold chose grid `max_depth` 10.

### duplication

n = 16,717 (13,117 positive / 3,600 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 16717 | 0.925 | [0.921, 0.930] | 0.851 | 0.931 | 0.960 | not reachable |
| logistic f1 | 16717 | 0.971 | [0.967, 0.975] | 0.930 | 0.964 | 0.899 | 0.947 |
| logistic f1+context | 16717 | 0.973 | [0.969, 0.976] | 0.933 | 0.963 | 0.901 | 0.947 |
| xgboost f1 | 16717 | 0.990 | [0.987, 0.993] | 0.974 | 0.986 | 0.901 | 0.977 |
| xgboost f1+context | 16717 | 0.992 | [0.987, 0.994] | 0.978 | 0.990 | 0.903 | 0.982 |

Auto choice: `xgboost f1+context` — highest recall at precision ≥ 0.95 (0.982). Every outer fold chose grid `max_depth` 10.

### complex

n = 9,625 (4,276 positive / 5,349 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 9625 | 0.680 | [0.669, 0.692] | 0.800 | 0.444 | 1.000 | not reachable |
| logistic f1 | 9625 | 0.897 | [0.889, 0.906] | 0.912 | 0.723 | 0.901 | 0.537 |
| logistic f1+context | 9625 | 0.898 | [0.890, 0.907] | 0.913 | 0.732 | 0.900 | 0.550 |
| xgboost f1 | 9625 | 0.955 | [0.951, 0.960] | 0.955 | 0.854 | 0.905 | 0.831 |
| xgboost f1+context | 9625 | 0.959 | [0.955, 0.963] | 0.957 | 0.872 | 0.904 | 0.843 |

Auto choice: `xgboost f1` — its recall at precision ≥ 0.95 (0.831) sits inside the `xgboost f1+context` interval [0.830, 0.854], and f1 is the simpler of the two. Fold depths were 6 or 10 from the grid.

At the recall ≥ 0.90 threshold the best precision is 0.872. Precision ≥ 0.95 is reachable, at recall 0.843 (`xgboost f1+context`) or 0.831 (`xgboost f1`).

### inversion

n = 241 (145 positive / 96 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 241 | 0.660 | [0.611, 0.712] | 0.610 | 0.643 | 0.917 | not reachable |
| logistic f1 | 241 | 0.970 | [0.949, 0.989] | 0.953 | 0.930 | 0.917 | 0.841 |
| logistic f1+context | 241 | 0.985 | [0.974, 0.994] | 0.974 | 0.971 | 0.938 | 0.952 |
| xgboost f1 | 241 | 0.971 | [0.929, 0.999] | 0.975 | 0.986 | 0.938 | 0.972 |
| xgboost f1+context | 241 | 0.975 | [0.941, 0.998] | 0.972 | 0.979 | 0.972 | 0.972 |
| rf f1 | 241 | 0.993 | [0.985, 0.999] | 0.986 | 0.993 | 0.931 | 0.972 |
| rf f1+context | 241 | 0.993 | [0.985, 1.000] | 0.986 | 0.985 | 0.924 | 0.972 |

`sr>0` AP interval width is 0.101. The other AP interval widths are 0.041 (logistic f1), 0.020 (logistic f1+context), 0.069 (xgboost f1), 0.057 (xgboost f1+context), 0.014 (rf f1), and 0.014 (rf f1+context).

Every fitted model's AP interval lies entirely above the `sr>0` interval. On this split the models beat the rule on average precision. The rule does not reach precision 0.95. The models do, with recall from 0.841 (logistic f1) to 0.972 (xgboost and random forest).

Auto choice: `logistic f1+context`. Random forest f1 has the highest point estimate of recall at precision ≥ 0.95 (0.972), and its interval is [0.938, 1.000]. Logistic f1+context (0.952) sits inside that interval, so the simpler model is kept. With 241 rows that interval still runs from 0.94 to 1.00, so the ranking among the models that clear precision 0.95 is not resolved.

## PCA ablation

Same folds, four PCA columns removed (`rearr_pc1`, `rearr_pc2`, `ref_pc1`, `ref_pc2`). The difference is no-pca minus with-pca.

| Class | Model | AP difference | AP difference CI |
|:--|:-|--:|:--|
| deletion | logistic f1+context | 0.000 | [0.000, 0.000] |
| deletion | xgboost f1+context | 0.000 | [0.000, 0.000] |
| duplication | logistic f1+context | 0.000 | [−0.001, 0.000] |
| duplication | xgboost f1+context | 0.000 | [−0.001, 0.001] |
| complex | logistic f1+context | −0.001 | [−0.002, −0.000] |
| complex | xgboost f1+context | 0.000 | [−0.002, 0.001] |
| inversion | logistic f1+context | 0.002 | [−0.001, 0.006] |
| inversion | xgboost f1+context | 0.010 | [0.000, 0.028] |

For XGBoost the AP change is inside the bootstrap interval on every class, including a complex interval of about 0.003. Logistic complex loses about 0.001 AP, and that small drop sits just outside zero. Removing PCA does not change the auto choice. It would drop the SVD step for a change at the size of these intervals.

## Precision at the recall ≥ 0.90 threshold, by superpopulation

### deletion

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.955 | 0.955 | 0.948 | 0.953 | 0.953 |
| logistic f1 | 0.971 | 0.970 | 0.966 | 0.969 | 0.970 |
| logistic f1+context | 0.971 | 0.970 | 0.967 | 0.969 | 0.970 |
| xgboost f1 | 0.990 | 0.993 | 0.991 | 0.992 | 0.994 |
| xgboost f1+context | 0.992 | 0.996 | 0.992 | 0.995 | 0.996 |

### duplication

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.930 | 0.929 | 0.926 | 0.939 | 0.932 |
| logistic f1 | 0.960 | 0.966 | 0.959 | 0.971 | 0.963 |
| logistic f1+context | 0.961 | 0.963 | 0.959 | 0.971 | 0.961 |
| xgboost f1 | 0.981 | 0.988 | 0.985 | 0.990 | 0.988 |
| xgboost f1+context | 0.985 | 0.990 | 0.989 | 0.993 | 0.993 |

### inversion

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.636 | 0.605 | 0.735 | 0.585 | 0.667 |
| logistic f1 | 0.926 | 0.893 | 0.963 | 0.889 | 0.971 |
| logistic f1+context | 0.962 | 0.963 | 1.000 | 0.926 | 1.000 |
| xgboost f1 | 0.964 | 1.000 | 1.000 | 0.962 | 1.000 |
| xgboost f1+context | 0.933 | 1.000 | 1.000 | 0.962 | 1.000 |
| rf f1 | 1.000 | 1.000 | 1.000 | 0.962 | 1.000 |
| rf f1+context | 0.960 | 1.000 | 1.000 | 0.962 | 1.000 |

### complex

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.468 | 0.438 | 0.442 | 0.434 | 0.438 |
| logistic f1 | 0.724 | 0.711 | 0.718 | 0.716 | 0.744 |
| logistic f1+context | 0.730 | 0.723 | 0.727 | 0.727 | 0.755 |
| xgboost f1 | 0.853 | 0.842 | 0.863 | 0.851 | 0.861 |
| xgboost f1+context | 0.876 | 0.870 | 0.880 | 0.860 | 0.875 |

## Top XGBoost gains

Mean gain over the outer training fits.

Complex, `xgboost f1`: split_support_min (109), n_bp_no_split_support (37), n_duplicated_blocks (13), split_support_mean (8.25), num_sv (7.77), n_inverted_blocks (5.3), is_het (5.05), af (4.77), sv_bp_uncertainty_max (4.69), n_novel_adjacencies (4.49), split_support_max (4.12), n_bp_no_pair_support (3.9), n_tokens_alt (3.67), log10_span (3.35), bp_uncertainty_sum (3.11).

Complex, `xgboost f1+context`: split_support_min (76.9), n_bp_no_split_support (59.7), n_duplicated_blocks (11.6), split_support_mean (10.4), af (5.68), num_sv (5.54), split_support_max (4.8), n_novel_adjacencies (4.76), n_inverted_blocks (4.55), sv_bp_uncertainty_max (4.41), is_het (4.14), n_bp_no_pair_support (3.54), n_tokens_alt (2.88), log10_len_affected (2.71), span (2.61). No context column is in this list.

Inversion, `xgboost f1`: split_support_min (15.9), n_bp_no_split_support (8.9), len_affected (5.52), log10_span (5.44), log10_len_affected (4.34), len_affected_over_span (2.7), split_support_sum (2.49), is_het (1.91), score_vs_ref_per_bp (1.76), score_vs_ref (1.48), rearr_pc2 (1.05), bp_uncertainty_max (1.02), bp_uncertainty_mean (1.01), inner_block_len_mean (0.997), split_support_mean (0.975).

Inversion, `xgboost f1+context`: split_support_min (16), n_bp_no_split_support (7.08), len_affected (5.77), log10_len_affected (5.7), split_support_mean (4.37), len_affected_over_span (3.83), span (3.16), split_support_sum (2.48), pair_support_mean (1.73), score_vs_ref_per_bp (1.58), score_vs_ref (1.32), min_bp_total_support (1.07), bp_uncertainty_mean (1.03), log10_span (1.03), frac_simple_repeat (1.02).

## What would move complex and inversion

Complex is the class that still misses the precision target at high recall: the best precision at recall ≥ 0.90 is 0.872, and precision ≥ 0.95 only appears once recall falls to 0.843. The gain list is split-read support and a few rearrangement counts; the six context columns do not enter the top 15, `coverage_ratio` is missing on all 210,141 rows, and `frac_blacklist` is constant. Another tree depth will not close that gap — depths 6 and 10 won, and depth 143 never did. More labeled complex calls, and features that are not a restatement of `sr_support`, are what the numbers point to. Inversion already beats `sr>0` on average precision, with the rule's interval 0.101 wide and the tightest model interval 0.014 wide, but 241 rows leave the recall-at-precision-0.95 interval at [0.938, 1.000]. That is too wide to pick a model, which is why auto stays with logistic f1+context. More labeled inversions, across samples, are what would shrink it.
