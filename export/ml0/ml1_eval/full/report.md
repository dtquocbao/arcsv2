# ML-1 model comparison

Outer `StratifiedGroupKFold` by sample, stratified by superpopulation. Hyperparameters and thresholds are chosen with inner `GroupKFold` on the training samples only. Intervals are 95% bootstrap-over-samples percentiles.

constant 0.0 on every ML-0 row; not used as a feature.
NA on every ML-0 row; kept, so logistic/RF see a missing indicator and XGBoost sees NaN.
Included as extra XGBoost/RF candidates. max_depth 143 (complex) and 84 (deletion) are far above the grid {3, 6, 10} and are not defaults.

## Per class

### deletion

n = 103426 (89009 positive / 14417 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 103426 | 0.949 | [0.947, 0.951] | 0.833 | 0.953 | 0.955 | 0.955 |
| logistic f1 | 103426 | 0.978 | [0.977, 0.979] | 0.909 | 0.969 | 0.900 | 0.958 |
| logistic f1+context | 103426 | 0.978 | [0.976, 0.979] | 0.909 | 0.969 | 0.900 | 0.958 |
| logistic f1+context no-pca | 103426 | 0.978 | [0.977, 0.979] | 0.909 | 0.969 | 0.900 | 0.958 |
| xgboost f1 | 103426 | 0.995 | [0.994, 0.996] | 0.974 | 0.992 | 0.899 | 0.985 |
| xgboost f1+context | 103426 | 0.996 | [0.995, 0.997] | 0.980 | 0.994 | 0.901 | 0.991 |
| xgboost f1+context no-pca | 103426 | 0.996 | [0.995, 0.997] | 0.980 | 0.994 | 0.900 | 0.991 |

Auto choice: `xgboost f1+context` — highest recall_at_precision95 (0.9906).

`logistic f1+context no-pca` minus the same model with PCA: AP 0.000 (CI [0.000, 0.000]), precision @ recall≥0.90 -0.000 (CI [-0.000, 0.000]).
`xgboost f1+context no-pca` minus the same model with PCA: AP -0.000 (CI [-0.000, 0.000]), precision @ recall≥0.90 0.000 (CI [-0.000, 0.000]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (123), split_support_mean (57.2), n_bp_no_split_support (41.1), sv_bp_uncertainty_max (3.8), bp_uncertainty_sum (3.6), n_blocks_ref (2.69), len_affected_over_span (2.33), n_tokens_alt (2.28), score_vs_next (2.13), ref_pc1 (2.11), af (2.02), log10_span (1.99), is_het (1.98), split_support_max (1.94), span (1.91)
`xgboost f1+context`: split_support_min (123), split_support_mean (99.5), n_bp_no_split_support (48.8), split_support_max (4.31), sv_bp_uncertainty_max (4.21), bp_uncertainty_sum (3.56), split_support_sum (3.34), n_tokens_alt (3.25), ref_pc1 (3.01), frac_segdup (2.73), n_blocks_ref (2.51), len_affected_over_span (2.32), bp_uncertainty_max (2.28), n_bp_no_pair_support (2.27), n_novel_adjacencies (2.22)

### duplication

n = 16717 (13117 positive / 3600 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 16717 | 0.925 | [0.921, 0.930] | 0.851 | 0.931 | 0.960 | not reachable |
| logistic f1 | 16717 | 0.971 | [0.967, 0.975] | 0.930 | 0.964 | 0.899 | 0.947 |
| logistic f1+context | 16717 | 0.973 | [0.969, 0.976] | 0.933 | 0.963 | 0.901 | 0.947 |
| logistic f1+context no-pca | 16717 | 0.972 | [0.968, 0.976] | 0.932 | 0.963 | 0.900 | 0.947 |
| xgboost f1 | 16717 | 0.990 | [0.987, 0.993] | 0.974 | 0.986 | 0.901 | 0.977 |
| xgboost f1+context | 16717 | 0.992 | [0.987, 0.994] | 0.978 | 0.990 | 0.903 | 0.982 |
| xgboost f1+context no-pca | 16717 | 0.992 | [0.988, 0.994] | 0.978 | 0.990 | 0.905 | 0.982 |

Auto choice: `xgboost f1+context` — highest recall_at_precision95 (0.9817).

`logistic f1+context no-pca` minus the same model with PCA: AP -0.000 (CI [-0.001, 0.000]), precision @ recall≥0.90 -0.000 (CI [-0.001, 0.001]).
`xgboost f1+context no-pca` minus the same model with PCA: AP 0.000 (CI [-0.001, 0.001]), precision @ recall≥0.90 -0.000 (CI [-0.001, 0.001]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (55.1), n_duplicated_blocks (44.3), n_bp_no_split_support (37.9), n_novel_adjacencies (14.2), split_support_mean (8.35), n_blocks_ref (2.76), rearr_pc2 (1.89), log10_span (1.65), split_support_sum (1.47), next_is_ref_het (1.38), is_het (1.34), len_affected_over_span (1.21), len_affected (1.14), rearr_pc1 (1.14), span (1.12)
`xgboost f1+context`: n_bp_no_split_support (100), split_support_min (52.9), n_duplicated_blocks (44.7), split_support_mean (14), split_support_max (11.6), n_novel_adjacencies (6.53), n_blocks_ref (3.06), rearr_pc2 (2.8), frac_simple_repeat (1.77), log10_len_affected (1.63), split_support_sum (1.54), span (1.49), n_tokens_alt (1.41), rearr_pc1 (1.31), bp_uncertainty_max (1.27)

### inversion

n = 241 (145 positive / 96 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 241 | 0.660 | [0.611, 0.712] | 0.610 | 0.643 | 0.917 | not reachable |
| logistic f1 | 241 | 0.970 | [0.949, 0.989] | 0.953 | 0.930 | 0.917 | 0.841 |
| logistic f1+context | 241 | 0.985 | [0.974, 0.994] | 0.974 | 0.971 | 0.938 | 0.952 |
| logistic f1+context no-pca | 241 | 0.987 | [0.977, 0.996] | 0.978 | 0.979 | 0.945 | 0.952 |
| xgboost f1 | 241 | 0.971 | [0.929, 0.999] | 0.975 | 0.986 | 0.938 | 0.972 |
| xgboost f1+context | 241 | 0.975 | [0.941, 0.998] | 0.972 | 0.979 | 0.972 | 0.972 |
| xgboost f1+context no-pca | 241 | 0.985 | [0.966, 0.998] | 0.975 | 0.993 | 0.966 | 0.972 |
| rf f1 | 241 | 0.993 | [0.985, 0.999] | 0.986 | 0.993 | 0.931 | 0.972 |
| rf f1+context | 241 | 0.993 | [0.985, 1.000] | 0.986 | 0.985 | 0.924 | 0.972 |

Auto choice: `logistic f1+context` — logistic f1+context (0.9517) is inside the rf f1 recall_at_precision95 CI [0.9375, 1.0000] and is simpler.

`logistic f1+context no-pca` minus the same model with PCA: AP 0.002 (CI [-0.001, 0.006]), precision @ recall≥0.90 0.007 (CI [0.000, 0.021]).
`xgboost f1+context no-pca` minus the same model with PCA: AP 0.010 (CI [-0.000, 0.028]), precision @ recall≥0.90 0.014 (CI [-0.000, 0.034]).

Beats `sr>0` on AP outside the bootstrap CI: logistic f1, logistic f1+context, logistic f1+context no-pca, xgboost f1, xgboost f1+context, xgboost f1+context no-pca, rf f1, rf f1+context.
`sr>0` AP CI width is 0.101 ([0.611, 0.712]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (15.9), n_bp_no_split_support (8.9), len_affected (5.52), log10_span (5.44), log10_len_affected (4.34), len_affected_over_span (2.7), split_support_sum (2.49), is_het (1.91), score_vs_ref_per_bp (1.76), score_vs_ref (1.48), rearr_pc2 (1.05), bp_uncertainty_max (1.02), bp_uncertainty_mean (1.01), inner_block_len_mean (0.997), split_support_mean (0.975)
`xgboost f1+context`: split_support_min (16), n_bp_no_split_support (7.08), len_affected (5.77), log10_len_affected (5.7), split_support_mean (4.37), len_affected_over_span (3.83), span (3.16), split_support_sum (2.48), pair_support_mean (1.73), score_vs_ref_per_bp (1.58), score_vs_ref (1.32), min_bp_total_support (1.07), bp_uncertainty_mean (1.03), log10_span (1.03), frac_simple_repeat (1.02)

### complex

n = 9625 (4276 positive / 5349 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 9625 | 0.680 | [0.669, 0.692] | 0.800 | 0.444 | 1.000 | not reachable |
| logistic f1 | 9625 | 0.897 | [0.889, 0.906] | 0.912 | 0.723 | 0.901 | 0.537 |
| logistic f1+context | 9625 | 0.898 | [0.890, 0.907] | 0.913 | 0.732 | 0.900 | 0.550 |
| logistic f1+context no-pca | 9625 | 0.897 | [0.889, 0.906] | 0.912 | 0.738 | 0.899 | 0.545 |
| xgboost f1 | 9625 | 0.955 | [0.951, 0.960] | 0.955 | 0.854 | 0.905 | 0.831 |
| xgboost f1+context | 9625 | 0.959 | [0.955, 0.963] | 0.957 | 0.872 | 0.904 | 0.843 |
| xgboost f1+context no-pca | 9625 | 0.959 | [0.954, 0.963] | 0.957 | 0.869 | 0.902 | 0.840 |

Auto choice: `xgboost f1` — xgboost f1 (0.8314) is inside the xgboost f1+context recall_at_precision95 CI [0.8299, 0.8544] and is simpler.

`logistic f1+context no-pca` minus the same model with PCA: AP -0.001 (CI [-0.002, -0.000]), precision @ recall≥0.90 0.006 (CI [0.002, 0.010]).
`xgboost f1+context no-pca` minus the same model with PCA: AP -0.000 (CI [-0.002, 0.001]), precision @ recall≥0.90 -0.003 (CI [-0.008, 0.002]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (109), n_bp_no_split_support (37), n_duplicated_blocks (13), split_support_mean (8.25), num_sv (7.77), n_inverted_blocks (5.3), is_het (5.05), af (4.77), sv_bp_uncertainty_max (4.69), n_novel_adjacencies (4.49), split_support_max (4.12), n_bp_no_pair_support (3.9), n_tokens_alt (3.67), log10_span (3.35), bp_uncertainty_sum (3.11)
`xgboost f1+context`: split_support_min (76.9), n_bp_no_split_support (59.7), n_duplicated_blocks (11.6), split_support_mean (10.4), af (5.68), num_sv (5.54), split_support_max (4.8), n_novel_adjacencies (4.76), n_inverted_blocks (4.55), sv_bp_uncertainty_max (4.41), is_het (4.14), n_bp_no_pair_support (3.54), n_tokens_alt (2.88), log10_len_affected (2.71), span (2.61)

## Precision at the recall ≥ 0.90 threshold, by superpopulation

### deletion

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.955 | 0.955 | 0.948 | 0.953 | 0.953 |
| logistic f1 | 0.971 | 0.970 | 0.966 | 0.969 | 0.970 |
| logistic f1+context | 0.971 | 0.970 | 0.967 | 0.969 | 0.970 |
| logistic f1+context no-pca | 0.971 | 0.970 | 0.967 | 0.969 | 0.970 |
| xgboost f1 | 0.990 | 0.993 | 0.991 | 0.992 | 0.994 |
| xgboost f1+context | 0.992 | 0.996 | 0.992 | 0.995 | 0.996 |
| xgboost f1+context no-pca | 0.992 | 0.996 | 0.993 | 0.995 | 0.995 |

### duplication

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.930 | 0.929 | 0.926 | 0.939 | 0.932 |
| logistic f1 | 0.960 | 0.966 | 0.959 | 0.971 | 0.963 |
| logistic f1+context | 0.961 | 0.963 | 0.959 | 0.971 | 0.961 |
| logistic f1+context no-pca | 0.961 | 0.963 | 0.960 | 0.970 | 0.961 |
| xgboost f1 | 0.981 | 0.988 | 0.985 | 0.990 | 0.988 |
| xgboost f1+context | 0.985 | 0.990 | 0.989 | 0.993 | 0.993 |
| xgboost f1+context no-pca | 0.984 | 0.991 | 0.988 | 0.995 | 0.992 |

### inversion

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.636 | 0.605 | 0.735 | 0.585 | 0.667 |
| logistic f1 | 0.926 | 0.893 | 0.963 | 0.889 | 0.971 |
| logistic f1+context | 0.962 | 0.963 | 1.000 | 0.926 | 1.000 |
| logistic f1+context no-pca | 1.000 | 0.964 | 1.000 | 0.926 | 1.000 |
| xgboost f1 | 0.964 | 1.000 | 1.000 | 0.962 | 1.000 |
| xgboost f1+context | 0.933 | 1.000 | 1.000 | 0.962 | 1.000 |
| xgboost f1+context no-pca | 0.964 | 1.000 | 1.000 | 1.000 | 1.000 |
| rf f1 | 1.000 | 1.000 | 1.000 | 0.962 | 1.000 |
| rf f1+context | 0.960 | 1.000 | 1.000 | 0.962 | 1.000 |

### complex

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.468 | 0.438 | 0.442 | 0.434 | 0.438 |
| logistic f1 | 0.724 | 0.711 | 0.718 | 0.716 | 0.744 |
| logistic f1+context | 0.730 | 0.723 | 0.727 | 0.727 | 0.755 |
| logistic f1+context no-pca | 0.737 | 0.724 | 0.733 | 0.732 | 0.767 |
| xgboost f1 | 0.853 | 0.842 | 0.863 | 0.851 | 0.861 |
| xgboost f1+context | 0.876 | 0.870 | 0.880 | 0.860 | 0.875 |
| xgboost f1+context no-pca | 0.863 | 0.860 | 0.884 | 0.862 | 0.878 |

