# ML-1 model comparison

Outer `StratifiedGroupKFold` by sample, stratified by superpopulation. Hyperparameters and thresholds are chosen with inner `GroupKFold` on the training samples only. Intervals are 95% bootstrap-over-samples percentiles.

constant 0.0 on every ML-0 row; not used as a feature.
NA on every ML-0 row; kept, so logistic/RF see a missing indicator and XGBoost sees NaN.
Included as extra XGBoost/RF candidates. max_depth 143 (complex) and 84 (deletion) are far above the grid {3, 6, 10} and are not defaults.

## Per class

### deletion

n = 14904 (12846 positive / 2058 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 14904 | 0.950 | [0.949, 0.952] | 0.835 | 0.954 | 0.957 | 0.957 |
| logistic f1 | 14904 | 0.977 | [0.976, 0.979] | 0.908 | 0.970 | 0.900 | 0.960 |
| logistic f1+context | 14904 | 0.976 | [0.976, 0.978] | 0.908 | 0.970 | 0.901 | 0.960 |
| logistic f1+context no-pca | 14904 | 0.977 | [0.976, 0.979] | 0.909 | 0.970 | 0.901 | 0.960 |
| xgboost f1 | 14904 | 0.987 | [0.986, 0.990] | 0.943 | 0.980 | 0.901 | 0.972 |
| xgboost f1+context | 14904 | 0.989 | [0.988, 0.990] | 0.949 | 0.983 | 0.901 | 0.976 |
| xgboost f1+context no-pca | 14904 | 0.990 | [0.988, 0.991] | 0.951 | 0.982 | 0.899 | 0.977 |

Auto choice: `xgboost f1+context` — highest recall_at_precision95 (0.9759).

`logistic f1+context no-pca` minus the same model with PCA: AP 0.000 (CI [0.000, 0.000]), precision @ recall≥0.90 -0.000 (CI [-0.001, 0.000]).
`xgboost f1+context no-pca` minus the same model with PCA: AP 0.001 (CI [-0.000, 0.002]), precision @ recall≥0.90 -0.000 (CI [-0.002, 0.000]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (54.1), split_support_mean (27), n_bp_no_split_support (20.3), split_support_sum (2.87), sv_bp_uncertainty_max (2.43), bp_uncertainty_sum (1.99), split_support_max (1.91), score_vs_next (1.78), n_blocks_ref (1.71), n_tokens_alt (1.47), is_het (1.4), ref_pc1 (1.31), len_affected_over_span (1.27), bp_uncertainty_mean (1.19), score_vs_ref (1.13)
`xgboost f1+context`: split_support_min (43.3), split_support_mean (31), n_bp_no_split_support (26.4), split_support_max (6.99), split_support_sum (2.51), sv_bp_uncertainty_max (2.01), ref_pc1 (1.58), n_tokens_alt (1.52), bp_uncertainty_sum (1.45), n_blocks_ref (1.43), score_vs_next (1.27), score_vs_ref_per_bp (1.01), len_affected_over_span (0.968), is_het (0.95), af (0.916)

### duplication

n = 2327 (1847 positive / 480 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 2327 | 0.920 | [0.905, 0.934] | 0.831 | 0.925 | 0.961 | not reachable |
| logistic f1 | 2327 | 0.957 | [0.944, 0.976] | 0.903 | 0.949 | 0.918 | 0.899 |
| logistic f1+context | 2327 | 0.958 | [0.944, 0.976] | 0.905 | 0.951 | 0.917 | 0.927 |
| logistic f1+context no-pca | 2327 | 0.958 | [0.945, 0.975] | 0.904 | 0.952 | 0.917 | 0.931 |
| xgboost f1 | 2327 | 0.978 | [0.972, 0.984] | 0.938 | 0.962 | 0.904 | 0.946 |
| xgboost f1+context | 2327 | 0.979 | [0.972, 0.986] | 0.943 | 0.968 | 0.906 | 0.952 |
| xgboost f1+context no-pca | 2327 | 0.975 | [0.970, 0.983] | 0.940 | 0.970 | 0.914 | 0.958 |

Auto choice: `xgboost f1` — xgboost f1 (0.9459) is inside the xgboost f1+context recall_at_precision95 CI [0.9404, 0.9618] and is simpler.

`logistic f1+context no-pca` minus the same model with PCA: AP -0.000 (CI [-0.003, 0.002]), precision @ recall≥0.90 0.001 (CI [-0.001, 0.003]).
`xgboost f1+context no-pca` minus the same model with PCA: AP -0.004 (CI [-0.007, 0.000]), precision @ recall≥0.90 0.002 (CI [-0.000, 0.004]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (21), n_bp_no_split_support (16.2), n_duplicated_blocks (8.81), split_support_mean (8.67), n_novel_adjacencies (8.31), rearr_pc2 (2.5), split_support_sum (2.27), n_blocks_ref (1.16), rearr_pc1 (0.971), af (0.946), split_support_max (0.877), log10_span (0.866), is_het (0.853), score_vs_ref (0.847), score_vs_ref_per_bp (0.794)
`xgboost f1+context`: n_bp_no_split_support (23.3), split_support_min (22.8), split_support_mean (15.5), n_duplicated_blocks (10.6), split_support_max (4.72), rearr_pc2 (3), split_support_sum (2.26), frac_simple_repeat (2.14), n_blocks_ref (1.93), ref_pc1 (1.29), is_het (1.26), score_vs_ref_per_bp (1.15), rearr_pc1 (0.934), min_bp_total_support (0.912), af (0.912)

### inversion

n = 30 (17 positive / 13 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 30 | 0.630 | [0.517, 0.731] | 0.615 | 0.630 | 1.000 | not reachable |
| logistic f1 | 30 | 0.982 | [0.962, 1.000] | 0.973 | 0.941 | 0.941 | 0.824 |
| logistic f1+context | 30 | 0.980 | [0.962, 1.000] | 0.968 | 0.941 | 0.941 | 0.765 |
| logistic f1+context no-pca | 30 | 0.982 | [0.962, 1.000] | 0.973 | 0.941 | 0.941 | 0.824 |
| xgboost f1 | 30 | 0.990 | [0.980, 1.000] | 0.986 | 0.938 | 0.882 | 0.882 |
| xgboost f1+context | 30 | 0.990 | [0.980, 1.000] | 0.986 | 0.938 | 0.882 | 0.882 |
| xgboost f1+context no-pca | 30 | 0.990 | [0.980, 1.000] | 0.986 | 0.938 | 0.882 | 0.882 |
| rf f1 | 30 | 0.993 | [0.981, 1.000] | 0.991 | 1.000 | 0.882 | 0.882 |
| rf f1+context | 30 | 0.997 | [0.986, 1.000] | 0.995 | 1.000 | 0.882 | 0.882 |

Auto choice: `logistic f1` — logistic f1 (0.8235) is inside the rf f1 recall_at_precision95 CI [0.7500, 1.0000] and is simpler.

`logistic f1+context no-pca` minus the same model with PCA: AP 0.003 (CI [0.000, 0.009]), precision @ recall≥0.90 0.000 (CI [0.000, 0.000]).
`xgboost f1+context no-pca` minus the same model with PCA: AP 0.000 (CI [-0.000, 0.000]), precision @ recall≥0.90 0.000 (CI [0.000, 0.000]).

Beats `sr>0` on AP outside the bootstrap CI: logistic f1, logistic f1+context, logistic f1+context no-pca, xgboost f1, xgboost f1+context, xgboost f1+context no-pca, rf f1, rf f1+context.
`sr>0` AP CI width is 0.214 ([0.517, 0.731]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_sum (4.21), split_support_min (2.73), split_support_mean (1.39), n_bp_no_split_support (0.588), split_support_max (0.44), min_bp_total_support (0.224), inner_block_len_mean (0.207), len_affected (0.151), score_vs_ref (0.136), score_vs_next (0.131), inner_block_len_min (0.117), flank_len_min (0.0816), log10_len_affected (0.0703), len_affected_over_span (0.0652), score_ratio_next_ref (0.063)
`xgboost f1+context`: split_support_mean (4.01), split_support_sum (3.75), split_support_min (3.37), n_bp_no_split_support (1.68), split_support_max (0.478), score_vs_ref (0.345), log10_len_affected (0.208), min_bp_total_support (0.206), inner_block_len_min (0.141), len_affected (0.137), score_vs_next (0.0788), pair_support_min (0.0628), gc_flank_500 (0.0447), score_vs_ref_per_bp (0.0385), flank_len_min (0.0351)

### complex

n = 1244 (575 positive / 669 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 1244 | 0.686 | [0.665, 0.710] | 0.793 | 0.647 | 0.917 | not reachable |
| logistic f1 | 1244 | 0.868 | [0.849, 0.886] | 0.899 | 0.726 | 0.894 | not reachable |
| logistic f1+context | 1244 | 0.867 | [0.847, 0.886] | 0.898 | 0.719 | 0.896 | not reachable |
| logistic f1+context no-pca | 1244 | 0.869 | [0.848, 0.887] | 0.900 | 0.724 | 0.901 | not reachable |
| xgboost f1 | 1244 | 0.911 | [0.900, 0.920] | 0.923 | 0.767 | 0.903 | 0.652 |
| xgboost f1+context | 1244 | 0.905 | [0.894, 0.916] | 0.918 | 0.781 | 0.896 | 0.678 |
| xgboost f1+context no-pca | 1244 | 0.914 | [0.899, 0.922] | 0.924 | 0.786 | 0.894 | 0.670 |

Auto choice: `xgboost f1` — xgboost f1 (0.6522) is inside the xgboost f1+context recall_at_precision95 CI [0.6334, 0.7208] and is simpler.

`logistic f1+context no-pca` minus the same model with PCA: AP 0.002 (CI [-0.001, 0.004]), precision @ recall≥0.90 0.005 (CI [-0.009, 0.022]).
`xgboost f1+context no-pca` minus the same model with PCA: AP 0.010 (CI [0.003, 0.013]), precision @ recall≥0.90 0.004 (CI [-0.014, 0.019]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (44.1), n_bp_no_split_support (20.4), split_support_mean (9.15), n_duplicated_blocks (7.47), split_support_max (6.82), split_support_sum (6.37), inner_block_len_mean (3.43), n_tokens_alt (3.4), sv_bp_uncertainty_max (3.19), n_inverted_blocks (3.14), n_blocks_ref (2.84), bp_uncertainty_sum (2.75), log10_span (2.73), bp_uncertainty_mean (2.66), min_bp_total_support (2.6)
`xgboost f1+context`: split_support_min (39.7), n_bp_no_split_support (38.9), split_support_mean (15.5), split_support_max (9.29), n_duplicated_blocks (7.75), split_support_sum (6.56), n_inverted_blocks (3.75), sv_bp_uncertainty_max (3.74), inner_block_len_mean (3.57), n_tokens_alt (3.45), ref_pc1 (3.14), bp_uncertainty_mean (3.14), min_bp_total_support (3.07), is_het (2.82), inner_block_len_max (2.75)

## Precision at the recall ≥ 0.90 threshold, by superpopulation

### deletion

| Model | AFR |
|:--|---:|
| sr>0 | 0.954 |
| logistic f1 | 0.970 |
| logistic f1+context | 0.970 |
| logistic f1+context no-pca | 0.970 |
| xgboost f1 | 0.980 |
| xgboost f1+context | 0.983 |
| xgboost f1+context no-pca | 0.982 |

### duplication

| Model | AFR |
|:--|---:|
| sr>0 | 0.925 |
| logistic f1 | 0.949 |
| logistic f1+context | 0.951 |
| logistic f1+context no-pca | 0.952 |
| xgboost f1 | 0.962 |
| xgboost f1+context | 0.968 |
| xgboost f1+context no-pca | 0.970 |

### inversion

| Model | AFR |
|:--|---:|
| sr>0 | 0.630 |
| logistic f1 | 0.941 |
| logistic f1+context | 0.941 |
| logistic f1+context no-pca | 0.941 |
| xgboost f1 | 0.938 |
| xgboost f1+context | 0.938 |
| xgboost f1+context no-pca | 0.938 |
| rf f1 | 1.000 |
| rf f1+context | 1.000 |

### complex

| Model | AFR |
|:--|---:|
| sr>0 | 0.647 |
| logistic f1 | 0.726 |
| logistic f1+context | 0.719 |
| logistic f1+context no-pca | 0.724 |
| xgboost f1 | 0.767 |
| xgboost f1+context | 0.781 |
| xgboost f1+context no-pca | 0.786 |

