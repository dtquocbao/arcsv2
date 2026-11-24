# ML-1 model comparison

Outer `GroupKFold` by chromosome, balanced by row count. Hyperparameters and thresholds use inner `GroupKFold` on the training chromosomes. Intervals are 95% bootstrap-over-samples percentiles.

constant 0.0 on every ML-0 row; not used as a feature.
NA on every ML-0 row; kept, so logistic/RF see a missing indicator and XGBoost sees NaN.
Included as extra XGBoost/RF candidates. max_depth 143 (complex) and 84 (deletion) are far above the grid {3, 6, 10} and are not defaults.

## Per class

### deletion

n = 103426 (89009 positive / 14417 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 103426 | 0.949 | [0.947, 0.951] | 0.833 | 0.953 | 0.955 | 0.955 |
| logistic f1 | 103426 | 0.976 | [0.975, 0.978] | 0.906 | 0.969 | 0.900 | 0.957 |
| logistic f1+context | 103426 | 0.976 | [0.975, 0.977] | 0.905 | 0.969 | 0.900 | 0.957 |
| xgboost f1 | 103426 | 0.981 | [0.980, 0.982] | 0.922 | 0.973 | 0.898 | 0.966 |
| xgboost f1+context | 103426 | 0.982 | [0.981, 0.983] | 0.925 | 0.974 | 0.898 | 0.969 |

Auto choice: `xgboost f1+context` — highest recall_at_precision95 (0.9687).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (753), split_support_mean (211), n_bp_no_split_support (158), n_tokens_alt (39.2), bp_uncertainty_sum (33.1), score_vs_next (31.9), split_support_sum (27.9), sv_bp_uncertainty_max (26.1), bp_uncertainty_max (20), n_blocks_ref (19.5), len_affected_over_span (17.8), af (16.3), bp_uncertainty_mean (15.6), ref_pc2 (15.4), score_vs_ref (14.5)
`xgboost f1+context`: split_support_min (598), split_support_mean (371), n_bp_no_split_support (158), split_support_max (28.5), n_tokens_alt (27.1), split_support_sum (24.9), score_vs_next (23.6), bp_uncertainty_sum (22.6), bp_uncertainty_max (22.3), sv_bp_uncertainty_max (18.8), ref_pc1 (16.5), score_vs_ref_per_bp (14.2), n_blocks_ref (14), ref_pc2 (12.3), len_affected_over_span (12.2)

### duplication

n = 16717 (13117 positive / 3600 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 16717 | 0.925 | [0.921, 0.930] | 0.851 | 0.931 | 0.960 | not reachable |
| logistic f1 | 16717 | 0.968 | [0.964, 0.972] | 0.923 | 0.960 | 0.903 | 0.945 |
| logistic f1+context | 16717 | 0.968 | [0.964, 0.972] | 0.923 | 0.961 | 0.903 | 0.945 |
| xgboost f1 | 16717 | 0.967 | [0.963, 0.970] | 0.924 | 0.961 | 0.899 | 0.941 |
| xgboost f1+context | 16717 | 0.972 | [0.968, 0.975] | 0.932 | 0.963 | 0.904 | 0.945 |

Auto choice: `logistic f1` — highest recall_at_precision95 (0.9450).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (85), n_duplicated_blocks (45.7), n_bp_no_split_support (43.5), split_support_mean (24.2), n_novel_adjacencies (17.8), rearr_pc2 (3.76), n_blocks_ref (2.77), af (2.08), next_is_ref_het (2.03), rearr_pc1 (1.96), len_affected (1.75), len_affected_over_span (1.74), log10_span (1.71), inner_block_len_min (1.69), span (1.66)
`xgboost f1+context`: n_bp_no_split_support (197), split_support_min (190), split_support_max (78.9), n_duplicated_blocks (52.5), split_support_mean (31.4), rearr_pc2 (26.2), n_novel_adjacencies (23.6), n_blocks_ref (8.43), split_support_sum (8.01), n_tokens_alt (6.36), frac_simple_repeat (6.26), bp_uncertainty_max (4.96), gc_flank_500 (4.82), ref_pc1 (4.62), span (4.47)

### inversion

n = 241 (145 positive / 96 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 241 | 0.660 | [0.611, 0.712] | 0.610 | 0.601 | 0.924 | not reachable |
| logistic f1 | 241 | 0.843 | [0.808, 0.878] | 0.770 | 0.686 | 0.917 | 0.131 |
| logistic f1+context | 241 | 0.856 | [0.823, 0.888] | 0.798 | 0.654 | 0.924 | 0.366 |
| xgboost f1 | 241 | 0.775 | [0.725, 0.825] | 0.735 | 0.684 | 0.924 | 0.379 |
| xgboost f1+context | 241 | 0.800 | [0.752, 0.846] | 0.749 | 0.700 | 0.869 | not reachable |

Auto choice: `logistic f1+context` — logistic f1+context (0.3655) is inside the xgboost f1 recall_at_precision95 CI [0.3333, 0.4286] and is simpler.

Beats `sr>0` on AP outside the bootstrap CI: logistic f1, logistic f1+context, xgboost f1, xgboost f1+context.
`sr>0` AP CI width is 0.101 ([0.611, 0.712]).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (22.5), n_bp_no_split_support (8.82), len_affected (5.26), span (4.62), log10_len_affected (3.73), len_affected_over_span (3.57), split_support_sum (2.98), is_het (2.15), sv_bp_uncertainty_max (2.14), score_vs_ref (2.03), score_vs_ref_per_bp (1.78), inner_block_len_mean (1.68), log10_span (1.56), split_support_mean (1.54), bp_uncertainty_max (1.42)
`xgboost f1+context`: split_support_min (21), n_bp_no_split_support (10.1), split_support_mean (6.99), len_affected_over_span (5.31), len_affected (4.57), log10_len_affected (4.36), span (2.5), bp_uncertainty_sum (2.22), is_het (2.21), split_support_sum (2.18), gc_flank_500 (2.01), score_vs_ref (1.78), split_support_max (1.71), score_vs_ref_per_bp (1.61), score_vs_ref_per_kb (1.59)

### complex

n = 9625 (4276 positive / 5349 negative).

| Model | n | AP | AP CI | ROC AUC | Precision @ recall≥0.90 | Recall @ recall≥0.90 | Recall @ precision≥0.95 |
|:--|---:|---:|:--|---:|---:|---:|---:|
| sr>0 | 9625 | 0.680 | [0.669, 0.692] | 0.800 | 0.515 | 0.945 | not reachable |
| logistic f1 | 9625 | 0.873 | [0.863, 0.884] | 0.903 | 0.704 | 0.906 | not reachable |
| logistic f1+context | 9625 | 0.873 | [0.865, 0.884] | 0.902 | 0.705 | 0.901 | not reachable |
| xgboost f1 | 9625 | 0.914 | [0.907, 0.920] | 0.921 | 0.738 | 0.900 | 0.564 |
| xgboost f1+context | 9625 | 0.914 | [0.907, 0.921] | 0.922 | 0.739 | 0.903 | 0.598 |

Auto choice: `xgboost f1+context` — highest recall_at_precision95 (0.5975).

XGBoost gain (mean over outer training fits), top 15:

`xgboost f1`: split_support_min (284), n_bp_no_split_support (71.5), split_support_mean (34.2), n_duplicated_blocks (24.4), split_support_max (18.1), num_sv (17.4), sv_bp_uncertainty_max (14.3), bp_uncertainty_sum (10.6), n_inverted_blocks (9.71), n_novel_adjacencies (9.41), n_tokens_alt (9.41), bp_uncertainty_mean (8.62), score_vs_ref_per_kb (8.1), log10_span (8.03), bp_uncertainty_max (7.88)
`xgboost f1+context`: split_support_min (197), n_bp_no_split_support (129), split_support_mean (32.3), n_duplicated_blocks (24.5), split_support_max (19.1), sv_bp_uncertainty_max (11.4), num_sv (10), n_novel_adjacencies (9.94), n_inverted_blocks (9.36), af (8.45), is_het (7.69), n_tokens_alt (6.75), n_bp_no_pair_support (6.54), bp_uncertainty_sum (6.23), log10_len_affected (6.14)

## Precision at the recall ≥ 0.90 threshold, by superpopulation

### deletion

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.955 | 0.955 | 0.948 | 0.953 | 0.953 |
| logistic f1 | 0.970 | 0.969 | 0.966 | 0.968 | 0.969 |
| logistic f1+context | 0.970 | 0.969 | 0.966 | 0.968 | 0.969 |
| xgboost f1 | 0.975 | 0.974 | 0.970 | 0.973 | 0.974 |
| xgboost f1+context | 0.976 | 0.976 | 0.971 | 0.973 | 0.974 |

### duplication

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.930 | 0.929 | 0.926 | 0.939 | 0.932 |
| logistic f1 | 0.958 | 0.961 | 0.954 | 0.969 | 0.958 |
| logistic f1+context | 0.960 | 0.961 | 0.957 | 0.969 | 0.959 |
| xgboost f1 | 0.962 | 0.961 | 0.955 | 0.968 | 0.960 |
| xgboost f1+context | 0.964 | 0.962 | 0.958 | 0.969 | 0.961 |

### inversion

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.609 | 0.529 | 0.676 | 0.571 | 0.638 |
| logistic f1 | 0.659 | 0.659 | 0.781 | 0.632 | 0.714 |
| logistic f1+context | 0.659 | 0.600 | 0.735 | 0.571 | 0.721 |
| xgboost f1 | 0.659 | 0.659 | 0.750 | 0.600 | 0.769 |
| xgboost f1+context | 0.730 | 0.667 | 0.786 | 0.579 | 0.771 |

### complex

| Model | AFR | AMR | EAS | EUR | SAS |
|:--|---:|---:|---:|---:|---:|
| sr>0 | 0.539 | 0.507 | 0.514 | 0.497 | 0.517 |
| logistic f1 | 0.710 | 0.693 | 0.691 | 0.699 | 0.723 |
| logistic f1+context | 0.704 | 0.694 | 0.703 | 0.701 | 0.721 |
| xgboost f1 | 0.753 | 0.723 | 0.734 | 0.720 | 0.758 |
| xgboost f1+context | 0.754 | 0.722 | 0.744 | 0.714 | 0.761 |

