# Analysis scripts

## `compare_orig_vs_new_chr20.py`

Head-to-head precision on chr20 for HG01884 + NA19700 using the same haplotype
labeler. Inputs are directories of `*.labeled.tsv` files.

```bash
python scripts/analysis/compare_orig_vs_new_chr20.py \
  --orig-labels export/orig/labels \
  --new-labels export/v2fix11/labels
```

Reproduces (pre-ML filter):

| Caller | Validated / labeled | Precision |
|:--|--:|--:|
| Original ARC-SV | 243 / 449 | 54.1% |
| ARC-SV2 (`v2fix11`) | 222 / 376 | 59.0% |

Overall Fisher *p* ≈ 0.16; duplication *p* ≈ 0.011; 221/243 original true
positives matched by ARC-SV2 (±100 bp, same sample/class).

## `eval_drbo_xgb.py`

Applies Dr. Bo's unpublished XGBoost filter to preprocessed matrices. The model
package is **not** in git — see `docs/ML_DATASET.md` and
`third_party/drbo_model/` (gitignored).

After XGB keep (chr20): orig 193/215 (89.8%) vs new 165/201 (82.1%), *p* ≈ 0.033.

Requires: `pandas`, `numpy`, `scipy`, `scikit-learn`, `joblib`.
