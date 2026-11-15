#!/usr/bin/env python3
"""Evaluate Dr. Bo's XGBoost filter on preprocessed feature matrices.

This script does **not** ship Dr. Bo's unpublished ``model_replication`` package.
Place that material under a gitignored path (e.g. ``third_party/drbo_model/``) and
pass paths explicitly.

Expected headline (chr20 HG01884+NA19700 after XGB keep):
  orig 193/215 (89.8%) vs new 165/201 (82.1%), Fisher p ≈ 0.033

Example (paths depend on how the model package lays out CSV/PKL files):

  python scripts/analysis/eval_drbo_xgb.py \\
    --tag orig \\
    --meta third_party/drbo_model/metaorig.csv \\
    --features third_party/drbo_model/inorig_preprocess.csv \\
    --raw-input third_party/drbo_model/inorig.csv \\
    --predictions third_party/drbo_model/outorig.csv \\
    --model third_party/drbo_model/model.pkl

  python scripts/analysis/eval_drbo_xgb.py --tag new ...
"""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def min_sr(value: object) -> float:
    try:
        return min(float(v) for v in str(value).split(",") if v not in ("", "NA"))
    except (TypeError, ValueError):
        return float("nan")


def report(name: str, y: np.ndarray, mask: np.ndarray) -> None:
    keep = int(mask.sum())
    prec = float(y[mask].mean()) if keep else float("nan")
    rec = float(y[mask].sum() / y.sum()) if y.sum() else float("nan")
    print(f"  {name:22s} keep={keep:6d} ({keep / len(y):.2f})  precision={prec:.3f}  true kept={rec:.2f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True, help="Short label printed in the header (e.g. orig/new)")
    ap.add_argument("--meta", type=Path, required=True, help="CSV with at least a label column")
    ap.add_argument("--features", type=Path, required=True, help="Preprocessed feature matrix CSV")
    ap.add_argument("--raw-input", type=Path, required=True, help="Raw call table CSV (for sr_support col)")
    ap.add_argument("--predictions", type=Path, required=True, help="Model output TSV/CSV (last col keep flag)")
    ap.add_argument("--model", type=Path, required=True, help="joblib/pickle model with predict_proba")
    ap.add_argument("--sr-col", type=int, default=18, help="0-based column index of sr_support in raw input")
    args = ap.parse_args()

    for path in (args.meta, args.features, args.raw-input, args.predictions, args.model):
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing {path}. Obtain Dr. Bo's model_replication package and place it "
                "under third_party/drbo_model/ (gitignored); see docs/ML_DATASET.md."
            )

    meta = pd.read_csv(args.meta)
    features = pd.read_csv(args.features)
    model = joblib.load(args.model)
    probs = model.predict_proba(features)[:, 1]
    pred = (
        pd.read_csv(args.predictions, sep="\t", header=None)
        .iloc[:, -1]
        .astype(str)
        .isin(["True", "1"])
        .to_numpy()
    )
    y = meta["label"].astype(float).to_numpy()
    raw = pd.read_csv(args.raw_input, sep="\t", header=None)
    sr = raw.iloc[:, args.sr_col].map(min_sr).to_numpy()

    print(
        f"== {args.tag}: n={len(y)} AUC={roc_auc_score(y, probs):.3f} "
        f"AP={average_precision_score(y, probs):.3f}"
    )
    report("no filter", y, np.ones(len(y), dtype=bool))
    report("sr>0 rule", y, sr > 0)
    report("Dr.Bo XGB predict", y, pred)
    for thr in (0.7, 0.9):
        report(f"XGB prob>={thr}", y, probs >= thr)
    report("XGB & sr>0", y, pred & (sr > 0))

    if "superpopulation" in meta.columns:
        for sp, _group in meta.groupby("superpopulation"):
            mk = (meta["superpopulation"] == sp).to_numpy()
            kept = mk & pred
            prec = float(y[kept].mean()) if kept.any() else float("nan")
            denom = float(y[mk].sum()) or float("nan")
            frac = float(y[kept].sum() / denom) if denom == denom else float("nan")
            print(f"     {sp}: XGB precision {prec:.3f}  kept {frac:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
