#!/usr/bin/env python3
"""Compare original ARC-SV vs ARC-SV2 chr20 labels (HG01884 + NA19700).

Expected headline numbers (pre-ML filter, labeled rows only):
  orig 243/449 (54.1%) vs new 222/376 (59.0%), overall Fisher p ≈ 0.16
  duplication Fisher p ≈ 0.011
  221/243 orig true calls matched by new (±100 bp, same sample/class)

Example:
  python scripts/analysis/compare_orig_vs_new_chr20.py \\
    --orig-labels export/orig/labels \\
    --new-labels export/v2fix11/labels
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

SAMPLES = ("HG01884", "NA19700")


def load_labels(directory: Path) -> pd.DataFrame:
    frames = []
    for sample in SAMPLES:
        path = directory / f"{sample}.labeled.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing labeled table: {path}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        if "sample" not in frame.columns:
            frame["sample"] = sample
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def add_label_numeric(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["y"] = pd.to_numeric(out["label"], errors="coerce")
    return out


def match_calls(query: pd.DataFrame, haystack: pd.DataFrame, tol: int = 100) -> np.ndarray:
    hits = []
    for _, row in query.iterrows():
        same = haystack[
            (haystack["sample"] == row["sample"])
            & (haystack["sv_class"] == row["sv_class"])
            & ((haystack["minbp"] - row["minbp"]).abs() <= tol)
            & ((haystack["maxbp"] - row["maxbp"]).abs() <= tol)
        ]
        hits.append(len(same) > 0)
    return np.asarray(hits, dtype=bool)


def fisher_prec(a: pd.Series, b: pd.Series) -> float:
    """Two-sided Fisher exact on precision contingency (b vs a)."""
    a = a.dropna()
    b = b.dropna()
    table = [
        [int(b.sum()), int(len(b) - b.sum())],
        [int(a.sum()), int(len(a) - a.sum())],
    ]
    return float(fisher_exact(table)[1])


def summarize(name: str, frame: pd.DataFrame) -> None:
    labeled = frame["y"].notna().sum()
    validated = int(frame["y"].sum())
    prec = float(frame["y"].mean()) if labeled else float("nan")
    print(f"{name}: rows={len(frame)} labeled={labeled} validated={validated} prec={prec:.3f}")
    by = (
        frame.groupby("sv_class")
        .agg(calls=("y", "size"), labeled=("y", "count"), validated=("y", "sum"), prec=("y", "mean"))
        .round(3)
    )
    print(by.to_string())
    na = frame[frame["y"].isna()]
    if len(na) and "reason" in na.columns:
        print(" NA reasons:", na["reason"].value_counts().to_dict())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--orig-labels",
        type=Path,
        default=Path("export/orig/labels"),
        help="Directory with HG01884.labeled.tsv and NA19700.labeled.tsv (original)",
    )
    ap.add_argument(
        "--new-labels",
        type=Path,
        default=Path("export/v2fix11/labels"),
        help="Directory with HG01884.labeled.tsv and NA19700.labeled.tsv (ARC-SV2)",
    )
    ap.add_argument("--tol", type=int, default=100, help="Match tolerance in bp")
    args = ap.parse_args()

    orig = add_label_numeric(load_labels(args.orig_labels))
    new = add_label_numeric(load_labels(args.new_labels))

    summarize("orig", orig)
    summarize("new", new)

    for sv_class in ("deletion", "duplication", "complex", "inversion"):
        a = orig.loc[orig["sv_class"] == sv_class, "y"]
        b = new.loc[new["sv_class"] == sv_class, "y"]
        if a.notna().any() and b.notna().any():
            p = fisher_prec(a, b)
            print(
                f"{sv_class}: new {int(b.sum())}/{int(b.notna().sum())} "
                f"orig {int(a.sum())}/{int(a.notna().sum())} p={p:.3g}"
            )

    overall_p = fisher_prec(orig["y"], new["y"])
    a = orig["y"].dropna()
    b = new["y"].dropna()
    print(
        f"overall: orig {int(a.sum())}/{len(a)} ({100 * a.mean():.1f}%) "
        f"vs new {int(b.sum())}/{len(b)} ({100 * b.mean():.1f}%) p={overall_p:.3g}"
    )

    ot = orig[orig["y"] == 1]
    nt = new[new["y"] == 1]
    mo = match_calls(ot, new, tol=args.tol)
    mn = match_calls(nt, orig, tol=args.tol)
    print(f"orig TPs found by new (±{args.tol}bp same class): {int(mo.sum())}/{len(ot)}")
    print(f"new TPs found by orig: {int(mn.sum())}/{len(nt)}")
    if (~mo).any():
        print(
            "orig-only TPs by class:",
            ot.loc[~mo, "sv_class"].value_counts().to_dict(),
            "median size",
            float(ot.loc[~mo, "len_affected"].median()),
        )
    if (~mn).any():
        print("new-only TPs by class:", nt.loc[~mn, "sv_class"].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
