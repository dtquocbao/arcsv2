"""Outer-fold comparison for ML-1. Also available as ``arcsv2 ml-eval``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from arcsv2.ml1 import write_evaluation  # noqa: E402


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Compare ML-1 models on outer sample folds.")
    parser.add_argument("--labels", type=Path, required=True, help="dataset_v0 TSV or TSV.GZ")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-samples", type=int, default=None, help="One superpopulation, this many samples")
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--outer-group", choices=("sample", "chrom"), default="sample")
    parser.add_argument("--include-rf", choices=("inversion", "all", "none"), default="inversion")
    parser.add_argument("--no-ablation", action="store_true", help="Skip the PCA-dropped models")
    args = parser.parse_args(argv)
    write_evaluation(
        args.labels,
        args.out,
        max_samples=args.max_samples,
        n_boot=args.n_boot,
        outer_group=args.outer_group,
        include_rf=args.include_rf,
        ablation=not args.no_ablation,
    )


if __name__ == "__main__":
    main()
