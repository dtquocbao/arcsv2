#!/usr/bin/env python3
"""Assemble ML-0 S5 dataset_v0 from labeled.tsv + sample metadata + context.

Keeps NA labels (unlike ``arcsv2 ml-dataset``). Writes parquet + tsv.gz under
``/data/ml/`` and a JSON summary for docs.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

from arcsv2.beds import read_bed_by_chrom
from arcsv2.merge_filter import _frac, inslen_is_lower_bound

FLANK = 500

IDENTITY_OUT = [
    "sample",
    "superpopulation",
    "population",
    "sex",
    "family_id",
    "chrom",
    "start",
    "end",
    "call_id",
]

LABEL_OUT = [
    "label",
    "reason",
    "gt_concordant",
    "ref_validated",
    "identity",
    "mapq",
    "max_indel",
    "strand",
    "label_bp_shift",
    "best_hap",
]

CONTEXT_OUT = [
    "frac_simple_repeat",
    "frac_segdup",
    "frac_blacklist",
    "frac_gap",
    "gc_flank_500",
    "coverage_ratio",
    "length_is_lower_bound",
]


def load_meta(path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        lines = [ln for ln in handle if ln.strip() and not ln.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        sample = (row.get("sample") or "").strip()
        if sample:
            meta[sample] = row
    return meta


def gc_fraction(fasta: pysam.FastaFile, chrom: str, start: int, end: int, flank: int) -> float:
    lo = max(0, start - flank)
    hi = end + flank
    try:
        seq = fasta.fetch(chrom, lo, hi).upper()
    except (ValueError, KeyError):
        return float("nan")
    n_gc = n = 0
    for base in seq:
        if base in "ACGT":
            n += 1
            if base in "GC":
                n_gc += 1
    return (n_gc / n) if n else float("nan")


def summarize(df: pd.DataFrame) -> dict:
    reasons = Counter(df["reason"].fillna("").astype(str))
    by_class: dict[str, dict[str, int]] = {}
    prec_class: dict[str, dict] = {}
    for sv, sub in df.groupby("sv_class", dropna=False):
        counts = Counter(sub["label"].astype(str))
        by_class[str(sv)] = dict(counts)
        pos, neg = counts.get("1", 0), counts.get("0", 0)
        den = pos + neg
        prec_class[str(sv)] = {
            "precision": (pos / den) if den else None,
            "n1": pos,
            "n0": neg,
            "nNA": counts.get("NA", 0),
        }

    def slice_prec(col: str) -> dict[str, dict]:
        out: dict[str, dict] = {}
        for key, sub in df.groupby(col, dropna=False):
            counts = Counter(sub["label"].astype(str))
            pos, neg = counts.get("1", 0), counts.get("0", 0)
            den = pos + neg
            out[str(key)] = {
                "precision": (pos / den) if den else None,
                "n1": pos,
                "n0": neg,
                "nNA": counts.get("NA", 0),
            }
        return out

    label_tot = Counter(df["label"].astype(str))
    return {
        "n_rows": int(len(df)),
        "n_samples": int(df["sample"].nunique()),
        "label_totals": dict(label_tot),
        "reasons": dict(reasons),
        "by_sv_class": by_class,
        "precision_by_sv_class": prec_class,
        "precision_by_superpopulation": slice_prec("superpopulation"),
        "precision_by_chrom": slice_prec("chrom"),
        "feature_schema": "arcsv2-f1",
        "label_version": "arcsv2-label-1",
        "dataset_version": "ml0-dataset-v0",
        "notes": {
            "na_kept": True,
            "blacklist_bed": "omitted (HTTP 403 on 10x svblacklist); frac_blacklist=0",
            "coverage_ratio": "NA (BAMs deleted after S4 labeling)",
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels-dir", type=Path, default=Path("/data/labels_ml0_s3"))
    ap.add_argument("--samples", type=Path, default=Path("/work/workflows/config/samples.ml.tsv"))
    ap.add_argument("--reference", type=Path, default=Path("/data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"))
    ap.add_argument("--gaps", type=Path, default=Path("/data/ref/gaps.bed"))
    ap.add_argument("--segdups", type=Path, default=Path("/data/ref/segdups.bed"))
    ap.add_argument("--simple-repeats", type=Path, default=Path("/data/ref/simple_repeats.bed"))
    ap.add_argument("--blacklist", type=Path, default=Path("/data/ref/blacklist.bed"))
    ap.add_argument("--out-dir", type=Path, default=Path("/data/ml"))
    args = ap.parse_args()

    meta = load_meta(args.samples)
    gaps = read_bed_by_chrom(args.gaps) if args.gaps.is_file() else {}
    seg = read_bed_by_chrom(args.segdups) if args.segdups.is_file() else {}
    simp = read_bed_by_chrom(args.simple_repeats) if args.simple_repeats.is_file() else {}
    black = read_bed_by_chrom(args.blacklist) if args.blacklist.is_file() else {}

    tables = sorted(args.labels_dir.glob("*.labeled.tsv"))
    if not tables:
        print(f"no labeled.tsv under {args.labels_dir}", file=sys.stderr)
        return 1

    fasta = pysam.FastaFile(str(args.reference))
    rows: list[dict] = []
    try:
        for path in tables:
            with path.open(encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    sample = row.get("sample") or path.name.split(".")[0]
                    m = meta.get(sample, {})
                    chrom = row["chrom"]
                    start = int(row["minbp"])
                    end = int(row["maxbp"])
                    out = {
                        "sample": sample,
                        "superpopulation": m.get("superpopulation", ""),
                        "population": m.get("population", ""),
                        "sex": m.get("sex", ""),
                        "family_id": m.get("family_id", ""),
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "call_id": row.get("id", ""),
                    }
                    # raw caller + label columns
                    for key, val in row.items():
                        out[key] = val
                    out["frac_simple_repeat"] = _frac(start, end, simp.get(chrom, []))
                    out["frac_segdup"] = _frac(start, end, seg.get(chrom, []))
                    out["frac_blacklist"] = _frac(start, end, black.get(chrom, []))
                    out["frac_gap"] = _frac(start, end, gaps.get(chrom, []))
                    out["gc_flank_500"] = gc_fraction(fasta, chrom, start, end, FLANK)
                    out["coverage_ratio"] = float("nan")
                    out["length_is_lower_bound"] = inslen_is_lower_bound(row.get("inslen"))
                    rows.append(out)
    finally:
        fasta.close()

    df = pd.DataFrame(rows)
    if "exclude" not in df.columns:
        df["exclude"] = 0
    # Column order: identity, labels, remaining call cols, context, exclude
    call_cols = [
        c
        for c in df.columns
        if c not in set(IDENTITY_OUT + LABEL_OUT + CONTEXT_OUT + ["exclude"])
    ]
    ordered = [
        c
        for c in IDENTITY_OUT + LABEL_OUT + call_cols + CONTEXT_OUT + ["exclude"]
        if c in df.columns
    ]
    df = df[ordered]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = args.out_dir / "dataset_v0.parquet"
    tsv_path = args.out_dir / "dataset_v0.tsv.gz"
    summary_path = args.out_dir / "dataset_v0.summary.json"

    df.to_parquet(parquet_path, index=False)
    df.to_csv(tsv_path, sep="\t", index=False, compression="gzip", na_rep="NA")

    summary = summarize(df)
    summary["paths"] = {
        "parquet": str(parquet_path),
        "tsv_gz": str(tsv_path),
        "n_input_tables": len(tables),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    export = Path("/work/export/ml0")
    export.mkdir(parents=True, exist_ok=True)
    shutil.copy2(parquet_path, export / "dataset_v0.parquet")
    shutil.copy2(tsv_path, export / "dataset_v0.tsv.gz")
    shutil.copy2(summary_path, export / "dataset_v0.summary.json")

    print(
        f"wrote {len(df)} rows -> {parquet_path} + {tsv_path}; "
        f"labels 1/0/NA = {summary['label_totals']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
