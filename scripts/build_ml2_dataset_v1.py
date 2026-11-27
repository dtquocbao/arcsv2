#!/usr/bin/env python3
"""Build dataset_v1 from dataset_v0 plus the ML2 labeled tables.

ML-0 rows are copied from dataset_v0. ML2 rows come from ``*.labeled.tsv``
(NA labels kept). Context columns that are already present are kept. Rows that
lack ``frac_simple_repeat`` are annotated with the dataset_v0 helpers, which
import pysam, so that import happens only then.

Coverage is joined on ``(sample, call_id)`` and is never imputed. ``frac_blacklist``
is not recomputed: ML-0 and ML2 both labeled with an empty blacklist.

Run inside the 1.6.1 image. The login node is Python 3.6.

    singularity exec -B /scratch "$ARCSV2_SIF" env PYTHONPATH=python \\
      python scripts/build_ml2_dataset_v1.py \\
      --data-root "$SCRATCH/arcsv2" \\
      --out "$SCRATCH/arcsv2/repo/export/ml2"
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import os
from collections import Counter
from pathlib import Path

import pandas as pd

from arcsv2.train import class_of_row

EXTRA = ["coverage_ratio_chrom", "split", "coverage_source"]
SOURCE_RANK = {"missing": 0, "backfill": 1, "computed": 2}
_MISSING = {"", "NA", "na", "NaN", "nan", "None", "."}


def default_data_root() -> Path:
    scratch = os.environ.get("SCRATCH")
    if scratch:
        return Path(scratch) / "arcsv2"
    return Path("/scratch/arcsv2")


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    lines: list[str] = []
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            lines.append(raw)
    if not lines:
        return []
    return list(csv.DictReader(lines, delimiter="\t"))


def header_of(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            return raw.rstrip("\n").split("\t")
    return []


def load_meta(path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    if not path.is_file():
        return meta
    for row in read_tsv(path):
        sample = (row.get("sample") or "").strip()
        if sample:
            meta[sample] = row
    return meta


def as_coverage(value: object) -> str:
    """Coverage text for the table. Missing stays ``NA``; nothing is imputed."""
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return format(value, ".12g")
    if value is None:
        return "NA"
    text = str(value).strip()
    if text in _MISSING:
        return "NA"
    return text


def label_text(value: object) -> str:
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    if value is None:
        return "NA"
    text = str(value).strip()
    if text in _MISSING or text.lower() == "none":
        return "NA"
    return text


def _has_context(row: dict) -> bool:
    text = row.get("frac_simple_repeat")
    if text is None:
        return False
    return str(text).strip() not in ("",)


def normalize_call(row: dict, meta: dict[str, dict[str, str]]) -> dict:
    out = dict(row)
    sample = (out.get("sample") or "").strip()
    info = meta.get(sample, {})
    out["sample"] = sample
    for key in ("superpopulation", "population", "sex", "family_id"):
        if not str(out.get(key) or "").strip():
            out[key] = info.get(key, "")
    if not str(out.get("start") or "").strip() and out.get("minbp"):
        out["start"] = out["minbp"]
    if not str(out.get("end") or "").strip() and out.get("maxbp"):
        out["end"] = out["maxbp"]
    if not str(out.get("call_id") or "").strip():
        out["call_id"] = out.get("id") or ""
    if not str(out.get("id") or "").strip():
        out["id"] = out.get("call_id") or ""
    return out


def call_key(row: dict) -> tuple[str, str]:
    return (str(row.get("sample") or ""), str(row.get("call_id") or row.get("id") or ""))


def _v0_module():
    path = Path(__file__).resolve().with_name("build_ml0_dataset_v0.py")
    spec = importlib.util.spec_from_file_location("build_ml0_dataset_v0", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def annotate_context(
    rows: list[dict],
    reference: Path | None,
    gaps: Path | None,
    segdups: Path | None,
    simple_repeats: Path | None,
) -> None:
    """Fill context with the dataset_v0 helpers. Blacklist stays 0."""
    missing = [row for row in rows if not _has_context(row)]
    if not missing:
        return
    if reference is None or not reference.is_file():
        raise SystemExit(
            "labeled rows lack frac_simple_repeat; pass --reference, --gaps, "
            "--segdups, and --simple-repeats"
        )
    v0 = _v0_module()
    import pysam

    from arcsv2.beds import read_bed_by_chrom

    gap_bed = read_bed_by_chrom(gaps) if gaps is not None and gaps.is_file() else {}
    seg_bed = read_bed_by_chrom(segdups) if segdups is not None and segdups.is_file() else {}
    simp_bed = (
        read_bed_by_chrom(simple_repeats)
        if simple_repeats is not None and simple_repeats.is_file()
        else {}
    )
    fasta = pysam.FastaFile(str(reference))
    try:
        for row in missing:
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            row["frac_simple_repeat"] = v0._frac(start, end, simp_bed.get(chrom, []))
            row["frac_segdup"] = v0._frac(start, end, seg_bed.get(chrom, []))
            row["frac_blacklist"] = 0.0
            row["frac_gap"] = v0._frac(start, end, gap_bed.get(chrom, []))
            row["gc_flank_500"] = v0.gc_fraction(fasta, chrom, start, end, v0.FLANK)
            if not str(row.get("length_is_lower_bound") or "").strip():
                row["length_is_lower_bound"] = v0.inslen_is_lower_bound(row.get("inslen"))
            if "coverage_ratio" not in row:
                row["coverage_ratio"] = float("nan")
    finally:
        fasta.close()


def sample_from_coverage_path(path: Path) -> str:
    name = path.name
    if name.endswith(".coverage.tsv"):
        return name[: -len(".coverage.tsv")]
    if path.parent.name.lower().startswith("chr"):
        return path.parent.parent.name
    return path.parent.name


def index_coverage_file(index: dict[tuple[str, str], dict], path: Path, sample: str) -> None:
    for row in read_tsv(path):
        sid = (row.get("sample") or sample or "").strip()
        cid = (row.get("id") or row.get("call_id") or "").strip()
        if sid and cid:
            index[(sid, cid)] = row


def index_results(directories: list[Path]) -> dict[tuple[str, str], dict]:
    index: dict[tuple[str, str], dict] = {}
    for directory in directories:
        if not directory.is_dir():
            continue
        for path in sorted(directory.rglob("coverage.tsv")):
            index_coverage_file(index, path, sample_from_coverage_path(path))
        for path in sorted(directory.rglob("*.coverage.tsv")):
            index_coverage_file(index, path, sample_from_coverage_path(path))
    return index


def index_folder_datasets(folders: list[Path]) -> dict[tuple[str, str], dict]:
    """Coverage columns on a folder ``dataset.tsv`` win over per-chromosome files."""
    index: dict[tuple[str, str], dict] = {}
    for folder in folders:
        path = folder / "dataset.tsv"
        if not path.is_file():
            continue
        header = header_of(path)
        if "coverage_ratio" not in header or "coverage_ratio_chrom" not in header:
            continue
        for row in read_tsv(path):
            sample = (row.get("sample") or "").strip()
            cid = (row.get("id") or row.get("call_id") or "").strip()
            if sample and cid:
                index[(sample, cid)] = row
    return index


def apply_coverage(
    row: dict,
    computed: dict[tuple[str, str], dict],
    backfill: dict[tuple[str, str], dict],
) -> None:
    key = call_key(row)
    hit = computed.get(key)
    source = "computed"
    if hit is None:
        hit = backfill.get(key)
        source = "backfill"
    if hit is None:
        ratio = as_coverage(row.get("coverage_ratio"))
        chrom_ratio = as_coverage(row.get("coverage_ratio_chrom"))
        row["coverage_ratio"] = ratio
        row["coverage_ratio_chrom"] = chrom_ratio
        row["coverage_source"] = "missing" if ratio == "NA" and chrom_ratio == "NA" else "computed"
        return
    row["coverage_ratio"] = as_coverage(hit.get("coverage_ratio"))
    row["coverage_ratio_chrom"] = as_coverage(hit.get("coverage_ratio_chrom"))
    row["coverage_source"] = source


def discover_label_dirs(data_root: Path, given: list[Path] | None) -> list[Path]:
    if given:
        return given
    data = data_root / "data"
    names = [f"labels_ml2_batch{i:02d}" for i in range(1, 11)]
    names.extend(["labels_ml2_test", "labels_ml2_pilot"])
    found = [data / name for name in names if (data / name).is_dir()]
    if not found and data.is_dir():
        found = sorted(path for path in data.glob("labels_ml2_*") if path.is_dir())
    return found


def discover_results(data_root: Path, given: list[Path] | None) -> list[Path]:
    if given:
        return given
    data = data_root / "data"
    if not data.is_dir():
        return []
    return sorted(path for path in data.glob("results_ml2_*") if path.is_dir())


def find_ml0(data_root: Path, given: Path | None) -> Path:
    if given is not None:
        if not given.is_file():
            raise SystemExit(f"missing ML-0 dataset {given}")
        return given
    candidates = [
        data_root / "ml" / "dataset_v0.tsv.gz",
        data_root / "export" / "ml0" / "dataset_v0.tsv.gz",
        data_root / "repo" / "export" / "ml0" / "dataset_v0.tsv.gz",
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit("missing dataset_v0; pass --ml0 or place it under the data root")


def row_class(row: dict) -> str:
    sv = str(row.get("sv_class") or "").strip()
    if sv and sv not in _MISSING:
        return sv
    return class_of_row({key: "" if value is None else str(value) for key, value in row.items()}) or "NA"


def column_order(rows: list[dict]) -> list[str]:
    seen: list[str] = []
    for row in rows:
        for key in row:
            if key not in seen and key not in EXTRA:
                seen.append(key)
    for key in EXTRA:
        if key not in seen:
            seen.append(key)
    return seen


def summarize(rows: list[dict], inputs: dict) -> dict:
    labels: dict[str, dict[str, dict[str, int]]] = {}
    supers: dict[str, dict[str, str]] = {}
    sources: dict[str, dict[str, int]] = {}
    for row in rows:
        split = row["split"]
        kind = row_class(row)
        label = label_text(row.get("label"))
        bucket = labels.setdefault(split, {}).setdefault(kind, {"0": 0, "1": 0, "NA": 0})
        bucket[label] = bucket.get(label, 0) + 1
        supers.setdefault(split, {})[row["sample"]] = str(row.get("superpopulation") or "NA") or "NA"
        src = sources.setdefault(split, {"computed": 0, "backfill": 0, "missing": 0})
        src[row["coverage_source"]] = src.get(row["coverage_source"], 0) + 1
    return {
        "dataset_version": "ml2-dataset-v1",
        "n_rows": len(rows),
        "n_samples": len({row["sample"] for row in rows}),
        "label_counts": labels,
        "superpopulation_counts": {
            split: dict(Counter(mapping.values())) for split, mapping in supers.items()
        },
        "coverage_source_counts": sources,
        "blacklist": (
            'ML-0 and ML2 used blacklist: "". frac_blacklist was not recomputed '
            "(0 where the v0 builder had no bed)."
        ),
        "inputs": inputs,
        "notes": {"na_kept": True, "coverage_imputed": False},
    }


def _dedupe(rows: list[dict]) -> list[dict]:
    kept: dict[tuple[str, str], dict] = {}
    order: list[tuple[str, str]] = []
    for row in rows:
        key = call_key(row)
        if not key[0] or not key[1]:
            raise SystemExit(f"row missing sample or call id: {key}")
        previous = kept.get(key)
        if previous is None:
            kept[key] = row
            order.append(key)
            continue
        left = label_text(previous.get("label"))
        right = label_text(row.get("label"))
        if left != right:
            raise SystemExit(
                f"conflicting labels for {key[0]} {key[1]}: {left} vs {right}"
            )
        if SOURCE_RANK.get(row["coverage_source"], 0) > SOURCE_RANK.get(previous["coverage_source"], 0):
            kept[key] = row
    return [kept[key] for key in order]


def build(args: argparse.Namespace) -> int:
    data_root = args.data_root
    config_dir = args.config_dir or (Path(__file__).resolve().parents[1] / "workflows" / "config")
    test_path = args.test_samples or (config_dir / "samples.ml2.test.tsv")
    if not test_path.is_file():
        raise SystemExit(f"missing test sample list {test_path}")
    ml0_path = find_ml0(data_root, args.ml0)
    label_dirs = discover_label_dirs(data_root, args.labels_dir)
    if not label_dirs:
        raise SystemExit(f"no ML2 label folders under {data_root / 'data'}")
    results_dirs = discover_results(data_root, args.results_dir)
    backfill_dir = args.backfill_dir or (data_root / "ml" / "coverage_backfill")

    meta: dict[str, dict[str, str]] = {}
    if config_dir.is_dir():
        for path in sorted(config_dir.glob("samples.ml*.tsv")):
            meta.update(load_meta(path))
    meta.update(load_meta(test_path))

    test_samples = set(load_meta(test_path))
    ml0_rows = [normalize_call(row, meta) for row in read_tsv(ml0_path)]
    ml0_samples = {row["sample"] for row in ml0_rows}
    all40 = config_dir / "samples.ml.all40.tsv"
    if all40.is_file():
        ml0_samples |= set(load_meta(all40))
    overlap = sorted(test_samples & ml0_samples)
    if overlap:
        raise SystemExit(
            "samples are in both ML-0 and samples.ml2.test.tsv: " + ", ".join(overlap)
        )

    labeled_tables: list[Path] = []
    ml2_rows: list[dict] = []
    for folder in label_dirs:
        tables = sorted(folder.glob("*.labeled.tsv"))
        if not tables:
            raise SystemExit(f"no labeled.tsv under {folder}")
        labeled_tables.extend(tables)
        for path in tables:
            for raw in read_tsv(path):
                row = normalize_call(raw, meta)
                if not row["sample"]:
                    row["sample"] = path.name.split(".")[0]
                ml2_rows.append(row)
    annotate_context(ml2_rows, args.reference, args.gaps, args.segdups, args.simple_repeats)

    computed = index_results(results_dirs)
    computed.update(index_folder_datasets(label_dirs))
    backfill: dict[tuple[str, str], dict] = {}
    if backfill_dir.is_dir():
        for path in sorted(backfill_dir.glob("*.coverage.tsv")):
            index_coverage_file(backfill, path, sample_from_coverage_path(path))

    def assign(row: dict) -> None:
        sample = row["sample"]
        if sample in test_samples:
            row["split"] = "test"
        elif sample in ml0_samples:
            row["split"] = "ml0"
        else:
            row["split"] = "train_add"
        apply_coverage(row, computed, backfill)

    for row in ml0_rows:
        assign(row)
    for row in ml2_rows:
        assign(row)
    rows = _dedupe(ml0_rows + ml2_rows)
    rows.sort(
        key=lambda row: (
            row["sample"],
            str(row.get("chrom") or ""),
            str(row.get("start") or ""),
            row["call_id"],
        )
    )

    frame = pd.DataFrame(rows)
    order = [key for key in column_order(rows) if key in frame.columns]
    frame = frame[order]
    args.out.mkdir(parents=True, exist_ok=True)
    tsv_path = args.out / "dataset_v1.tsv.gz"
    summary_path = args.out / "dataset_v1.summary.json"
    frame.to_csv(tsv_path, sep="\t", index=False, compression="gzip", na_rep="NA")

    parquet_path = args.out / "dataset_v1.parquet"
    wrote_parquet = False
    try:
        import pyarrow  # noqa: F401

        frame.to_parquet(parquet_path, index=False)
        wrote_parquet = True
    except ImportError:
        parquet_path = None

    inputs = {
        "data_root": str(data_root),
        "ml0": str(ml0_path),
        "label_dirs": [str(path) for path in label_dirs],
        "labeled_tables": [str(path) for path in labeled_tables],
        "results_dirs": [str(path) for path in results_dirs],
        "backfill_dir": str(backfill_dir),
        "test_samples": str(test_path),
        "config_dir": str(config_dir),
        "tsv_gz": str(tsv_path),
        "parquet": str(parquet_path) if wrote_parquet else None,
    }
    summary = summarize(rows, inputs)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    totals = Counter(label_text(row.get("label")) for row in rows)
    print(
        f"wrote {len(rows)} rows -> {tsv_path}"
        + (f" + {parquet_path}" if wrote_parquet else "")
        + f"; labels 1/0/NA = {totals.get('1', 0)}/{totals.get('0', 0)}/{totals.get('NA', 0)}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--ml0", type=Path, default=None, help="dataset_v0 TSV or TSV.GZ")
    parser.add_argument("--config-dir", type=Path, default=None)
    parser.add_argument("--test-samples", type=Path, default=None)
    parser.add_argument("--labels-dir", type=Path, action="append", default=None)
    parser.add_argument("--results-dir", type=Path, action="append", default=None)
    parser.add_argument("--backfill-dir", type=Path, default=None)
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--gaps", type=Path, default=None)
    parser.add_argument("--segdups", type=Path, default=None)
    parser.add_argument("--simple-repeats", type=Path, default=None)
    args = parser.parse_args(argv)
    args.data_root = args.data_root or default_data_root()
    if args.out is None:
        args.out = args.data_root / "repo" / "export" / "ml2"
    if args.reference is None:
        args.reference = args.data_root / "ref" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    if args.gaps is None:
        args.gaps = args.data_root / "ref" / "gaps.bed"
    if args.segdups is None:
        args.segdups = args.data_root / "ref" / "segdups.bed"
    if args.simple_repeats is None:
        args.simple_repeats = args.data_root / "ref" / "simple_repeats.bed"
    return build(args)


if __name__ == "__main__":
    raise SystemExit(main())
