#!/usr/bin/env python3
"""Choose ML2 samples from catalogs already on disk. Does not download anything.

Candidates are HPRC R2 samples with both haplotypes and a 1KGP high-coverage
CRAM, minus the ML-0 40. At most one sample is kept per family, and a family
already used by an ML-0 sample is dropped. The rest are split with a fixed
seed, balanced by superpopulation: 20 frozen test samples, then train_add in
batches of 10, and never more than 120 new samples.

Reads only:

  <data>/catalog/1000G_2504_high_coverage.sequence.index
  <data>/catalog/1000G_698_related_high_coverage.sequence.index
  <data>/catalog/assemblies_release2_v1.0.<commit12>.csv
  <data>/catalog/pedigree/20130606_g1k.ped
  <data>/catalog/pedigree/20130606_sample_info.txt   (optional; sex)

Grace::

    cd $SCRATCH/arcsv2/repo
    singularity exec -B /scratch $SCRATCH/arcsv2/sif/arcsv2_1.6.1.sif \\
      python scripts/select_ml2_samples.py --data $SCRATCH/arcsv2/data
"""

from __future__ import annotations

import argparse
import os
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    from fetch_data import HPRC_COMMIT, intersect_candidates, parse_1kgp_index, parse_hprc_csv
except ImportError:  # pytest imports this as scripts.select_ml2_samples
    from scripts.fetch_data import (
        HPRC_COMMIT,
        intersect_candidates,
        parse_1kgp_index,
        parse_hprc_csv,
    )

SEED = 204178949
TEST_N = 20
MAX_NEW = 120
BATCH_N = 10
SUPERS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

POP_TO_SUPER = {
    "ACB": "AFR",
    "ASW": "AFR",
    "ESN": "AFR",
    "GWD": "AFR",
    "LWK": "AFR",
    "MSL": "AFR",
    "YRI": "AFR",
    "CLM": "AMR",
    "MXL": "AMR",
    "PEL": "AMR",
    "PUR": "AMR",
    "CDX": "EAS",
    "CHB": "EAS",
    "CHS": "EAS",
    "JPT": "EAS",
    "KHV": "EAS",
    "CEU": "EUR",
    "FIN": "EUR",
    "GBR": "EUR",
    "IBS": "EUR",
    "TSI": "EUR",
    "BEB": "SAS",
    "GIH": "SAS",
    "ITU": "SAS",
    "PJL": "SAS",
    "STU": "SAS",
}

TSV_COLUMNS = (
    "sample",
    "population",
    "superpopulation",
    "sex",
    "family_id",
    "split",
    "cram_url",
    "crai_url",
    "hap1_url",
    "hap2_url",
    "cram_checksum",
    "crai_checksum",
    "hap1_checksum",
    "hap2_checksum",
)


def parse_ped(text: str) -> dict[str, str]:
    """Individual id to family id. PED is whitespace-separated, optional # header."""
    family_of: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        fid, iid = parts[0], parts[1]
        if iid.lower() in {"individual_id", "iid", "sample"}:
            continue
        family_of[iid] = fid
    return family_of


def parse_sex(text: str) -> dict[str, str]:
    sex_map: dict[str, str] = {}
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        return sex_map
    header = lines[0].split("\t")
    for line in lines[1:]:
        cells = line.split("\t")
        row = {header[i]: cells[i] if i < len(cells) else "" for i in range(len(header))}
        sid = (row.get("Sample") or row.get("sample") or "").strip()
        gender = (row.get("Gender") or row.get("gender") or row.get("Sex") or "").strip().lower()
        if not sid or not gender:
            continue
        if gender.startswith("f"):
            sex_map[sid] = "female"
        elif gender.startswith("m"):
            sex_map[sid] = "male"
    return sex_map


def read_ml0(path: Path) -> tuple[set[str], set[str]]:
    ids: set[str] = set()
    families: set[str] = set()
    lines = path.read_text(encoding="utf-8").splitlines()
    header: list[str] | None = None
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cells = line.split("\t")
        if header is None:
            header = cells
            continue
        row = {header[i]: cells[i] if i < len(cells) else "" for i in range(len(header))}
        sample = row.get("sample", "").strip()
        if sample:
            ids.add(sample)
        family = row.get("family_id", "").strip()
        if family:
            families.add(family)
    return ids, families


def annotate(
    candidates: list[dict[str, str]],
    family_of: dict[str, str],
    sex_map: dict[str, str],
) -> None:
    for row in candidates:
        row["superpopulation"] = POP_TO_SUPER.get(row.get("population", ""), "UNK")
        row["family_id"] = family_of.get(row["sample"], row["sample"])
        row["sex"] = sex_map.get(row["sample"], row.get("sex", ""))


def apply_family_filter(
    candidates: list[dict[str, str]],
    ml0_ids: set[str],
    ml0_families: set[str],
) -> tuple[list[dict[str, str]], list[tuple[str, str, str]], int]:
    """Drop ML-0 ids, their families, and extra members of a new family.

    The kept member of a new family is the lowest sample id. Returns
    ``(kept, dropped, n_excluded_ml0)``. ``dropped`` rows are
    ``(sample, reason, family_id)``.
    """
    excluded = 0
    dropped: list[tuple[str, str, str]] = []
    by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        if row["sample"] in ml0_ids:
            excluded += 1
            continue
        family = row.get("family_id") or row["sample"]
        if family in ml0_families:
            dropped.append((row["sample"], "ml0_family", family))
            continue
        by_family[family].append(row)
    kept: list[dict[str, str]] = []
    for family, rows in by_family.items():
        rows.sort(key=lambda item: item["sample"])
        kept.append(rows[0])
        for extra in rows[1:]:
            dropped.append((extra["sample"], "duplicate_family", family))
    kept.sort(key=lambda item: item["sample"])
    return kept, dropped, excluded


def interleave(rows: list[dict[str, str]], seed: int, limit: int) -> list[dict[str, str]]:
    """Round-robin across superpopulations after a seeded shuffle within each."""
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row.get("superpopulation") or "UNK"].append(row)
    rng = random.Random(seed)
    queues: dict[str, list[dict[str, str]]] = {}
    for superpop in SUPERS + ["UNK"]:
        bucket = sorted(groups.get(superpop, []), key=lambda item: item["sample"])
        rng.shuffle(bucket)
        queues[superpop] = bucket
    order = [superpop for superpop in SUPERS + ["UNK"] if queues[superpop]]
    picked: list[dict[str, str]] = []
    while order and len(picked) < limit:
        progressed = False
        for superpop in order:
            if not queues[superpop]:
                continue
            picked.append(queues[superpop].pop(0))
            progressed = True
            if len(picked) >= limit:
                break
        if not progressed:
            break
    return picked


def choose(
    candidates: list[dict[str, str]],
    *,
    ml0_ids: set[str],
    ml0_families: set[str],
    seed: int = SEED,
    test_n: int = TEST_N,
    max_new: int = MAX_NEW,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[tuple[str, str, str]], int]:
    """Return ``(test, train_add, dropped, n_excluded_ml0)``.

    ``test`` has ``test_n`` samples when at least that many remain after the
    family filter. ``train_add`` is the rest, capped so test + train_add
    never exceeds ``max_new``.
    """
    kept, dropped, excluded = apply_family_filter(candidates, ml0_ids, ml0_families)
    picked = interleave(kept, seed, max_new)
    if len(picked) < test_n:
        return [], [], dropped, excluded
    test = picked[:test_n]
    train = picked[test_n:]
    for row in test:
        row["split"] = "test"
    for row in train:
        row["split"] = "train_add"
    return test, train, dropped, excluded


def super_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = Counter(row.get("superpopulation") or "UNK" for row in rows)
    return {superpop: counts.get(superpop, 0) for superpop in SUPERS + ["UNK"]}


def batches(rows: list[dict[str, str]], size: int = BATCH_N) -> list[list[dict[str, str]]]:
    return [rows[i : i + size] for i in range(0, len(rows), size)]


def write_tsv(path: Path, rows: list[dict[str, str]], comment: str) -> None:
    lines = [comment.rstrip(), "\t".join(TSV_COLUMNS)]
    for row in rows:
        lines.append("\t".join(row.get(col, "") for col in TSV_COLUMNS))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_label_yaml(
    path: Path,
    rows: list[dict[str, str]],
    *,
    scratch: str,
    results_dir: str,
    labels_dir: str,
    comment: str,
) -> None:
    """Absolute Grace paths. CRAMs stay as URLs (streamed, no prefetch)."""
    root = scratch.rstrip("/") + "/arcsv2"
    data = root + "/data"
    ref = root + "/ref"
    lines = [
        comment.rstrip(),
        'container: ""',
        f"data_root: {data}",
        f"results_dir: {results_dir}",
        f"labels_dir: {labels_dir}",
        f"reference: {ref}/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        f"gaps: {ref}/gaps.bed",
        f"unreliable: {ref}/unreliable.bed",
        f"simple_repeats: {ref}/simple_repeats.bed",
        f"segdups: {ref}/segdups.bed",
        f"blacklist: {ref}/blacklist.bed",
        "chroms:",
    ]
    for chrom in CHROMS:
        lines.append(f"  - {chrom}")
    lines.append("samples:")
    for row in rows:
        lines.append(f"  - {row['sample']}")
    lines.append("cram:")
    for row in rows:
        lines.append(f"  {row['sample']}: {row['cram_url']}")
    lines.extend(
        [
            "sequential: false",
            "threads: 4",
            "label_threads: 4",
            "label_mem_mb: 30000",
            "mem_mb_floor: 28000",
            "call_min_per_mb: 2.0",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read(path: Path) -> str:
    if not path.is_file():
        raise SystemExit(
            f"missing catalog {path}\n"
            "This script only reads files already under the data directory. "
            "It does not download catalogs or sample data."
        )
    return path.read_text(encoding="utf-8", errors="replace")


def load_candidates(data: Path, ped: Path | None) -> tuple[list[dict[str, str]], str]:
    catalog = data / "catalog"
    commit12 = HPRC_COMMIT[:12]
    kg_u = _read(catalog / "1000G_2504_high_coverage.sequence.index")
    kg_r = _read(catalog / "1000G_698_related_high_coverage.sequence.index")
    hprc = _read(catalog / f"assemblies_release2_v1.0.{commit12}.csv")
    ped_path = ped or (catalog / "pedigree" / "20130606_g1k.ped")
    ped_text = _read(ped_path)
    sex_path = catalog / "pedigree" / "20130606_sample_info.txt"
    sex_map = parse_sex(sex_path.read_text(encoding="utf-8", errors="replace")) if sex_path.is_file() else {}
    kg = parse_1kgp_index(kg_u, related=False)
    kg.update(parse_1kgp_index(kg_r, related=True))
    candidates = intersect_candidates(kg, parse_hprc_csv(hprc))
    annotate(candidates, parse_ped(ped_text), sex_map)
    return candidates, ped_path.name


def _print_counts(label: str, rows: list[dict[str, str]]) -> None:
    counts = super_counts(rows)
    parts = " ".join(f"{key}={value}" for key, value in counts.items())
    print(f"{label}\tn={len(rows)}\t{parts}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="Catalog root, e.g. $SCRATCH/arcsv2/data")
    parser.add_argument(
        "--ml0",
        type=Path,
        default=Path("workflows/config/samples.ml.all40.tsv"),
        help="ML-0 sample table (ids and family_id)",
    )
    parser.add_argument("--ped", type=Path, default=None, help="PED file (default: <data>/catalog/pedigree/20130606_g1k.ped)")
    parser.add_argument("--out-dir", type=Path, default=Path("workflows/config"))
    parser.add_argument("--scratch", default="", help="Grace $SCRATCH (default: $SCRATCH or /scratch/user/dtquocbao)")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)

    scratch = args.scratch or os.environ.get("SCRATCH") or "/scratch/user/dtquocbao"
    ml0_ids, ml0_families = read_ml0(args.ml0)
    # A relative is also blocked when the pedigree family of an ML-0 sample
    # differs from the family_id stored on the ML-0 table.
    candidates, ped_name = load_candidates(args.data, args.ped)
    family_of = {row["sample"]: row["family_id"] for row in candidates}
    for sample in ml0_ids:
        if sample in family_of:
            ml0_families.add(family_of[sample])

    test, train, dropped, excluded = choose(
        candidates,
        ml0_ids=ml0_ids,
        ml0_families=ml0_families,
        seed=args.seed,
    )
    print(f"pedigree={ped_name}")
    print(f"candidates={len(candidates)} excluded_ml0={excluded} dropped_family={len(dropped)}")
    reasons = Counter(reason for _sample, reason, _family in dropped)
    for reason, count in sorted(reasons.items()):
        print(f"  dropped_{reason}={count}")
    if len(test) < TEST_N:
        remaining = len(candidates) - excluded - len(dropped)
        print(
            f"only {remaining} samples remain after the family filter; "
            f"need {TEST_N} for the frozen test split. Wrote nothing.",
            file=sys.stderr,
        )
        return 2

    comment = (
        f"# ML2 selection seed={args.seed} ped={ped_name} "
        f"hprc={HPRC_COMMIT[:12]} dropped_family={len(dropped)}"
    )
    out = args.out_dir
    write_tsv(out / "samples.ml2.test.tsv", test, comment + " split=test (frozen; not for training)")
    print(f"wrote {out / 'samples.ml2.test.tsv'}")
    _print_counts("test", test)
    write_label_yaml(
        out / "label.ml2.test.yaml",
        test,
        scratch=scratch,
        results_dir="results_ml2_test",
        labels_dir="labels_ml2_test",
        comment=(
            "# Frozen ML2 test samples. Do not use for training, model selection, or thresholds.\n"
            "# CRAMs are streamed (no prefetch). Generated by scripts/select_ml2_samples.py"
        ),
    )
    print(f"wrote {out / 'label.ml2.test.yaml'}")

    for index, chunk in enumerate(batches(train), start=1):
        name = f"batch{index:02d}"
        tsv = out / f"samples.ml2.{name}.tsv"
        yaml_path = out / f"label.ml2.{name}.yaml"
        write_tsv(tsv, chunk, comment + f" split=train_add {name}")
        write_label_yaml(
            yaml_path,
            chunk,
            scratch=scratch,
            results_dir=f"results_ml2_{name}",
            labels_dir=f"labels_ml2_{name}",
            comment=(
                f"# ML2 train_add {name} ({len(chunk)} samples). CRAMs are streamed (no prefetch).\n"
                "# Generated by scripts/select_ml2_samples.py. Do not start until the pilot is reviewed."
            ),
        )
        print(f"wrote {tsv}")
        print(f"wrote {yaml_path}")
    _print_counts("train_add", train)
    print(f"total_new={len(test) + len(train)} test={len(test)} train_add={len(train)} batches={len(batches(train))}")
    if len(test) + len(train) < MAX_NEW:
        print(f"stopped_early kept={len(test) + len(train)} cap={MAX_NEW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
