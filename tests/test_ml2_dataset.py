"""dataset_v1 builder, test-split exclusion, and bedcov retry."""

from __future__ import annotations

import csv
import gzip
import json
import os
import sys
from pathlib import Path

import pytest
from scripts.build_ml2_dataset_v1 import main

from arcsv2.coverage import bedcov
from arcsv2.ml1 import _outer_chrom_folds, _outer_sample_folds, apply_exclude_splits

CONTEXT = [
    "frac_simple_repeat",
    "frac_segdup",
    "frac_blacklist",
    "frac_gap",
    "gc_flank_500",
    "coverage_ratio",
    "length_is_lower_bound",
]


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = list(rows[0])
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _call(sample: str, call_id: str, label: str, superpopulation: str = "AFR") -> dict[str, str]:
    row = {
        "sample": sample,
        "superpopulation": superpopulation,
        "population": "ACB",
        "sex": "female",
        "family_id": "F1",
        "chrom": "chr1",
        "minbp": "100",
        "maxbp": "200",
        "id": call_id,
        "svtype": "DEL",
        "label": label,
        "sv_class": "deletion",
    }
    for name in CONTEXT:
        row[name] = "0" if name != "coverage_ratio" else "NA"
    row["length_is_lower_bound"] = "False"
    return row


def _ml0(sample: str, call_id: str, label: str) -> dict[str, str]:
    row = _call(sample, call_id, label)
    row["start"] = row.pop("minbp")
    row["end"] = row.pop("maxbp")
    row["call_id"] = row["id"]
    return row


def _read_out(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _build(root: Path, **extra: str) -> list[dict[str, str]]:
    out = root / "out"
    argv = [
        "--data-root",
        str(root),
        "--out",
        str(out),
        "--ml0",
        str(root / "ml0" / "dataset_v0.tsv"),
        "--config-dir",
        str(root / "config"),
        "--test-samples",
        str(root / "config" / "samples.ml2.test.tsv"),
        "--labels-dir",
        str(root / "labels" / "batch"),
        "--labels-dir",
        str(root / "labels" / "pilot"),
        "--labels-dir",
        str(root / "labels" / "test"),
        "--results-dir",
        str(root / "results"),
        "--backfill-dir",
        str(root / "backfill"),
    ]
    for key, value in extra.items():
        argv.extend([key, value])
    assert main(argv) == 0
    return _read_out(out / "dataset_v1.tsv.gz")


def _fixture(root: Path, pilot_label: str = "1") -> None:
    _write_tsv(
        root / "ml0" / "dataset_v0.tsv",
        [_ml0("M0A", "c1", "1"), _ml0("M0B", "c2", "0")],
    )
    _write_tsv(root / "labels" / "batch" / "A.labeled.tsv", [_call("ADD1", "c4", "1")])
    _write_tsv(root / "labels" / "batch" / "B.labeled.tsv", [_call("ADD2", "c3", "1", "EUR")])
    _write_tsv(
        root / "labels" / "batch" / "dataset.tsv",
        [{"sample": "ADD2", "id": "c3", "label": "1"}],
    )
    _write_tsv(root / "labels" / "pilot" / "A.labeled.tsv", [_call("ADD1", "c4", pilot_label)])
    _write_tsv(root / "labels" / "test" / "T.labeled.tsv", [_call("HOLD", "c5", "0", "EAS")])
    _write_tsv(
        root / "labels" / "test" / "dataset.tsv",
        [
            {
                "sample": "HOLD",
                "id": "c5",
                "coverage_ratio": "1.5",
                "coverage_ratio_chrom": "0.8",
            }
        ],
    )
    _write_tsv(
        root / "results" / "run1" / "ADD2" / "chr1" / "coverage.tsv",
        [
            {
                "id": "c3",
                "chrom": "chr1",
                "cov_inside": "4",
                "cov_flank": "4",
                "cov_chrom_median": "4",
                "coverage_ratio": "1.25",
                "coverage_ratio_chrom": "1.0",
            }
        ],
    )
    _write_tsv(
        root / "backfill" / "M0A.coverage.tsv",
        [{"id": "c1", "coverage_ratio": "2.0", "coverage_ratio_chrom": "1.1"}],
    )
    (root / "config").mkdir(parents=True)
    (root / "config" / "samples.ml2.test.tsv").write_text(
        "# frozen test\nsample\tsuperpopulation\nHOLD\tEAS\n",
        encoding="utf-8",
        newline="\n",
    )


def test_dataset_v1_assigns_splits_and_coverage(tmp_path: Path):
    _fixture(tmp_path)
    rows = _build(tmp_path)
    by_sample = {row["sample"]: row for row in rows}
    assert set(by_sample) == {"M0A", "M0B", "ADD1", "ADD2", "HOLD"}
    assert [row["split"] for row in rows if row["sample"] == "ADD1"] == ["train_add"]
    assert by_sample["M0A"]["split"] == "ml0"
    assert by_sample["M0B"]["split"] == "ml0"
    assert by_sample["ADD2"]["split"] == "train_add"
    assert by_sample["HOLD"]["split"] == "test"
    assert by_sample["M0A"]["coverage_source"] == "backfill"
    assert by_sample["M0A"]["coverage_ratio"] == "2.0"
    assert by_sample["M0A"]["coverage_ratio_chrom"] == "1.1"
    assert by_sample["M0B"]["coverage_source"] == "missing"
    assert by_sample["M0B"]["coverage_ratio"] == "NA"
    assert by_sample["M0B"]["coverage_ratio_chrom"] == "NA"
    assert by_sample["ADD2"]["coverage_source"] == "computed"
    assert by_sample["ADD2"]["coverage_ratio"] == "1.25"
    assert by_sample["HOLD"]["coverage_source"] == "computed"
    assert by_sample["HOLD"]["coverage_ratio"] == "1.5"
    assert by_sample["ADD1"]["coverage_source"] == "missing"
    summary = json.loads((tmp_path / "out" / "dataset_v1.summary.json").read_text(encoding="utf-8"))
    assert summary["label_counts"]["ml0"]["deletion"] == {"0": 1, "1": 1, "NA": 0}
    assert summary["superpopulation_counts"]["test"] == {"EAS": 1}
    assert summary["coverage_source_counts"]["ml0"] == {"computed": 0, "backfill": 1, "missing": 1}
    assert 'blacklist: ""' in summary["blacklist"]
    assert summary["inputs"]["ml0"].endswith("dataset_v0.tsv")
    assert summary["notes"]["coverage_imputed"] is False


def test_dataset_v1_rejects_conflicting_labels(tmp_path: Path):
    _fixture(tmp_path, pilot_label="0")
    with pytest.raises(SystemExit, match="conflicting labels for ADD1 c4"):
        _build(tmp_path)


def test_dataset_v1_rejects_test_sample_in_ml0(tmp_path: Path):
    _fixture(tmp_path)
    text = (tmp_path / "config" / "samples.ml2.test.tsv").read_text(encoding="utf-8")
    (tmp_path / "config" / "samples.ml2.test.tsv").write_text(text + "M0A\tAFR\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="both ML-0 and samples.ml2.test.tsv"):
        _build(tmp_path)


def test_exclude_split_drops_samples_before_folds():
    rows = []
    for index, sample in enumerate(["S1", "S2", "S3", "S4"]):
        chroms = ("chr1", "chr2") if index % 2 == 0 else ("chr2", "chr1")
        for chrom in chroms:
            rows.append(
                {
                    "sample": sample,
                    "superpopulation": "AFR",
                    "chrom": chrom,
                    "label": "1",
                    "svtype": "DEL",
                    "split": "ml0",
                }
            )
    rows.append(
        {
            "sample": "HOLD",
            "superpopulation": "AFR",
            "chrom": "chr1",
            "label": "0",
            "svtype": "DEL",
            "split": "test",
        }
    )
    rows.append(
        {
            "sample": "HOLD",
            "superpopulation": "AFR",
            "chrom": "chr2",
            "label": "1",
            "svtype": "DEL",
            "split": "test",
        }
    )
    kept, info = apply_exclude_splits(rows, ["test"])
    assert info["samples"] == 1
    assert info["rows"] == 2
    assert info["sample_ids"] == ["HOLD"]
    assert {row["sample"] for row in kept} == {"S1", "S2", "S3", "S4"}
    for train, test in _outer_sample_folds(kept, 2):
        assert "HOLD" not in train
        assert "HOLD" not in test
        assert set(train).isdisjoint(test)
    seen: set[str] = set()
    for train, test in _outer_sample_folds(kept, 2):
        seen.update(train)
        seen.update(test)
    assert seen == {"S1", "S2", "S3", "S4"}
    chroms: set[str] = set()
    for train, test in _outer_chrom_folds(kept, 2):
        assert set(train).isdisjoint(test)
        chroms.update(train)
        chroms.update(test)
    assert chroms == {"chr1", "chr2"}
    untouched, empty = apply_exclude_splits([{"sample": "S1"}], [])
    assert empty["rows"] == 0
    assert untouched == [{"sample": "S1"}]
    with pytest.raises(ValueError, match="no split column"):
        apply_exclude_splits([{"sample": "S1", "label": "1"}], ["test"])


_FAKE_SAMTOOLS = """
import os
import sys
from pathlib import Path

args = sys.argv[1:]
if not args or args[0] != "bedcov":
    sys.stderr.write("expected bedcov\\n")
    sys.exit(2)
index = 1
bed = None
while index < len(args):
    token = args[index]
    if token in ("-Q", "--reference"):
        index += 2
        continue
    if token.startswith("-"):
        index += 1
        continue
    bed = token
    break
state = Path(os.environ["FAKE_SAMTOOLS_STATE"])
count = int(state.read_text(encoding="utf-8") or "0") + 1
state.write_text(str(count), encoding="utf-8")
if count <= int(os.environ.get("FAKE_SAMTOOLS_FAIL", "0")):
    sys.stderr.write("retry_reconnect samtools bedcov: error reading from input file\\n")
    sys.exit(1)
for line in Path(bed).read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    chrom, start, end, name = line.split("\\t")[:4]
    length = int(end) - int(start)
    sys.stdout.write(f"{chrom}\\t{start}\\t{end}\\t{name}\\t{length * 3}\\n")
"""


def _install_fake_samtools(bindir: Path, state: Path) -> None:
    bindir.mkdir(parents=True)
    script = bindir / "fake_samtools.py"
    script.write_text(_FAKE_SAMTOOLS, encoding="utf-8", newline="\n")
    if os.name == "nt":
        launcher = bindir / "samtools.bat"
        launcher.write_text(
            f'@echo off\r\n"{sys.executable}" "{script}" %*\r\nexit /b %ERRORLEVEL%\r\n',
            encoding="utf-8",
        )
    else:
        launcher = bindir / "samtools"
        launcher.write_text(
            f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n',
            encoding="utf-8",
            newline="\n",
        )
        launcher.chmod(0o755)
    os.environ["FAKE_SAMTOOLS_STATE"] = str(state)
    state.write_text("0", encoding="utf-8")


def test_bedcov_retry_matches_one_shot(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    bindir = tmp_path / "bin"
    state = tmp_path / "state.txt"
    _install_fake_samtools(bindir, state)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    monkeypatch.setenv("FAKE_SAMTOOLS_STATE", str(state))
    regions = [("chr1", 0, 100, "w0"), ("chr1", 100, 250, "w1")]
    expected = {"w0": 300, "w1": 450}

    def run(fail: int, chunk: int, retries: int = 5, cache: Path | None = None):
        state.write_text("0", encoding="utf-8")
        monkeypatch.setenv("FAKE_SAMTOOLS_FAIL", str(fail))
        delays: list[float] = []
        got = bedcov(
            "reads.bam",
            regions,
            chunk_size=chunk,
            retries=retries,
            cache_dir=cache,
            sleeper=delays.append,
        )
        return got, delays

    once, once_delays = run(0, 500)
    retried, delays = run(2, 500)
    chunked, chunk_delays = run(2, 1)
    assert once == retried == chunked == expected
    assert once_delays == []
    assert delays == [5.0, 15.0]
    assert chunk_delays == [5.0, 15.0]

    cache = tmp_path / "chunks"
    cached, _ = run(0, 1, cache=cache)
    assert cached == expected
    assert (cache / "chunk_00000.tsv").is_file()
    state.write_text("0", encoding="utf-8")
    monkeypatch.setenv("FAKE_SAMTOOLS_FAIL", "100")
    resumed = bedcov(
        "reads.bam",
        regions,
        chunk_size=1,
        retries=1,
        cache_dir=cache,
        sleeper=lambda _delay: (_ for _ in ()).throw(AssertionError("cached chunk was rerun")),
    )
    assert resumed == expected
    assert state.read_text(encoding="utf-8") == "0"

    state.write_text("0", encoding="utf-8")
    monkeypatch.setenv("FAKE_SAMTOOLS_FAIL", "100")
    exhausted: list[float] = []
    with pytest.raises(RuntimeError, match="retry_reconnect"):
        bedcov(
            "reads.bam",
            regions,
            chunk_size=500,
            retries=5,
            sleeper=exhausted.append,
        )
    assert exhausted == [5.0, 15.0, 45.0, 135.0]
