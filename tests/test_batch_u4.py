"""U4 batch, mmi validation, and atomic I/O tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from arcsv2.atomic_io import atomic_write_text
from arcsv2.batch import (
    load_status_table,
    parse_batch_manifest,
    write_failed_txt,
    write_status_table,
)
from arcsv2.input_checks import InputCheckError
from arcsv2.mmi_validate import (
    check_no_hit_fraction,
    mmi_smoke_aligns,
    no_hit_fraction,
    validate_mmi_size,
)
from arcsv2.qc_report import build_qc_payload, write_qc_report


def test_atomic_write(tmp_path: Path):
    path = tmp_path / "sub" / "out.txt"
    atomic_write_text(path, "hello")
    assert path.read_text(encoding="utf-8") == "hello"


def test_batch_manifest(tmp_path: Path):
    manifest = tmp_path / "batch.tsv"
    manifest.write_text(
        "# comment\nsample\tinput\treference\nS1\t/path/a.cram\tGRCh38\nS2\t/path/b.cram\n",
        encoding="utf-8",
    )
    rows = parse_batch_manifest(manifest)
    assert len(rows) == 2
    assert rows[0]["sample"] == "S1"
    assert rows[1]["reference"] == ""


def test_status_table_roundtrip(tmp_path: Path):
    path = tmp_path / "batch_status.tsv"
    table = {
        "A": {"sample": "A", "status": "done", "message": "ok", "updated_utc": "t"},
    }
    write_status_table(path, table)
    loaded = load_status_table(path)
    assert loaded["A"]["status"] == "done"


def test_failed_txt(tmp_path: Path):
    p = write_failed_txt(tmp_path, "something broke")
    assert "something broke" in p.read_text(encoding="utf-8")


def test_mmi_size_too_small(tmp_path: Path):
    fa = tmp_path / "hap1.fa.gz"
    fa.write_bytes(b"x" * 1000)
    mmi = tmp_path / "hap1.mmi"
    mmi.write_bytes(b"y" * 100)
    with pytest.raises(InputCheckError, match="too small"):
        validate_mmi_size(mmi, fa, hap_name="hap1")


def test_mmi_validate_cli_rejects_small_index(tmp_path: Path):
    from typer.testing import CliRunner

    from arcsv2.cli import app

    fa = tmp_path / "hap1.fa.gz"
    fa.write_bytes(b"x" * 1000)
    mmi = tmp_path / "hap1.mmi"
    mmi.write_bytes(b"y" * 100)
    result = CliRunner().invoke(app, ["mmi-validate", str(mmi), "--fasta", str(fa)])
    assert result.exit_code != 0
    assert "too small" in result.output


def test_index_haps_runs_mmi_validate():
    text = Path(__file__).resolve().parents[1].joinpath("workflows", "label.smk").read_text(encoding="utf-8")
    assert "arcsv2 mmi-validate {output:q} --fasta {input:q}" in text


def test_mmi_smoke_fails_closed_without_aligner(tmp_path: Path, monkeypatch):
    real_import = __import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if name in {"mappy", "pysam"}:
            raise ImportError(name)
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", blocked)
    mmi = tmp_path / "hap.mmi"
    mmi.write_bytes(b"x")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    with pytest.raises(InputCheckError, match="mappy and pysam"):
        mmi_smoke_aligns(mmi, ref)


def test_mmi_smoke_fails_closed_when_chrom_or_read_fails(tmp_path: Path, monkeypatch):
    class Fasta:
        def __init__(self, path, references, fetch_error=None):
            self.references = references
            self._error = fetch_error

        def fetch(self, chrom, start, end):
            if self._error is not None:
                raise self._error
            return "A" * 1000

        def close(self):
            return None

    class Pysam:
        def __init__(self, references, fetch_error=None):
            self._references = references
            self._error = fetch_error

        def FastaFile(self, path):
            return Fasta(path, self._references, self._error)

    class Hit:
        is_primary = True

    class Aligner:
        def __init__(self, *args, **kwargs):
            return None

        def __bool__(self):
            return True

        def map(self, seq):
            yield Hit()

    class Mappy:
        pass

    Mappy.Aligner = Aligner

    mmi = tmp_path / "hap.mmi"
    mmi.write_bytes(b"x" * 100)
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + ("A" * 200) + "\n", encoding="utf-8")

    monkeypatch.setitem(sys.modules, "pysam", Pysam(["chr1"]))
    monkeypatch.setitem(sys.modules, "mappy", Mappy())
    with pytest.raises(InputCheckError, match="no chr20"):
        mmi_smoke_aligns(mmi, ref)

    monkeypatch.setitem(sys.modules, "pysam", Pysam(["chr1", "chr20"], OSError("unreadable")))
    with pytest.raises(InputCheckError, match="Could not fetch"):
        mmi_smoke_aligns(mmi, ref)


def test_mmi_smoke_probe_is_outside_chr20_centromere():
    from arcsv2.mmi_validate import SMOKE_PROBES

    assert ("chr20", 40_000_000, 40_001_000) in SMOKE_PROBES
    assert ("chr1", 40_000_000, 40_001_000) in SMOKE_PROBES
    assert all(not (chrom == "chr20" and start < 31_000_000) for chrom, start, _end in SMOKE_PROBES)


def test_no_hit_fraction_gate():
    rows = [{"label": "0", "reason": "no_hit"}] * 3 + [{"label": "1", "reason": "validated"}] * 97
    frac, nh, tot = no_hit_fraction(rows)
    assert nh == 3 and tot == 100
    assert frac == pytest.approx(0.03)
    with pytest.raises(InputCheckError, match="no_hit"):
        check_no_hit_fraction(rows)
    rows_ok = [{"label": "1", "reason": "validated"}] * 50
    check_no_hit_fraction(rows_ok)


def test_qc_report_files(tmp_path: Path):
    payload = build_qc_payload(
        {
            "sample": "S1",
            "runtime_sec": 1.0,
            "median_read_length": 150,
            "median_insert_size": 400,
            "coverage_x": 12.0,
            "merged_calls": 10,
            "stats": {"total": 12, "by_type": {"Del": 5}, "size_bins": {}},
        }
    )
    assert "low_coverage" in payload["flags"][0]
    j, h = write_qc_report(tmp_path, payload)
    assert j.is_file() and h.is_file()
