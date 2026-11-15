"""U3 input validation tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from arcsv2.input_checks import (
    InputCheckError,
    check_coverage,
    check_disk_space,
    check_index_present,
    check_mapping_stats,
    check_out_writable,
    check_read_length,
    check_reads_readable,
    check_reference_files,
    estimate_coverage_x,
    estimate_run_disk_bytes,
    index_path,
    validate_run_inputs_early,
)


def test_reads_missing(tmp_path: Path):
    missing = tmp_path / "nope.cram"
    with pytest.raises(InputCheckError, match="does not exist"):
        check_reads_readable(missing)


def test_reads_bad_extension(tmp_path: Path):
    bad = tmp_path / "reads.txt"
    bad.write_text("x", encoding="utf-8")
    with pytest.raises(InputCheckError, match="Unrecognized alignment"):
        check_reads_readable(bad)


def test_reads_ok(tmp_path: Path):
    cram = tmp_path / "s.cram"
    cram.write_bytes(b"\x00")
    check_reads_readable(cram)


def test_index_path_and_missing_bam(tmp_path: Path):
    bam = tmp_path / "s.bam"
    bam.write_bytes(b"")
    assert index_path(bam) == tmp_path / "s.bam.bai"
    with pytest.raises(InputCheckError, match="index missing"):
        check_index_present(bam, auto_build=False)


def test_index_auto_build_without_samtools(tmp_path: Path, monkeypatch):
    bam = tmp_path / "s.bam"
    bam.write_bytes(b"")
    monkeypatch.setattr("shutil.which", lambda _name: None)
    with pytest.raises(InputCheckError, match="samtools is not on PATH"):
        check_index_present(bam, auto_build=True)


def test_reference_and_gaps(tmp_path: Path):
    with pytest.raises(InputCheckError, match="Reference FASTA not found"):
        check_reference_files(tmp_path / "ref.fa")
    ref = tmp_path / "ref.fa"
    ref.write_text(">a\nA\n", encoding="utf-8")
    check_reference_files(ref)
    gaps = tmp_path / "gaps.bed"
    gaps.write_text("chr1\t0\t1\n", encoding="utf-8")
    check_reference_files(ref, gaps)


def test_out_writable(tmp_path: Path):
    check_out_writable(tmp_path / "out")


def test_read_length_bounds():
    check_read_length(150)
    with pytest.raises(InputCheckError, match="below"):
        check_read_length(30)
    with pytest.raises(InputCheckError, match="exceeds"):
        check_read_length(900)


def test_coverage_gate():
    check_coverage(15.0, "chr20")
    with pytest.raises(InputCheckError, match="Estimated coverage"):
        check_coverage(3.0, "chr20")


def test_estimate_coverage_from_idxstats():
    rows = [("chr20", 1000, 100)]
    cov = estimate_coverage_x(
        Path("dummy.cram"),
        "chr20",
        10_000,
        150,
        idxstats_rows=rows,
    )
    assert cov == pytest.approx(1.5)
    assert estimate_coverage_x(Path("x"), "chr21", 1000, 150, idxstats_rows=rows) == 0.0


def test_disk_space(tmp_path: Path, monkeypatch):
    class Usage:
        free = 1000

    monkeypatch.setattr("shutil.disk_usage", lambda _p: Usage())
    with pytest.raises(InputCheckError, match="Not enough disk space"):
        check_disk_space(tmp_path, 10_000_000_000)


def test_estimate_run_disk_grows_with_chroms(tmp_path: Path):
    cram = tmp_path / "s.cram"
    cram.write_bytes(b"x" * 1000)
    small = estimate_run_disk_bytes(cram, n_chroms=1, need_ref_download=False)
    large = estimate_run_disk_bytes(cram, n_chroms=22, need_ref_download=True)
    assert large > small


def test_mapping_stats():
    check_mapping_stats([("chr1", 1000, 100)])
    with pytest.raises(InputCheckError, match="No mapped reads"):
        check_mapping_stats([("chr1", 1000, 0)])
    with pytest.raises(InputCheckError, match="Too many unmapped"):
        check_mapping_stats([("chr1", 1000, 10), ("*", 0, 900)])


def test_validate_run_early(tmp_path: Path):
    cram = tmp_path / "s.cram"
    cram.write_bytes(b"x")
    validate_run_inputs_early(cram, "GRCh38", tmp_path / "out", ("chr20",))
