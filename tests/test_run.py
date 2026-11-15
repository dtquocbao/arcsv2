"""Tests for ``arcsv2 run`` helpers (U2)."""

from __future__ import annotations

from pathlib import Path

import pytest

from arcsv2.input_checks import InputCheckError
from arcsv2.run import (
    PRIMARY_CHROMS,
    _gaps_from_fasta,
    _summarize_calls,
    detect_chrom_style,
    merge_chrom_tables,
    parse_chroms,
    resolve_reference,
    resolve_threads,
)


def test_parse_chroms_default():
    assert parse_chroms(None) == PRIMARY_CHROMS
    assert parse_chroms("") == PRIMARY_CHROMS
    assert parse_chroms("chr20, chr21") == ("chr20", "chr21")


def test_resolve_threads_auto():
    assert resolve_threads("auto") >= 1
    assert resolve_threads("4") == 4


def test_gaps_and_local_reference(tmp_path: Path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGTNNNAC\n>chr2\nA\n", encoding="utf-8")
    gaps = tmp_path / "gaps.bed"
    _gaps_from_fasta(fa, gaps)
    text = gaps.read_text(encoding="utf-8")
    assert "chr1\t4\t7\n" in text

    fasta, gaps_out = resolve_reference(str(fa), cache_dir=tmp_path / "cache")
    assert fasta == fa.resolve()
    assert gaps_out.is_file()


def test_detect_chrom_style(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr20\t1000\t0\t60\t61\nchr21\t1000\t0\t60\t61\n", encoding="utf-8")
    cram = tmp_path / "x.cram"
    cram.write_text("", encoding="utf-8")
    assert detect_chrom_style(cram, fai, ("chr20",)) == ("chr20",)
    assert detect_chrom_style(cram, fai, ("20",)) == ("chr20",)


def test_summarize_and_merge(tmp_path: Path):
    header = "chrom\tminbp\tmaxbp\tlen_affected\tsvtype\tcomplextype\tinslen\tgt"
    tab = tmp_path / "chr20" / "arcsv_out.tab"
    tab.parent.mkdir()
    tab.write_text(
        header + "\n"
        "chr20\t100\t200\t150\tDel\tNA\tNA\tHET\n"
        "chr20\t1000\t2000\t800\tDup\tNA\tNA\tHET\n",
        encoding="utf-8",
    )
    stats = _summarize_calls([tab])
    assert stats["total"] == 2
    assert stats["by_type"]["Del"] == 1

    merged = tmp_path / "merged.tab"
    n = merge_chrom_tables([tab], merged, sample="S1", min_size=50)
    assert n == 2
    body = merged.read_text(encoding="utf-8")
    assert body.count("\n") >= 2


def test_detect_chrom_style_missing(tmp_path: Path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t100\n", encoding="utf-8")
    cram = tmp_path / "x.cram"
    cram.write_text("", encoding="utf-8")
    with pytest.raises(InputCheckError, match="not in reference"):
        detect_chrom_style(cram, fai, ("chr99",))
