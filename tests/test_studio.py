"""ARC-SV Studio helpers (U6)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from arcsv2.dashboard import data
from arcsv2.dashboard.studio import checks, jobs


def test_rearrangement_plain():
    row = {"reference": "ABC", "rearrangement": "ABBC", "svtype": "Dup", "complextype": "NA"}
    text = data.rearrangement_plain(row)
    assert "segment" in text.lower() or "ABBC" in text


def test_batch_status_loader(tmp_path: Path):
    path = tmp_path / "batch_status.tsv"
    path.write_text("sample\tstatus\tmessage\tupdated_utc\nS1\tdone\tok\tt\n", encoding="utf-8")
    frame = data.load_batch_status(path)
    assert len(frame) == 1
    assert frame.iloc[0]["status"] == "done"


def test_checks_cpu():
    r = checks.check_cpu()
    assert r.name == "CPU cores"


def test_manifest_roundtrip(tmp_path: Path):
    cram = tmp_path / "s.cram"
    cram.write_bytes(b"x")
    manifest = tmp_path / "m.tsv"
    jobs.write_batch_manifest([("S1", cram)], manifest)
    assert "S1" in manifest.read_text(encoding="utf-8")


def test_install_example(tmp_path: Path):
    dest = jobs.install_example_run(tmp_path)
    assert dest.is_dir()
    runs = data.discover_runs(dest)
    assert runs


def test_example_package_has_bam():
    from arcsv2.dashboard.studio.example import package_dir

    pkg = package_dir()
    assert (pkg / "example.bam").is_file()
    assert (pkg / "ref.fa").is_file()


def test_infer_progress(tmp_path: Path):
    from arcsv2.dashboard.studio import progress

    sample = tmp_path / "S1"
    (sample / "chr20").mkdir(parents=True)
    (sample / "chr20" / "arcsv_out.tab").write_text("chrom\n", encoding="utf-8")
    prog = progress.infer_progress(sample, ["chr20"])
    assert prog["stage"] in ("merge", "done", "call", "qc")
    assert prog["percent"] >= 10
