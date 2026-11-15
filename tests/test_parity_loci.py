"""Parity loci from tests/golden (caller-v1 / v2fix11 freeze)."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
GOLDEN = REPO / "tests" / "golden"
MANIFEST = json.loads((GOLDEN / "manifest.json").read_text(encoding="utf-8"))


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _entry(center: int) -> dict:
    for item in MANIFEST:
        if int(item["center"]) == center:
            return item
    raise KeyError(center)


def test_manifest_covers_required_loci():
    centers = {int(e["center"]) for e in MANIFEST}
    assert centers == {2659058, 5962348, 52764531, 60792115, 13358681, 23857951, 63458085}


def test_dup_abbc_2659058():
    rows = _rows(GOLDEN / _entry(2659058)["file"])
    assert len(rows) == 1
    assert rows[0]["svtype"] == "DUP"
    assert rows[0]["rearrangement"] == "ABBC"


def test_del_5962348():
    rows = _rows(GOLDEN / _entry(5962348)["file"])
    assert len(rows) == 1
    assert rows[0]["svtype"] == "DEL"


def test_dup_abbc_60792115():
    rows = _rows(GOLDEN / _entry(60792115)["file"])
    assert len(rows) == 1
    assert rows[0]["svtype"] == "DUP"
    assert rows[0]["rearrangement"] == "ABBC"


@pytest.mark.parametrize("center", [13358681, 23857951, 63458085])
def test_weak_del_not_called(center: int):
    rows = _rows(GOLDEN / _entry(center)["file"])
    assert rows == []


def test_same_block_support_52764531():
    rows = _rows(GOLDEN / _entry(52764531)["file"])
    assert len(rows) == 1
    assert rows[0]["svtype"] == "DUP"
    assert rows[0]["rearrangement"] == "ABBC"
    sr = int(rows[0]["sr_support"])
    pe = int(rows[0]["pe_support"])
    assert sr + pe == 2


def test_check_golden_fast_vs_export():
    """Instant: window slices match export/v2fix11 full merged.tab."""
    export_results = REPO / "export" / "v2fix11" / "results"
    if not export_results.is_dir():
        pytest.skip("export/v2fix11/results not present")
    proc = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / "check_golden.py"),
            "--fast",
            "--results",
            str(export_results),
        ],
        cwd=REPO,
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


@pytest.mark.slow
@pytest.mark.io
def test_window_recall_loci_against_bam(tmp_path: Path):
    """Optional BAM re-call of recall windows; skipped without /data BAMs."""
    ref = Path("/data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    gaps_bed = Path("/data/ref/gaps.bed")
    if not ref.is_file() or not gaps_bed.is_file():
        pytest.skip("/data ref not mounted")

    from arcsv2.beds import read_bed
    from arcsv2.config import CallConfig
    from arcsv2.pipeline import call_sam

    cfg = CallConfig()
    out_root = tmp_path / "results"
    # One merged.tab per sample accumulating all window calls for check_golden --fast.
    per_sample: dict[str, list[str]] = {}
    header: str | None = None
    for entry in MANIFEST:
        sample = entry["sample"]
        center = int(entry["center"])
        bam = Path(f"/data/work/{sample}/chr20.bam")
        if not bam.is_file():
            pytest.skip(f"missing {bam}")
        start0 = max(0, center - 20_000)
        end0 = center + 20_000
        tab = tmp_path / f"{sample}_{center}.tab"
        call_sam(
            str(bam),
            str(tab),
            cfg,
            chrom="chr20",
            start=start0,
            end=end0,
            reference=str(ref),
            gaps=read_bed(gaps_bed, "chr20"),
        )
        lines = tab.read_text(encoding="utf-8").splitlines()
        if not lines:
            continue
        header = lines[0]
        per_sample.setdefault(sample, []).extend(lines[1:])

    assert header is not None
    for sample, rows in per_sample.items():
        dest = out_root / sample
        dest.mkdir(parents=True, exist_ok=True)
        (dest / "merged.tab").write_text(
            header + "\n" + "\n".join(rows) + ("\n" if rows else ""),
            encoding="utf-8",
        )

    proc = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / "check_golden.py"),
            "--fast",
            "--results",
            str(out_root),
        ],
        cwd=REPO,
        check=False,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
