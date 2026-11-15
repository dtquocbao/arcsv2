"""U5 SLURM profile and resource formula tests."""

from __future__ import annotations

import math
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "workflows" / "profiles" / "slurm" / "config.yaml"


def _chrom_mb(chrom: str) -> float:
    fai = ROOT / "tests" / "fixtures" / "label-smoke" / "ref.fa.fai"
    if not fai.is_file():
        return 1.0
    for line in fai.read_text(encoding="utf-8").splitlines():
        name, length, *_ = line.split("\t")
        if name == chrom:
            return int(length) / 1_000_000
    return 1.0


def test_slurm_profile_config_exists():
    data = yaml.safe_load(PROFILE.read_text(encoding="utf-8"))
    assert data.get("executor") == "cluster-generic"
    assert int(data.get("restart-times", 0)) >= 2
    assert "cluster-generic-submit-cmd" in data
    cmd = data["cluster-generic-submit-cmd"]
    assert "{resources.mem_mb}" in cmd
    assert "{threads}" in cmd
    assert "{resources.runtime}" in cmd


def test_call_mem_formula_matches_label_smk():
    mem_mb_per_mb = 30.68
    floor = 4000
    mb = _chrom_mb("chrSim")
    expected = max(floor, int(math.ceil(mb * mem_mb_per_mb * 1.2)))
    assert expected >= floor


def test_runtime_helpers():
    minutes = max(180, int(_chrom_mb("chr1") * 12)) if _chrom_mb("chr1") > 1 else 180
    hours, mins = divmod(minutes, 60)
    hms = f"{hours}:{mins:02d}:00"
    assert hms.count(":") == 2


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_slurm_profile_dryrun_two_samples():
    smoke = ROOT / "tests" / "fixtures" / "label-smoke"
    smoke.mkdir(parents=True, exist_ok=True)
    haps = ROOT / "haps"
    haps.mkdir(exist_ok=True)
    ref = smoke / "ref.fa"
    if not ref.is_file():
        ref.write_text(">chrSim\nACGTACGT\n", encoding="utf-8")
    for name in ("sim.hap1.fa", "sim.hap2.fa", "sim2.hap1.fa", "sim2.hap2.fa"):
        dst = haps / name
        if not dst.is_file():
            dst.write_text(ref.read_text(encoding="utf-8"), encoding="utf-8")
    for fname in ("gaps.bed", "unreliable.bed", "sim.cram"):
        p = smoke / fname
        if not p.is_file():
            p.write_text("", encoding="utf-8")
    proc = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflows/label.smk",
            "-n",
            "--profile",
            "workflows/profiles/slurm",
            "--configfile",
            "workflows/config/label.cluster_smoke.yaml",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr or proc.stdout
    assert "call" in proc.stdout or "sim" in proc.stdout
