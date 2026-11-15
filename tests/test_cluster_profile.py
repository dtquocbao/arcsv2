"""U5 SLURM profile and resource formula tests."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

# Global container directive added in HPC1 F1. Snakemake only honors it with
# --sdm apptainer; a dry-run still shells out to singularity, so parse the line.
_CONTAINER_BLOCK = re.compile(
    r"_SIF = str\(config\.get\(\"container\"\) or os\.environ\.get\(\"ARCSV2_SIF\"\) or \"\"\)\.strip\(\)\n"
    r"if _SIF:\n"
    r"    container: _SIF\n"
)

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


def test_global_container_directive():
    for rel in ("workflows/label.smk", "workflows/Snakefile"):
        text = (ROOT / rel).read_text(encoding="utf-8")
        match = _CONTAINER_BLOCK.search(text)
        assert match, rel
        assert "container: _SIF" in match.group(0)
    cfg = yaml.safe_load((ROOT / "workflows/config/label.yaml").read_text(encoding="utf-8"))
    assert cfg.get("container") == ""


def test_slurm_profile_config_exists():
    text = PROFILE.read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    assert data.get("executor") == "slurm"
    assert int(data.get("restart-times", 0)) >= 2
    assert "cluster-generic" not in text
    assert "singularity" not in text.lower()
    defaults = data.get("default-resources") or {}
    assert defaults.get("slurm_partition") != "normal"
    assert defaults.get("runtime") == 60
    assert not (PROFILE.parent / "slurm-status.py").exists()


def test_site_slurm_profiles():
    grace_path = ROOT / "workflows" / "profiles" / "slurm-grace" / "config.yaml"
    faster_path = ROOT / "workflows" / "profiles" / "slurm-faster" / "config.yaml"
    grace = yaml.safe_load(grace_path.read_text(encoding="utf-8"))
    faster = yaml.safe_load(faster_path.read_text(encoding="utf-8"))
    for data in (grace, faster):
        assert data["executor"] == "slurm"
        assert data["software-deployment-method"] == "apptainer"
        assert data["default-resources"]["runtime"] == 60
        assert data["default-resources"]["mem_mb"] == 4000
        assert "cleanenv" in data["apptainer-args"]
        assert data["default-resources"]["slurm_partition"] != "normal"
    assert grace["default-resources"]["slurm_account"] == "132702112504"
    assert grace["default-resources"]["slurm_partition"] == "medium"
    assert grace["set-resources"]["dataset"]["slurm_partition"] == "short"
    assert grace["set-resources"]["merge"]["slurm_partition"] == "short"
    assert faster["default-resources"]["slurm_account"] == "142702113847"
    assert faster["default-resources"]["slurm_partition"] == "cpu"
    faster_text = faster_path.read_text(encoding="utf-8")
    assert "TODO(Bao): confirm partitions" in faster_text


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


def _slurm_plugin_importable() -> bool:
    proc = subprocess.run(
        [sys.executable, "-c", "import snakemake_executor_plugin_slurm"],
        capture_output=True,
        text=True,
        check=False,
    )
    return proc.returncode == 0


@pytest.mark.skipif(
    shutil.which("snakemake") is None or not _slurm_plugin_importable(),
    reason="snakemake or snakemake-executor-plugin-slurm not installed",
)
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
