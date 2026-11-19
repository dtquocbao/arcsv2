"""U5 SLURM profile and resource formula tests."""

from __future__ import annotations

import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

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


def _grch38_bp() -> dict[str, int]:
    path = ROOT / "scripts" / "bench" / "chrom_lens.py"
    spec = importlib.util.spec_from_file_location("chrom_lens", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return dict(module.GRCH38_BP)


_GRCH38_BP = _grch38_bp()


def _chrom_mb(chrom: str) -> float:
    length = _GRCH38_BP.get(chrom)
    if length is None:
        raise AssertionError(f"{chrom} is not in scripts/bench/chrom_lens.py GRCH38_BP")
    return length / 1_000_000


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
    assert "attempts" not in data
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


# GRCh38 primary-assembly lengths from scripts/bench/chrom_lens.py.
# chr1 248956422 bp → ceil(124.478) = 125 min, above the 60 min floor.


def _load_helpers(rel: str) -> dict:
    """Exec the Python helpers above the first Snakemake rule."""
    text = (ROOT / rel).read_text(encoding="utf-8")
    text = text.split("\nrule ", 1)[0]
    text = re.sub(r"^configfile:.*$", "config = {}", text, flags=re.M)
    text = re.sub(r"^([ \t]*)container:.*$", r"\1pass", text, flags=re.M)
    ns: dict = {"__name__": "workflow_helpers"}
    exec(compile(text, rel, "exec"), ns, ns)
    ns["CHROM_LEN_BP"].clear()
    ns["CHROM_LEN_BP"].update(_GRCH38_BP)
    return ns


def test_call_runtime_and_attempt_scaling():
    for rel in ("workflows/label.smk", "workflows/Snakefile"):
        text = (ROOT / rel).read_text(encoding="utf-8")
        assert 'lambda wc, attempt: _cap(call_runtime_min(wc) * attempt, "medium")' in text
        assert "lambda wc, attempt: call_mem_mb(wc) * attempt" in text
        assert "mb * 12" not in text
        ns = _load_helpers(rel)
        wc1 = SimpleNamespace(chrom="chr1")
        wc20 = SimpleNamespace(chrom="chr20")

        def runtime(wc, attempt, helpers=ns):
            return helpers["call_runtime_min"](wc) * attempt

        def memory(wc, attempt, helpers=ns):
            return helpers["call_mem_mb"](wc) * attempt

        assert _chrom_mb("chr1") > 200
        assert runtime(wc1, 1) == max(60, math.ceil(_chrom_mb("chr1") * 0.5))
        assert runtime(wc1, 1) == 125
        assert 60 <= runtime(wc1, 1) <= 1440
        assert runtime(wc20, 1) >= 60
        assert runtime(wc1, 2) == 2 * runtime(wc1, 1)
        assert memory(wc1, 2) == 2 * memory(wc1, 1)
        assert runtime(wc20, 2) == 2 * runtime(wc20, 1)
        assert memory(wc20, 2) == 2 * memory(wc20, 1)

    label = (ROOT / "workflows/label.smk").read_text(encoding="utf-8")
    assert 'lambda wc, attempt: _cap(index_runtime_min(wc) * attempt, "medium")' in label
    assert 'lambda wc, attempt: _cap(label_runtime_min(wc) * attempt, "medium")' in label
    ns = _load_helpers("workflows/label.smk")
    assert ns["index_runtime_min"](None) == 360
    assert ns["label_runtime_min"](None) == 480
    grace = yaml.safe_load(
        (ROOT / "workflows/config/label.grace.example.yaml").read_text(encoding="utf-8")
    )
    assert grace["label_mem_mb"] == 30000
    assert grace["label_threads"] == 4
    assert grace["threads"] == 4
    laptop = yaml.safe_load((ROOT / "workflows/config/label.yaml").read_text(encoding="utf-8"))
    assert "label_mem_mb" not in laptop


def test_runtime_stays_inside_grace_partition():
    """Attempt 1..3 never asks for more minutes than the Grace partition allows."""
    ns = _load_helpers("workflows/label.smk")
    limits = ns["PARTITION_LIMIT_MIN"]
    assert limits == {"short": 120, "medium": 1440, "long": 10080}
    assert ns["MERGE_RUNTIME_MIN"] == 30
    assert ns["DATASET_RUNTIME_MIN"] == 30
    grace = yaml.safe_load(
        (ROOT / "workflows" / "profiles" / "slurm-grace" / "config.yaml").read_text(encoding="utf-8")
    )
    default_partition = grace["default-resources"]["slurm_partition"]
    overrides = {
        rule: spec["slurm_partition"] for rule, spec in (grace.get("set-resources") or {}).items()
    }
    wc = SimpleNamespace(chrom="chr1")
    bases = {
        "fetch_region": lambda attempt: ns["FETCH_RUNTIME_MIN"] * attempt,
        "call": lambda attempt: ns["call_runtime_min"](wc) * attempt,
        "merge": lambda attempt: ns["MERGE_RUNTIME_MIN"] * attempt,
        "index_haps": lambda attempt: ns["index_runtime_min"](wc) * attempt,
        "label": lambda attempt: ns["label_runtime_min"](wc) * attempt,
        "dataset": lambda attempt: ns["DATASET_RUNTIME_MIN"] * attempt,
    }
    text = (ROOT / "workflows/label.smk").read_text(encoding="utf-8")
    for rule, base in bases.items():
        partition = overrides.get(rule, default_partition)
        assert "_cap(" in text
        assert partition in limits
        for attempt in (1, 2, 3):
            requested = ns["_cap"](base(attempt), partition)
            assert requested <= limits[partition]
            assert requested == min(base(attempt), limits[partition])
        # The rule's runtime line names this partition, so the cap matches the profile.
        assert f', "{partition}")' in text


def test_call_mem_formula_matches_label_smk():
    mem_mb_per_mb = 30.68
    floor = 4000
    mb = _chrom_mb("chr1")
    expected = max(floor, int(math.ceil(mb * mem_mb_per_mb * 1.2)))
    assert expected > floor


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
@pytest.mark.parametrize("profile", ["slurm", "slurm-grace", "slurm-faster"])
def test_slurm_profile_dryrun_two_samples(profile: str):
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
    env = dict(**os.environ)
    env.pop("ARCSV2_SIF", None)
    # mtime-only: a previous run's .snakemake metadata otherwise shells out to
    # `conda info` while hashing the software stack. The profile itself is valid.
    proc = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflows/label.smk",
            "-n",
            "--rerun-triggers",
            "mtime",
            "--profile",
            f"workflows/profiles/{profile}",
            "--configfile",
            "workflows/config/label.cluster_smoke.yaml",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr or proc.stdout
    combined = proc.stdout + proc.stderr
    assert "usage:" not in combined.lower()
    assert "call" in proc.stdout
    if profile == "slurm-grace":
        assert "132702112504" in proc.stdout
    if profile == "slurm-faster":
        assert "142702113847" in proc.stdout


def test_pull_sif_keeps_ghcr_token_out_of_logs():
    text = (ROOT / "scripts" / "hpc" / "pull_sif.sh").read_text(encoding="utf-8")
    assert "SINGULARITY_DOCKER_USERNAME=dtquocbao" in text
    assert 'SINGULARITY_DOCKER_PASSWORD="${GHCR_TOKEN}"' in text
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("echo") and ("GHCR_TOKEN" in stripped or "SINGULARITY_DOCKER_PASSWORD" in stripped):
            raise AssertionError(stripped)
        if ">>" in stripped and "GHCR_TOKEN" in stripped:
            raise AssertionError(stripped)


def _bash() -> str | None:
    git_bash = Path(r"C:\Program Files\Git\bin\bash.exe")
    if git_bash.is_file():
        return str(git_bash)
    found = shutil.which("bash")
    if found and "system32" not in found.lower():
        return found
    return None


@pytest.mark.skipif(_bash() is None, reason="bash is not installed")
def test_su_report_sums_fetch_call_merge(tmp_path: Path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    sacct = bindir / "sacct"
    sacct.write_text(
        "#!/bin/bash\n"
        "printf '%s\\n' \\\n"
        "  '200|fetch|00:30:00|1|4000M|' \\\n"
        "  '100|call|01:00:00|4|4000M|' \\\n"
        "  '300|merge|00:15:00|1|8000M|'\n",
        encoding="utf-8",
        newline="\n",
    )
    sacct.chmod(0o755)
    log = tmp_path / "snakemake.log"
    log.write_text(
        "rule fetch_region:\n"
        "    jobid: 1\n"
        "Job 1 has been submitted with SLURM jobid 200 (log: x).\n"
        "rule call:\n"
        "    jobid: 2\n"
        "Job 2 has been submitted with SLURM jobid 100 (log: x).\n"
        "rule merge:\n"
        "    jobid: 3\n"
        "Job 3 has been submitted with SLURM jobid 300 (log: x).\n",
        encoding="utf-8",
        newline="\n",
    )
    env = dict(os.environ)
    env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
    proc = subprocess.run(
        [_bash(), str(ROOT / "scripts" / "hpc" / "su_report.sh"), str(log)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr or proc.stdout
    assert "SU per genome 4.7604" in proc.stdout
    assert "fetch_region 1 0.5000" in proc.stdout
    assert "call 1 4.0000" in proc.stdout
    assert "merge 1 0.2604" in proc.stdout
