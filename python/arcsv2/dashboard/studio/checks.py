"""Setup checks for ARC-SV Studio (U6)."""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str
    fix: str = ""


def check_docker() -> CheckResult:
    if Path("/.dockerenv").is_file():
        return CheckResult("Docker / container", True, "Running inside the ARC-SV2 container.")
    exe = shutil.which("docker")
    if not exe:
        return CheckResult(
            "Docker Desktop",
            False,
            "Docker was not found on PATH.",
            "Install Docker Desktop: https://docs.docker.com/desktop/",
        )
    try:
        proc = subprocess.run([exe, "info"], capture_output=True, text=True, timeout=15, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return CheckResult(
            "Docker Desktop",
            False,
            "Docker did not respond.",
            "Start Docker Desktop and wait until it shows “Running”.",
        )
    if proc.returncode != 0:
        return CheckResult(
            "Docker Desktop",
            False,
            "Docker is installed but not running.",
            "Start Docker Desktop: https://docs.docker.com/desktop/",
        )
    return CheckResult("Docker Desktop", True, "Docker daemon is running.")


def check_cpu() -> CheckResult:
    n = os.cpu_count() or 1
    ok = n >= 4
    return CheckResult(
        "CPU cores",
        ok,
        f"{n} logical cores detected.",
        "" if ok else "Whole-genome calling needs at least 4 cores for reasonable runtime.",
    )


def check_ram() -> CheckResult:
    gb: float | None = None
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    gb = int(line.split()[1]) / (1024 * 1024)
                    break
    except OSError:
        pass
    if gb is None and platform.system() == "Windows":
        try:
            import ctypes

            class MemStatus(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            stat = MemStatus()
            stat.dwLength = ctypes.sizeof(MemStatus)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
                gb = stat.ullTotalPhys / (1024**3)
        except (OSError, AttributeError):
            gb = None
    if gb is None:
        return CheckResult("RAM", True, "Could not measure RAM (check skipped).")
    ok = gb >= 16
    return CheckResult(
        "RAM",
        ok,
        f"{gb:.1f} GiB total.",
        "" if ok else "Labeling with two haplotype indexes needs 16+ GiB; use Quick test (one chromosome) on smaller laptops.",
    )


def check_disk(data_root: Path) -> CheckResult:
    try:
        usage = shutil.disk_usage(data_root if data_root.is_dir() else data_root.parent)
    except OSError as exc:
        return CheckResult("Free disk", False, str(exc), "Choose a data folder on a drive with more space.")
    free_gb = usage.free / (1024**3)
    ok = free_gb >= 50
    return CheckResult(
        "Free disk",
        ok,
        f"{free_gb:.1f} GiB free at {data_root}.",
        "" if ok else "Reserve at least 50 GiB for reference + one whole-genome CRAM run.",
    )


def check_reference(data_root: Path) -> CheckResult:
    from arcsv2.run import GRCH38_MAX_FA_BYTES, GRCH38_MIN_FA_BYTES, GRCH38_REF_NAME, default_cache_dir

    candidates = [
        data_root / "ref" / GRCH38_REF_NAME,
        default_cache_dir() / "GRCh38" / GRCH38_REF_NAME,
        Path.home() / ".cache" / "arcsv2" / "GRCh38" / GRCH38_REF_NAME,
    ]

    for path in candidates:
        if not path.is_file():
            continue
        size = path.stat().st_size
        if GRCH38_MIN_FA_BYTES <= size <= GRCH38_MAX_FA_BYTES:
            return CheckResult("Reference genome", True, f"Found {path.name}.")
        if size > GRCH38_MAX_FA_BYTES:
            return CheckResult(
                "Reference genome",
                False,
                f"Corrupt oversized cache at {path} ({size} bytes).",
                "Delete that file (or use Setup → Download reference) and try again.",
            )
    return CheckResult(
        "Reference genome",
        False,
        "GRCh38 reference not cached yet.",
        "It downloads automatically on the first `arcsv2 run` (needs network once).",
    )


def download_grch38(*, cache_dir: Path, on_progress=None) -> tuple[Path, Path]:
    """Download/cache GRCh38 under ``cache_dir`` with optional progress callback."""
    from arcsv2.run import (
        GRCH38_MAX_FA_BYTES,
        GRCH38_MIN_FA_BYTES,
        GRCH38_REF_NAME,
        GRCH38_REF_URL,
        _gaps_from_fasta,
        download_with_progress,
        ensure_fai,
    )

    root = Path(cache_dir) / "GRCh38"
    fasta = root / GRCH38_REF_NAME
    gaps = root / "gaps.bed"
    size = fasta.stat().st_size if fasta.is_file() else 0
    if size > GRCH38_MAX_FA_BYTES:
        fasta.unlink(missing_ok=True)
        fasta.with_suffix(fasta.suffix + ".fai").unlink(missing_ok=True)
        fasta.with_suffix(fasta.suffix + ".partial").unlink(missing_ok=True)
        size = 0
    if size < GRCH38_MIN_FA_BYTES:
        download_with_progress(
            GRCH38_REF_URL, fasta, on_progress=on_progress, max_bytes=GRCH38_MAX_FA_BYTES
        )
    ensure_fai(fasta)
    if not gaps.is_file():
        _gaps_from_fasta(fasta, gaps)
    return fasta, gaps


def run_all_checks(data_root: Path) -> list[CheckResult]:
    return [
        check_docker(),
        check_cpu(),
        check_ram(),
        check_disk(data_root),
        check_reference(data_root),
    ]
