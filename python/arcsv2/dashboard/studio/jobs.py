"""Background batch jobs for Studio (U6). Uses U4 batch status files."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from arcsv2.atomic_io import atomic_write_text


def data_root() -> Path:
    raw = os.environ.get("ARCSV2_RESULTS", "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    return Path("/data" if Path("/data").is_dir() else Path.cwd())


def jobs_registry(root: Path | None = None) -> Path:
    root = root or data_root()
    return root / ".arcsv2" / "studio_jobs.json"


def load_jobs(root: Path | None = None) -> list[dict]:
    path = jobs_registry(root)
    if not path.is_file():
        return []
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []


def save_jobs(jobs: list[dict], root: Path | None = None) -> None:
    path = jobs_registry(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(path, json.dumps(jobs, indent=2) + "\n")


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def list_read_files(root: Path) -> list[Path]:
    out: list[Path] = []
    for pattern in ("*.cram", "*.bam", "*.CRAM", "*.BAM"):
        out.extend(sorted(root.rglob(pattern)))
    return out


def write_batch_manifest(rows: list[tuple[str, Path]], path: Path, reference: str = "GRCh38") -> None:
    lines = ["sample\tinput\treference"]
    for sample, reads in rows:
        lines.append(f"{sample}\t{reads}\t{reference}")
    atomic_write_text(path, "\n".join(lines) + "\n")


def start_batch_job(
    *,
    manifest: Path,
    out_dir: Path,
    chroms: str,
    threads: str = "auto",
    reference: str = "GRCh38",
) -> dict:
    """Start ``arcsv2 batch`` detached (survives browser close)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "studio_batch.log"
    log_handle = log_path.open("a", encoding="utf-8")
    log_handle.write(f"\n=== start {_utc()} ===\n")
    cmd = [
        sys.executable,
        "-m",
        "arcsv2",
        "batch",
        "--manifest",
        str(manifest),
        "--out",
        str(out_dir),
        "--reference",
        reference,
        "--threads",
        threads,
    ]
    if chroms:
        cmd.extend(["--chroms", chroms])
    proc = subprocess.Popen(
        cmd,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    job = {
        "id": f"batch-{int(time.time())}",
        "pid": proc.pid,
        "manifest": str(manifest),
        "out_dir": str(out_dir),
        "chroms": chroms,
        "started_utc": _utc(),
        "log": str(log_path),
    }
    jobs = load_jobs(out_dir.parent if out_dir.name != "data" else out_dir)
    jobs.append(job)
    save_jobs(jobs, data_root())
    return job


def install_example_run(root: Path | None = None) -> Path:
    """Copy bundled simulation demo into the data folder as ``studio-example``."""
    from shutil import copytree

    root = root or data_root()
    dest = root / "studio-example"
    if (dest / "simulated_2Mb_30x" / "arcsv_out.tab").is_file() or (dest / "merged.tab").is_file():
        return dest
    src = Path(__file__).resolve().parent.parent / "demo" / "simulated_2Mb_30x"
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / "simulated_2Mb_30x"
    if src.is_dir():
        copytree(src, target, dirs_exist_ok=True)
    meta = {"sample": "studio-example", "coverage": 30, "note": "Bundled 2 Mb simulation demo"}
    atomic_write_text(dest / "run.json", json.dumps(meta, indent=2) + "\n")
    return dest


def estimate_run(chroms: str, n_samples: int, input_bytes: int) -> tuple[str, str]:
    from arcsv2.input_checks import estimate_run_disk_bytes

    n_chrom = len([c for c in chroms.replace(",", " ").split() if c.strip()]) if chroms.strip() else 23

    class _FakePath:
        def stat(self):
            class _S:
                st_size = max(int(input_bytes), 500_000_000)

            return _S()

    need = estimate_run_disk_bytes(_FakePath(), n_chroms=n_chrom, need_ref_download=True) * max(1, n_samples)
    disk_gb = need / (1024**3)
    if chroms.strip():
        time_est = "roughly 15–45 minutes per sample (one chromosome)"
    else:
        time_est = "roughly 4–12 hours per sample (whole genome, depends on CPU)"
    return time_est, f"about {disk_gb:.0f} GiB for {n_samples} sample(s)"
