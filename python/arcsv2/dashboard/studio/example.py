"""Bundled example BAM runner for Studio (U6)."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

from arcsv2.atomic_io import atomic_write_text
from arcsv2.dashboard.studio import jobs, progress

_PACKAGE_DEMO = Path(__file__).resolve().parent.parent / "demo" / "example_bam"


def package_dir() -> Path:
    """Return the first directory that contains ``example.bam`` + ``ref.fa``."""
    candidates: list[Path] = []
    env = os.environ.get("ARCSV2_EXAMPLE_DIR", "").strip()
    if env:
        candidates.append(Path(env))
    candidates.append(_PACKAGE_DEMO)
    # Seeded by the Windows/macOS launcher into the mounted data folder.
    try:
        candidates.append(jobs.data_root() / "_bundled_example")
    except Exception:
        pass
    # Live checkout mounted at /work (developer --build).
    candidates.append(Path("/work/python/arcsv2/dashboard/demo/example_bam"))
    for path in candidates:
        if (path / "example.bam").is_file() and (path / "ref.fa").is_file():
            return path
    return _PACKAGE_DEMO


def _repo_root() -> Path | None:
    here = Path(__file__).resolve()
    for folder in [here, *here.parents]:
        script = folder / "scripts" / "simulate_haplotypes.py"
        if script.is_file():
            return folder
    work = Path("/work/scripts/simulate_haplotypes.py")
    if work.is_file():
        return Path("/work")
    return None


def ensure_example_assets(dest: Path) -> dict[str, Path]:
    """Copy (or generate) example BAM + mini reference into ``dest``."""
    dest.mkdir(parents=True, exist_ok=True)
    bam = dest / "example.bam"
    ref = dest / "ref.fa"
    gaps = dest / "gaps.bed"
    src = package_dir()
    if (src / "example.bam").is_file() and (src / "ref.fa").is_file():
        for name in ("example.bam", "example.bam.bai", "ref.fa", "ref.fa.fai", "gaps.bed", "truth.tsv"):
            s = src / name
            if s.is_file():
                shutil.copy2(s, dest / name)
        if not gaps.is_file():
            gaps.write_text("", encoding="utf-8")
        return {"bam": bam, "reference": ref, "gaps": gaps}

    root = _repo_root()
    if root is None:
        raise FileNotFoundError(
            "Bundled example BAM is missing from this install. "
            "From a git checkout run:  Start ARC-SV.bat --build   "
            "(or bash scripts/build_studio_example.sh and rebuild the image). "
            "The launcher also copies demo/example_bam into your data folder as _bundled_example."
        )
    script = root / "scripts" / "simulate_haplotypes.py"
    subprocess.check_call(
        [
            sys.executable,
            str(script),
            "--out",
            str(dest),
            "--length",
            "100000",
            "--events",
            "8",
            "--coverage",
            "20",
            "--seed",
            "42",
        ],
        cwd=str(root),
    )
    reads = dest / "reads.bam"
    if not reads.is_file():
        raise RuntimeError("Example BAM generation failed (samtools/pysam needed once).")
    shutil.copy2(reads, bam)
    bai = dest / "reads.bam.bai"
    if bai.is_file():
        shutil.copy2(bai, dest / "example.bam.bai")
    if not gaps.is_file():
        gaps.write_text("", encoding="utf-8")
    return {"bam": bam, "reference": ref, "gaps": gaps}


def start_example_run(root: Path | None = None) -> dict:
    """Start ``arcsv2 run`` on the bundled example BAM (about 2–5 minutes)."""
    root = root or jobs.data_root()
    work = root / "studio-example"
    assets = ensure_example_assets(work / "input")
    out_dir = work / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    sample = "sim"
    sample_dir = out_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    progress.write_progress(
        sample_dir,
        stage="pending",
        percent=0,
        started_utc=progress._utc(),
        chroms=["chrSim"],
    )
    status = out_dir / "batch_status.tsv"
    atomic_write_text(
        status,
        "sample\tstatus\tmessage\tupdated_utc\n"
        f"{sample}\trunning\texample run\t{progress._utc()}\n",
    )
    manifest = out_dir / "manifest.tsv"
    jobs.write_batch_manifest([(sample, assets["bam"])], manifest, reference=str(assets["reference"]))
    log_path = out_dir / "studio_example.log"
    log_handle = log_path.open("a", encoding="utf-8")
    log_handle.write(f"\n=== example start {progress._utc()} ===\n")
    cmd = [
        sys.executable,
        "-m",
        "arcsv2",
        "run",
        "--input",
        str(assets["bam"]),
        "--reference",
        str(assets["reference"]),
        "--out",
        str(out_dir),
        "--chroms",
        "chrSim",
        "--threads",
        "auto",
    ]
    proc = subprocess.Popen(
        cmd,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    job = {
        "id": f"example-{int(time.time())}",
        "pid": proc.pid,
        "manifest": str(manifest),
        "out_dir": str(out_dir),
        "chroms": "chrSim",
        "started_utc": progress._utc(),
        "log": str(log_path),
        "kind": "example",
    }
    registry = jobs.load_jobs(root)
    registry.append(job)
    jobs.save_jobs(registry, root)

    def _watch() -> None:
        code = proc.wait()
        tab = sample_dir / "chrSim" / "arcsv_out.tab"
        merged = sample_dir / "merged.tab"
        ok = code == 0 and (merged.is_file() or (tab.is_file() and tab.stat().st_size > 0))
        if ok:
            atomic_write_text(
                status,
                "sample\tstatus\tmessage\tupdated_utc\n"
                f"{sample}\tdone\texample complete\t{progress._utc()}\n",
            )
            progress.write_progress(sample_dir, stage="done", percent=100)
        else:
            atomic_write_text(
                status,
                "sample\tstatus\tmessage\tupdated_utc\n"
                f"{sample}\tfailed\tsee studio_example.log\t{progress._utc()}\n",
            )
            progress.write_progress(sample_dir, stage="failed", percent=0)

    threading.Thread(target=_watch, daemon=True).start()
    return job


def example_results_ready(root: Path | None = None) -> Path | None:
    root = root or jobs.data_root()
    candidates = [
        root / "studio-example" / "results" / "sim",
        root / "studio-example" / "results" / "example",
        root / "studio-example" / "results",
    ]
    for path in candidates:
        if (path / "merged.tab").is_file() or (path / "chrSim" / "arcsv_out.tab").is_file():
            return path
    return None
