"""Per-sample progress, cancel, and resume helpers for Studio (U6)."""

from __future__ import annotations

import json
import os
import signal
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from arcsv2.atomic_io import atomic_write_json, atomic_write_text
from arcsv2.dashboard.studio import jobs
from arcsv2.dashboard import data as dash_data


STAGES = ("pending", "extract", "call", "merge", "qc", "done", "failed", "cancelled")


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def progress_path(sample_dir: Path) -> Path:
    return sample_dir / "progress.json"


def write_progress(sample_dir: Path, **fields: Any) -> None:
    path = progress_path(sample_dir)
    payload = {"updated_utc": _utc()}
    if path.is_file():
        try:
            payload.update(json.loads(path.read_text(encoding="utf-8")))
        except json.JSONDecodeError:
            pass
    payload.update(fields)
    atomic_write_json(path, payload)


def read_progress(sample_dir: Path) -> dict[str, Any]:
    path = progress_path(sample_dir)
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def infer_progress(sample_dir: Path, chroms: list[str] | None = None) -> dict[str, Any]:
    """Derive stage / percent from on-disk artifacts (resumable)."""
    stored = read_progress(sample_dir)
    status_file = sample_dir.parent / "batch_status.tsv"
    status = "pending"
    message = ""
    if status_file.is_file():
        table = dash_data.load_batch_status(status_file)
        name = sample_dir.name
        if not table.empty and "sample" in table.columns:
            hit = table[table["sample"] == name]
            if not hit.empty:
                status = str(hit.iloc[0].get("status") or "pending")
                message = str(hit.iloc[0].get("message") or "")

    chrom_list = chroms or stored.get("chroms") or []
    if isinstance(chrom_list, str):
        chrom_list = [c for c in chrom_list.replace(",", " ").split() if c]
    n = max(1, len(chrom_list) if chrom_list else 1)
    done_chroms = 0
    for chrom in chrom_list or ["chr20"]:
        tab = sample_dir / str(chrom) / "arcsv_out.tab"
        if tab.is_file() and tab.stat().st_size > 0:
            done_chroms += 1
    merged = sample_dir / "merged.tab"
    has_merged = merged.is_file() and merged.stat().st_size > 0
    failed = (sample_dir / "FAILED.txt").is_file()

    if status == "cancelled" or stored.get("stage") == "cancelled":
        stage = "cancelled"
        pct = float(stored.get("percent") or 0)
    elif failed or status == "failed":
        stage = "failed"
        pct = float(stored.get("percent") or (100.0 * done_chroms / n))
    elif status == "done" or has_merged:
        stage = "done"
        pct = 100.0
    elif done_chroms == 0 and (sample_dir / "_work").is_dir():
        stage = "extract"
        pct = 5.0
    elif done_chroms < n:
        stage = "call"
        pct = 10.0 + 70.0 * (done_chroms / n)
    elif not has_merged:
        stage = "merge"
        pct = 90.0
    else:
        stage = "qc"
        pct = 95.0

    started = stored.get("started_utc") or ""
    eta = stored.get("eta_utc") or ""
    if stage not in ("done", "failed", "cancelled") and started and pct > 5:
        try:
            t0 = datetime.strptime(started, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
            elapsed = max(1.0, (datetime.now(timezone.utc) - t0).total_seconds())
            if pct < 100:
                remaining = elapsed * (100.0 - pct) / pct
                eta_dt = datetime.now(timezone.utc).timestamp() + remaining
                eta = datetime.fromtimestamp(eta_dt, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        except ValueError:
            pass

    return {
        "sample": sample_dir.name,
        "stage": stage,
        "percent": round(pct, 1),
        "done_chroms": done_chroms,
        "n_chroms": n,
        "eta_utc": eta,
        "started_utc": started,
        "message": message,
        "status": status,
    }


def cancel_job(job_id: str, root: Path | None = None) -> str:
    """Signal a background batch PID to stop; mark samples cancelled."""
    root = root or jobs.data_root()
    registry = jobs.load_jobs(root)
    target = next((j for j in registry if j.get("id") == job_id), None)
    if target is None:
        return f"No job {job_id} in the registry."
    pid = int(target.get("pid") or 0)
    out_dir = Path(target["out_dir"])
    if pid > 0:
        try:
            os.kill(pid, signal.SIGTERM)
        except (OSError, ProcessLookupError):
            pass
    target["cancelled"] = True
    target["cancelled_utc"] = _utc()
    jobs.save_jobs(registry, root)
    status_path = out_dir / "batch_status.tsv"
    if status_path.is_file():
        table = dash_data.load_batch_status(status_path)
        for _, row in table.iterrows():
            if str(row.get("status")) in ("running", "pending"):
                sample = str(row["sample"])
                write_progress(out_dir / sample, stage="cancelled", percent=0)
                failed = out_dir / sample / "FAILED.txt"
                atomic_write_text(failed, "Cancelled by user from ARC-SV Studio.\n")
    return f"Sent cancel to job {job_id} (pid {pid})."


def resume_batch(out_dir: Path, *, chroms: str = "", threads: str = "auto") -> dict:
    """Re-launch batch on the same folder (skips finished samples / chrom tabs)."""
    manifest = out_dir / "manifest.tsv"
    if not manifest.is_file():
        raise FileNotFoundError(f"No manifest.tsv in {out_dir}")
    # Clear FAILED / cancelled markers for pending samples so they retry
    status = dash_data.load_batch_status(out_dir / "batch_status.tsv")
    for _, row in status.iterrows():
        if str(row.get("status")) in ("failed", "cancelled"):
            sample = str(row["sample"])
            (out_dir / sample / "FAILED.txt").unlink(missing_ok=True)
            prog = progress_path(out_dir / sample)
            if prog.is_file():
                prog.unlink(missing_ok=True)
    return jobs.start_batch_job(
        manifest=manifest,
        out_dir=out_dir,
        chroms=chroms,
        threads=threads,
    )


def list_analysis_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for status_file in dash_data.find_batch_status_files(root):
        batch_dir = status_file.parent
        table = dash_data.load_batch_status(status_file)
        chroms_hint: list[str] = []
        job_match = next(
            (j for j in jobs.load_jobs(root) if Path(j.get("out_dir", "")) == batch_dir),
            None,
        )
        if job_match and job_match.get("chroms"):
            chroms_hint = [c for c in str(job_match["chroms"]).replace(",", " ").split() if c]
        for _, row in table.iterrows():
            sample = str(row["sample"])
            sample_dir = batch_dir / sample
            prog = infer_progress(sample_dir, chroms_hint or None)
            prog["batch_dir"] = str(batch_dir)
            prog["job_id"] = job_match["id"] if job_match else ""
            rows.append(prog)
    return rows
