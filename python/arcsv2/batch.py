"""U4 batch runner: isolated samples, retries, status table, FAILED.txt."""

from __future__ import annotations

import csv
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from arcsv2 import __version__
from arcsv2.atomic_io import atomic_write_json, atomic_write_text
from arcsv2.input_checks import InputCheckError
from arcsv2.qc_report import build_qc_payload, write_qc_report
from arcsv2.run import run_sample

STATUS_COLUMNS = ("sample", "status", "message", "updated_utc")

TRANSIENT_HINTS = (
    "timed out",
    "timeout",
    "connection reset",
    "connection refused",
    "temporary failure",
    "resource temporarily unavailable",
    "503",
    "502",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_batch_manifest(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.lstrip("\ufeff").strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0].strip().lower() == "sample":
            continue
        if len(parts) < 2:
            raise InputCheckError(
                f"Batch manifest row needs at least sample and input: {line!r}",
                "Use columns: sample\\tinput\\t[reference]",
            )
        row = {
            "sample": parts[0].strip(),
            "input": parts[1].strip(),
            "reference": parts[2].strip() if len(parts) > 2 and parts[2].strip() else "",
        }
        if not row["sample"] or not row["input"]:
            raise InputCheckError("Empty sample or input in manifest row.", "Fill sample and input paths.")
        rows.append(row)
    if not rows:
        raise InputCheckError(f"No samples in batch manifest: {path}", "Add one row per sample.")
    return rows


def _is_transient(exc: BaseException) -> bool:
    msg = str(exc).lower()
    return any(h in msg for h in TRANSIENT_HINTS)


def _git_describe() -> str:
    try:
        proc = subprocess.run(
            ["git", "describe", "--tags", "--always", "--dirty"],
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
        if proc.returncode == 0:
            return proc.stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        pass
    return "unknown"


def _input_fingerprint(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {"path": str(path), "size": st.st_size, "mtime": int(st.st_mtime)}


def write_run_info(sample_dir: Path, *, sample: str, params: dict[str, Any], reads: Path) -> Path:
    info = {
        "arcsv2_version": __version__,
        "git": _git_describe(),
        "sample": sample,
        "parameters": params,
        "input": _input_fingerprint(reads),
        "time_utc": _utc_now(),
    }
    path = sample_dir / "run_info.json"
    atomic_write_json(path, info)
    return path


def write_failed_txt(sample_dir: Path, reason: str) -> Path:
    path = sample_dir / "FAILED.txt"
    atomic_write_text(path, reason.strip() + "\n")
    return path


def load_status_table(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    out: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["sample"]] = row
    return out


def write_status_table(path: Path, rows: dict[str, dict[str, str]]) -> None:
    lines = ["\t".join(STATUS_COLUMNS)]
    for sample in sorted(rows):
        r = rows[sample]
        lines.append("\t".join(r.get(c, "") for c in STATUS_COLUMNS))
    atomic_write_text(path, "\n".join(lines) + "\n")


def _set_status(
    table: dict[str, dict[str, str]],
    sample: str,
    status: str,
    message: str = "",
) -> None:
    table[sample] = {
        "sample": sample,
        "status": status,
        "message": message.replace("\t", " ").replace("\n", " ")[:500],
        "updated_utc": _utc_now(),
    }


def run_batch(
    manifest: Path,
    out_dir: Path,
    *,
    reference_default: str = "GRCh38",
    chroms: str | None = None,
    threads: str | int = "auto",
    seed: int = 204178949,
    engine: str = "rust",
    max_retries: int = 3,
) -> dict[str, Any]:
    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    entries = parse_batch_manifest(manifest)
    status_path = out_dir / "batch_status.tsv"
    status = load_status_table(status_path)
    for entry in entries:
        if entry["sample"] not in status:
            _set_status(status, entry["sample"], "pending")
    write_status_table(status_path, status)

    results: dict[str, Any] = {"done": [], "failed": [], "out_dir": str(out_dir)}
    params_base = {
        "chroms": chroms or "",
        "threads": str(threads),
        "seed": seed,
        "engine": engine,
    }

    for entry in entries:
        sample = entry["sample"]
        reads = Path(entry["input"]).expanduser().resolve()
        ref = entry["reference"] or reference_default
        sample_root = out_dir / sample
        sample_root.mkdir(parents=True, exist_ok=True)

        _set_status(status, sample, "running")
        write_status_table(status_path, status)

        last_err = ""
        for attempt in range(1, max(1, max_retries) + 1):
            try:
                meta = run_sample(
                    reads=reads,
                    reference=ref,
                    out_dir=sample_root,
                    chroms=chroms,
                    threads=threads,
                    seed=seed,
                    engine=engine,
                    nest_sample_dir=False,
                )
                run_info_path = write_run_info(
                    sample_root,
                    sample=sample,
                    params={**params_base, "reference": ref},
                    reads=reads,
                )
                payload = build_qc_payload(meta, run_info={"path": str(run_info_path)})
                write_qc_report(sample_root, payload)
                _set_status(status, sample, "done", f"merged={meta.get('merged_calls')}")
                write_status_table(status_path, status)
                results["done"].append(sample)
                last_err = ""
                break
            except InputCheckError as exc:
                last_err = str(exc)
                break
            except Exception as exc:
                last_err = str(exc)
                if attempt < max_retries and _is_transient(exc):
                    time.sleep(min(60, 2**attempt))
                    continue
                break

        if last_err:
            write_failed_txt(sample_root, last_err)
            _set_status(status, sample, "failed", last_err)
            write_status_table(status_path, status)
            results["failed"].append({"sample": sample, "error": last_err})

    atomic_write_json(
        out_dir / "batch_summary.json",
        {
            "time_utc": _utc_now(),
            "manifest": str(manifest.resolve()),
            "done": results["done"],
            "failed": results["failed"],
        },
    )
    return results
