"""Per-sample QC report JSON + HTML (U4)."""

from __future__ import annotations

import html
from pathlib import Path
from typing import Any

from arcsv2.atomic_io import atomic_write_json, atomic_write_text


def _flags(meta: dict[str, Any]) -> list[str]:
    flags: list[str] = []
    cov = meta.get("coverage_x")
    if cov is not None and float(cov) < 15.0:
        flags.append(f"low_coverage ({float(cov):.1f}×)")
    isize = meta.get("median_insert_size")
    if isize is not None and (isize < 200 or isize > 800):
        flags.append(f"unusual_insert_size ({isize} bp)")
    rlen = meta.get("median_read_length")
    if rlen is not None and (rlen < 100 or rlen > 200):
        flags.append(f"unusual_read_length ({rlen} bp)")
    return flags


def build_qc_payload(meta: dict[str, Any], *, run_info: dict[str, Any] | None = None) -> dict[str, Any]:
    stats = meta.get("stats") or {}
    return {
        "sample": meta.get("sample"),
        "runtime_sec": meta.get("runtime_sec"),
        "median_read_length": meta.get("median_read_length"),
        "median_insert_size": meta.get("median_insert_size"),
        "coverage_x": meta.get("coverage_x"),
        "calls_total": stats.get("total"),
        "calls_by_type": stats.get("by_type"),
        "calls_by_size": stats.get("size_bins"),
        "merged_calls": meta.get("merged_calls"),
        "flags": _flags(meta),
        "paths": meta.get("paths"),
        "run_info": run_info or {},
    }


def write_qc_report(out_dir: Path, payload: dict[str, Any]) -> tuple[Path, Path]:
    qc_dir = out_dir / "qc"
    json_path = qc_dir / "report.json"
    html_path = qc_dir / "report.html"
    atomic_write_json(json_path, payload)

    sample = html.escape(str(payload.get("sample") or ""))
    flags = payload.get("flags") or []
    flag_html = (
        "<ul>" + "".join(f"<li>{html.escape(str(f))}</li>" for f in flags) + "</ul>"
        if flags
        else "<p>None</p>"
    )
    by_type = payload.get("calls_by_type") or {}
    type_rows = "".join(
        f"<tr><td>{html.escape(k)}</td><td>{v}</td></tr>" for k, v in sorted(by_type.items())
    )
    body = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>QC {sample}</title>
<style>body{{font-family:system-ui;max-width:720px;margin:2rem auto;}}
table{{border-collapse:collapse}} td,th{{border:1px solid #ccc;padding:4px 8px}}</style>
</head><body>
<h1>ARC-SV2 QC — {sample}</h1>
<p>Runtime: {payload.get('runtime_sec', '—')} s</p>
<p>Median read length: {payload.get('median_read_length', '—')} bp</p>
<p>Median insert size: {payload.get('median_insert_size', '—')} bp</p>
<p>Merged calls (≥50 bp): {payload.get('merged_calls', '—')}</p>
<h2>Flags</h2>{flag_html}
<h2>Calls by type</h2>
<table><tr><th>Type</th><th>Count</th></tr>{type_rows or '<tr><td colspan=2>—</td></tr>'}</table>
</body></html>
"""
    atomic_write_text(html_path, body)
    return json_path, html_path
