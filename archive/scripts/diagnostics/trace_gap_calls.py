#!/usr/bin/env python3
"""Trace gap-analysis call candidates against both callers' output tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_tab(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_int(value: str | None, default: int = 0) -> int:
    if value is None or value in ("", "NA"):
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def nearby(rows: list[dict], chrom: str, start: int, end: int, slop: int = 200) -> list[dict]:
    hits = []
    for r in rows:
        if r.get("chrom") != chrom:
            continue
        lo, hi = as_int(r.get("minbp")), as_int(r.get("maxbp"))
        if not (hi + slop < start or end + slop < lo):
            hits.append(r)
    return hits


def fmt(row: dict | None) -> str:
    if row is None:
        return "ABSENT"
    return (
        f"{row.get('svtype')} {row.get('minbp')}-{row.get('maxbp')} "
        f"gt={row.get('gt')} sr={row.get('sr_support')} pe={row.get('pe_support')} "
        f"score={row.get('score_vs_ref')} rearr={row.get('rearrangement')}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sample", default="HG01884")
    args = parser.parse_args()
    root = args.root
    candidates = root / "scratch" / "gap_analysis" / f"{args.sample}_call_trace_candidates.tsv"
    orig_calls = load_tab(root / "export" / "orig" / args.sample / "arcsv_out.tab")
    v2_calls = load_tab(root / "export" / "results" / args.sample / "chr20" / "arcsv_out.tab")
    out = root / "scratch" / "gap_analysis" / f"{args.sample}_call_traces.txt"
    lines = [f"sample={args.sample}", ""]
    with candidates.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = row["chrom"]
            start, end = as_int(row["minbp"]), as_int(row["maxbp"])
            o_hits = nearby(orig_calls, chrom, start, end)
            v_hits = nearby(v2_calls, chrom, start, end)
            lines.append(f"## {row['set']} {row['svtype']} {start}-{end} labeled_sr={row['sr_support']} pe={row['pe_support']}")
            lines.append(f"   candidate bp={row.get('bp')}")
            lines.append(f"   original nearby ({len(o_hits)}): " + "; ".join(fmt(h) for h in o_hits[:3]))
            lines.append(f"   arcsv2 nearby ({len(v_hits)}): " + "; ".join(fmt(h) for h in v_hits[:3] or [None]))
            lines.append("")
    out.write_text("\n".join(lines), encoding="utf-8")
    print(out.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
