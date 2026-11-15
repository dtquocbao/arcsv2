#!/usr/bin/env python3
"""Stage-ish counts for chr20:63-66Mb from existing full-chr exports + optional window call dirs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def in_window(row: dict, start: int, end: int) -> bool:
    lo = int(float(row["minbp"]))
    hi = int(float(row["maxbp"]))
    return not (hi < start or lo > end)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--start", type=int, default=63_000_000)
    parser.add_argument("--end", type=int, default=66_000_000)
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output path (default: stdout only)",
    )
    args = parser.parse_args()
    root = args.root
    o = [r for r in load(root / "export/orig/HG01884/arcsv_out.tab") if in_window(r, args.start, args.end)]
    v = [r for r in load(root / "export/results/HG01884/chr20/arcsv_out.tab") if in_window(r, args.start, args.end)]

    def stats(rows: list[dict]) -> dict:
        sr0 = sum(1 for r in rows if str(r.get("sr_support", "0")).split(",")[0] in ("0", "NA", ""))
        by = {}
        for r in rows:
            by[r.get("svtype", "?")] = by.get(r.get("svtype", "?"), 0) + 1
        return {"n": len(rows), "sr0": sr0, "by_type": by}

    text = "\n".join(
        [
            f"window chr20:{args.start}-{args.end}",
            f"original calls: {stats(o)}",
            f"arcsv2 calls: {stats(v)}",
            "",
            "Note: full stage counts (clips/DRPs/clusters) require windowed re-runs.",
        ]
    ) + "\n"
    print(text, end="")
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
