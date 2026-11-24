#!/usr/bin/env python3
"""Write coverage.tsv for one ML-0 sample from its CRAM URL.

Does not download the CRAM and does not rewrite merged.tab. The sbatch array
runs this once per sample inside the SIF.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from arcsv2.coverage import write_coverage


def _rows(path: Path) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cells = line.split("\t")
        if header is None:
            header = cells
            continue
        rows.append({header[i]: cells[i] if i < len(cells) else "" for i in range(len(header))})
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=int, required=True, help="1-based sample index (SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--merged-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--gaps", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    rows = _rows(args.samples)
    if args.index < 1 or args.index > len(rows):
        print(f"index {args.index} outside 1..{len(rows)}", file=sys.stderr)
        return 2
    row = rows[args.index - 1]
    sample = row["sample"]
    calls = args.merged_dir / f"{sample}.merged.tab"
    if not calls.is_file():
        matches = sorted(args.merged_dir.glob(f"{sample}*.merged.tab"))
        if len(matches) != 1:
            print(f"missing merged table for {sample} under {args.merged_dir}", file=sys.stderr)
            return 1
        calls = matches[0]
    out = args.out_dir / f"{sample}.coverage.tsv"
    n = write_coverage(
        calls,
        row["cram_url"],
        out,
        reference=args.reference,
        gaps=args.gaps if args.gaps.is_file() else None,
    )
    print(f"wrote {n} rows for {sample} to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
