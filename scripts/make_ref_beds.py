"""Build gaps.bed and unreliable.bed for labeling.

    python scripts/make_ref_beds.py --fasta /data/ref/GRCh38....fa \\
        --centromeres /data/ref/centromeres.txt --out-dir /data/ref

``gaps.bed`` is every run of Ns in the FASTA. ``unreliable.bed`` is the
1KGP centromere/other-regions file, optionally plus UCSC segmental duplications
(``--segdups``).
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def n_runs(fasta: Path) -> list[tuple[str, int, int]]:
    """Half-open intervals of consecutive Ns (case-insensitive)."""
    intervals: list[tuple[str, int, int]] = []
    chrom = ""
    pos = 0
    run_start: int | None = None
    with fasta.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if chrom and run_start is not None:
                    intervals.append((chrom, run_start, pos))
                    run_start = None
                chrom = line[1:].split()[0]
                pos = 0
                continue
            for base in line.strip():
                if base in "Nn":
                    if run_start is None:
                        run_start = pos
                elif run_start is not None:
                    intervals.append((chrom, run_start, pos))
                    run_start = None
                pos += 1
        if chrom and run_start is not None:
            intervals.append((chrom, run_start, pos))
    return intervals


def parse_centromeres(path: Path) -> list[tuple[str, int, int]]:
    """1KGP ``20150713_location_of_centromeres_and_other_regions.txt`` rows."""
    rows: list[tuple[str, int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        chrom, start, end = parts[0], parts[1], parts[2]
        if not start.isdigit() or not end.isdigit():
            continue
        rows.append((chrom, int(start), int(end)))
    return rows


def parse_segdups(path: Path) -> list[tuple[str, int, int]]:
    """UCSC ``genomicSuperDups.txt.gz`` (chrom, chromStart, chromEnd)."""
    open_fn = gzip.open if path.name.endswith(".gz") else open
    rows: list[tuple[str, int, int]] = []
    with open_fn(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[arg-type]
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            # bin, chrom, chromStart, chromEnd, ...
            chrom, start, end = parts[1], parts[2], parts[3]
            if not start.isdigit() or not end.isdigit():
                continue
            rows.append((chrom, int(start), int(end)))
    return rows


def write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{chrom}\t{start}\t{end}" for chrom, start, end in intervals]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--centromeres", type=Path, default=None)
    parser.add_argument("--segdups", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    if not args.fasta.is_file():
        raise SystemExit(f"missing FASTA: {args.fasta}")
    gaps = n_runs(args.fasta)
    write_bed(args.out_dir / "gaps.bed", gaps)
    unreliable: list[tuple[str, int, int]] = []
    if args.centromeres and args.centromeres.is_file():
        unreliable.extend(parse_centromeres(args.centromeres))
    if args.segdups and args.segdups.is_file():
        unreliable.extend(parse_segdups(args.segdups))
    write_bed(args.out_dir / "unreliable.bed", unreliable)
    print(f"wrote {args.out_dir / 'gaps.bed'} ({len(gaps)} gaps)")
    print(f"wrote {args.out_dir / 'unreliable.bed'} ({len(unreliable)} intervals)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
