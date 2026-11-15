"""BED intervals. Coordinates stay 0-based and half-open."""

from __future__ import annotations

from pathlib import Path


def read_bed_by_chrom(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Read ``chrom start end`` rows, grouped and sorted per chromosome."""
    grouped: dict[str, list[tuple[int, int]]] = {}
    text = Path(path).read_text(encoding="utf-8")
    for line in text.splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        name, start_s, end_s = parts[0], parts[1], parts[2]
        start, end = int(start_s), int(end_s)
        if end > start:
            grouped.setdefault(name, []).append((start, end))
    for intervals in grouped.values():
        intervals.sort()
    return grouped


def read_bed(path: str | Path, chrom: str | None = None) -> list[tuple[int, int]]:
    """Read ``chrom start end`` rows. Header and track lines are skipped.

    With ``chrom`` set, only that chromosome is returned. Without it, intervals
    from every chromosome are concatenated.
    """
    grouped = read_bed_by_chrom(path)
    if chrom:
        return list(grouped.get(chrom, []))
    intervals = [iv for rows in grouped.values() for iv in rows]
    intervals.sort()
    return intervals
