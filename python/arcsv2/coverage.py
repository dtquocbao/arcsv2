"""Per-call depth from a BAM or a CRAM URL.

``samtools bedcov -Q 20`` is the only depth source. MAPQ >= 20 everywhere.
A span of 50 kb or less is one window. A longer span is 20 evenly spaced
1 kb windows. Flanks are the 1 kb left of minbp and the 1 kb right of maxbp.
The chromosome median is 200 random 1 kb windows, seed ``204178949`` plus
the chromosome name, skipping gaps. ``coverage_ratio`` is NA when the flank
depth is 0.
"""

from __future__ import annotations

import random
import shutil
import subprocess
import tempfile
from pathlib import Path

SEED = 204178949
MAPQ = 20
LONG_SPAN_BP = 50_000
WINDOW_BP = 1_000
FLANK_BP = 1_000
N_SAMPLE_WINDOWS = 20
N_CHROM_WINDOWS = 200

COLUMNS = (
    "id",
    "chrom",
    "cov_inside",
    "cov_flank",
    "cov_chrom_median",
    "coverage_ratio",
    "coverage_ratio_chrom",
)


def inside_windows(start: int, end: int) -> list[tuple[int, int]]:
    """Windows whose mean depth is ``cov_inside``."""
    span = end - start
    if span <= 0:
        return []
    if span <= LONG_SPAN_BP:
        return [(start, end)]
    width = WINDOW_BP
    if span <= width:
        return [(start, end)]
    room = span - width
    windows: list[tuple[int, int]] = []
    for i in range(N_SAMPLE_WINDOWS):
        pos = start + int(round(i * room / (N_SAMPLE_WINDOWS - 1)))
        windows.append((pos, min(end, pos + width)))
    return windows


def flank_windows(start: int, end: int, chrom_len: int) -> list[tuple[int, int]]:
    """1 kb left of ``start`` and 1 kb right of ``end``, clipped to the chromosome."""
    windows: list[tuple[int, int]] = []
    left = max(0, start - FLANK_BP)
    if start > left:
        windows.append((left, start))
    right = min(chrom_len, end + FLANK_BP)
    if right > end:
        windows.append((end, right))
    return windows


def _overlaps(start: int, end: int, gaps: list[tuple[int, int]]) -> bool:
    for gap_start, gap_end in gaps:
        if gap_start >= end:
            break
        if gap_end > start and gap_start < end:
            return True
    return False


def chrom_median_windows(
    chrom: str,
    chrom_len: int,
    gaps: list[tuple[int, int]],
    *,
    seed: int = SEED,
    n: int = N_CHROM_WINDOWS,
) -> list[tuple[int, int]]:
    """Fixed-seed 1 kb windows on ``chrom`` that do not overlap ``gaps``."""
    if chrom_len <= 0:
        return []
    if chrom_len <= WINDOW_BP:
        return [(0, chrom_len)]
    blocked = sorted(gaps)
    rng = random.Random(f"{seed}:{chrom}")
    hi = chrom_len - WINDOW_BP
    chosen: list[tuple[int, int]] = []
    tries = 0
    limit = n * 50
    while len(chosen) < n and tries < limit:
        tries += 1
        start = rng.randrange(0, hi + 1)
        end = start + WINDOW_BP
        if _overlaps(start, end, blocked):
            continue
        chosen.append((start, end))
    return chosen


def _fmt(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def _mean(lengths: list[int], totals: list[int]) -> float | None:
    bases = sum(lengths)
    if bases <= 0:
        return None
    return sum(totals) / bases


def _ratio(inside: float | None, denom: float | None) -> float | None:
    if inside is None or denom is None or denom == 0.0:
        return None
    return inside / denom


def read_gaps(path: Path | None, chrom: str) -> list[tuple[int, int]]:
    if path is None or not path.is_file():
        return []
    gaps: list[tuple[int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        parts = line.split()
        if len(parts) < 3 or parts[0] != chrom:
            continue
        if not parts[1].isdigit() or not parts[2].isdigit():
            continue
        start, end = int(parts[1]), int(parts[2])
        if end > start:
            gaps.append((start, end))
    gaps.sort()
    return gaps


def fai_length(reference: Path | None, chrom: str) -> int | None:
    if reference is None:
        return None
    fai = Path(str(reference) + ".fai")
    if not fai.is_file():
        return None
    for line in fai.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] == chrom and parts[1].isdigit():
            return int(parts[1])
    return None


def read_calls(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def plan_regions(
    calls: list[dict[str, str]],
    lengths: dict[str, int],
    gaps_for: dict[str, list[tuple[int, int]]],
    *,
    seed: int = SEED,
) -> tuple[list[tuple[str, int, int, str]], dict[str, tuple[str, str, int]]]:
    """BED regions and ``name -> (kind, key, length)``.

    ``key`` is the call id, or the chromosome name when ``kind`` is ``median``.
    """
    regions: list[tuple[str, int, int, str]] = []
    meta: dict[str, tuple[str, str, int]] = {}
    seen_median: set[str] = set()
    n = 0

    def add(chrom: str, start: int, end: int, kind: str, key: str) -> None:
        nonlocal n
        if end <= start:
            return
        name = f"w{n}"
        n += 1
        regions.append((chrom, start, end, name))
        meta[name] = (kind, key, end - start)

    for row in calls:
        chrom = row["chrom"]
        start = int(row["minbp"])
        end = int(row["maxbp"])
        call_id = row["id"]
        length = lengths.get(chrom)
        if length is None:
            raise ValueError(f"no length for {chrom}")
        for win_start, win_end in inside_windows(start, end):
            add(chrom, win_start, win_end, "inside", call_id)
        for win_start, win_end in flank_windows(start, end, length):
            add(chrom, win_start, win_end, "flank", call_id)
        if chrom not in seen_median:
            seen_median.add(chrom)
            for win_start, win_end in chrom_median_windows(
                chrom, length, gaps_for.get(chrom, []), seed=seed
            ):
                add(chrom, win_start, win_end, "median", chrom)
    return regions, meta


def rows_from_totals(
    calls: list[dict[str, str]],
    meta: dict[str, tuple[str, str, int]],
    totals: dict[str, int],
) -> list[dict[str, str]]:
    """Join bedcov sums back onto calls. A zero flank depth yields NA."""
    grouped: dict[str, dict[str, list[tuple[int, int]]]] = {}
    median: dict[str, list[tuple[int, int]]] = {}
    for name, (kind, key, length) in meta.items():
        pair = (length, totals.get(name, 0))
        if kind == "median":
            median.setdefault(key, []).append(pair)
        else:
            grouped.setdefault(key, {}).setdefault(kind, []).append(pair)

    median_mean = {
        chrom: _mean([length for length, _total in pairs], [total for _length, total in pairs])
        for chrom, pairs in median.items()
    }
    out: list[dict[str, str]] = []
    for row in calls:
        call_id = row["id"]
        chrom = row["chrom"]
        parts = grouped.get(call_id, {})
        inside_pairs = parts.get("inside", [])
        flank_pairs = parts.get("flank", [])
        cov_inside = _mean(
            [length for length, _total in inside_pairs],
            [total for _length, total in inside_pairs],
        )
        cov_flank = _mean(
            [length for length, _total in flank_pairs],
            [total for _length, total in flank_pairs],
        )
        cov_chrom = median_mean.get(chrom)
        out.append(
            {
                "id": call_id,
                "chrom": chrom,
                "cov_inside": _fmt(cov_inside),
                "cov_flank": _fmt(cov_flank),
                "cov_chrom_median": _fmt(cov_chrom),
                "coverage_ratio": _fmt(_ratio(cov_inside, cov_flank)),
                "coverage_ratio_chrom": _fmt(_ratio(cov_inside, cov_chrom)),
            }
        )
    return out


def write_coverage_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(COLUMNS)]
    for row in rows:
        lines.append("\t".join(row[col] for col in COLUMNS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def bedcov(
    reads: str,
    regions: list[tuple[str, int, int, str]],
    *,
    reference: Path | None = None,
    mapq: int = MAPQ,
) -> dict[str, int]:
    """Sum of per-base depths for each named region. MAPQ threshold is ``mapq``."""
    if not regions:
        return {}
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools is not on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        bed = Path(tmp) / "windows.bed"
        bed.write_text(
            "".join(f"{chrom}\t{start}\t{end}\t{name}\n" for chrom, start, end, name in regions),
            encoding="utf-8",
        )
        cmd = ["samtools", "bedcov", "-Q", str(mapq)]
        if reference is not None:
            cmd.extend(["--reference", str(reference)])
        cmd.extend([str(bed), reads])
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "samtools bedcov failed").strip()
        raise RuntimeError(detail)
    totals: dict[str, int] = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        totals[parts[3]] = int(parts[-1])
    return totals


def idxstats_length(reads: str, chrom: str, reference: Path | None = None) -> int | None:
    cmd = ["samtools", "idxstats"]
    if reference is not None:
        cmd.extend(["--reference", str(reference)])
    cmd.append(reads)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        return None
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] == chrom and parts[1].isdigit():
            return int(parts[1])
    return None


def chrom_lengths(
    calls: list[dict[str, str]],
    reads: str,
    reference: Path | None,
) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for row in calls:
        chrom = row["chrom"]
        if chrom in lengths:
            continue
        length = fai_length(reference, chrom)
        if length is None:
            length = idxstats_length(reads, chrom, reference)
        if length is None:
            raise RuntimeError(f"no length for {chrom} (need a reference .fai or an indexed BAM)")
        lengths[chrom] = length
    return lengths


def coverage_rows(
    calls: list[dict[str, str]],
    reads: str,
    *,
    reference: Path | None = None,
    gaps: Path | None = None,
    seed: int = SEED,
    totals: dict[str, int] | None = None,
    lengths: dict[str, int] | None = None,
) -> list[dict[str, str]]:
    """Depth columns for ``calls``.

    Pass ``totals`` to skip samtools (tests). Otherwise ``reads`` is a BAM
    path or a CRAM URL; a CRAM needs ``reference``.
    """
    if not calls:
        return []
    if lengths is None:
        lengths = chrom_lengths(calls, reads, reference)
    gaps_for = {chrom: read_gaps(gaps, chrom) for chrom in lengths}
    regions, meta = plan_regions(calls, lengths, gaps_for, seed=seed)
    if totals is None:
        totals = bedcov(reads, regions, reference=reference)
    return rows_from_totals(calls, meta, totals)


def write_coverage(
    calls_path: Path,
    reads: str,
    out: Path,
    *,
    reference: Path | None = None,
    gaps: Path | None = None,
    seed: int = SEED,
) -> int:
    calls = read_calls(calls_path)
    rows = coverage_rows(calls, reads, reference=reference, gaps=gaps, seed=seed)
    write_coverage_tsv(out, rows)
    return len(rows)
