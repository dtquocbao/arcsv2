"""Download and build GRCh38 paper merge-filter annotation BEDs.

Paper Methods p. 27 / Key Resources:

- simple repeats + low complexity + satellite + STR (UCSC tracks, union)
- segmental duplications (UCSC genomicSuperDups)
- SV blacklist (10x Genomics)
- reference gaps (already built by ``make_ref_beds.py`` from N-runs)

Writes under ``--dest`` (default ``/data/ref``):

- ``simple_repeats.bed``  — union of simpleRepeat, rmsk Low_complexity/Satellite/Simple_repeat, microsat
- ``segdups.bed``
- ``blacklist.bed``

Usage (inside the container):

    python scripts/fetch_filter_beds.py --dest /data/ref
    python scripts/fetch_filter_beds.py --dest /data/ref --dry-run

Paths are meant for ``workflows/config/label.yaml`` keys ``simple_repeats``, ``segdups``,
``blacklist`` (plus existing ``gaps``).
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

USER_AGENT = "arcsv2-fetch-filter-beds/1"

# Documented in docs/reviews/29092026-round13.md and paper Key Resources.
SIMPLE_REPEAT_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"
)
RMSK_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
MICROSAT_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/microsat.txt.gz"
SEGDUPS_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"
)
BLACKLIST_URL = "https://cf.10xgenomics.com/supp/genome/GRCh38/svblacklist.bed"

# RepeatMasker classes that match the paper's >85% simple-repeat/low-complexity/satellite/STR cut.
RMSK_REPEAT_CLASSES = frozenset({"Low_complexity", "Satellite", "Simple_repeat"})


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 0:
        print(f"keep {dest}")
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"download {url} -> {dest}")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=600) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(dest)
    return dest


def _open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _parse_ucsc_chrom_start_end(
    path: Path,
    *,
    chrom_i: int,
    start_i: int,
    end_i: int,
    class_i: int | None = None,
    class_allow: frozenset[str] | None = None,
) -> list[tuple[str, int, int]]:
    rows: list[tuple[str, int, int]] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            need = max(chrom_i, start_i, end_i, class_i or 0) + 1
            if len(parts) < need:
                continue
            if class_allow is not None and class_i is not None:
                if parts[class_i] not in class_allow:
                    continue
            chrom, start_s, end_s = parts[chrom_i], parts[start_i], parts[end_i]
            if not start_s.isdigit() or not end_s.isdigit():
                continue
            start, end = int(start_s), int(end_s)
            if end > start:
                rows.append((chrom, start, end))
    return rows


def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge overlapping/adjacent half-open intervals per chromosome."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))
    out: list[tuple[str, int, int]] = []
    for chrom in sorted(by_chrom):
        ivs = sorted(by_chrom[chrom])
        cur_a, cur_b = ivs[0]
        for a, b in ivs[1:]:
            if a <= cur_b:
                cur_b = max(cur_b, b)
            else:
                out.append((chrom, cur_a, cur_b))
                cur_a, cur_b = a, b
        out.append((chrom, cur_a, cur_b))
    return out


def write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{chrom}\t{start}\t{end}" for chrom, start, end in intervals]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    print(f"wrote {path} ({len(intervals)} intervals)")


def build_filter_beds(dest: Path, *, dry_run: bool = False) -> dict[str, Path]:
    """Download sources and write paper merge-filter BEDs under ``dest``."""
    dest = dest.resolve()
    raw = dest / "raw"
    plan = {
        "simple_repeats": dest / "simple_repeats.bed",
        "segdups": dest / "segdups.bed",
        "blacklist": dest / "blacklist.bed",
    }
    sources = {
        "simpleRepeat": (SIMPLE_REPEAT_URL, raw / "simpleRepeat.txt.gz"),
        "rmsk": (RMSK_URL, raw / "rmsk.txt.gz"),
        "microsat": (MICROSAT_URL, raw / "microsat.txt.gz"),
        "segdups": (SEGDUPS_URL, raw / "genomicSuperDups.txt.gz"),
        "blacklist": (BLACKLIST_URL, raw / "svblacklist.bed"),
    }
    print("paper merge-filter annotation URLs:")
    for name, (url, path) in sources.items():
        print(f"  {name}: {url} -> {path}")
    for key, path in plan.items():
        print(f"  out {key}: {path}")
    print(f"  gaps (existing): {dest / 'gaps.bed'}  # from make_ref_beds.py / FASTA N-runs")

    if dry_run:
        print("dry-run: no downloads")
        return plan

    simple_path = _download(*sources["simpleRepeat"])
    rmsk_path = _download(*sources["rmsk"])
    microsat_path = _download(*sources["microsat"])
    segdups_src = _download(*sources["segdups"])
    blacklist_src = _download(*sources["blacklist"])

    # UCSC simpleRepeat / microsat / genomicSuperDups: bin, chrom, chromStart, chromEnd, ...
    repeats: list[tuple[str, int, int]] = []
    repeats.extend(
        _parse_ucsc_chrom_start_end(simple_path, chrom_i=1, start_i=2, end_i=3)
    )
    repeats.extend(
        _parse_ucsc_chrom_start_end(microsat_path, chrom_i=1, start_i=2, end_i=3)
    )
    # rmsk: genoName=5, genoStart=6, genoEnd=7, repClass=11
    repeats.extend(
        _parse_ucsc_chrom_start_end(
            rmsk_path,
            chrom_i=5,
            start_i=6,
            end_i=7,
            class_i=11,
            class_allow=RMSK_REPEAT_CLASSES,
        )
    )
    write_bed(plan["simple_repeats"], merge_intervals(repeats))

    segdups = _parse_ucsc_chrom_start_end(segdups_src, chrom_i=1, start_i=2, end_i=3)
    write_bed(plan["segdups"], merge_intervals(segdups))

    # Blacklist is already BED; normalize and merge.
    black_rows: list[tuple[str, int, int]] = []
    with blacklist_src.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3 or not parts[1].isdigit() or not parts[2].isdigit():
                continue
            start, end = int(parts[1]), int(parts[2])
            if end > start:
                black_rows.append((parts[0], start, end))
    write_bed(plan["blacklist"], merge_intervals(black_rows))

    gaps = dest / "gaps.bed"
    if not gaps.is_file():
        print(
            f"note: {gaps} missing — build it with scripts/make_ref_beds.py "
            "(FASTA N-runs) before paper merge-filter",
            file=sys.stderr,
        )
    return plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path("/data/ref"),
        help="Directory for BED outputs and raw/ downloads (default: /data/ref)",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    build_filter_beds(args.dest, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
