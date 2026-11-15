"""U3: validate inputs before heavy work; plain-language errors."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

# Paper cohort is ~30× short WGS; reject obvious mismatches early.
MIN_COVERAGE_X = 10.0
MIN_READ_LENGTH = 50
MAX_READ_LENGTH = 600
MAX_UNMAPPED_FRACTION = 0.85

GRCH38_REF_BYTES = 3_200_000_000


class InputCheckError(Exception):
    """User-facing validation error (what + how to fix)."""

    def __init__(self, message: str, fix: str = "") -> None:
        self.message = message
        self.fix = fix
        text = message if not fix else f"{message} Fix: {fix}"
        super().__init__(text)


def index_path(reads: Path) -> Path | None:
    if reads.suffix == ".cram":
        return Path(str(reads) + ".crai")
    if reads.suffix == ".bam":
        return Path(str(reads) + ".bai")
    return None


def check_reads_readable(reads: Path) -> None:
    path = reads.expanduser()
    if not path.exists():
        raise InputCheckError(
            f"Input file does not exist: {path}",
            "Pass a valid CRAM, BAM, or SAM path to --input.",
        )
    if not path.is_file():
        raise InputCheckError(f"Input is not a file: {path}", "Use a CRAM/BAM/SAM file path.")
    try:
        with path.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise InputCheckError(
            f"Cannot read input file: {path} ({exc})",
            "Check permissions and that the path is not locked by another program.",
        ) from exc
    if path.suffix.lower() not in (".cram", ".bam", ".sam"):
        raise InputCheckError(
            f"Unrecognized alignment format (expected .cram, .bam, or .sam): {path.name}",
            "Rename the file or convert with samtools.",
        )


def check_out_writable(out_dir: Path) -> None:
    root = out_dir.expanduser()
    try:
        root.mkdir(parents=True, exist_ok=True)
        probe = root / ".arcsv2_write_probe"
        probe.write_text("", encoding="utf-8")
        probe.unlink(missing_ok=True)
    except OSError as exc:
        raise InputCheckError(
            f"Output directory is not writable: {root} ({exc})",
            "Choose another --out path or fix directory permissions.",
        ) from exc


def check_reference_files(reference: Path, gaps: Path | None = None) -> None:
    if not reference.is_file():
        raise InputCheckError(
            f"Reference FASTA not found: {reference}",
            "Pass -R/--reference to a GRCh38 FASTA or use --reference GRCh38 with network/cache.",
        )
    try:
        with reference.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise InputCheckError(f"Cannot read reference FASTA: {reference} ({exc})", "Check permissions.") from exc
    if gaps is not None and gaps != Path("") and gaps.is_file():
        try:
            with gaps.open("rb") as handle:
                handle.read(1)
        except OSError as exc:
            raise InputCheckError(f"Cannot read gaps BED: {gaps} ({exc})", "Regenerate gaps or fix the path.") from exc


def check_index_present(reads: Path, *, auto_build: bool) -> None:
    if reads.suffix.lower() == ".sam":
        raise InputCheckError(
            "SAM input is not supported for whole-genome run without conversion.",
            "samtools view -bS in.sam | samtools sort -o sorted.bam && samtools index sorted.bam",
        )
    idx = index_path(reads)
    if idx is None:
        return
    if idx.is_file():
        return
    if auto_build:
        if not shutil.which("samtools"):
            raise InputCheckError(
                f"Index missing ({idx.name}) and samtools is not on PATH.",
                f"Run: samtools index {reads}",
            )
        return
    raise InputCheckError(
        f"Alignment index missing: {idx}",
        f"Run: samtools index {reads}",
    )


def _samtools_idxstats(reads: Path) -> list[tuple[str, int, int]]:
    exe = shutil.which("samtools")
    if not exe:
        return []
    proc = subprocess.run(
        [exe, "idxstats", str(reads)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return []
    rows: list[tuple[str, int, int]] = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        chrom, length, mapped = parts[0], parts[1], parts[2]
        if not length.isdigit() or not mapped.isdigit():
            continue
        if chrom == "*":
            continue
        rows.append((chrom, int(length), int(mapped)))
    return rows


def check_sorted_and_mapped(
    reads: Path,
    reference: Path | None,
    *,
    idxstats_rows: list[tuple[str, int, int]] | None = None,
) -> None:
    try:
        import pysam
    except ImportError:
        return
    kwargs: dict = {}
    if str(reads).lower().endswith(".cram") and reference and reference.is_file():
        kwargs["reference_filename"] = str(reference)
    try:
        with pysam.AlignmentFile(str(reads), "rb", **kwargs) as handle:
            hd = handle.header.to_dict().get("HD", {})
            if isinstance(hd, list):
                so = (hd[0].get("SO", "") if hd else "") or ""
            elif isinstance(hd, dict):
                so = str(hd.get("SO", "") or "")
            else:
                so = ""
            if so and so != "coordinate":
                if so == "unsorted":
                    raise InputCheckError(
                        "BAM/CRAM is not coordinate-sorted (SO:unsorted).",
                        f"samtools sort -o sorted.bam {reads} && samtools index sorted.bam",
                    )
                raise InputCheckError(
                    f"BAM/CRAM sort order {so!r} is not supported (need coordinate-sorted data).",
                    f"samtools sort -o sorted.bam {reads} && samtools index sorted.bam",
                )
    except InputCheckError:
        raise
    except (OSError, ValueError) as exc:
        msg = str(exc).lower()
        if "could not retrieve alignment file" in msg or "missing reference" in msg:
            raise InputCheckError(
                "CRAM cannot be read with this reference (encoding mismatch).",
                "Use the same GRCh38 FASTA the CRAM was aligned to (-R / --reference).",
            ) from exc
        raise InputCheckError(
            f"Cannot open alignment file: {exc}",
            "Verify the file is intact and indexed.",
        ) from exc

    rows = idxstats_rows if idxstats_rows is not None else _samtools_idxstats(reads)
    check_mapping_stats(rows)


def check_mapping_stats(rows: list[tuple[str, int, int]]) -> None:
    if not rows:
        return
    total_mapped = sum(m for c, _l, m in rows if c != "*")
    if total_mapped == 0:
        raise InputCheckError(
            "No mapped reads found (idxstats shows zero mapped alignments).",
            "Use a BAM/CRAM aligned to the reference, not unmapped reads.",
        )
    unmapped = sum(m for c, _l, m in rows if c == "*")
    denom = total_mapped + unmapped
    if denom > 0 and unmapped / denom > MAX_UNMAPPED_FRACTION:
        raise InputCheckError(
            "Too many unmapped reads (>85% unmapped).",
            "Re-align to the reference FASTA or pass a properly aligned CRAM/BAM.",
        )


def check_read_length(median_length: int | None) -> None:
    if median_length is None:
        return
    if median_length < MIN_READ_LENGTH:
        raise InputCheckError(
            f"Median read length {median_length} bp is below the supported minimum ({MIN_READ_LENGTH} bp).",
            "Use standard Illumina WGS (≈150 bp) or contact support for exotic platforms.",
        )
    if median_length > MAX_READ_LENGTH:
        raise InputCheckError(
            f"Median read length {median_length} bp exceeds the supported maximum ({MAX_READ_LENGTH} bp).",
            "Subsample or trim reads, or use arcsv2 call on a converted BAM after QC.",
        )


def estimate_coverage_x(
    reads: Path,
    chrom: str,
    chrom_length: int,
    median_read_length: int | None,
    *,
    idxstats_rows: list[tuple[str, int, int]] | None = None,
) -> float | None:
    if chrom_length <= 0:
        return None
    rows = idxstats_rows if idxstats_rows is not None else _samtools_idxstats(reads)
    mapped = 0
    for name, _length, count in rows:
        if name == chrom:
            mapped = count
            break
    if mapped <= 0:
        return 0.0
    rlen = median_read_length or 150
    return (mapped * rlen) / chrom_length


def check_coverage(coverage_x: float | None, chrom: str) -> None:
    if coverage_x is None:
        return
    if coverage_x < MIN_COVERAGE_X:
        raise InputCheckError(
            f"Estimated coverage on {chrom} is {coverage_x:.1f}× (need at least {MIN_COVERAGE_X:.0f}×).",
            "Use deeper WGS, fewer --chroms for a pilot, or confirm the CRAM matches this reference.",
        )


def estimate_run_disk_bytes(reads: Path, *, n_chroms: int, need_ref_download: bool) -> int:
    input_size = reads.stat().st_size
    extract = min(input_size * max(1, n_chroms) / 24.0, input_size * 2)
    ref = GRCH38_REF_BYTES if need_ref_download else 0
    return int(input_size + extract + ref + 2_000_000_000)


def check_disk_space(path: Path, needed_bytes: int) -> None:
    try:
        usage = shutil.disk_usage(path if path.is_dir() else path.parent)
    except OSError:
        return
    if usage.free < needed_bytes:
        need_gb = needed_bytes / (1024**3)
        free_gb = usage.free / (1024**3)
        raise InputCheckError(
            f"Not enough disk space at {path}: need about {need_gb:.1f} GiB free, have {free_gb:.1f} GiB.",
            "Free space, use --chroms for a smaller pilot, or set --out on a larger volume.",
        )


def check_memory_available(required_mb: int) -> None:
    avail_mb: int | None = None
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    avail_mb = int(line.split()[1]) // 1024
                    break
    except OSError:
        pass
    if avail_mb is None:
        return
    if avail_mb < required_mb:
        raise InputCheckError(
            f"Not enough RAM for the largest chromosome job (need about {required_mb} MiB, "
            f"{avail_mb} MiB available).",
            "Use fewer parallel jobs, a machine with more memory, or --chroms chr20 for a pilot.",
        )


def needs_grch38_download(reference: str) -> bool:
    from arcsv2.run import GRCH38_ALIASES, GRCH38_MIN_FA_BYTES, GRCH38_REF_NAME, default_cache_dir

    ref = reference.strip()
    if ref.lower() in GRCH38_ALIASES or ref.upper() == "GRCH38":
        fasta = default_cache_dir() / "GRCh38" / GRCH38_REF_NAME
        return not fasta.is_file() or fasta.stat().st_size < GRCH38_MIN_FA_BYTES
    return False


def validate_call_inputs(
    reads: Path,
    reference: Path,
    gaps: Path,
    *,
    region: str = "",
) -> None:
    """Pre-flight checks for ``arcsv2 call``."""
    check_reads_readable(reads)
    check_reference_files(reference, gaps)
    check_index_present(reads, auto_build=False)
    check_sorted_and_mapped(reads, reference)


def validate_run_inputs(
    reads: Path,
    reference: str,
    out_dir: Path,
    chroms: tuple[str, ...],
    fasta: Path,
    *,
    median_read_length: int | None = None,
    idxstats_rows: list[tuple[str, int, int]] | None = None,
) -> None:
    """Pre-flight checks after reference resolve, before per-chrom extraction."""
    check_reads_readable(reads)
    check_out_writable(out_dir)
    check_reference_files(fasta)
    check_index_present(reads, auto_build=True)
    check_sorted_and_mapped(reads, fasta, idxstats_rows=idxstats_rows)
    check_read_length(median_read_length)

    from arcsv2.pipeline import fasta_length

    probe_chrom = chroms[0] if chroms else "chr20"
    try:
        chrom_len = fasta_length(str(fasta), probe_chrom)
    except (OSError, ValueError):
        chrom_len = 0
    cov = estimate_coverage_x(
        reads,
        probe_chrom,
        chrom_len,
        median_read_length,
        idxstats_rows=idxstats_rows,
    )
    check_coverage(cov, probe_chrom)

    need_ref = needs_grch38_download(reference)
    needed_disk = estimate_run_disk_bytes(reads, n_chroms=len(chroms), need_ref_download=need_ref)
    check_disk_space(out_dir, needed_disk)

    chrom_mbs = []
    for c in chroms:
        try:
            chrom_mbs.append(fasta_length(str(fasta), c) / 1_000_000)
        except (OSError, ValueError):
            chrom_mbs.append(100.0)
    if chrom_mbs:
        largest = max(chrom_mbs)
        required_mb = max(4000, int(largest * 30.68 * 1.2))
        check_memory_available(required_mb)


def validate_run_inputs_early(reads: Path, reference: str, out_dir: Path, chroms: tuple[str, ...]) -> None:
    """Checks that do not need the resolved FASTA (before download)."""
    check_reads_readable(reads)
    check_out_writable(out_dir)
    need_ref = needs_grch38_download(reference)
    needed_disk = estimate_run_disk_bytes(reads, n_chroms=len(chroms), need_ref_download=need_ref)
    check_disk_space(out_dir, needed_disk)
