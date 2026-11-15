"""U2: end-to-end ``arcsv2 run`` orchestration (CRAM → per-chrom calls → merged.tab)."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from collections.abc import Callable
from pathlib import Path

from arcsv2.beds import read_bed, read_bed_by_chrom
from arcsv2.config import CallConfig, OUTPUT_COLUMNS
from arcsv2.merge_filter import inslen_is_lower_bound, merge_calls, passes_filter
from arcsv2.pipeline import call_sam, fasta_length
from arcsv2.atomic_io import atomic_write_text
from arcsv2.input_checks import InputCheckError, validate_run_inputs, validate_run_inputs_early
from arcsv2.vcf import sample_from_bam

PRIMARY_CHROMS: tuple[str, ...] = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX",)

GRCH38_ALIASES = frozenset(
    {"grch38", "hg38", "grch38_full_analysis_set_plus_decoy_hla", "grch38_full_analysis_set_plus_decoy_hla.fa"}
)

GRCH38_REF_URL = (
    "s3://1000genomes/technical/reference/GRCh38_reference_genome/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)
GRCH38_REF_NAME = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
GRCH38_MIN_FA_BYTES = 2_800_000_000
# Uncompressed GRCh38_full_analysis_set_plus_decoy_hla.fa is ~3.04 GiB; reject
# doubled/corrupt caches from concurrent resume downloads (seen at ~6 GiB).
GRCH38_MAX_FA_BYTES = 3_500_000_000

S3_HTTPS = {
    "1000genomes": "https://s3.amazonaws.com/1000genomes",
}

USER_AGENT = "arcsv2-run/1"


def _resolve_url(url: str) -> str:
    if not url.startswith("s3://"):
        return url
    rest = url[5:]
    bucket, _, key = rest.partition("/")
    base = S3_HTTPS.get(bucket)
    if base is None:
        raise RuntimeError(
            f"Unknown S3 bucket {bucket!r}. Pass an https:// URL or a local FASTA path for --reference."
        )
    return f"{base}/{key}"


def default_cache_dir() -> Path:
    base = os.environ.get("ARCSV2_CACHE", "").strip()
    if base:
        return Path(base).expanduser()
    return Path.home() / ".cache" / "arcsv2"


def parse_chroms(chroms: str | None) -> tuple[str, ...]:
    if not chroms or not chroms.strip():
        return PRIMARY_CHROMS
    return tuple(c.strip() for c in chroms.replace(",", " ").split() if c.strip())


def resolve_threads(threads: str | int) -> int:
    if isinstance(threads, str) and threads.lower() == "auto":
        n = os.cpu_count() or 4
        return max(1, min(n, 16))
    return max(1, int(threads))


def resolve_jobs(threads: int, chrom_lengths_mb: list[float]) -> tuple[int, int]:
    """Return (parallel chrom jobs, threads per job)."""
    if not chrom_lengths_mb:
        return 1, threads
    largest = max(chrom_lengths_mb)
    mem_mb_per_job = max(4000, int(largest * 30.68 * 1.2))
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    avail_kb = int(line.split()[1])
                    avail_mb = avail_kb // 1024
                    jobs = max(1, min(threads // 2 or 1, avail_mb // mem_mb_per_job))
                    jobs = min(jobs, 4)
                    if jobs <= 1:
                        return 1, threads
                    per = max(1, threads // jobs)
                    return jobs, per
    except OSError:
        pass
    return 1, threads


def _digest_file(path: Path, algo: str) -> str:
    h = hashlib.md5() if algo == "md5" else hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def _parse_checksum(text: str) -> tuple[str, str] | None:
    text = (text or "").strip()
    if not text or text.lower() == "none":
        return None
    if text.startswith("md5:"):
        return "md5", text[4:].lower()
    if text.startswith("sha256:"):
        return "sha256", text[7:].lower()
    if len(text) == 32 and all(c in "0123456789abcdef" for c in text.lower()):
        return "md5", text.lower()
    return "sha256", text.lower()


def download_with_progress(
    url: str,
    dest: Path,
    *,
    expected_checksum: str = "",
    on_progress: Callable[[int, int | None], None] | None = None,
    max_bytes: int | None = None,
) -> None:
    """Download ``url`` to ``dest`` with optional resume and checksum verify."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    lock = dest.with_suffix(dest.suffix + ".lock")
    lock_fh = None
    try:
        lock_fh = lock.open("a+b")
        try:
            if sys.platform == "win32":
                import msvcrt

                msvcrt.locking(lock_fh.fileno(), msvcrt.LK_LOCK, 1)
            else:
                import fcntl

                fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
        except OSError:
            pass
        parsed = _parse_checksum(expected_checksum)
        if parsed and dest.is_file():
            algo, expected = parsed
            if _digest_file(dest, algo) == expected:
                return
        if (
            dest.is_file()
            and dest.stat().st_size >= GRCH38_MIN_FA_BYTES
            and (max_bytes is None or dest.stat().st_size <= max_bytes)
            and not parsed
        ):
            return
        resolved = _resolve_url(url)
        tmp = dest.with_suffix(dest.suffix + ".partial")
        existing = tmp.stat().st_size if tmp.is_file() else 0
        if max_bytes is not None and existing > max_bytes:
            tmp.unlink(missing_ok=True)
            existing = 0
        headers = {"User-Agent": USER_AGENT}
        if existing:
            headers["Range"] = f"bytes={existing}-"
        request = urllib.request.Request(resolved, headers=headers)
        with urllib.request.urlopen(request, timeout=120) as response:
            total_raw = response.headers.get("Content-Length")
            total = int(total_raw) + existing if total_raw and existing else (int(total_raw) if total_raw else None)
            status = getattr(response, "status", 200)
            if existing and status == 200:
                tmp.unlink(missing_ok=True)
                existing = 0
                if total_raw:
                    total = int(total_raw)
            if max_bytes is not None and total is not None and total > max_bytes:
                raise RuntimeError(
                    f"Download size {total} exceeds max_bytes={max_bytes} for {dest.name}; "
                    "refusing to write (possible corrupt resume)."
                )
            mode = "ab" if existing else "wb"
            written = existing
            with tmp.open(mode) as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    written += len(chunk)
                    if max_bytes is not None and written > max_bytes:
                        handle.close()
                        tmp.unlink(missing_ok=True)
                        raise RuntimeError(
                            f"Download exceeded max_bytes={max_bytes} for {dest.name}; "
                            "removed partial (concurrent resume corruption)."
                        )
                    if on_progress:
                        on_progress(written, total)
        if parsed:
            algo, expected = parsed
            digest = _digest_file(tmp, algo)
            if digest != expected:
                tmp.unlink(missing_ok=True)
                raise RuntimeError(
                    f"Checksum mismatch for {dest.name}: got {algo}:{digest}, expected {algo}:{expected}"
                )
        tmp.replace(dest)
    finally:
        if lock_fh is not None:
            try:
                if sys.platform == "win32":
                    import msvcrt

                    lock_fh.seek(0)
                    msvcrt.locking(lock_fh.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl

                    fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
            lock_fh.close()
            lock.unlink(missing_ok=True)


def _gaps_from_fasta(fasta: Path, out_bed: Path) -> None:
    """Write half-open N-run intervals (same rule as scripts/make_ref_beds.py)."""
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
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with out_bed.open("w", encoding="utf-8", newline="\n") as bed:
        for chrom_name, start, end in intervals:
            bed.write(f"{chrom_name}\t{start}\t{end}\n")


def _write_fai_from_fasta(fasta: Path, fai: Path) -> None:
    """Minimal .fai for small FASTAs when samtools is unavailable (tests only)."""
    text = fasta.read_text(encoding="utf-8", errors="replace")
    lines_out: list[str] = []
    name = ""
    length = 0
    line_bases = 0
    offset = 0
    pos = 0
    while pos < len(text):
        line_end = text.find("\n", pos)
        if line_end < 0:
            line_end = len(text)
        raw = text[pos:line_end]
        pos = line_end + 1
        if raw.startswith(">"):
            if name:
                lines_out.append(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_bases + 1}\n")
            name = raw[1:].split()[0]
            length = 0
            line_bases = 0
            offset = pos
            continue
        strip = raw.strip()
        if strip:
            length += len(strip)
            line_bases = len(strip)
    if name:
        lines_out.append(f"{name}\t{length}\t{offset}\t{line_bases}\t{line_bases + 1}\n")
    fai.write_text("".join(lines_out), encoding="utf-8")


def ensure_fai(fasta: Path) -> Path:
    fai = fasta.with_suffix(fasta.suffix + ".fai") if fasta.suffix else Path(str(fasta) + ".fai")
    if not fai.is_file():
        if shutil.which("samtools"):
            _run_samtools(["faidx", str(fasta)])
        else:
            _write_fai_from_fasta(fasta, fai)
    return fai


def resolve_reference(reference: str, *, cache_dir: Path | None = None) -> tuple[Path, Path]:
    """Return (fasta, gaps.bed). ``reference`` may be GRCh38 alias or a local path."""
    ref = reference.strip()
    if ref.lower() in GRCH38_ALIASES or ref.upper() == "GRCH38":
        root = (cache_dir or default_cache_dir()) / "GRCh38"
        fasta = root / GRCH38_REF_NAME
        gaps = root / "gaps.bed"
        size = fasta.stat().st_size if fasta.is_file() else 0
        if size > GRCH38_MAX_FA_BYTES:
            print(
                f"Removing oversized/corrupt reference cache ({size} bytes) at {fasta}",
                file=sys.stderr,
                flush=True,
            )
            fasta.unlink(missing_ok=True)
            fasta.with_suffix(fasta.suffix + ".fai").unlink(missing_ok=True)
            fasta.with_suffix(fasta.suffix + ".partial").unlink(missing_ok=True)
            size = 0
        if size < GRCH38_MIN_FA_BYTES:
            last_pct = [-1]

            def progress(done: int, total: int | None) -> None:
                if total:
                    pct = min(100, int(100 * done / total))
                    if pct != last_pct[0] and pct % 5 == 0:
                        mb = done // (1024 * 1024)
                        tot_mb = total // (1024 * 1024)
                        print(f"  reference download: {pct}% ({mb}/{tot_mb} MiB)", file=sys.stderr, flush=True)
                        last_pct[0] = pct

            print(f"Downloading GRCh38 reference to {fasta} …", file=sys.stderr, flush=True)
            download_with_progress(
                GRCH38_REF_URL, fasta, on_progress=progress, max_bytes=GRCH38_MAX_FA_BYTES
            )
        ensure_fai(fasta)
        if not gaps.is_file():
            print(f"Building {gaps.name} from reference N-runs …", file=sys.stderr, flush=True)
            _gaps_from_fasta(fasta, gaps)
        return fasta, gaps
    fasta = Path(ref).expanduser().resolve()
    if not fasta.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta}")
    ensure_fai(fasta)
    gaps = fasta.parent / "gaps.bed"
    if not gaps.is_file():
        _gaps_from_fasta(fasta, gaps)
    return fasta, gaps


def _run_samtools(args: list[str]) -> None:
    exe = shutil.which("samtools")
    if not exe:
        raise RuntimeError("samtools is required on PATH for arcsv2 run (indexing and CRAM extraction).")
    proc = subprocess.run([exe, *args], capture_output=True, text=True)
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"samtools {' '.join(args)} failed: {msg}")


def ensure_index(reads: Path) -> None:
    if reads.suffix == ".cram":
        idx = Path(str(reads) + ".crai")
        if not idx.is_file():
            _run_samtools(["index", str(reads)])
        return
    if reads.suffix in (".bam", ".sam"):
        idx = reads.with_suffix(reads.suffix + ".bai") if reads.suffix == ".bam" else None
        if idx and not idx.is_file():
            _run_samtools(["index", str(reads)])


def fai_chroms(fai: Path) -> list[str]:
    names: list[str] = []
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                names.append(line.split("\t", 1)[0])
    return names


def detect_chrom_style(reads: Path, fai: Path, chroms: tuple[str, ...]) -> tuple[str, ...]:
    """Map requested chrom names to names present in the reference index."""
    del reads
    ref_names = set(fai_chroms(fai))
    if all(c in ref_names for c in chroms):
        return chroms
    alt = tuple(c[3:] if c.startswith("chr") else f"chr{c}" for c in chroms)
    if all(c in ref_names for c in alt):
        return alt
    missing = [c for c in chroms if c not in ref_names][:5]
    raise InputCheckError(
        f"Chromosome names not in reference index: {', '.join(missing)}.",
        "Set --chroms to names that match the FASTA (e.g. chr20 vs 20).",
    )


def cram_reference_warning(reads: Path, fasta: Path) -> str | None:
    try:
        import pysam
    except ImportError:
        return None
    try:
        with pysam.AlignmentFile(str(reads), "rb") as handle:
            cram_ref = handle.header.get("HD", {}).get("SQ", [])
            if not cram_ref:
                return None
            # Best-effort: compare first SQ length to fai
            fai = fasta.with_suffix(fasta.suffix + ".fai") if fasta.suffix else Path(str(fasta) + ".fai")
            if not fai.is_file():
                return None
            sq0 = handle.header.references[0] if handle.header.references else None
            if sq0 and sq0 not in fai_chroms(fai):
                return (
                    f"CRAM header references contig {sq0!r} which is not in {fasta.name}; "
                    "ensure -T/--reference matches the CRAM encoding."
                )
    except (OSError, ValueError):
        return None
    return None


def estimate_read_stats(reads: Path, reference: Path, *, max_reads: int = 20_000) -> tuple[int | None, int | None]:
    """Median read length and insert size from a small header scan."""
    try:
        import pysam
    except ImportError:
        return None, None
    lengths: list[int] = []
    isizes: list[int] = []
    try:
        kwargs = {}
        if str(reads).lower().endswith(".cram") and reference.is_file():
            kwargs["reference_filename"] = str(reference)
        with pysam.AlignmentFile(str(reads), "rb", **kwargs) as handle:
            for i, read in enumerate(handle):
                if i >= max_reads:
                    break
                if read.query_length:
                    lengths.append(int(read.query_length))
                if read.is_paired and read.is_proper_pair and not read.is_unmapped and read.template_length:
                    isizes.append(abs(int(read.template_length)))
    except (OSError, ValueError, TypeError):
        return None, None
    if not lengths:
        return None, None
    lengths.sort()
    med_len = lengths[len(lengths) // 2]
    med_isize = None
    if isizes:
        isizes.sort()
        med_isize = isizes[len(isizes) // 2]
    return med_len, med_isize


def extract_chrom_bam(reads: Path, reference: Path, chrom: str, dest_bam: Path) -> None:
    dest_bam.parent.mkdir(parents=True, exist_ok=True)
    bai = dest_bam.with_suffix(dest_bam.suffix + ".bai")
    if dest_bam.is_file() and dest_bam.stat().st_size > 0 and bai.is_file():
        return
    samtools = shutil.which("samtools")
    if not samtools:
        raise RuntimeError("samtools is required on PATH for arcsv2 run.")
    tmp = dest_bam.with_suffix(".bam.tmp")
    view = subprocess.Popen(
        [samtools, "view", "-u", "-T", str(reference), str(reads), chrom],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert view.stdout is not None
    sort = subprocess.run(
        [samtools, "sort", "-o", str(tmp), "-"],
        stdin=view.stdout,
        capture_output=True,
        text=True,
    )
    view.wait()
    if view.returncode != 0 or sort.returncode != 0:
        tmp.unlink(missing_ok=True)
        err = (sort.stderr or "") + (view.stderr.read().decode() if view.stderr else "")
        raise RuntimeError(f"Failed to extract {chrom} from {reads.name}: {err.strip()}")
    tmp.replace(dest_bam)
    _run_samtools(["index", str(dest_bam), str(bai)])


def merge_chrom_tables(
    tables: list[Path],
    out_merged: Path,
    *,
    sample: str = "SAMPLE",
    gaps_bed: Path | None = None,
    min_size: int = 50,
) -> int:
    gap_iv = read_bed_by_chrom(gaps_bed) if gaps_bed and gaps_bed.is_file() else {}
    records = []
    for table in tables:
        if not table.is_file():
            continue
        lines = table.read_text(encoding="utf-8").splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        for line in lines[1:]:
            row = dict(zip(header, line.split("\t")))
            start = int(row["minbp"])
            end = int(row["maxbp"])
            affected = int(float(row["len_affected"]))
            chrom = row["chrom"]
            ok, _reason = passes_filter(
                start,
                end,
                affected,
                gaps=gap_iv.get(chrom, []),
                min_affected=min_size,
                length_is_lower_bound=inslen_is_lower_bound(row.get("inslen")),
            )
            if not ok:
                continue
            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "sv_class": row.get("complextype")
                    if row.get("complextype") not in ("", "NA")
                    else row.get("svtype"),
                    "gt": row.get("gt", "HET"),
                    "sample": sample,
                    "row": row,
                }
            )
    merged = merge_calls(records, n_samples=1)
    lines = ["\t".join(OUTPUT_COLUMNS)]
    for rec in merged:
        row = dict(rec["row"])
        row["af"] = f"{rec['allele_frequency']:.6f}"
        lines.append("\t".join(row.get(col, "NA") for col in OUTPUT_COLUMNS))
    atomic_write_text(out_merged, "\n".join(lines) + "\n")
    return len(merged)


def _summarize_calls(tables: list[Path]) -> dict[str, object]:
    by_type: dict[str, int] = {}
    size_bins = {"50-499": 0, "500-4999": 0, "5000-49999": 0, "50000+": 0}
    total = 0
    for table in tables:
        if not table.is_file():
            continue
        lines = table.read_text(encoding="utf-8").splitlines()
        if len(lines) < 2:
            continue
        header = lines[0].split("\t")
        for line in lines[1:]:
            if not line.strip():
                continue
            row = dict(zip(header, line.split("\t")))
            total += 1
            st = row.get("svtype") or "UNK"
            by_type[st] = by_type.get(st, 0) + 1
            try:
                length = int(float(row.get("len_affected", 0)))
            except ValueError:
                length = 0
            if length < 500:
                size_bins["50-499"] += 1
            elif length < 5000:
                size_bins["500-4999"] += 1
            elif length < 50000:
                size_bins["5000-49999"] += 1
            else:
                size_bins["50000+"] += 1
    return {"total": total, "by_type": by_type, "size_bins": size_bins}


def _call_one_chrom(
    chrom: str,
    bam: Path,
    reference: Path,
    gaps_path: Path,
    out_tab: Path,
    cfg: CallConfig,
) -> None:
    gap_list = read_bed(gaps_path, chrom) if gaps_path.is_file() else []
    chrom_len = fasta_length(str(reference), chrom)
    call_sam(
        str(bam),
        str(out_tab),
        cfg,
        chrom=chrom,
        chrom_len=chrom_len,
        gaps=gap_list,
        reference=str(reference),
    )


def run_sample(
    *,
    reads: Path,
    reference: str,
    out_dir: Path,
    chroms: str | None = None,
    threads: str | int = "auto",
    seed: int = 204178949,
    engine: str = "rust",
    min_size: int = 50,
    nest_sample_dir: bool = True,
) -> dict[str, object]:
    """Run caller on all requested chromosomes; return summary metadata."""
    t0 = time.perf_counter()
    reads = reads.expanduser().resolve()
    out_dir = out_dir.expanduser().resolve()
    requested = parse_chroms(chroms)
    validate_run_inputs_early(reads, reference, out_dir, requested)

    fasta, gaps_path = resolve_reference(reference)
    ensure_index(reads)
    sample = sample_from_bam(reads, default=reads.stem)
    warn = cram_reference_warning(reads, fasta)
    if warn:
        print(f"Warning: {warn}", file=sys.stderr)

    fai = ensure_fai(fasta)
    chrom_list = detect_chrom_style(reads, fai, requested)
    med_rlen, med_isize = estimate_read_stats(reads, fasta)
    validate_run_inputs(
        reads,
        reference,
        out_dir,
        chrom_list,
        fasta,
        median_read_length=med_rlen,
    )

    thread_n = resolve_threads(threads)
    chrom_mbs = []
    for c in chrom_list:
        try:
            chrom_mbs.append(fasta_length(str(fasta), c) / 1_000_000)
        except (OSError, ValueError):
            chrom_mbs.append(100.0)
    jobs, threads_per_job = resolve_jobs(thread_n, chrom_mbs)
    cfg = CallConfig(seed=seed, threads=threads_per_job, engine=engine)

    work = out_dir / "_work"
    sample_root = (out_dir / sample) if nest_sample_dir else out_dir
    sample_root.mkdir(parents=True, exist_ok=True)

    pending: list[tuple[str, Path, Path]] = []
    done_tabs: list[Path] = []
    for chrom in chrom_list:
        out_tab = sample_root / chrom / "arcsv_out.tab"
        if out_tab.is_file() and out_tab.stat().st_size > 0:
            done_tabs.append(out_tab)
            continue
        bam = work / f"{chrom}.bam"
        pending.append((chrom, bam, out_tab))

    def _progress(**fields: object) -> None:
        try:
            from arcsv2.atomic_io import atomic_write_json

            path = sample_root / "progress.json"
            payload: dict[str, object] = {}
            if path.is_file():
                import json as _json

                payload.update(_json.loads(path.read_text(encoding="utf-8")))
            payload.update(fields)
            payload["updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            atomic_write_json(path, payload)
        except Exception:
            pass

    _progress(
        stage="extract" if pending else "call",
        percent=5.0 if pending else 50.0,
        started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        chroms=list(chrom_list),
        n_chroms=len(chrom_list),
    )

    total_pending = len(pending)
    for i, (chrom, bam, out_tab) in enumerate(pending):
        extract_chrom_bam(reads, fasta, chrom, bam)
        _progress(
            stage="call",
            percent=10.0 + 70.0 * (i / max(1, total_pending)),
            done_chroms=i,
            n_chroms=len(chrom_list),
        )

    def run_chrom(item: tuple[str, Path, Path]) -> Path:
        chrom, bam, out_tab = item
        _call_one_chrom(chrom, bam, fasta, gaps_path, out_tab, cfg)
        return out_tab

    if jobs <= 1 or len(pending) <= 1:
        for i, item in enumerate(pending):
            done_tabs.append(run_chrom(item))
            _progress(
                stage="call",
                percent=10.0 + 70.0 * ((i + 1) / max(1, total_pending or 1)),
                done_chroms=len(done_tabs),
                n_chroms=len(chrom_list),
            )
    else:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=jobs) as pool:
            done_tabs.extend(pool.map(run_chrom, pending))

    _progress(stage="merge", percent=90.0)

    merged = sample_root / "merged.tab"
    n_merged = merge_chrom_tables(
        sorted(done_tabs), merged, sample=sample, gaps_bed=gaps_path, min_size=min_size
    )
    stats = _summarize_calls(done_tabs)
    elapsed = time.perf_counter() - t0
    _progress(stage="done", percent=100.0)
    coverage_x = None
    if chrom_list:
        try:
            from arcsv2.input_checks import estimate_coverage_x

            clen = fasta_length(str(fasta), chrom_list[0])
            coverage_x = estimate_coverage_x(reads, chrom_list[0], clen, med_rlen)
        except (OSError, ValueError):
            coverage_x = None

    return {
        "sample": sample,
        "coverage_x": coverage_x,
        "reference": str(fasta),
        "gaps": str(gaps_path),
        "chroms": chrom_list,
        "threads": thread_n,
        "jobs": jobs,
        "threads_per_job": threads_per_job,
        "median_read_length": med_rlen,
        "median_insert_size": med_isize,
        "merged_calls": n_merged,
        "stats": stats,
        "paths": {
            "out_dir": str(out_dir),
            "sample_dir": str(sample_root),
            "merged_tab": str(merged),
            "chrom_tabs": [str(p) for p in sorted(done_tabs)],
        },
        "runtime_sec": elapsed,
    }


def print_run_summary(meta: dict[str, object]) -> None:
    stats = meta.get("stats") or {}
    by_type = stats.get("by_type") or {}
    size_bins = stats.get("size_bins") or {}
    paths = meta.get("paths") or {}
    print("")
    print("=== arcsv2 run summary ===")
    print(f"sample: {meta.get('sample')}")
    print(f"runtime: {float(meta.get('runtime_sec', 0)):.1f}s")
    print(f"threads: {meta.get('threads')} ({meta.get('jobs')} chrom job(s), {meta.get('threads_per_job')} per job)")
    if meta.get("median_read_length"):
        print(f"median read length: {meta.get('median_read_length')} bp")
    if meta.get("median_insert_size"):
        print(f"median insert size: {meta.get('median_insert_size')} bp")
    print(f"raw calls (all chroms): {stats.get('total', 0)}")
    print(f"merged/filtered (≥{50} bp): {meta.get('merged_calls')}")
    if by_type:
        print("by type:")
        for key in sorted(by_type):
            print(f"  {key}: {by_type[key]}")
    if size_bins:
        print("by affected length:")
        for key in ("50-499", "500-4999", "5000-49999", "50000+"):
            print(f"  {key} bp: {size_bins.get(key, 0)}")
    print("outputs:")
    print(f"  {paths.get('merged_tab')}")
    for p in paths.get("chrom_tabs") or []:
        vcf = Path(str(p)).with_name("arcsv_out.vcf.gz")
        if vcf.is_file():
            print(f"  {p}")
            print(f"  {vcf}")
        else:
            print(f"  {p}")
