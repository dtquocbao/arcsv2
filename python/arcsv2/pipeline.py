"""Orchestrate a region call and a chromosome-sharded genome call."""

from __future__ import annotations

from pathlib import Path

from arcsv2.config import CallConfig
from arcsv2.engine.caller import CallRow, call_evidence, write_calls
from arcsv2.progress import Progress
from arcsv2.samio import load_region


def fasta_length(path: str, chrom: str | None) -> int | None:
    """Length of ``chrom`` in a FASTA file. The first record is used when ``chrom`` is empty."""
    length = 0
    seen = False
    current = ""
    target = chrom or ""
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if seen and (not target or current == target):
                    return length
                current = line[1:].split()[0]
                seen = True
                length = 0
                continue
            if seen and (not target or current == target):
                length += len(line.strip())
    if seen and (not target or current == target):
        return length
    return None


def parse_region(region: str) -> tuple[str | None, int | None, int | None]:
    """``chr`` or ``chr:start-end``.

    Coordinates follow samtools: 1-based and inclusive. The returned window is
    0-based and half-open.
    """
    region = region.strip()
    if not region:
        return None, None, None
    if ":" not in region:
        return region, None, None
    chrom, span = region.split(":", 1)
    if "-" not in span:
        return chrom, None, None
    left, right = span.split("-", 1)
    start = int(left) - 1 if left else None
    end = int(right) if right else None
    return chrom, start, end


def call_sam(
    sam_path: str,
    out_tab: str,
    cfg: CallConfig | None = None,
    chrom: str | None = None,
    chrom_len: int | None = None,
    gaps: list[tuple[int, int]] | None = None,
    start: int | None = None,
    end: int | None = None,
    reference: str | None = None,
    progress: Progress | None = None,
) -> list[CallRow]:
    cfg = cfg or CallConfig()
    prog = progress or Progress(label=f"call:{chrom or 'region'}")
    window = ""
    if chrom:
        if start is not None and end is not None:
            window = f"{chrom}:{start + 1}-{end}"
        else:
            window = chrom
    prog.stage("start", f"{sam_path}" + (f" ({window})" if window else ""))
    pairs, extras, splits, chrom_name, inferred_len, softclips = load_region(
        sam_path, chrom, cfg, start, end, reference, progress=prog
    )
    prog.stage(
        "loaded",
        f"{len(pairs)} pairs, {len(extras)} clip/split bps, {len(splits)} splits",
    )
    rows, logs = call_evidence(
        pairs,
        chrom=chrom or chrom_name,
        chrom_len=chrom_len or inferred_len,
        cfg=cfg,
        gaps=gaps,
        extra_breakpoints=extras,
        splits=splits,
        softclips=softclips,
        progress=prog,
    )
    dest = Path(out_tab)
    dest.parent.mkdir(parents=True, exist_ok=True)
    write_calls(rows, out_tab)
    dest.with_name("graph_log.txt").write_text("\n".join(logs) + ("\n" if logs else ""), encoding="utf-8")
    # U1: VCF beside the tab (same calls; .tab stays the golden artifact).
    try:
        from arcsv2.vcf import (
            contigs_from_fai,
            rows_from_call_rows,
            sample_from_bam,
            write_vcf_gz,
        )

        sample = sample_from_bam(sam_path)
        chrom_name = chrom or (rows[0].values.get("chrom") if rows else None) or "unknown"
        contigs = contigs_from_fai(
            reference,
            fallback=[(chrom_name, int(chrom_len or 0))],
        )
        vcf_gz = dest.with_name("arcsv_out.vcf.gz")
        write_vcf_gz(rows_from_call_rows(rows), vcf_gz, sample=sample, contigs=contigs)
        prog.stage("vcf", f"{vcf_gz.name} (+.tbi)")
    except Exception as exc:  # noqa: BLE001 — tab already written; surface VCF issues
        prog.stage("vcf", f"skipped: {exc}")
    prog.stage("done", f"{len(rows)} calls -> {dest}")
    return rows


def _call_one(chrom: str, path: str, dest: str, cfg: CallConfig) -> tuple[str, str | None]:
    try:
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        call_sam(path, dest, cfg, chrom=chrom)
        return dest, None
    except Exception as exc:  # recorded, not swallowed silently
        return dest, f"{chrom}\t{exc}"


def call_genome(
    sam_paths: dict[str, str],
    out_dir: str,
    cfg: CallConfig | None = None,
    jobs: int = 1,
) -> list[str]:
    """One job per chromosome, largest file first. Skip a chromosome whose tab exists."""
    cfg = cfg or CallConfig()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ordered = sorted(sam_paths.items(), key=lambda item: Path(item[1]).stat().st_size, reverse=True)
    pending: list[tuple[str, str, str]] = []
    written: list[str] = []
    for chrom, path in ordered:
        dest = out / chrom / "arcsv_out.tab"
        if dest.exists() and dest.stat().st_size > 0:
            written.append(str(dest))
            continue
        pending.append((chrom, path, str(dest)))
    failures: list[str] = []
    if jobs <= 1 or len(pending) <= 1:
        results = [_call_one(chrom, path, dest, cfg) for chrom, path, dest in pending]
    else:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = [pool.submit(_call_one, chrom, path, dest, cfg) for chrom, path, dest in pending]
            results = [fut.result() for fut in futures]
    for dest, failure in results:
        if failure:
            failures.append(failure)
        else:
            written.append(dest)
    if failures:
        (out / "failures.tsv").write_text("\n".join(failures) + "\n", encoding="utf-8")
    return written
