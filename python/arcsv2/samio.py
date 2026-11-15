"""SAM → pair records, soft clips, and split reads.

BAM/CRAM goes through pysam when it is installed. Filters follow §4.1: not
secondary, not duplicate, not QC-fail. MAPQ cuts are applied later, per use
(model vs soft-clip discovery). Insert length follows orig:process_insert_len
(signed outer span + soft-clip / trim adjustments), not abs(TLEN).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace

from arcsv2.config import CallConfig
from arcsv2.engine.breakpoints import Breakpoint
from arcsv2.engine.caller import PairRecord
from arcsv2.engine.pairs import Pairs
from arcsv2.engine.evidence import (
    SplitHit,
    _softclip_counts,
    aligned_length,
    cigar_ops,
    clip_breakpoints,
    parse_split,
    softclips_from_alignment,
    split_breakpoints,
)
from arcsv2.progress import Progress

FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_REVERSE = 0x10
FLAG_MATE_REVERSE = 0x20
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUP = 0x400
FLAG_SUPPLEMENTARY = 0x800
# Drop before building SamRead (same mask as _primary).
_SKIP_FLAGS = FLAG_SECONDARY | FLAG_DUP | FLAG_QCFAIL | FLAG_SUPPLEMENTARY


@dataclass
class SamRead:
    qname: str
    flag: int
    rname: str
    start: int  # 0-based
    mapq: int
    cigar: str
    seq: str
    qual: str
    tlen: int
    tags: dict[str, str]


def _parse_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields:
        if ":" not in field:
            continue
        key, typ, value = field.split(":", 2)
        tags[key] = value
    return tags


def iter_sam(path: str):
    """Yield primary-alignment records one line at a time."""
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            yield SamRead(
                qname=parts[0],
                flag=int(parts[1]),
                rname=parts[2],
                start=int(parts[3]) - 1,
                mapq=int(parts[4]),
                cigar=parts[5],
                seq=parts[9],
                qual=parts[10],
                tlen=int(parts[8]),
                tags=_parse_tags(parts[11:]),
            )


def read_sam(path: str) -> list[SamRead]:
    return list(iter_sam(path))


def _aligned_end(read: SamRead) -> int:
    alen = aligned_length(read.cigar)
    if read.seq and read.seq != "*":
        return read.start + max(alen, len(read.seq))
    return read.start + alen


def _in_window(read: SamRead, start: int | None, end: int | None) -> bool:
    """True when the alignment overlaps half-open ``[start, end)``."""
    if start is None and end is None:
        return True
    read_end = _aligned_end(read)
    if start is not None and read_end <= start:
        return False
    if end is not None and read.start >= end:
        return False
    return True


def _primary(read: SamRead) -> bool:
    return not (read.flag & (FLAG_SECONDARY | FLAG_DUP | FLAG_QCFAIL | FLAG_SUPPLEMENTARY))


def _on_target(read: SamRead, chrom: str | None, start: int | None, end: int | None) -> bool:
    """Mapped read on the requested chromosome and inside the half-open window."""
    if read.flag & FLAG_UNMAPPED:
        return False
    if chrom and read.rname != chrom:
        return False
    return _in_window(read, start, end)


def _query_len(read: SamRead) -> int:
    if read.seq and read.seq != "*":
        return len(read.seq)
    # CIGAR query-consuming ops when SEQ is absent (BAM often omits it).
    return sum(n for n, op in cigar_ops(read.cigar) if op in "MIS=X")


def _hanging_fields(anchor: SamRead, mate: SamRead, distant: bool) -> dict:
    mate_unmapped = bool(mate.flag & FLAG_UNMAPPED)
    insert = 0
    if not mate_unmapped:
        computed = insert_length(anchor, mate)
        insert = 0 if computed is None else computed
    return dict(
        start1=anchor.start,
        end1=anchor.start + aligned_length(anchor.cigar),
        start2=0 if mate_unmapped else mate.start,
        end2=0 if mate_unmapped else mate.start + aligned_length(mate.cigar),
        reverse1=bool(anchor.flag & FLAG_REVERSE),
        reverse2=bool(mate.flag & FLAG_REVERSE),
        mapq1=anchor.mapq,
        mapq2=0 if mate_unmapped else mate.mapq,
        insert=insert,
        lib_id=_lib_id(anchor.tags),
        hanging=True,
        qname_hash=_qhash(anchor.qname),
        distant=distant,
        mate_len=_query_len(mate),
    )


def _hanging_pair(anchor: SamRead, mate: SamRead, distant: bool) -> PairRecord:
    return PairRecord(**_hanging_fields(anchor, mate, distant))


def _lib_id(tags: dict[str, str]) -> int:
    rg = tags.get("LB") or tags.get("RG") or "lib1"
    return int(hashlib.sha256(rg.encode()).hexdigest()[:4], 16) % 1000


def _qhash(qname: str) -> int:
    return int(hashlib.sha256(qname.encode()).hexdigest()[:8], 16)


def insert_length(
    a: SamRead,
    b: SamRead,
    *,
    read_len: int = 0,
    lib_is_rf: bool = False,
    lib_insert_is_inner: bool = False,
) -> int | None:
    """Signed insert length. orig:bamparser_streaming.process_insert_len.

    FR (default): first = forward, last = reverse; outer insert is
    ``last.reference_end − first.reference_start``. Outward Dup pairs yield a
    small or negative value. Same-orientation pairs return ``None``.
    """
    a_rev = bool(a.flag & FLAG_REVERSE)
    b_rev = bool(b.flag & FLAG_REVERSE)
    if a_rev == b_rev:
        return None
    which_minus = 0 if a_rev else 1
    which_first = which_minus if lib_is_rf else (1 - which_minus)
    which_last = 1 - which_first
    pair = (a, b)
    first, last = pair[which_first], pair[which_last]
    first_end = first.start + aligned_length(first.cigar)
    last_end = last.start + aligned_length(last.cigar)
    if lib_insert_is_inner:
        ilen = last.start - first_end
    else:
        ilen = last_end - first.start
    qlen1 = _query_len(a)
    qlen2 = _query_len(b)
    if read_len != 0:
        ilen += 2 * read_len - qlen1 - qlen2
    # Soft-clip of first's 5' and last's 3' (pysam query_alignment_start /
    # query_length - query_alignment_end).
    first_left, _ = _softclip_counts(first.cigar)
    _, last_right = _softclip_counts(last.cigar)
    ilen += first_left + last_right
    return int(ilen)


def _cigar_has_softclip(cigar: str | None) -> bool:
    return bool(cigar) and "S" in cigar


def _from_pysam(aln, *, want_qual: bool = False) -> SamRead:
    """One pysam record in the same fields as a SAM line.

    Qualities are materialised only when ``want_qual`` is set (soft-clipped
    reads). Sequence is left as ``*``; query length comes from the CIGAR.
    """
    tags: dict[str, str] = {}
    for key in ("SA", "RG", "LB"):
        if aln.has_tag(key):
            tags[key] = str(aln.get_tag(key))
    if want_qual and aln.query_qualities is not None:
        qual = "".join(chr(int(base) + 33) for base in aln.query_qualities)
    else:
        qual = ""
    unmapped = bool(aln.flag & FLAG_UNMAPPED)
    start = 0 if unmapped or aln.reference_start is None else int(aln.reference_start)
    return SamRead(
        qname=aln.query_name or "",
        flag=int(aln.flag),
        rname="*" if unmapped or not aln.reference_name else aln.reference_name,
        start=start,
        mapq=int(aln.mapping_quality or 0),
        cigar=aln.cigarstring or "*",
        seq="*",
        qual=qual,
        tlen=int(aln.template_length or 0),
        tags=tags,
    )


def iter_bam(
    path: str,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    reference: str | None = None,
    threads: int = 1,
    max_nm: int | None = None,
):
    """Stream a BAM or CRAM. ``start`` and ``end`` are 0-based and half-open.

    A CRAM needs ``reference``. The fetch uses the index, so mates that do not
    overlap the window are not returned. Secondary / duplicate / QC-fail /
    supplementary alignments are dropped before building a ``SamRead``.
    When ``max_nm`` is set (C5 ``filter_high_nm``), alignments with NM above
    that threshold are skipped. Default callers leave ``max_nm`` as None.
    """
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError(
            "BAM/CRAM input needs pysam. There is no native Windows build; "
            "use the Linux container in docs/CROSS_PLATFORM.md (pip install arcsv2[io])."
        ) from exc
    lower = path.lower()
    cram = lower.endswith(".cram")
    if cram and not reference:
        raise ValueError("CRAM input requires a reference FASTA")
    kwargs: dict = {}
    if cram:
        kwargs["reference_filename"] = reference
    if threads and threads > 1:
        kwargs["threads"] = int(threads)
    with pysam.AlignmentFile(path, "rc" if cram else "rb", **kwargs) as handle:
        if chrom:
            fetched = handle.fetch(chrom, start or 0, end)
        else:
            fetched = handle.fetch(until_eof=True)
        for aln in fetched:
            flag = int(aln.flag)
            if flag & _SKIP_FLAGS:
                continue
            if max_nm is not None and aln.has_tag("NM") and int(aln.get_tag("NM")) > max_nm:
                continue
            cigar = aln.cigarstring or "*"
            # Qualities only when a soft clip may need the Phred window check.
            want_qual = _cigar_has_softclip(cigar)
            yield _from_pysam(aln, want_qual=want_qual)


def _iter_alignments(
    path: str,
    chrom: str | None,
    start: int | None,
    end: int | None,
    reference: str | None,
    threads: int = 1,
    max_nm: int | None = None,
):
    lower = path.lower()
    if lower.endswith(".bam") or lower.endswith(".cram"):
        yield from iter_bam(
            path, chrom, start, end, reference, threads=threads, max_nm=max_nm
        )
    else:
        yield from iter_sam(path)


def _load_region_rust(
    path: str,
    chrom: str | None,
    cfg: CallConfig,
    start: int | None,
    end: int | None,
    reference: str | None,
    progress: Progress | None,
) -> tuple[Pairs, list[Breakpoint], list[SplitHit], str, int, list]:
    """BAM/CRAM via arcsv2_native; softclip/split discovery stays in evidence.py."""
    from arcsv2 import native as native_mod

    if progress is not None:
        progress.stage("load_bam", f"{path} (engine=rust)")
    raw = native_mod.load_bam_pairs(
        path,
        chrom=chrom,
        start=start,
        end=end,
        reference=reference,
        max_pair_distance=int(cfg.max_pair_distance),
        threads=max(1, int(cfg.threads)),
        filter_high_nm=bool(cfg.filter_high_nm),
        max_nm=int(cfg.max_nm),
        filter_both_mates_clipped=bool(cfg.filter_both_mates_clipped),
    )
    pairs = Pairs.from_arrays(
        start1=raw["start1"],
        end1=raw["end1"],
        start2=raw["start2"],
        end2=raw["end2"],
        insert=raw["insert"],
        mate_len=raw["mate_len"],
        qname_hash=raw["qname_hash"],
        mapq1=raw["mapq1"],
        mapq2=raw["mapq2"],
        lib_id=raw["lib_id"],
        flags=raw["flags"],
    )
    splits: list[SplitHit] = []
    clips = []
    # M4: columnar takes (preferred) or legacy list-of-dicts.
    n_takes = int(raw.get("n_takes", 0))
    if n_takes > 0 and "take_start" in raw:
        t_start = raw.pop("take_start")
        t_mapq = raw.pop("take_mapq")
        t_reverse = raw.pop("take_reverse")
        t_qname = raw.pop("take_qname")
        t_rname = raw.pop("take_rname")
        t_cigar = raw.pop("take_cigar")
        t_sa = raw.pop("take_sa")
        t_qual = raw.pop("take_qual")
        for i in range(n_takes):
            sa = t_sa[i] or ""
            cigar = t_cigar[i]
            has_clip = _cigar_has_softclip(cigar)
            if not sa and not has_clip:
                continue
            hit = parse_split(
                chrom=t_rname[i],
                start=int(t_start[i]),
                mapq=int(t_mapq[i]),
                cigar=cigar,
                reverse=bool(int(t_reverse[i])),
                sa=sa,
                cfg=cfg,
            )
            if hit is not None:
                splits.append(replace(hit, qname_hash=_qhash(t_qname[i])))
                continue
            if not has_clip:
                continue
            clips.extend(
                softclips_from_alignment(
                    start=int(t_start[i]),
                    mapq=int(t_mapq[i]),
                    cigar=cigar,
                    qual=t_qual[i] or "",
                    cfg=cfg,
                )
            )
        del t_start, t_mapq, t_reverse, t_qname, t_rname, t_cigar, t_sa, t_qual
    else:
        takes = raw.pop("takes", [])
        for take in takes:
            sa = take.get("sa") or ""
            has_clip = _cigar_has_softclip(take.get("cigar"))
            if not sa and not has_clip:
                continue
            hit = parse_split(
                chrom=take["rname"],
                start=int(take["start"]),
                mapq=int(take["mapq"]),
                cigar=take["cigar"],
                reverse=bool(take["reverse"]),
                sa=sa,
                cfg=cfg,
            )
            if hit is not None:
                splits.append(replace(hit, qname_hash=_qhash(take["qname"])))
                continue
            if not has_clip:
                continue
            clips.extend(
                softclips_from_alignment(
                    start=int(take["start"]),
                    mapq=int(take["mapq"]),
                    cigar=take["cigar"],
                    qual=take.get("qual") or "",
                    cfg=cfg,
                )
            )
        del takes
    extras = clip_breakpoints(clips, cfg) + split_breakpoints(splits)
    chrom_name = str(raw.get("chrom_name") or chrom or "")
    chrom_end = int(raw.get("chrom_end") or 0)
    n_aln = int(raw.get("n_aln", 0))
    n_insert_stat = int(raw.get("n_insert_stat", 0))
    # M4: keep only softclip positions for stamp/block-prep; drop Softclip objects.
    import gc

    import numpy as np

    softclip_pos = np.asarray([c.pos for c in clips], dtype=np.int64)
    del clips
    del raw
    gc.collect()
    try:
        import ctypes

        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass
    if progress is not None:
        progress.stage(
            "load_bam done",
            f"{n_aln:,} alignments -> {len(pairs):,} pairs, "
            f"{len(extras)} breakpoints, {len(splits)} splits "
            f"(insert_stats={n_insert_stat:,}, engine=rust)",
        )
    return pairs, extras, splits, chrom_name, chrom_end, softclip_pos


def load_region(
    path: str,
    chrom: str | None = None,
    cfg: CallConfig | None = None,
    start: int | None = None,
    end: int | None = None,
    reference: str | None = None,
    progress: Progress | None = None,
) -> tuple[list[PairRecord] | Pairs, list[Breakpoint], list[SplitHit], str, int, list]:
    """Pairs, clip/split breakpoints, split hits, chromosome, and max aligned end.

    SAM text is read directly. BAM and CRAM go through pysam and the file index.
    Records are consumed as the file is read. Only the mate that has not yet
    arrived is kept, so a coordinate-sorted SAM does not hold the whole file.
    ``start`` and ``end`` are a 0-based half-open window. A mapped read whose
    mate is unmapped is a hanging pair. A mate on another chromosome, or
    farther than ``max_pair_distance``, is a distant hanging pair
    (orig:sv_parse_reads.py). Soft clips are skipped on reads that already
    produced a split (orig:process_softclip).

    When ``cfg.columnar_pairs`` is True, pairs are stored as a numpy SoA
    (``Pairs``) instead of a list of ``PairRecord`` objects.

    When ``cfg.engine == "rust"`` and the input is BAM/CRAM, the Rust loader
    fills the SoA; softclip/split discovery reuses the Python evidence helpers
    on the same take-record stream. Missing ``arcsv2_native`` raises — only
    ``--engine python`` uses the pysam path.
    """
    cfg = cfg or CallConfig()
    lower = path.lower()
    if (
        cfg.engine == "rust"
        and (lower.endswith(".bam") or lower.endswith(".cram"))
        and cfg.columnar_pairs
    ):
        from arcsv2 import native as native_mod

        if not native_mod.rust_available():
            raise RuntimeError(
                "engine=rust requires arcsv2_native (rebuild arcsv2:dev or pass --engine python)"
            )
        return _load_region_rust(path, chrom, cfg, start, end, reference, progress)

    use_soa = bool(cfg.columnar_pairs)
    pending: dict[str, SamRead] = {}
    pairs_list: list[PairRecord] = []
    pairs_soa: Pairs | None = Pairs(capacity=1_000_000) if use_soa else None
    splits: list[SplitHit] = []
    clips = []
    chrom_name = chrom or ""
    chrom_end = 0
    n_aln = 0
    # Insert lengths counted for library stats even when a pair is otherwise
    # uninteresting (proper, no clip, no SA). PairRecords are still emitted —
    # build_graph needs concordant pairs for adjacency support.
    n_insert_stat = 0
    last_report = 0
    if progress is not None:
        progress.stage("load_bam", path)

    def n_pairs() -> int:
        return len(pairs_soa) if pairs_soa is not None else len(pairs_list)

    def take(read: SamRead) -> None:
        nonlocal chrom_name, chrom_end
        chrom_name = chrom_name or read.rname
        chrom_end = max(chrom_end, _aligned_end(read))
        sa = read.tags.get("SA", "")
        has_clip = _cigar_has_softclip(read.cigar)
        # Cheap path: no SA and no soft clip → nothing for split/clip discovery.
        if not sa and not has_clip:
            return
        hit = parse_split(
            chrom=read.rname,
            start=read.start,
            mapq=read.mapq,
            cigar=read.cigar,
            reverse=bool(read.flag & FLAG_REVERSE),
            sa=sa,
            cfg=cfg,
        )
        if hit is not None:
            splits.append(replace(hit, qname_hash=_qhash(read.qname)))
            return
        if not has_clip:
            return
        clips.extend(
            softclips_from_alignment(
                start=read.start,
                mapq=read.mapq,
                cigar=read.cigar,
                qual=read.qual,
                cfg=cfg,
            )
        )

    def push_pair(**fields) -> None:
        if pairs_soa is not None:
            pairs_soa.append(**fields)
        else:
            pairs_list.append(PairRecord(**fields))

    def emit_pair(a: SamRead, b: SamRead) -> None:
        nonlocal n_insert_stat
        if a.flag & FLAG_REVERSE and not (b.flag & FLAG_REVERSE):
            a, b = b, a
        ilen = insert_length(a, b)
        # Always count usable inserts for library / coverage stats.
        if ilen is not None and ilen >= 0:
            n_insert_stat += 1
        # C5 optional smoove-style filter (default OFF): drop both-clipped pairs
        # from evidence, but keep the insert in the library sample above.
        if cfg.filter_both_mates_clipped and _cigar_has_softclip(a.cigar) and _cigar_has_softclip(
            b.cigar
        ):
            return
        push_pair(
            start1=a.start,
            end1=a.start + aligned_length(a.cigar),
            start2=b.start,
            end2=b.start + aligned_length(b.cigar),
            reverse1=bool(a.flag & FLAG_REVERSE),
            reverse2=bool(b.flag & FLAG_REVERSE),
            mapq1=a.mapq,
            mapq2=b.mapq,
            # Same-orient pairs have no signed insert; keep 0 for Inv DRPs.
            insert=0 if ilen is None else ilen,
            lib_id=_lib_id(a.tags),
            hanging=False,
            qname_hash=_qhash(a.qname),
        )
        take(a)
        take(b)

    def emit_joined(a: SamRead, b: SamRead) -> None:
        a_ok = _on_target(a, chrom, start, end)
        b_ok = _on_target(b, chrom, start, end)
        if not a_ok and not b_ok:
            return
        a_unmapped = bool(a.flag & FLAG_UNMAPPED)
        b_unmapped = bool(b.flag & FLAG_UNMAPPED)
        if a_unmapped and b_unmapped:
            return
        if a_unmapped or b_unmapped:
            anchor, mate = (b, a) if a_unmapped else (a, b)
            if _on_target(anchor, chrom, start, end):
                push_pair(**_hanging_fields(anchor, mate, distant=False))
                take(anchor)
            return
        distant = a.rname != b.rname or abs(a.start - b.start) > cfg.max_pair_distance
        if distant:
            if a_ok:
                push_pair(**_hanging_fields(a, b, distant=True))
                take(a)
            if b_ok:
                push_pair(**_hanging_fields(b, a, distant=True))
                take(b)
            return
        emit_pair(a, b)

    bam_threads = max(1, int(cfg.threads))
    max_nm = int(cfg.max_nm) if cfg.filter_high_nm else None
    for read in _iter_alignments(
        path, chrom, start, end, reference, threads=bam_threads, max_nm=max_nm
    ):
        # BAM path already dropped skip-flags; SAM text still needs the check.
        if not _primary(read):
            continue
        if max_nm is not None:
            nm_tag = read.tags.get("NM")
            if nm_tag is not None and int(float(nm_tag)) > max_nm:
                continue
        n_aln += 1
        if progress is not None and n_aln - last_report >= 500_000:
            progress.stage(
                "load_bam",
                f"{n_aln:,} alignments, {n_pairs():,} pairs "
                f"(insert_stats={n_insert_stat:,}, pos ~{read.start:,})",
            )
            last_report = n_aln
        earlier = pending.pop(read.qname, None)
        if earlier is None:
            pending[read.qname] = read
            continue
        emit_joined(earlier, read)
    for read in list(pending.values()):
        if not _on_target(read, chrom, start, end):
            continue
        take(read)
        if not (read.flag & FLAG_MATE_UNMAPPED):
            continue
        push_pair(
            start1=read.start,
            end1=read.start + aligned_length(read.cigar),
            start2=0,
            end2=0,
            reverse1=bool(read.flag & FLAG_REVERSE),
            reverse2=False,
            mapq1=read.mapq,
            mapq2=0,
            insert=0,
            lib_id=_lib_id(read.tags),
            hanging=True,
            qname_hash=_qhash(read.qname),
        )
    # Free BAM-side mate buffer after pairs are complete.
    pending.clear()
    del pending

    if pairs_soa is not None:
        pairs_soa.trim()
        pairs_out: list[PairRecord] | Pairs = pairs_soa
    else:
        pairs_out = pairs_list

    extras = clip_breakpoints(clips, cfg) + split_breakpoints(splits)
    import numpy as np

    softclip_pos = np.asarray([c.pos for c in clips], dtype=np.int64)
    del clips
    if progress is not None:
        progress.stage(
            "load_bam done",
            f"{n_aln:,} alignments -> {n_pairs():,} pairs, "
            f"{len(extras)} breakpoints, {len(splits)} splits "
            f"(insert_stats={n_insert_stat:,})",
        )
    return pairs_out, extras, splits, chrom_name, chrom_end, softclip_pos


def pairs_from_sam(path: str, chrom: str | None = None) -> tuple[list[PairRecord] | Pairs, str, int]:
    """Pair primary alignments by query name. Returns pairs, chrom, max end."""
    pairs, _extras, _splits, chrom_name, chrom_end, _clips = load_region(path, chrom)
    return pairs, chrom_name, chrom_end
