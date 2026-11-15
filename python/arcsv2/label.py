"""Assembly validation labels (§4.10).

A call is validated when one minimap2 alignment (default parameters, no preset)
covers the whole query, has MAPQ ≥ 60, BLAST identity ``mlen/blen`` > 0.90,
and no CIGAR indel of 49 bp or more. ``_`` and ``=`` tokens have no sequence
to validate and are labeled missing.

``label_version`` is ``arcsv2-label-1``. Sequences use ``blocks_point`` when
the caller wrote it, so blocks meet at the clip or the interval midpoint, and
the outer blocks are cut to 150 bp next to the event.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

LABEL_VERSION = "arcsv2-label-1"

# Haplotypes of megabase-scale multi-block loci can build multi-Mb ref queries.
# Check the span before fetching/aligning; mappy otherwise hangs for hours.
MAX_LABEL_QUERY_LEN = 1_000_000
# Wall-clock cap per mappy map of one query against one haplotype (SIGALRM).
LABEL_MAP_TIMEOUT_S = 60.0

LABEL_COLUMNS = (
    "sample",
    "label",
    "reason",
    "ref_validated",
    "gt_concordant",
    "best_hap",
    "identity",
    "mapq",
    "max_indel",
    "strand",
    "label_bp_shift",
    "label_version",
    "feature_schema",
    "sv_class",
    "group",
)

TOKEN_RE = re.compile(r"([A-Za-z]\d*)('?)|(_)|(=)")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class LabelResult:
    label: float  # 1, 0, or nan
    reason: str
    ref_validated: bool = False
    gt_concordant: str = "NA"
    best_hap: str = "NA"
    identity: str = "NA"
    mapq: str = "NA"
    max_indel: str = "NA"
    strand: str = "NA"
    label_bp_shift: str = "0"


def largest_indel(cigar: str) -> int:
    biggest = 0
    for count, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar or ""):
        if op in ("I", "D"):
            biggest = max(biggest, int(count))
    return biggest


def is_validated(
    *,
    q_start: int,
    q_end: int,
    q_len: int,
    mapq: int,
    mlen: int,
    blen: int,
    cigar: str,
    mapq_min: int = 60,
    identity_min: float = 0.90,
    max_gap: int = 48,
) -> tuple[bool, str]:
    """Paper acceptance test for one alignment (§4.10)."""
    if q_len <= 0 or q_end <= q_start:
        return False, "no_hit"
    if q_start != 0 or q_end != q_len:
        return False, "clipped"
    if mapq < mapq_min:
        return False, "low_mapq"
    if blen <= 0 or (mlen / blen) <= identity_min:
        return False, "low_identity"
    if largest_indel(cigar) >= max_gap + 1:
        return False, "large_gap"
    return True, "validated"


def sequence_from_call(
    reference: str,
    rearrangement: str,
    block_seqs: dict[str, str],
    flank_left: str,
    flank_right: str,
) -> str | None:
    """Build the predicted haplotype plus flanks. None for ``_`` or ``=``."""
    pieces = [flank_left]
    for match in TOKEN_RE.finditer(rearrangement or ""):
        if match.group(3) or match.group(4):
            return None
        name = match.group(1)
        seq = block_seqs.get(name)
        if seq is None:
            return None
        if match.group(2):
            seq = seq.translate(COMP)[::-1]
        pieces.append(seq)
    pieces.append(flank_right)
    return "".join(pieces)


def label_alignment_hit(**kwargs) -> LabelResult:
    ok, reason = is_validated(**kwargs)
    if reason == "no_hit":
        return LabelResult(float("nan") if kwargs.get("novel") else 0.0, reason)
    return LabelResult(1.0 if ok else 0.0, reason)


def read_fasta(path: str | Path) -> dict[str, str]:
    """Contig name to sequence. Wrapped FASTA lines are joined."""
    seqs: dict[str, list[str]] = {}
    name = ""
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {key: "".join(parts) for key, parts in seqs.items()}


def blocks_from_row(row: dict[str, str], column: str = "blocks") -> dict[str, tuple[int, int]]:
    """Map reference block names to 0-based half-open intervals.

    The ``blocks`` column is the caller's locus list. ``blocks_point`` is the
    same list closed at each breakpoint. Names come from ``reference`` in that
    same order. ``_`` entries are not reference sequence.
    """
    from arcsv2.engine.names import block_idx_to_name

    intervals: list[tuple[int, int]] = []
    for part in (row.get(column) or "").split(","):
        part = part.strip()
        if not part or part == "_":
            continue
        left, right = part.split("-", 1)
        intervals.append((int(left), int(right)))
    names = [match.group(1) for match in TOKEN_RE.finditer(row.get("reference") or "") if match.group(1)]
    if len(names) != len(intervals):
        names = [block_idx_to_name(i) for i in range(len(intervals))]
    return {name: interval for name, interval in zip(names, intervals)}


class ReferenceView:
    """One open reference. An indexed FASTA is fetched by interval.

    ``pysam.FastaFile`` reads ``chrom:start-end`` from the ``.fai`` index.
    Without pysam, or without an index, the file is loaded once.
    """

    def __init__(self, fasta: dict[str, str] | str | Path | ReferenceView) -> None:
        self._dict: dict[str, str] | None = None
        self._file = None
        if isinstance(fasta, ReferenceView):
            self._dict = fasta._dict
            self._file = fasta._file
            self._owned = False
            return
        self._owned = False
        if isinstance(fasta, dict):
            self._dict = fasta
            return
        path = Path(fasta)
        try:
            import pysam
        except ImportError:
            pysam = None
        if pysam is not None:
            try:
                self._file = pysam.FastaFile(str(path))
                self._owned = True
                return
            except (OSError, ValueError):
                self._file = None
        self._dict = read_fasta(path)

    def contigs(self) -> list[str]:
        if self._file is not None:
            return list(self._file.references)
        return list(self._dict or ())

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open slice. Coordinates are clipped to the contig."""
        start = max(0, start)
        if self._file is not None:
            length = int(self._file.get_reference_length(chrom))
            end = min(end, length)
            if start >= end:
                return ""
            return self._file.fetch(chrom, start, end)
        seq = (self._dict or {})[chrom]
        return seq[start:end]

    def close(self) -> None:
        if self._owned and self._file is not None:
            self._file.close()
            self._owned = False


def _reference_block_spans(column: str) -> list[tuple[int, int]]:
    spans: list[tuple[int, int]] = []
    for part in column.split(","):
        part = part.strip()
        if not part or part == "_":
            continue
        left, right = part.split("-", 1)
        spans.append((int(left), int(right)))
    return spans


def _internal_junctions(row: dict[str, str]) -> list[tuple[int, int, int, str]]:
    """Each internal junction: uncovered ``(lo, hi)``, current ``point``, ``clip`` or ``mid``."""
    raw = _reference_block_spans(row.get("blocks") or "")
    point = _reference_block_spans(row.get("blocks_point") or "")
    if len(raw) != len(point) or len(raw) < 2:
        return []
    sources = [part.strip() for part in (row.get("bp_point_source") or "NA").split(",")]
    out: list[tuple[int, int, int, str]] = []
    for index in range(len(raw) - 1):
        lo, hi = raw[index][1], raw[index + 1][0]
        junction = point[index][1]
        source = sources[index] if index < len(sources) else "mid"
        if source not in ("clip", "mid"):
            source = "mid"
        out.append((lo, hi, junction, source))
    return out


def _midpoint_candidates(lo: int, hi: int) -> list[int]:
    if hi <= lo:
        return [lo]
    cands = {lo, hi - 1, (lo + hi) // 2}
    for pos in range(lo, hi, 10):
        cands.add(pos)
    return sorted(cands)


def _blocks_point_from_junctions(row: dict[str, str], junctions: list[int]) -> str:
    tokens = [part.strip() for part in (row.get("blocks_point") or "").split(",") if part.strip()]
    parsed: list[tuple[int, int] | str] = []
    for token in tokens:
        if token == "_":
            parsed.append("_")
        else:
            left, right = token.split("-", 1)
            parsed.append((int(left), int(right)))
    ref = [item for item in parsed if item != "_"]
    if len(ref) < 2 or len(junctions) != len(ref) - 1:
        return row.get("blocks_point") or ""
    rebuilt: list[tuple[int, int]] = []
    for index, (left, right) in enumerate(ref):
        edge_left = left if index == 0 else junctions[index - 1]
        edge_right = right if index == len(ref) - 1 else junctions[index]
        if edge_right < edge_left:
            edge_right = edge_left
        rebuilt.append((edge_left, edge_right))
    parts: list[str] = []
    ref_index = 0
    for item in parsed:
        if item == "_":
            parts.append("_")
        else:
            left, right = rebuilt[ref_index]
            parts.append(f"{left}-{right}")
            ref_index += 1
    return ",".join(parts)


def _prepared_blocks(
    row: dict[str, str], flank: int
) -> tuple[dict[str, tuple[int, int]], int, int]:
    """Point-estimate blocks, with outer blocks cut to ``flank`` bp.

    A shorter outer block is kept whole. The returned pads are the extra
    reference bases that bring that side up to ``flank``.
    """
    point = (row.get("blocks_point") or "").strip()
    column = "blocks_point" if point and point != "NA" else "blocks"
    named = blocks_from_row(row, column)
    ordered = list(named)
    if len(ordered) < 2 or flank <= 0:
        return named, 0, 0
    trimmed = dict(named)
    left_name = ordered[0]
    left_start, left_end = trimmed[left_name]
    if left_end - left_start > flank:
        trimmed[left_name] = (left_end - flank, left_end)
        pad_left = 0
    else:
        pad_left = flank - (left_end - left_start)
    right_name = ordered[-1]
    right_start, right_end = trimmed[right_name]
    if right_end - right_start > flank:
        trimmed[right_name] = (right_start, right_start + flank)
        pad_right = 0
    else:
        pad_right = flank - (right_end - right_start)
    return trimmed, pad_left, pad_right


def _estimated_allele_lens(row: dict[str, str], flank: int) -> tuple[int | None, int]:
    """``(alt_len or None, ref_len)`` from block spans without reading the FASTA.

    ``None`` alt means the rearrangement has no sequence (``_`` / ``=``).
    """
    named, pad_left, pad_right = _prepared_blocks(row, flank)
    if not named:
        return None, 0
    lengths = {name: end - start for name, (start, end) in named.items()}
    ordered = list(named)
    ref_len = pad_left + sum(lengths[name] for name in ordered) + pad_right
    rearrangement = row.get("rearrangement") or ""
    if "_" in rearrangement or "=" in rearrangement:
        return None, ref_len
    alt_len = pad_left + pad_right
    for match in TOKEN_RE.finditer(rearrangement):
        if match.group(3) or match.group(4):
            return None, ref_len
        name = match.group(1)
        if name not in lengths:
            return None, ref_len
        alt_len += lengths[name]
    return alt_len, ref_len


def _map_with_timeout(aligner: object, seq: str, seconds: float) -> list[object] | None:
    """Run ``aligner.map(seq)`` under ``SIGALRM``. ``None`` means the timer fired.

    Only available on Unix (the Linux label container). Without ``SIGALRM``,
    the call runs uncapped and the length guard is the primary protection.
    """
    try:
        import signal
    except ImportError:
        return list(aligner.map(seq))
    if not hasattr(signal, "SIGALRM") or not hasattr(signal, "setitimer"):
        return list(aligner.map(seq))

    class _TimedOut(Exception):
        pass

    def _handler(_signum: int, _frame: object) -> None:
        raise _TimedOut()

    previous = signal.signal(signal.SIGALRM, _handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return list(aligner.map(seq))
    except _TimedOut:
        return None
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, previous)


def haplotype_sequences(
    row: dict[str, str],
    fasta: dict[str, str] | str | Path | ReferenceView,
    flank: int = 150,
) -> tuple[str | None, str]:
    """``(alt, ref)`` for one call. ``alt`` is None when the call has ``_`` or ``=``.

    Sequences come from ``blocks_point`` and keep ``flank`` bp of each outer
    block next to the event. Only that window is read. Pass a ``ReferenceView``
    when labeling many rows so the reference is opened once.
    """
    view = fasta if isinstance(fasta, ReferenceView) else ReferenceView(fasta)
    owned = not isinstance(fasta, ReferenceView)
    try:
        named, pad_left, pad_right = _prepared_blocks(row, flank)
        if not named:
            return None, ""
        ordered = list(named)
        chrom = row.get("chrom") or view.contigs()[0]
        first = named[ordered[0]][0]
        last = named[ordered[-1]][1]
        window_start = max(0, first - pad_left)
        window = view.fetch(chrom, window_start, last + pad_right)

        def cut(start: int, end: int) -> str:
            return window[max(0, start - window_start) : max(0, end - window_start)]

        block_seqs = {name: cut(start, end) for name, (start, end) in named.items()}
        flank_left = cut(window_start, first)
        flank_right = cut(last, last + pad_right)
        alt = sequence_from_call(
            row.get("reference") or "",
            row.get("rearrangement") or "",
            block_seqs,
            flank_left,
            flank_right,
        )
        ref = flank_left + "".join(block_seqs[name] for name in ordered) + flank_right
        return alt, ref
    finally:
        if owned:
            view.close()


@dataclass
class Hit:
    hap: str
    q_start: int
    q_end: int
    q_len: int
    mapq: int
    mlen: int
    blen: int
    cigar: str
    strand: int


def _cigar_text(raw) -> str:
    text = getattr(raw, "cigar_str", None) or ""
    if text:
        return text
    ops = "MIDNSHP=X"
    parts: list[str] = []
    for op, length in getattr(raw, "cigar", ()) or ():
        parts.append(f"{int(length)}{ops[int(op)]}")
    return "".join(parts)


class HapAligner:
    """Two minimap2 indexes, default parameters, no preset.

    ``sequential=True`` keeps only one haplotype index in memory at a time.
    ``label_rows`` then maps every query against hap1, frees it, then hap2.
    """

    def __init__(
        self,
        hap1: str | Path,
        hap2: str | Path,
        hap1_index: str | Path | None = None,
        hap2_index: str | Path | None = None,
        threads: int = 1,
        sequential: bool = False,
    ) -> None:
        try:
            import mappy
        except ImportError as exc:
            raise RuntimeError(
                "Haplotype labeling needs mappy. There is no native Windows build; "
                "use the Linux container in docs/CROSS_PLATFORM.md (pip install arcsv2[label])."
            ) from exc
        self._mappy = mappy
        self._threads = threads
        self.sequential = sequential
        self._paths = {
            "hap1": str(hap1_index or hap1),
            "hap2": str(hap2_index or hap2),
        }
        self.aligners: dict[str, object] = {}
        if not sequential:
            for hap, path in self._paths.items():
                aligner = mappy.Aligner(path, n_threads=threads)
                if not aligner:
                    raise RuntimeError(f"failed to load haplotype index: {hap}")
                self.aligners[hap] = aligner

    def ensure(self, hap: str) -> None:
        """Load ``hap`` and drop any other haplotype when running sequentially."""
        if hap in self.aligners:
            return
        if self.sequential:
            self.aligners.clear()
        aligner = self._mappy.Aligner(self._paths[hap], n_threads=self._threads)
        if not aligner:
            raise RuntimeError(f"failed to load haplotype index: {hap}")
        self.aligners[hap] = aligner

    def map_loaded(self, hap: str, seq: str) -> list[Hit]:
        """Map ``seq`` against an already-loaded haplotype. Does not switch indexes.

        Raises ``TimeoutError`` when the map exceeds ``LABEL_MAP_TIMEOUT_S``.
        """
        aligner = self.aligners[hap]
        found: list[Hit] = []
        if not seq:
            return found
        raw_hits = _map_with_timeout(aligner, seq, LABEL_MAP_TIMEOUT_S)
        if raw_hits is None:
            raise TimeoutError(f"mappy map exceeded {LABEL_MAP_TIMEOUT_S:.0f}s on {hap}")
        for raw in raw_hits:
            if hasattr(raw, "is_primary") and not raw.is_primary:
                continue
            found.append(
                Hit(
                    hap=hap,
                    q_start=int(raw.q_st),
                    q_end=int(raw.q_en),
                    q_len=len(seq),
                    mapq=int(raw.mapq),
                    mlen=int(raw.mlen),
                    blen=int(raw.blen),
                    cigar=_cigar_text(raw),
                    strand=int(raw.strand),
                )
            )
        return found

    def hits_on(self, hap: str, seq: str) -> list[Hit]:
        """Primary alignments of ``seq`` to one haplotype. Loads that index if needed."""
        self.ensure(hap)
        return self.map_loaded(hap, seq)

    def hits(self, seq: str) -> list[Hit]:
        """Primary alignments of ``seq`` to each haplotype (joint mode)."""
        found: list[Hit] = []
        for hap in ("hap1", "hap2"):
            found.extend(self.hits_on(hap, seq))
        return found

    def best_hit(self, seq: str) -> Hit | None:
        """Highest-MAPQ primary hit of ``seq`` on either haplotype."""
        return _best(self.hits(seq))

    def close(self) -> None:
        self.aligners.clear()


def _passes(hit: Hit | None) -> tuple[bool, str]:
    """Paper criteria on either strand. Assembly contigs are not oriented to the reference."""
    if hit is None:
        return False, "no_hit"
    return is_validated(
        q_start=hit.q_start,
        q_end=hit.q_end,
        q_len=hit.q_len,
        mapq=hit.mapq,
        mlen=hit.mlen,
        blen=hit.blen,
        cigar=hit.cigar,
    )


def _best(hits: list[Hit]) -> Hit | None:
    if not hits:
        return None
    return max(hits, key=lambda hit: (hit.mapq, hit.mlen))


def _window_unreliable(
    row: dict[str, str],
    named: dict[str, tuple[int, int]],
    flank: int,
    unreliable: dict[str, list[tuple[int, int]]] | None,
) -> bool:
    if not unreliable:
        return False
    intervals = unreliable.get(row.get("chrom") or "", [])
    if not intervals or not named:
        return False
    start = max(0, min(left for left, _right in named.values()) - flank)
    end = max(right for _left, right in named.values()) + flank
    from arcsv2.merge_filter import _overlap_bp

    return _overlap_bp(start, end, intervals) > 0


def label_from_passes(
    alt_ok: dict[str, bool],
    ref_ok: dict[str, bool],
) -> tuple[float, str, bool]:
    """Label from per-haplotype pass flags.

    Ambiguous means alt and ref both pass on the same haplotype. Alt on one
    haplotype and ref on the other is a validated heterozygote.
    """
    haps = ("hap1", "hap2")
    ref_validated = any(bool(ref_ok.get(hap)) for hap in haps)
    if any(alt_ok.get(hap) and ref_ok.get(hap) for hap in haps):
        return float("nan"), "ambiguous", ref_validated
    if any(alt_ok.get(hap) for hap in haps):
        return 1.0, "validated", ref_validated
    return 0.0, "", ref_validated


def _group_hits(hits: list[Hit]) -> dict[str, list[Hit]]:
    grouped: dict[str, list[Hit]] = {"hap1": [], "hap2": []}
    for hit in hits:
        grouped.setdefault(hit.hap, []).append(hit)
    return grouped


def decide_from_hits(
    row: dict[str, str],
    alt_hits: list[Hit],
    ref_hits: list[Hit],
) -> LabelResult:
    """Label one call from precomputed alt and ref hits on both haplotypes."""
    alt_groups = _group_hits(alt_hits)
    ref_groups = _group_hits(ref_hits)
    alt_ok: dict[str, bool] = {}
    ref_ok: dict[str, bool] = {}
    failures: list[tuple[int, str]] = []
    passing_hits: list[Hit] = []
    for hap in ("hap1", "hap2"):
        best_alt = _best(alt_groups.get(hap, []))
        ok, reason = _passes(best_alt)
        alt_ok[hap] = ok
        if ok and best_alt is not None:
            passing_hits.append(best_alt)
        elif not ok:
            score = best_alt.mapq if best_alt is not None else -1
            failures.append((score, reason))
        best_ref = _best(ref_groups.get(hap, []))
        ref_hit_ok, _reason = _passes(best_ref)
        ref_ok[hap] = ref_hit_ok
    label, reason, ref_passed = label_from_passes(alt_ok, ref_ok)
    if label == 0.0:
        reason = max(failures)[1] if failures else "no_hit"
    best = _best(passing_hits) or _best([hit for group in alt_groups.values() for hit in group])
    identity = "NA"
    mapq = "NA"
    max_indel = "NA"
    best_hap = "NA"
    strand = "NA"
    if best is not None:
        best_hap = best.hap
        mapq = str(best.mapq)
        max_indel = str(largest_indel(best.cigar))
        strand = str(best.strand)
        if best.blen > 0:
            identity = f"{best.mlen / best.blen:.4f}"
    n_pass = sum(1 for hap in ("hap1", "hap2") if alt_ok.get(hap))
    gt = (row.get("gt") or "").upper()
    if gt == "HOM":
        concordant = n_pass == 2
    elif gt == "HET":
        concordant = n_pass == 1
    else:
        concordant = False
    gt_text = "1" if concordant else "0"
    return LabelResult(
        label, reason, ref_passed, gt_text, best_hap, identity, mapq, max_indel, strand
    )


def queries_for(
    row: dict[str, str],
    fasta: dict[str, str] | str | Path | ReferenceView,
    unreliable: dict[str, list[tuple[int, int]]] | None = None,
    flank: int = 150,
    *,
    include_shifts: bool = False,
) -> tuple[LabelResult | None, dict[str, str]]:
    """Early NA result, or ``{query_id: sequence}`` for alt/ref (and mid-shifts if asked).

    Query ids are ``alt``, ``ref``, and optionally ``shift:<j>:<pos>:alt`` / ``...:ref``.
    Shift queries are built only when ``include_shifts`` is true (after a failed base
    call); mapping every mid grid up front is expensive on real assemblies.
    """
    rearrangement = row.get("rearrangement") or ""
    if "_" in rearrangement or "=" in rearrangement:
        return LabelResult(float("nan"), "no_sequence"), {}
    named, _pad_left, _pad_right = _prepared_blocks(row, flank)
    if _window_unreliable(row, named, flank, unreliable):
        return LabelResult(float("nan"), "unreliable_assembly"), {}
    alt_est, ref_est = _estimated_allele_lens(row, flank)
    if alt_est is None:
        return LabelResult(float("nan"), "no_sequence"), {}
    if alt_est > MAX_LABEL_QUERY_LEN or ref_est > MAX_LABEL_QUERY_LEN:
        return LabelResult(float("nan"), "query_too_long"), {}
    alt, ref = haplotype_sequences(row, fasta, flank)
    if alt is None:
        return LabelResult(float("nan"), "no_sequence"), {}
    queries = {"alt": alt, "ref": ref}
    if include_shifts:
        queries.update(shift_queries_for(row, fasta, flank))
    return None, queries


def shift_queries_for(
    row: dict[str, str],
    fasta: dict[str, str] | str | Path | ReferenceView,
    flank: int = 150,
) -> dict[str, str]:
    """Mid-junction shift trials for a call that failed at the base point."""
    queries: dict[str, str] = {}
    junctions = _internal_junctions(row)
    points = [item[2] for item in junctions]
    mid_indexes = [index for index, item in enumerate(junctions) if item[3] == "mid"][:3]
    for j_idx in mid_indexes:
        lo, hi, original, _source = junctions[j_idx]
        for candidate in _midpoint_candidates(lo, hi):
            if candidate == original:
                continue
            trial_points = list(points)
            trial_points[j_idx] = candidate
            trial_row = dict(row)
            trial_row["blocks_point"] = _blocks_point_from_junctions(row, trial_points)
            alt_est, ref_est = _estimated_allele_lens(trial_row, flank)
            if (
                alt_est is None
                or alt_est > MAX_LABEL_QUERY_LEN
                or ref_est > MAX_LABEL_QUERY_LEN
            ):
                continue
            trial_alt, trial_ref = haplotype_sequences(trial_row, fasta, flank)
            if trial_alt is None:
                continue
            key = f"shift:{j_idx}:{candidate}"
            queries[f"{key}:alt"] = trial_alt
            queries[f"{key}:ref"] = trial_ref
    return queries


def row_needs_mid_shifts(row: dict[str, str]) -> bool:
    """True when the call has at least one ``mid`` internal junction to try."""
    return any(item[3] == "mid" for item in _internal_junctions(row)[:])


def decide_row_from_query_hits(
    row: dict[str, str],
    queries: dict[str, str],
    hits_by_query: dict[str, list[Hit]],
) -> LabelResult:
    """Pick the base label, or the first mid-shift that validates."""
    base = decide_from_hits(row, hits_by_query.get("alt", []), hits_by_query.get("ref", []))
    if base.label == 1.0 or math.isnan(base.label):
        base.label_bp_shift = "0"
        return base
    junctions = _internal_junctions(row)
    mid_indexes = [index for index, item in enumerate(junctions) if item[3] == "mid"][:3]
    for j_idx in mid_indexes:
        lo, hi, original, _source = junctions[j_idx]
        for candidate in _midpoint_candidates(lo, hi):
            if candidate == original:
                continue
            key = f"shift:{j_idx}:{candidate}"
            if f"{key}:alt" not in queries:
                continue
            trial = decide_from_hits(
                row,
                hits_by_query.get(f"{key}:alt", []),
                hits_by_query.get(f"{key}:ref", []),
            )
            if trial.label == 1.0:
                trial.label_bp_shift = str(candidate - original)
                return trial
    base.label_bp_shift = "0"
    return base


def label_call(
    row: dict[str, str],
    fasta: dict[str, str] | str | Path | ReferenceView,
    aligner: HapAligner | None,
    unreliable: dict[str, list[tuple[int, int]]] | None = None,
    flank: int = 150,
) -> LabelResult:
    """Section 3 of ``docs/LABELS_PLAN.md``. Genotype concordance is separate."""
    early, queries = queries_for(row, fasta, unreliable, flank, include_shifts=False)
    if early is not None:
        return early
    if aligner is None:
        raise ValueError("an aligner is required for a sequence call")
    try:
        hits_by_query = {qid: aligner.hits(seq) for qid, seq in queries.items()}
    except TimeoutError:
        return LabelResult(float("nan"), "timeout")
    base = decide_from_hits(row, hits_by_query.get("alt", []), hits_by_query.get("ref", []))
    if base.label == 1.0 or math.isnan(base.label) or not row_needs_mid_shifts(row):
        base.label_bp_shift = "0"
        return base
    shifts = shift_queries_for(row, fasta, flank)
    try:
        for qid, seq in shifts.items():
            hits_by_query[qid] = aligner.hits(seq)
    except TimeoutError:
        return LabelResult(float("nan"), "timeout")
    queries.update(shifts)
    return decide_row_from_query_hits(row, queries, hits_by_query)


def _label_fields(result: LabelResult, sample: str) -> dict[str, str]:
    from arcsv2.config import FEATURE_SCHEMA
    from arcsv2.train import class_of_row

    if math.isnan(result.label):
        label = "NA"
    else:
        label = str(int(result.label))
    return {
        "sample": sample,
        "label": label,
        "reason": result.reason,
        "ref_validated": "1" if result.ref_validated else "0",
        "gt_concordant": result.gt_concordant,
        "best_hap": result.best_hap,
        "identity": result.identity,
        "mapq": result.mapq,
        "max_indel": result.max_indel,
        "strand": result.strand,
        "label_bp_shift": result.label_bp_shift,
        "label_version": LABEL_VERSION,
        "feature_schema": FEATURE_SCHEMA,
        "sv_class": class_of_row({"svtype": "", "rearrangement": ""}) or "",
        "group": sample,
    }


def label_rows(
    rows: list[dict[str, str]],
    fasta: dict[str, str] | str | Path | ReferenceView,
    aligner: HapAligner | None,
    sample: str,
    unreliable: dict[str, list[tuple[int, int]]] | None = None,
    flank: int = 150,
) -> list[dict[str, str]]:
    """Every input column plus the label columns. The reference is opened once.

    Sequential mode: pass A maps base alt/ref for every call on hap1 then hap2
    (two loads). Pass B maps mid-shift queries only for rows still labeled 0
    (up to two more loads). Joint mode maps shifts only after a failed base call.
    """
    from arcsv2.train import class_of_row

    view = fasta if isinstance(fasta, ReferenceView) else ReferenceView(fasta)
    owned = not isinstance(fasta, ReferenceView)
    labeled: list[dict[str, str]] = []
    try:
        if aligner is not None and aligner.sequential:
            prepared: list[tuple[dict[str, str], LabelResult | None, dict[str, str]]] = []
            for row in rows:
                early, queries = queries_for(
                    row, view, unreliable, flank, include_shifts=False
                )
                prepared.append((row, early, queries))
            hits_by_row: list[dict[str, list[Hit]]] = [{} for _ in prepared]
            results: list[LabelResult | None] = [None] * len(prepared)
            for hap in ("hap1", "hap2"):
                aligner.ensure(hap)
                for index, (_row, early, queries) in enumerate(prepared):
                    if early is not None:
                        continue
                    if results[index] is not None and math.isnan(results[index].label):
                        continue
                    for qid, seq in queries.items():
                        try:
                            hits_by_row[index].setdefault(qid, []).extend(
                                aligner.map_loaded(hap, seq)
                            )
                        except TimeoutError:
                            results[index] = LabelResult(float("nan"), "timeout")
                            break
                aligner.close()
            need_shifts: list[tuple[int, dict[str, str]]] = []
            for index, (row, early, queries) in enumerate(prepared):
                if early is not None:
                    results[index] = early
                    continue
                if results[index] is not None:
                    continue
                base = decide_from_hits(
                    row,
                    hits_by_row[index].get("alt", []),
                    hits_by_row[index].get("ref", []),
                )
                if base.label == 1.0 or math.isnan(base.label) or not row_needs_mid_shifts(row):
                    base.label_bp_shift = "0"
                    results[index] = base
                    continue
                shifts = shift_queries_for(row, view, flank)
                if not shifts:
                    base.label_bp_shift = "0"
                    results[index] = base
                    continue
                queries.update(shifts)
                prepared[index] = (row, early, queries)
                need_shifts.append((index, shifts))
            if need_shifts:
                for hap in ("hap1", "hap2"):
                    aligner.ensure(hap)
                    for index, shifts in need_shifts:
                        if results[index] is not None and math.isnan(results[index].label):
                            continue
                        for qid, seq in shifts.items():
                            try:
                                hits_by_row[index].setdefault(qid, []).extend(
                                    aligner.map_loaded(hap, seq)
                                )
                            except TimeoutError:
                                results[index] = LabelResult(float("nan"), "timeout")
                                break
                    aligner.close()
                for index, _shifts in need_shifts:
                    if results[index] is not None:
                        continue
                    row, _early, queries = prepared[index]
                    results[index] = decide_row_from_query_hits(
                        row, queries, hits_by_row[index]
                    )
            for index, (row, _early, _queries) in enumerate(prepared):
                result = results[index]
                assert result is not None
                updated = dict(row)
                fields = _label_fields(result, sample)
                kind = class_of_row(row)
                fields["sv_class"] = kind or "insertion"
                updated.update(fields)
                labeled.append(updated)
        else:
            for row in rows:
                result = label_call(row, view, aligner, unreliable, flank)
                updated = dict(row)
                fields = _label_fields(result, sample)
                kind = class_of_row(row)
                fields["sv_class"] = kind or "insertion"
                updated.update(fields)
                labeled.append(updated)
    finally:
        if owned:
            view.close()
        if aligner is not None and aligner.sequential:
            aligner.close()
    return labeled


def dataset_tables(tables: list[list[dict[str, str]]]) -> tuple[list[dict[str, str]], dict]:
    """Concatenate labeled tables, drop NA labels, and summarize 1/0/NA.

    Tables must share ``label_version`` and ``feature_schema``.
    """
    versions = {row.get("label_version", "") for table in tables for row in table}
    schemas = {row.get("feature_schema", "") for table in tables for row in table}
    versions.discard("")
    schemas.discard("")
    if len(versions) > 1 or len(schemas) > 1:
        raise ValueError(f"mixed label tables: versions {sorted(versions)} schemas {sorted(schemas)}")
    summary: dict = {
        "label_version": next(iter(versions), LABEL_VERSION),
        "feature_schema": next(iter(schemas), ""),
        "classes": {},
        "samples": {},
    }
    kept: list[dict[str, str]] = []
    for table in tables:
        for row in table:
            kind = row.get("sv_class") or "NA"
            sample = row.get("sample") or row.get("group") or "NA"
            label = row.get("label", "NA")
            for bucket, key in ((summary["classes"], kind), (summary["samples"], sample)):
                slot = bucket.setdefault(key, {"1": 0, "0": 0, "NA": 0, "reasons": {}})
                slot[label if label in ("0", "1") else "NA"] += 1
                reason = row.get("reason") or "NA"
                slot["reasons"][reason] = slot["reasons"].get(reason, 0) + 1
            if label in ("0", "1"):
                out = dict(row)
                out["group"] = sample
                kept.append(out)
    return kept, summary


def write_dataset(rows: list[dict[str, str]], summary: dict, path: str | Path) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with dest.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(col, "NA") for col in columns) + "\n")
    dest.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
