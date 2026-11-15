"""Soft clips and split reads from one alignment (orig:softclip.py, splitreads.py).

Positions are 0-based. Split breakpoints follow orig:parse_splits, including
the forward-segment ``+ 1`` that turns the last aligned base into the same
breakpoint coordinate a discordant pair uses for ``reference_end``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from arcsv2.config import CallConfig
from arcsv2.engine.breakpoints import (
    Breakpoint,
    Softclip,
    classify_split_flags,
    cluster_softclips,
    merge_microhomology,
)

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_LOWQUAL = set("#!\"")


@dataclass(frozen=True)
class SplitHit:
    """One primary alignment with a single same-chromosome supplementary."""

    start1: int
    end1: int
    reverse1: bool
    start2: int
    end2: int
    reverse2: bool
    mapq: int
    split_type: str
    bp1: tuple[int, int]
    bp2: tuple[int, int]
    chrom: str
    # Primary alignment that carried the SA tag. -1 means the hit is not tied to a pair.
    primary_start: int = -1
    primary_reverse: bool = False
    qname_hash: int = -1


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in _CIGAR.findall(cigar or "")]


def aligned_length(cigar: str) -> int:
    return sum(n for n, op in cigar_ops(cigar) if op in "MDN=X")


def _softclip_counts(cigar: str) -> tuple[int, int]:
    """Leading and trailing soft-clip lengths. Hard clips are skipped."""
    ops = cigar_ops(cigar)
    i = 0
    while i < len(ops) and ops[i][1] == "H":
        i += 1
    left = ops[i][0] if i < len(ops) and ops[i][1] == "S" else 0
    j = len(ops) - 1
    while j >= 0 and ops[j][1] == "H":
        j -= 1
    right = ops[j][0] if j >= 0 and ops[j][1] == "S" and j != i else 0
    return left, right


def _lowqual_run(qual: str, extra: int) -> tuple[int, int]:
    """orig:helper.py::count_lowqual_bases. Phred 0–2 are ``!``, ``"``, ``#``."""
    left = 0
    for char in qual:
        if char in _LOWQUAL:
            left += 1
        else:
            break
    right = 0
    for char in reversed(qual):
        if char in _LOWQUAL:
            right += 1
        else:
            break
    if left:
        left += extra
    if right:
        right += extra
    return left, right


def _median(values: list[int]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def softclips_from_alignment(
    *,
    start: int,
    mapq: int,
    cigar: str,
    qual: str,
    cfg: CallConfig,
) -> list[Softclip]:
    """One alignment's left and right clips. orig:process_softclip for one read."""
    if mapq < cfg.min_mapq_softclip:
        return []
    left_n, right_n = _softclip_counts(cigar)
    if left_n == 0 and right_n == 0:
        return []
    low_l, low_r = _lowqual_run(qual or "", cfg.lowqual_trim_extra)
    quals = [ord(char) - 33 for char in qual] if qual and qual != "*" else []
    aligned_end = start + aligned_length(cigar)
    clips: list[Softclip] = []
    for side, nclip, low, pos in (
        ("L", left_n, low_l, start),
        ("R", right_n, low_r, aligned_end),
    ):
        if low > 0:
            nclip = max(0, nclip - low)
        if nclip < cfg.min_clipped_bases:
            continue
        if quals:
            if side == "L":
                window = quals[low : low + nclip]
            else:
                stop = -low if low else None
                window = quals[-(low + nclip) : stop]
            if window and _median(window) < cfg.min_clipped_qual:
                continue
        clips.append(
            Softclip(
                pos=pos,
                side=side,
                mapq=mapq,
                clipped_len=nclip,
                aligned_start=start,
                aligned_end=aligned_end,
            )
        )
    return clips


def clip_breakpoints(clips: list[Softclip], cfg: CallConfig) -> list[Breakpoint]:
    """Clusters at or above ``min_softclip_support``, then microhomology merge."""
    clusters = cluster_softclips(clips, cfg.softclip_merge_distance)
    clusters = [c for c in clusters if c.n_reads >= cfg.min_softclip_support]
    pieces: list[Breakpoint] = []
    for start, end, n_left, n_right in merge_microhomology(clusters, cfg.microhomology_distance):
        interval_end = end if end > start else start + 1
        # A 1-bp interval stores the clip at ``start``. A wider pair stores the right-hand clip.
        clip_pos = start if interval_end == start + 1 else end
        pieces.append(
            Breakpoint(
                start,
                interval_end,
                n_softclip=n_left + n_right,
                n_clip_left=n_left,
                n_clip_right=n_right,
                clip_pos=clip_pos,
            )
        )
    return pieces


def _query_coords(cigar: str, reverse: bool) -> tuple[int, int] | None:
    """1-based closed query interval of the aligned portion. orig:get_query_coords."""
    ops = cigar_ops(cigar)
    if reverse:
        ops = list(reversed(ops))
    start = 1
    end = 1
    first = True
    seen_m = False
    pending_i = 0
    for n, op in ops:
        if op in "M=X":
            end += n + pending_i
            pending_i = 0
            seen_m = True
        elif op in "HS":
            if first:
                start += n
                end += n
        elif op == "I":
            if seen_m:
                pending_i = n
            else:
                start += n
                end += n
        elif op in "DN":
            pass
        else:
            return None
        first = False
    if start == end:
        return None
    return start, end - 1


def _reference_positions(start: int, cigar: str) -> list[int | None]:
    """One entry per query base. Soft clips and insertions are None."""
    positions: list[int | None] = []
    ref = start
    for n, op in cigar_ops(cigar):
        if op in "M=X":
            positions.extend(range(ref, ref + n))
            ref += n
        elif op in "IS":
            positions.extend([None] * n)
        elif op in "DN":
            ref += n
        elif op == "H":
            continue
        else:
            return []
    return positions


def _min_not_none(values: list[int | None]) -> int:
    present = [i for i, value in enumerate(values) if value is not None]
    return -1 if not present else min(present)


def _max_not_none(values: list[int | None]) -> int:
    present = [i for i, value in enumerate(values) if value is not None]
    return -1 if not present else max(present)


def _segment_bp(positions: list[int | None], reverse: bool, end_idx: int, other_idx: int, *, is_first: bool) -> tuple[int, int] | None:
    if is_first:
        start_idx = min(end_idx, other_idx - 1)
        while start_idx < 0 or positions[start_idx] is None:
            start_idx += 1
            if start_idx >= len(positions):
                return None
        if reverse:
            pair = (positions[end_idx], positions[start_idx])
        else:
            pair = (positions[start_idx] + 1, positions[end_idx] + 1)
    else:
        last_idx = max(end_idx, other_idx + 1)
        while last_idx >= len(positions) or positions[last_idx] is None:
            last_idx -= 1
            if last_idx < 0:
                return None
        if reverse:
            pair = (positions[last_idx] + 1, positions[end_idx] + 1)
        else:
            pair = (positions[end_idx], positions[last_idx])
    if pair[0] is None or pair[1] is None:
        return None
    return tuple(sorted((int(pair[0]), int(pair[1]))))


def parse_split(
    *,
    chrom: str,
    start: int,
    mapq: int,
    cigar: str,
    reverse: bool,
    sa: str,
    cfg: CallConfig,
) -> SplitHit | None:
    """orig:parse_splits for one SA tag. More than one supplementary is ignored."""
    if not sa or mapq < cfg.min_mapq:
        return None
    records = [part for part in sa.split(";") if part]
    if len(records) != 1:
        return None
    fields = records[0].split(",")
    if len(fields) < 5:
        return None
    supp_chrom, pos_s, strand, supp_cigar, mapq_s = fields[:5]
    if supp_chrom != chrom:
        return None
    try:
        supp_mapq = int(mapq_s)
        supp_start = int(pos_s) - 1
    except ValueError:
        return None
    if supp_mapq < cfg.min_mapq:
        return None
    if abs(supp_start - start) > cfg.max_edge_span:
        return None
    supp_reverse = strand == "-"
    left_is_primary = start <= supp_start
    if left_is_primary:
        l_start, l_cigar, l_rev = start, cigar, reverse
        r_start, r_cigar, r_rev = supp_start, supp_cigar, supp_reverse
    else:
        l_start, l_cigar, l_rev = supp_start, supp_cigar, supp_reverse
        r_start, r_cigar, r_rev = start, cigar, reverse
    left_q = _query_coords(l_cigar, l_rev)
    right_q = _query_coords(r_cigar, r_rev)
    if left_q is None or right_q is None or left_q[0] == right_q[0]:
        return None
    left_first = left_q[0] < right_q[0]
    if left_first:
        first = (l_start, l_cigar, l_rev, True)
        second = (r_start, r_cigar, r_rev, False)
    else:
        first = (r_start, r_cigar, r_rev, False)
        second = (l_start, l_cigar, l_rev, True)
    first_plus = not first[2]
    second_plus = not second[2]
    first_end = first[0] + aligned_length(first[1])
    second_end = second[0] + aligned_length(second[1])
    if first_plus == second_plus:
        if first_plus:
            overlap = not (first_end < second[0])
        else:
            overlap = not (second_end < first[0])
        kind = classify_split_flags(first_plus, second_plus, overlap, False)
    else:
        kind = classify_split_flags(first_plus, second_plus, False, left_first)
    first_pos = _reference_positions(first[0], first[1])
    second_pos = _reference_positions(second[0], second[1])
    if not first_pos or not second_pos:
        return None
    if all(value is None for value in first_pos) or all(value is None for value in second_pos):
        return None
    if first[2]:
        first_pos = list(reversed(first_pos))
    if second[2]:
        second_pos = list(reversed(second_pos))
    first_end_idx = _max_not_none(first_pos)
    second_start_idx = _min_not_none(second_pos)
    bp_first = _segment_bp(first_pos, first[2], first_end_idx, second_start_idx, is_first=True)
    bp_second = _segment_bp(second_pos, second[2], second_start_idx, first_end_idx, is_first=False)
    if bp_first is None or bp_second is None:
        return None
    if left_first:
        bp_left, bp_right = bp_first, bp_second
    else:
        bp_left, bp_right = bp_second, bp_first
    mapped_left = sum(value is not None for value in _reference_positions(l_start, l_cigar))
    mapped_right = sum(value is not None for value in _reference_positions(r_start, r_cigar))
    if bp_left[1] - bp_left[0] >= mapped_left or bp_right[1] - bp_right[0] >= mapped_right:
        return None
    if bp_left[0] <= bp_right[0]:
        bp_a, bp_b = bp_left, bp_right
    else:
        bp_a, bp_b = bp_right, bp_left
    return SplitHit(
        start1=first[0],
        end1=first[0] + aligned_length(first[1]),
        reverse1=first[2],
        start2=second[0],
        end2=second[0] + aligned_length(second[1]),
        reverse2=second[2],
        mapq=min(mapq, supp_mapq),
        split_type=kind,
        bp1=bp_a,
        bp2=bp_b,
        chrom=chrom,
        primary_start=start,
        primary_reverse=reverse,
    )


def split_breakpoints(splits: list[SplitHit]) -> list[Breakpoint]:
    """Each closed split interval becomes a half-open breakpoint with ``n_split``."""
    pieces: list[Breakpoint] = []
    for hit in splits:
        for bp in (hit.bp1, hit.bp2):
            lo, hi = bp
            pieces.append(Breakpoint(lo, hi + 1 if hi >= lo else lo + 1, n_split=1))
    return pieces
