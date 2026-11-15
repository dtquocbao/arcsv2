"""Gap and insertion corrections to the haplotype normalizer.

Port of orig:sv_inference_insertions.py::genome_blocks_gaps,
get_gap_overlap_positions, get_insertion_overlap_positions, and
compute_normalizing_constant. Adjacent half-open pieces stay unmerged.
"""

from __future__ import annotations

from dataclasses import dataclass

from arcsv2.engine.intervals import IntervalSet, closed_iv, closedopen_iv, open_iv, overlaps
from arcsv2.engine.library import double_int


@dataclass
class HapBlock:
    """One block on a reference or a novel insertion. Length is ``end - start``."""

    start: int
    end: int
    is_de_novo: bool = False
    is_translocation: bool = False
    is_gap: bool = False

    def __len__(self) -> int:
        return self.end - self.start

    def is_insertion(self) -> bool:
        return self.is_de_novo or self.is_translocation


def block_gap(blocks: list[HapBlock], node: int) -> int:
    """Gap beside node ``node``. orig:helper.py::block_gap."""
    idx = node // 2
    block = blocks[idx]
    if block.is_insertion():
        return 0
    if node % 2 == 1:
        if idx + 1 < len(blocks) and not blocks[idx + 1].is_insertion():
            return blocks[idx + 1].start - blocks[idx].end
        return 0
    if idx > 0:
        return blocks[idx].start - blocks[idx - 1].end
    return 0


def genome_blocks_gaps(blocks: list[HapBlock], path: list[int] | tuple[int, ...]) -> list[HapBlock]:
    """Block sequence of ``path`` with the median of each crossed gap inserted.

    orig:sv_inference_insertions.py::genome_blocks_gaps. A zero-length gap is
    kept so later loops can skip it without shifting coordinates.
    Contiguous half-open block and gap pieces are joined by ``IntervalSet``.
    """
    out: list[HapBlock] = []
    out.append(HapBlock(0, block_gap(blocks, path[0]), is_gap=True))
    out.append(blocks[path[0] // 2])
    for i in range(1, len(path) - 1, 2):
        gap_size = (block_gap(blocks, path[i]) + block_gap(blocks, path[i + 1])) // 2
        out.append(HapBlock(0, gap_size, is_gap=True))
        out.append(blocks[path[i + 1] // 2])
    out.append(HapBlock(0, block_gap(blocks, path[-1]), is_gap=True))
    return out


def _to_open(result: IntervalSet) -> IntervalSet:
    """Drop empty pieces and rewrite every surviving piece as an open interval."""
    out = IntervalSet()
    for iv in result:
        lo = iv.lo - 1 if iv.lo_in else iv.lo
        hi = iv.hi + 1 if iv.hi_in else iv.hi
        if lo < hi - 1:
            out.add(open_iv(lo, hi))
    return out


def get_gap_overlap_positions(
    path: list[int] | tuple[int, ...],
    blocks: list[HapBlock],
    read_len: int,
    min_mappable: int = 20,
) -> IntervalSet:
    """Open intervals of read starts that overlap a gap. orig:get_gap_overlap_positions."""
    pieces = genome_blocks_gaps(blocks, path)
    gap_ref = IntervalSet()
    ref = IntervalSet()
    pos = 0
    for block in pieces:
        if len(block) == 0:
            continue
        if not block.is_insertion():
            gap_ref.add(closedopen_iv(pos, pos + len(block)))
            if not block.is_gap:
                ref.add(closedopen_iv(pos, pos + len(block)))
        pos += len(block)

    contained = IntervalSet()
    touches_ref = IntervalSet()
    for iv in gap_ref:
        if iv.lo <= iv.hi - min_mappable:
            contained.add(closed_iv(iv.lo, iv.hi - min_mappable))
    for iv in ref:
        touches_ref.add(closed_iv(iv.lo - min_mappable + 1, iv.hi - 1))
    both = contained.intersection(touches_ref)

    widened = IntervalSet()
    both_widened = IntervalSet()
    for iv in contained:
        widened.add(closed_iv(iv.lo - read_len + min_mappable, iv.hi))
    for iv in both:
        both_widened.add(closed_iv(iv.lo - read_len + min_mappable, iv.hi))
    return _to_open(widened.difference(both_widened))


def get_insertion_overlap_positions(
    path: list[int] | tuple[int, ...],
    blocks: list[HapBlock],
    read_len: int,
    min_mappable: int = 20,
) -> tuple[IntervalSet, list[bool], list[bool]]:
    """Open intervals of read starts that sit on a novel insertion.

    Returns the intervals plus, for each one, whether it overlaps a de novo
    insertion and whether it overlaps a translocation. orig:get_insertion_overlap_positions.
    """
    de_novo = IntervalSet()
    trans = IntervalSet()
    windows = IntervalSet()
    pos = 0
    for block in genome_blocks_gaps(blocks, path):
        span = len(block) - read_len + 2 * min_mappable
        if block.is_de_novo and span > 0:
            de_novo.add(open_iv(pos - min_mappable, pos + len(block) - read_len + min_mappable))
        elif block.is_translocation and span > 0:
            trans.add(open_iv(pos - min_mappable, pos + len(block) - read_len + min_mappable))
        if block.is_insertion():
            windows.add(open_iv(pos - min_mappable, pos + len(block)))
        pos += len(block)

    invalid = IntervalSet()
    for iv in windows:
        if iv.lo < iv.hi - (read_len - min_mappable):
            invalid.add(open_iv(iv.lo, iv.hi - (read_len - min_mappable)))

    on_de_novo: list[bool] = []
    on_trans: list[bool] = []
    for iv in invalid:
        on_de_novo.append(any(overlaps(iv, other) for other in de_novo))
        on_trans.append(any(overlaps(iv, other) for other in trans))
    return invalid, on_de_novo, on_trans


def compute_normalizing_constant(
    path: list[int] | tuple[int, ...],
    blocks: list[HapBlock],
    cdf_sum: list[float],
    class_prob: list[float],
    rlen1: int,
    rlen2: int,
) -> float:
    """Normalizer G_θ. orig:compute_normalizing_constant.

    ``class_prob`` is the mapping-outcome vector. The g[] coefficients follow
    the original even when several of the six placement classes are zero.
    """
    pieces = genome_blocks_gaps(blocks, path)
    length = sum(len(block) for block in pieces)
    if length <= 1:
        return 0.0

    def rect(a1: int, b1: int, a2: int, b2: int) -> float:
        return double_int(
            cdf_sum,
            a1,
            b1,
            a2,
            b2,
            xmin=0,
            xmax=length - 1,
            ymin=0,
            ymax=length - 1,
        )

    p_mapped = [
        1.0,
        class_prob[0] + class_prob[1] + class_prob[3],
        0.0,
        0.0,
        0.0,
        0.0,
    ]
    g = [
        p_mapped[0],
        p_mapped[1] - p_mapped[0],
        p_mapped[2] - 2 * p_mapped[1] + p_mapped[0],
        p_mapped[3] - p_mapped[0],
        p_mapped[4] - 2 * p_mapped[3] + p_mapped[0],
        p_mapped[5] - p_mapped[1] - p_mapped[3] + p_mapped[0],
    ]

    ins1, _, _ = get_insertion_overlap_positions(path, blocks, rlen1)
    ins2, _, _ = get_insertion_overlap_positions(path, blocks, rlen2)
    gap1 = list(get_gap_overlap_positions(path, blocks, rlen1))
    gap2 = list(get_gap_overlap_positions(path, blocks, rlen2))
    ins1_list = list(ins1)
    ins2_list = list(ins2)

    const = g[0] * rect(0, length - 1, 0, length - 1)
    for left in ins1_list:
        a, b = left.lo + 1, left.hi - 1
        const += g[1] * rect(a, b, 0, length - 1)
        for right in gap2:
            c, d = right.lo + 1, right.hi - 1
            const += g[5] * rect(a, b, c, d)
    for right in ins2_list:
        c, d = right.lo + 1, right.hi - 1
        const += g[1] * rect(0, length - 1, c, d)
        for left in ins1_list:
            a, b = left.lo + 1, left.hi - 1
            const += g[2] * rect(a, b, c, d)
        for left in gap1:
            a, b = left.lo + 1, left.hi - 1
            const += g[5] * rect(a, b, c, d)
    for left in gap1:
        a, b = left.lo + 1, left.hi - 1
        const += g[3] * rect(a, b, 0, length - 1)
    for right in gap2:
        c, d = right.lo + 1, right.hi - 1
        const += g[3] * rect(0, length - 1, c, d)
        for left in gap1:
            a, b = left.lo + 1, left.hi - 1
            const += g[4] * rect(a, b, c, d)
    return const
