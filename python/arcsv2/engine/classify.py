"""Classification and affected length (§4.8).

Complex type names are the lookup in orig:sv_classify.py::classify_complex
(the table near line 516). Affected length is Needleman–Wunsch on block tokens
(match 1000, mismatch -1, gap -1). Tokens match only when the block and the
orientation match. Every optimal traceback is kept, the strings are swapped,
and a block is affected if any of those alignments leaves it unmatched.
That is the union in orig:sv_affected_len.py::align_strings.
"""

from __future__ import annotations

# Paths are node sequences. The table is keyed by the simplified node path
# (orig:sv_classify.py::classify_complex).
COMPLEX_TYPES: dict[tuple[int, ...], str] = {
    (0, 1, 3, 2, 6, 7): "inv.del",
    (0, 1, 5, 4, 8, 9): "inv.2del",
    (0, 1, 4, 5, 2, 3, 6, 7): "trans",
    (0, 1, 5, 4, 2, 3, 6, 7): "trans.inv",
    (0, 1, 2, 3, 4, 5, 2, 3, 6, 7): "ddup",
    (0, 1, 2, 3, 4, 5, 3, 2, 6, 7): "invddup",
    (0, 1, 2, 3, 4, 5, 2, 3, 8, 9): "ddup.del",
    (0, 1, 2, 3, 4, 5, 3, 2, 8, 9): "invddup.del",
    (0, 1, 2, 3, 5, 4, 3, 2, 6, 7): "invddup.inv",
    (0, 1, 2, 3, 5, 4, 3, 2, 8, 9): "invddup.del.inv",
    (0, 1, 2, 3, 2, 3, 6, 7): "dup.del",
    (0, 1, 2, 3, 3, 2, 4, 5): "invdup",
    (0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7): "dup.dup",
}


def path_reversed(path: tuple[int, ...]) -> tuple[int, ...]:
    """orig:sv_classify.py::path_reversed."""
    if not path:
        return path
    big = max(path)
    return tuple(reversed([big - i for i in path]))


def complex_type(path: tuple[int, ...]) -> str | None:
    found = COMPLEX_TYPES.get(path)
    if found is None:
        found = COMPLEX_TYPES.get(path_reversed(path))
    return found


def align_tokens(
    s1: list[tuple],
    s2: list[tuple],
    match: int = 1000,
    mismatch: int = -1,
    gap: int = -1,
) -> tuple[set[int], set[int]]:
    """Needleman–Wunsch. Returns indices unaligned or mismatched in each string."""
    n, m = len(s1), len(s2)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    # 0 diagonal, 1 up (gap in s2), 2 left (gap in s1). Every co-optimal move is kept.
    tb: list[list[list[int]]] = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * gap
        tb[i][0] = [1]
    for j in range(1, m + 1):
        score[0][j] = j * gap
        tb[0][j] = [2]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            options = (
                (score[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch), 0),
                (score[i - 1][j] + gap, 1),
                (score[i][j - 1] + gap, 2),
            )
            best = max(value for value, _how in options)
            score[i][j] = best
            tb[i][j] = [how for value, how in options if value == best]

    aff1: set[int] = set()
    aff2: set[int] = set()

    def walk(i: int, j: int) -> None:
        if i == 0 and j == 0:
            return
        for how in tb[i][j]:
            if how == 0:
                if s1[i - 1] != s2[j - 1]:
                    aff1.add(i - 1)
                    aff2.add(j - 1)
                walk(i - 1, j - 1)
            elif how == 1:
                aff1.add(i - 1)
                walk(i - 1, j)
            else:
                aff2.add(j - 1)
                walk(i, j - 1)

    walk(n, m)
    return aff1, aff2


def affected_length(
    ref: list[tuple],
    alt: list[tuple],
    lengths: dict[int, int],
    insertion_lengths: dict[int, int] | None = None,
) -> tuple[int, set[int]]:
    """Union of the two orientations of the alignment (§4.8).

    ``ref`` / ``alt`` tokens are ``(block_index, strand)`` or ``(\"_\", ins_id)``.
    """
    aff_r, aff_a = align_tokens(ref, alt)
    aff_alt_as_query, aff_ref_swapped = align_tokens(alt, ref)
    ids: set[int] = set()
    for i in aff_r | aff_ref_swapped:
        token = ref[i]
        if isinstance(token[0], int):
            ids.add(token[0])
    for j in aff_a | aff_alt_as_query:
        token = alt[j]
        if isinstance(token[0], int):
            ids.add(token[0])
    total = sum(lengths.get(i, 0) for i in ids)
    insertion_lengths = insertion_lengths or {}
    for token in alt:
        if token[0] == "_":
            total += insertion_lengths.get(token[1], 0)
    return total, ids


def classify_simple(ref: list[tuple[int, int]], alt: list[tuple[int, int]]) -> str | None:
    """Return DEL, DUP, or INV when the alt is one simple event, else None."""
    ref_ids = [b for b, s in ref]
    alt_ids = [b for b, s in alt]
    alt_strands = {b: s for b, s in alt}
    alt_set = set(alt_ids)
    # Inversion: same multiset of blocks, exactly one reversed, order preserved.
    if alt_ids == ref_ids:
        flipped = [b for b, s in alt if s < 0]
        if len(flipped) == 1 and all(s > 0 or b in flipped for b, s in alt):
            return "INV"
        return None
    missing = [b for b in ref_ids if b not in alt_set]
    extra = [b for b in alt_ids if alt_ids.count(b) > ref_ids.count(b)]
    # Contiguous deletion of one or more blocks, flanks in order.
    if missing and not extra and alt_ids == [b for b in ref_ids if b not in set(missing)]:
        if all(alt_strands[b] > 0 for b in alt_ids):
            return "DEL"
    # Tandem duplication: one block repeated, everything else once and forward.
    if not missing and extra:
        duplicated = sorted(set(extra))
        if len(duplicated) == 1 and all(s > 0 for _, s in alt):
            return "DUP"
    return None


def simplify_blocks_diploid(
    blocks: list,
    path1: tuple[int, ...] | list[int],
    path2: tuple[int, ...] | list[int],
) -> tuple[list, tuple[int, ...], tuple[int, ...], list[list[int]]]:
    """Merge adjacent blocks whose mutual adjacency is unused as a novel cut.

    Port of ``orig:sv_validate.py::simplify_blocks_diploid``. When both
    haplotypes only traverse ``B_out → C_in`` (reference adjacency, degree 1
    on each side), the breakpoint between B and C is dropped and the blocks
    are fused. That turns a tandem dup of contiguous middle blocks
    (``ABCD`` / ``ABCBCD``) into ``ABC`` / ``ABBC`` for output and simple
    classification — the rule that removes ARC-SV2's extra cut at ABCBC sites.

    Returns ``(new_blocks, new_path1, new_path2, orig_groups)`` where
    ``orig_groups[i]`` lists the original block indexes fused into new block i.
    """
    from collections import defaultdict

    from arcsv2.engine.graph import Block

    def _block_nums(path: tuple[int, ...] | list[int]) -> list[int]:
        return [int(path[i] // 2) for i in range(1, len(path), 2)]

    block_nums = _block_nums(path1) + _block_nums(path2)
    neighbors: dict[int, set[int]] = defaultdict(set)
    for path in (path1, path2):
        if not path:
            continue
        neighbors[path[0]].add(-1)
        neighbors[path[-1]].add(-2)
        for i in range(1, len(path) - 1):
            if i % 2 == 0:
                neighbors[path[i]].add(path[i - 1])
            else:
                neighbors[path[i]].add(path[i + 1])

    genome = [b for b in block_nums if not blocks[b].is_insertion]
    if not genome:
        groups = [[i] for i in range(len(blocks))]
        return list(blocks), tuple(path1), tuple(path2), groups
    min_block = min(genome)
    max_block = max(genome)
    ins_blocks = [b for b in block_nums if blocks[b].is_insertion]
    new_blocks: list = []
    path_map: dict[int, int | None] = {}
    orig_groups: list[list[int]] = []
    idx = 0
    merging = False
    block_start = blocks[min_block].start
    pending_orig: list[int] = []
    for b in range(min_block, max_block + 1):
        if not merging:
            block_start = blocks[b].start
            pending_orig = []
        pending_orig.append(b)
        right_node = 2 * b + 1
        if (
            all(n == right_node + 1 for n in neighbors[right_node])
            and all(n == right_node for n in neighbors[right_node + 1])
            and b < max_block
        ):
            merging = True
            path_map[2 * b] = None
            path_map[2 * b + 1] = None
        else:
            src = blocks[b]
            new_blocks.append(
                Block(
                    block_start,
                    src.end,
                    is_insertion=False,
                    gap_left=blocks[pending_orig[0]].gap_left,
                    clip_pos=src.clip_pos,
                    clip_reads=src.clip_reads,
                )
            )
            path_map[2 * b] = 2 * idx
            path_map[2 * b + 1] = 2 * idx + 1
            orig_groups.append(list(pending_orig))
            merging = False
            idx += 1
            pending_orig = []
    for b in ins_blocks:
        src = blocks[b]
        new_blocks.append(
            Block(
                src.start,
                src.end,
                is_insertion=True,
                insertion_len=src.insertion_len or len(src),
                length_is_lower_bound=src.length_is_lower_bound,
            )
        )
        path_map[2 * b] = 2 * idx
        path_map[2 * b + 1] = 2 * idx + 1
        orig_groups.append([b])
        idx += 1

    def _remap(path: tuple[int, ...] | list[int]) -> tuple[int, ...]:
        return tuple(path_map[p] for p in path if path_map.get(p) is not None)

    return new_blocks, _remap(path1), _remap(path2), orig_groups
