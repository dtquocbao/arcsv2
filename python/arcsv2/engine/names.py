"""Block names. Port of orig:helper.py::block_idx_to_name."""

from __future__ import annotations

BLOCK_CHARS = 52
A_OFFSET = ord("A")
LOWER_OFFSET = 6  # skip the six non-letters between 'Z' and 'a'
INVERSION = "'"
INSERTION = "_"
TRANSLOCATION = "="


def block_idx_to_name(block_idx: int, start_idx: int = 0) -> str:
    """A..Z, a..z, A1..Z1, a1..z1, ... (orig:helper.py::block_idx_to_name)."""
    mult = (block_idx - start_idx) // BLOCK_CHARS
    mod = (block_idx - start_idx) % BLOCK_CHARS
    if mod > 25:
        mod += LOWER_OFFSET
    letter = chr(A_OFFSET + mod)
    if mult == 0:
        return letter
    return letter + str(mult)


def path_nodes_to_tokens(path: tuple[int, ...] | list[int]) -> list[tuple[int, int]]:
    """Odd node indexes name the block; an even node there is reverse orientation.

    Matches orig:helper.py::path_to_rearrangement, which reads path[1::2].
    Returns (block_index, strand) with strand +1 or -1.
    """
    tokens: list[tuple[int, int]] = []
    for node in path[1::2]:
        block = node // 2
        strand = -1 if node % 2 == 0 else 1
        tokens.append((block, strand))
    return tokens


def tokens_to_string(tokens: list[tuple[int, int]], start_idx: int = 0) -> str:
    parts: list[str] = []
    for block, strand in tokens:
        parts.append(block_idx_to_name(block, start_idx))
        if strand < 0:
            parts.append(INVERSION)
    return "".join(parts)


def reference_tokens(first: int, last_inclusive: int) -> list[tuple[int, int]]:
    return [(b, 1) for b in range(first, last_inclusive + 1)]
