"""Blocks, adjacency graph, cut points, and path enumeration (§4.4–§4.5).

Path search is the iterative DFS in orig:sv_inference.py::get_paths_iterative.
Neighbor lists are sorted by node id so the same graph always yields the same
paths. A locus whose traversal count exceeds ``max_paths`` is skipped.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Sequence

import numpy as np

# EdgeReadStore flag bits (match arcsv2_native / graph_native).
ER_HANGING = 1 << 0
ER_DISCORDANT = 1 << 1
ER_DISTANT = 1 << 2
ER_SUPPORTS_ADJ = 1 << 3
ER_IN_LIKELIHOOD = 1 << 4
KIND_PE = 0
KIND_SPLIT = 1


@dataclass
class Block:
    start: int
    end: int
    is_insertion: bool = False
    insertion_len: int = 0
    gap_left: bool = False
    # Soft-clip coordinate inside the uncovered interval to the right of this block.
    clip_pos: int | None = None
    clip_reads: int = 0
    # The tested insertion length is the largest size on the grid, or past the insert range.
    length_is_lower_bound: bool = False

    def __len__(self) -> int:
        if self.is_insertion:
            return self.insertion_len
        return max(0, self.end - self.start)


@dataclass
class EdgeRead:
    offset: int
    lib_id: int
    pmapped: float
    hanging: bool
    mapq: int
    kind: str  # "pe" or "split"
    # Template length outside the library quantile window. A concordant pair
    # that merely crosses the reference breakpoint does not span an insertion.
    discordant: bool = False
    # Distant mates use class_probs[3]. mate_len is the unmapped or distant query length.
    distant: bool = False
    mate_len: int = 0
    # True when this read incremented edge.support (PE adjacency gate or split).
    supports_adjacency: bool = False
    # Block path of each end, aligned to ``edge.u`` / ``edge.v`` (u < v).
    # None means a single-block end: every opposite-parity placement is allowed.
    # orig:add_pair adj1/adj2 via block_seq_to_path.
    adj1: tuple[int, ...] | None = None
    adj2: tuple[int, ...] | None = None
    # A split joined to its mate is still stored here for the sr column, but the
    # likelihood uses the pair edge and this read's block path as ``adj``.
    in_likelihood: bool = True


@dataclass
class EdgeReadStore:
    """Immutable CSR columns for edge-attached reads (R2 native / M4).

    Shared by all Edge shells on a Graph. Genotype workers fork-inherit the
    arrays read-only. Per-locus adj remaps live on Edge override dicts only.
    """

    offset: np.ndarray
    lib_id: np.ndarray
    pmapped: np.ndarray
    mapq: np.ndarray
    kind: np.ndarray
    flags: np.ndarray
    mate_len: np.ndarray
    adj1_ptr: np.ndarray
    adj1_len: np.ndarray
    adj2_ptr: np.ndarray
    adj2_len: np.ndarray
    adj_flat: np.ndarray

    def nbytes(self) -> int:
        return int(
            self.offset.nbytes
            + self.lib_id.nbytes
            + self.pmapped.nbytes
            + self.mapq.nbytes
            + self.kind.nbytes
            + self.flags.nbytes
            + self.mate_len.nbytes
            + self.adj1_ptr.nbytes
            + self.adj1_len.nbytes
            + self.adj2_ptr.nbytes
            + self.adj2_len.nbytes
            + self.adj_flat.nbytes
        )


class EdgeReadView:
    """Lightweight attribute view into an EdgeReadStore row.

    Optional adj overrides (per-locus remap) never mutate the shared store.
    """

    __slots__ = ("_store", "_i", "_adj1", "_adj2", "_has_adj_override")

    def __init__(
        self,
        store: EdgeReadStore,
        index: int,
        *,
        adj1: tuple[int, ...] | None = None,
        adj2: tuple[int, ...] | None = None,
        has_adj_override: bool = False,
    ) -> None:
        self._store = store
        self._i = int(index)
        self._adj1 = adj1
        self._adj2 = adj2
        self._has_adj_override = has_adj_override

    def _adj_from_store(self, ptr_arr, len_arr) -> tuple[int, ...] | None:
        ptr = int(ptr_arr[self._i])
        length = int(len_arr[self._i])
        if ptr < 0 or length <= 0:
            return None
        flat = self._store.adj_flat
        return tuple(int(x) for x in flat[ptr : ptr + length])

    @property
    def offset(self) -> int:
        return int(self._store.offset[self._i])

    @property
    def lib_id(self) -> int:
        return int(self._store.lib_id[self._i])

    @property
    def pmapped(self) -> float:
        return float(self._store.pmapped[self._i])

    @property
    def hanging(self) -> bool:
        return bool(int(self._store.flags[self._i]) & ER_HANGING)

    @property
    def mapq(self) -> int:
        return int(self._store.mapq[self._i])

    @property
    def kind(self) -> str:
        return "pe" if int(self._store.kind[self._i]) == KIND_PE else "split"

    @property
    def discordant(self) -> bool:
        return bool(int(self._store.flags[self._i]) & ER_DISCORDANT)

    @property
    def distant(self) -> bool:
        return bool(int(self._store.flags[self._i]) & ER_DISTANT)

    @property
    def mate_len(self) -> int:
        return int(self._store.mate_len[self._i])

    @property
    def supports_adjacency(self) -> bool:
        return bool(int(self._store.flags[self._i]) & ER_SUPPORTS_ADJ)

    @property
    def in_likelihood(self) -> bool:
        return bool(int(self._store.flags[self._i]) & ER_IN_LIKELIHOOD)

    @property
    def adj1(self) -> tuple[int, ...] | None:
        if self._has_adj_override:
            return self._adj1
        return self._adj_from_store(self._store.adj1_ptr, self._store.adj1_len)

    @property
    def adj2(self) -> tuple[int, ...] | None:
        if self._has_adj_override:
            return self._adj2
        return self._adj_from_store(self._store.adj2_ptr, self._store.adj2_len)


class _EdgeReadSeq(Sequence):
    """Sequence view over ``[read_start, read_end)`` in a shared store."""

    __slots__ = ("_edge",)

    def __init__(self, edge: "Edge") -> None:
        self._edge = edge

    def __len__(self) -> int:
        e = self._edge
        if e._reads is not None:
            return len(e._reads)
        return max(0, e.read_end - e.read_start)

    def __getitem__(self, index):
        e = self._edge
        if e._reads is not None:
            return e._reads[index]
        n = len(self)
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(n))]
        if index < 0:
            index += n
        if index < 0 or index >= n:
            raise IndexError(index)
        global_i = e.read_start + index
        assert e.store is not None
        ov = e._adj_override
        if ov is not None and index in ov:
            a1, a2 = ov[index]
            return EdgeReadView(e.store, global_i, adj1=a1, adj2=a2, has_adj_override=True)
        return EdgeReadView(e.store, global_i)

    def __iter__(self) -> Iterator:
        for i in range(len(self)):
            yield self[i]

    def append(self, read: EdgeRead) -> None:
        e = self._edge
        if e._reads is None:
            e._reads = []
            e.store = None
            e.read_start = 0
            e.read_end = 0
        e._reads.append(read)


class Edge:
    """Adjacency edge. Reads are either a Python list (engine=python) or a
    slice into a shared ``EdgeReadStore`` (engine=rust / M4).
    """

    __slots__ = (
        "u",
        "v",
        "support",
        "reference",
        "store",
        "read_start",
        "read_end",
        "_reads",
        "_adj_override",
    )

    def __init__(
        self,
        u: int,
        v: int,
        support: int = 0,
        reads: list[EdgeRead] | None = None,
        *,
        reference: bool = False,
        store: EdgeReadStore | None = None,
        read_start: int = 0,
        read_end: int = 0,
        _adj_override: dict
        | None = None,
        _reads: list[EdgeRead] | None = None,
    ) -> None:
        self.u = u
        self.v = v
        self.support = support
        self.reference = reference
        self.store = store
        self.read_start = read_start
        self.read_end = read_end
        self._adj_override = _adj_override
        if store is not None:
            self._reads = None
        elif _reads is not None:
            self._reads = _reads
        elif reads is not None:
            self._reads = list(reads)
        else:
            self._reads = []

    @property
    def reads(self) -> _EdgeReadSeq:
        return _EdgeReadSeq(self)

    @reads.setter
    def reads(self, value: list[EdgeRead]) -> None:
        self._reads = list(value)
        self.store = None
        self.read_start = 0
        self.read_end = 0
        self._adj_override = None

    def key(self) -> tuple[int, int]:
        return (min(self.u, self.v), max(self.u, self.v))

    def shell_copy(self) -> "Edge":
        """Copy mutable edge fields; share immutable CSR store (M3/M4 safe)."""
        if self.store is not None:
            return Edge(
                self.u,
                self.v,
                self.support,
                reference=self.reference,
                store=self.store,
                read_start=self.read_start,
                read_end=self.read_end,
                _adj_override=(
                    None if self._adj_override is None else dict(self._adj_override)
                ),
            )
        return Edge(
            self.u,
            self.v,
            self.support,
            reads=list(self._reads or ()),
            reference=self.reference,
        )


def create_blocks(
    intervals: list[tuple[int, int]],
    chrom_len: int,
    gaps: list[tuple[int, int]] | None = None,
) -> list[Block]:
    """Blocks between merged breakpoint intervals.

    orig:sv_parse_reads.py::create_blocks: a block ends at the breakpoint and
    the next begins there. Softclip/split merges are zero-width (point cuts),
    so blocks are contiguous. A half-open CI ``[start, end)`` with ``end > start``
    leaves that span uncovered (ARC-SV2 intervals). A point cut ``(p, p)``
    matches the original: no gap. A reference gap between two blocks marks the
    right-hand block ``gap_left`` so no reference edge is drawn across it.
    """
    cursor = 0
    blocks: list[Block] = []
    ordered = sorted((max(0, start), min(chrom_len, end)) for start, end in intervals)
    for start, end in ordered:
        if start >= chrom_len:
            continue
        if end < start:
            continue
        if start > cursor:
            blocks.append(Block(cursor, start))
        if end == start:
            # Point cut (orig): next block starts at the same coordinate.
            cursor = max(cursor, start)
        else:
            cursor = max(cursor, end)
    if cursor < chrom_len:
        blocks.append(Block(cursor, chrom_len))
    gap_list = gaps or []
    marked: list[Block] = []
    for i, block in enumerate(blocks):
        gap_left = False
        if i > 0:
            prev_end = blocks[i - 1].end
            for gs, ge in gap_list:
                if gs < block.start and ge > prev_end:
                    gap_left = True
                    break
        marked.append(Block(block.start, block.end, gap_left=gap_left))
    return [b for b in marked if len(b) > 0 or b.is_insertion]


def in_node(block: int) -> int:
    return 2 * block


def out_node(block: int) -> int:
    return 2 * block + 1


def find_block(blocks: list[Block], pos: int) -> int | None:
    """Smallest index with ``pos < block.end`` (orig:find_block_idx).

    Reference blocks are contiguous and ordered along the chromosome;
    novel-insertion blocks are only appended. Binary search over the
    reference prefix is O(log n) and matches the linear scan.
    """
    n = len(blocks)
    while n > 0 and blocks[n - 1].is_insertion:
        n -= 1
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if pos < blocks[mid].end:
            hi = mid
        else:
            lo = mid + 1
    if lo < n:
        return lo
    return None


def find_blocks_many(blocks: list[Block], positions):
    """Vectorised find_block for a 1-d int array of positions (T1).

    Returns int64 indexes; -1 means no block. Uses numpy.searchsorted so a
    thread pool can run over chromosome chunks without copying ``blocks``.
    """
    import numpy as np

    n = len(blocks)
    while n > 0 and blocks[n - 1].is_insertion:
        n -= 1
    pos = np.asarray(positions, dtype=np.int64)
    out = np.full(pos.shape, -1, dtype=np.int64)
    if n == 0 or pos.size == 0:
        return out
    ends = np.fromiter((blocks[i].end for i in range(n)), dtype=np.int64, count=n)
    idx = np.searchsorted(ends, pos, side="right")
    ok = idx < n
    out[ok] = idx[ok]
    return out


def read_vertex(block_idx: int, is_reverse: bool) -> int:
    """Forward reads sit on the out-vertex, reverse reads on the in-vertex.

    orig:sv_parse_reads.py::get_overlapping_blocks.
    """
    if is_reverse:
        return in_node(block_idx)
    return out_node(block_idx)


@dataclass
class Graph:
    blocks: list[Block]
    edges: dict[tuple[int, int], Edge] = field(default_factory=dict)
    _adj: dict[int, list[tuple[int, int]]] | None = field(default=None, repr=False)
    # Shared CSR for rust/M4 graphs. Fork workers inherit this read-only.
    edge_reads: EdgeReadStore | None = field(default=None, repr=False)

    @property
    def n_blocks(self) -> int:
        return len(self.blocks)

    def add_edge(
        self,
        u: int,
        v: int,
        *,
        support: int = 1,
        read: EdgeRead | None = None,
        reference: bool = False,
    ) -> Edge:
        self._adj = None
        key = (min(u, v), max(u, v))
        edge = self.edges.get(key)
        if edge is None:
            edge = Edge(key[0], key[1], support=0, reference=reference)
            self.edges[key] = edge
        if reference:
            edge.reference = True
            edge.support = max(edge.support, support)
        else:
            edge.support += support
        if read is not None:
            edge.reads.append(read)
        return edge

    def add_reference_edges(self, min_support: int) -> None:
        """orig:GenomeGraph.add_ref_path_support. No edge across a gap."""
        for i in range(self.n_blocks - 1):
            if self.blocks[i + 1].gap_left:
                continue
            if self.blocks[i].is_insertion or self.blocks[i + 1].is_insertion:
                continue
            self.add_edge(
                out_node(i),
                in_node(i + 1),
                support=min_support,
                reference=True,
            )

    def neighbors(self, node: int, min_support: int) -> list[int]:
        if self._adj is None:
            adj: dict[int, list[tuple[int, int]]] = {}
            for edge in self.edges.values():
                adj.setdefault(edge.u, []).append((edge.v, edge.support))
                if edge.v != edge.u:
                    adj.setdefault(edge.v, []).append((edge.u, edge.support))
            self._adj = adj
        found = [other for other, support in self._adj.get(node, ()) if support >= min_support]
        found.sort()
        return found

    def prune(self, min_support: int, max_span: int) -> None:
        """Drop adjacency edges below ``min_support`` or spanning more than ``max_span``.

        Reference edges inserted at exactly ``min_support`` are kept.
        """
        drop: list[tuple[int, int]] = []
        for key, edge in self.edges.items():
            if edge.support < min_support:
                # Reads are stored even when the edge is not a path adjacency
                # (out-of-range inserts, concordant pairs inside one block,
                # hanging reads). Path search ignores support below the minimum.
                if not edge.reads:
                    drop.append(key)
                    continue
                b1, b2 = edge.u // 2, edge.v // 2
                if b1 != b2:
                    left, right = sorted((b1, b2))
                    span = self.blocks[right].start - self.blocks[left].end
                    if span > max_span:
                        drop.append(key)
                continue
            b1, b2 = edge.u // 2, edge.v // 2
            if b1 == b2:
                continue
            left, right = sorted((b1, b2))
            if self.blocks[left].is_insertion or self.blocks[right].is_insertion:
                continue
            span = self.blocks[right].start - self.blocks[left].end
            if span > max_span:
                drop.append(key)
        for key in drop:
            del self.edges[key]
        self._adj = None


def _flip(node: int) -> int:
    return node - 1 if node % 2 else node + 1


def node_gap(blocks: list[Block], node: int) -> int:
    """Gap beside ``node``. orig:helper.py::block_gap."""
    idx = node // 2
    block = blocks[idx]
    if block.is_insertion:
        return 0
    if node % 2 == 1:
        if idx + 1 < len(blocks) and not blocks[idx + 1].is_insertion:
            return blocks[idx + 1].start - blocks[idx].end
        return 0
    if idx > 0 and not blocks[idx - 1].is_insertion:
        return blocks[idx].start - blocks[idx - 1].end
    return 0


def block_distance(path: list[int] | tuple[int, ...], blocks: list[Block], i: int, j: int) -> int:
    """Distance from path index ``i`` to ``j``. orig:helper.py::block_distance."""
    if i // 2 == j // 2:
        return 0
    dist = 0
    if path[i] % 2 == i % 2:
        dist += len(blocks[path[i] // 2])
    if path[j] % 2 != j % 2:
        dist += len(blocks[path[j] // 2])
    start = i + 2 if i % 2 == 0 else i + 1
    end = j - 1 if j % 2 == 1 else j
    for k in range(start, end, 2):
        dist += len(blocks[path[k] // 2])
    gap = sum(node_gap(blocks, path[k]) for k in range(start - 1, end + 1))
    dist += gap // 2
    if i % 2 == 0:
        dist = -dist
    return dist


class PathLayout:
    """Per-path indexes for bit-identical ``block_distance`` / placements.

    Prefix sums make each distance O(1). Node position lists avoid rescanning
    the path for every read. Used by genotype scoring; public helpers stay.
    """

    __slots__ = ("path", "blocks", "positions", "_block_len", "_even_cum", "_odd_cum", "_gap_cum")

    def __init__(self, path: tuple[int, ...] | list[int], blocks: list[Block]) -> None:
        self.path = path
        self.blocks = blocks
        positions: dict[int, list[int]] = {}
        n = len(path)
        block_len = [0] * n
        even_cum = [0] * (n + 1)
        odd_cum = [0] * (n + 1)
        gap_cum = [0] * (n + 1)
        for i, node in enumerate(path):
            positions.setdefault(node, []).append(i)
            block_len[i] = len(blocks[node // 2])
            even_cum[i + 1] = even_cum[i] + (block_len[i] if i % 2 == 0 else 0)
            odd_cum[i + 1] = odd_cum[i] + (block_len[i] if i % 2 == 1 else 0)
            gap_cum[i + 1] = gap_cum[i] + node_gap(blocks, node)
        self.positions = positions
        self._block_len = block_len
        self._even_cum = even_cum
        self._odd_cum = odd_cum
        self._gap_cum = gap_cum

    def distance(self, i: int, j: int) -> int:
        """Same integer result as ``block_distance(self.path, self.blocks, i, j)``."""
        if i // 2 == j // 2:
            return 0
        path = self.path
        block_len = self._block_len
        dist = 0
        if path[i] % 2 == i % 2:
            dist += block_len[i]
        if path[j] % 2 != j % 2:
            dist += block_len[j]
        start = i + 2 if i % 2 == 0 else i + 1
        end = j - 1 if j % 2 == 1 else j
        if start < end:
            if start % 2 == 0:
                dist += self._even_cum[end] - self._even_cum[start]
            else:
                dist += self._odd_cum[end] - self._odd_cum[start]
        gap_lo = start - 1
        gap_hi = end + 1
        if gap_lo < gap_hi:
            lo = 0 if gap_lo < 0 else gap_lo
            hi = gap_hi if gap_hi <= len(path) else len(path)
            if lo < hi:
                dist += (self._gap_cum[hi] - self._gap_cum[lo]) // 2
        if i % 2 == 0:
            dist = -dist
        return dist

    def placement_distances(
        self,
        v1: int,
        v2: int,
        adj1: tuple[int, ...] | None = None,
        adj2: tuple[int, ...] | None = None,
    ) -> list[int]:
        """Same order and values as module-level ``placement_distances``."""
        which1 = self.positions.get(v1, ())
        which2 = self.positions.get(v2, ())
        distances: list[int] = []
        path = self.path
        for i in which1:
            for j in which2:
                if (i - j) % 2 != 1:
                    continue
                if not is_adj_satisfied(adj1, path, i) or not is_adj_satisfied(adj2, path, j):
                    continue
                distances.append(self.distance(min(i, j), max(i, j)))
        return distances


def block_seq_to_path(seq: list[int]) -> tuple[int, ...]:
    """Expand a read's block vertices into the path ``is_adj_satisfied`` checks.

    orig:sv_parse_reads.py::block_seq_to_path. A one-vertex read is not an
    adjacency constraint; callers pass None instead of this result.
    """
    if not seq:
        return ()
    path = [seq[0]]
    for node in seq[1:-1]:
        flipped = node + 1 if node % 2 == 0 else node - 1
        path.extend((flipped, node))
    if len(seq) > 1:
        last = seq[-1]
        flipped = last + 1 if last % 2 == 0 else last - 1
        path.append(flipped)
    return tuple(path)


def is_adj_satisfied(adj: tuple[int, ...] | None, path: tuple[int, ...] | list[int], index: int) -> bool:
    """True when ``adj`` is a subpath of ``path`` starting at ``index``, either direction.

    orig:helper.py::is_adj_satisfied. ``None`` or an empty constraint always passes.
    """
    if not adj:
        return True
    n = len(adj)
    if index + n <= len(path) and tuple(path[index : index + n]) == adj:
        return True
    start = index + 1 - n
    if start >= 0 and tuple(path[start : index + 1]) == tuple(reversed(adj)):
        return True
    return False


def _in_one_bp_gap(blocks: list[Block], pos: int) -> tuple[int, int] | None:
    """If ``pos`` sits in a 1-bp uncovered gap between two ref blocks, return them.

    Half-open softclip CIs leave ``[prev.end, next.start)`` with
    ``next.start - prev.end == 1``. Point cuts close that gap; with gapped
    blocks a clip at ``prev.end`` must still resolve to a flanking block.
    """
    idx = find_block(blocks, pos)
    if idx is None or idx == 0:
        return None
    prev, nxt = blocks[idx - 1], blocks[idx]
    if prev.is_insertion or nxt.is_insertion:
        return None
    if nxt.start - prev.end == 1 and prev.end <= pos < nxt.start:
        return idx - 1, idx
    return None


def overlapping_vertices(blocks: list[Block], start: int, end: int, reverse: bool) -> list[int]:
    """Vertices a read's aligned interval touches, in the original's order.

    Forward reads are out-vertices left to right. Reverse reads are in-vertices
    right to left, so the anchor is the rightmost block.
    orig:get_overlapping_blocks, without CIGAR gaps.

    Lookup is O(log n + k): ``find_block`` then a short walk. A linear scan of
    every block per read made ``build_graph`` O(pairs × blocks) on chr20.

    A coordinate that falls in a 1-bp CI gap between blocks is treated as
    abutting the next block so Dup halves / multi-block PE mates still see
    the tandem block (NA19700 DUP 52764531 geometry).
    """
    if end <= start:
        return []
    first = find_block(blocks, start)
    if first is None:
        return []
    hits: list[int] = []
    for i in range(first, len(blocks)):
        block = blocks[i]
        if block.is_insertion:
            continue
        if block.start >= end:
            break
        if start < block.end and end > block.start:
            hits.append(i)
    # Read starts in a 1-bp gap just before the first hit: include that block.
    if not hits:
        gap = _in_one_bp_gap(blocks, start)
        if gap is not None:
            hits.append(gap[1])
    elif hits[0] > 0:
        gap = _in_one_bp_gap(blocks, start)
        if gap is not None and gap[1] == hits[0] and gap[0] not in hits:
            # start is in the gap; the next block is already first hit.
            pass
    if reverse:
        return [2 * block for block in reversed(hits)]
    return [2 * block + 1 for block in hits]


def placement_distances(
    path: tuple[int, ...] | list[int],
    blocks: list[Block],
    v1: int,
    v2: int,
    adj1: tuple[int, ...] | None = None,
    adj2: tuple[int, ...] | None = None,
) -> list[int]:
    """Block distances between placements of ``v1`` and ``v2`` whose adjacencies fit.

    orig:helper.py::get_block_distances_between_nodes. ``adj1`` belongs to ``v1``
    and ``adj2`` to ``v2``. A placement also needs opposite parity.
    """
    which1 = [i for i, node in enumerate(path) if node == v1]
    which2 = [i for i, node in enumerate(path) if node == v2]
    distances: list[int] = []
    for i in which1:
        for j in which2:
            if (i - j) % 2 != 1:
                continue
            if not is_adj_satisfied(adj1, path, i) or not is_adj_satisfied(adj2, path, j):
                continue
            distances.append(block_distance(path, blocks, min(i, j), max(i, j)))
    return distances


def _back_nodes(u: int, v: int, n_blocks: int) -> tuple[int, ...]:
    """orig:sv_inference.py::is_back_edge."""
    if u >= 2 * n_blocks or v >= 2 * n_blocks:
        return ()
    result: list[int] = []
    for n, m in ((u, v), (v, u)):
        is_out = n % 2 == 1
        if (is_out and m <= n) or ((not is_out) and n <= m):
            result.append(n)
    return tuple(result)


def enumerate_paths(
    graph: Graph,
    *,
    max_back: int = 2,
    max_paths: int = 1000,
    min_support: int = 2,
    start: int = 0,
    end: int | None = None,
) -> tuple[list[tuple[int, ...]], str]:
    """Iterative DFS. Status is ``ok``, ``skipped_npaths``, or ``skipped_backtrack``.

    orig:sv_inference.py::get_paths_iterative. The locus is skipped entirely
    when the number of traversals exceeds ``max_paths`` (paper: 1,000).
    """
    n_blocks = graph.n_blocks
    if end is None:
        end = 2 * n_blocks - 1
    start_in = start if start % 2 == 0 else start - 1
    end_out = end if end % 2 == 1 else end + 1
    if start_in > end_out:
        return [], "ok"

    paths: list[tuple[int, ...]] = []
    path: list[int] = []
    cur = start
    stacks: list[list[int]] = []
    back_count: dict[int, int] = {}
    npaths = 0
    backtracks = 0
    max_backtrack = 10 * max_paths

    while npaths <= max_paths and backtracks <= max_backtrack:
        nxt = _flip(cur)
        path.extend((cur, nxt))
        if nxt == end:
            npaths += 1
            paths.append(tuple(path))
            if npaths > max_paths:
                return [], "skipped_npaths"
        neigh: list[int] = []
        for n in graph.neighbors(nxt, min_support):
            backs = _back_nodes(nxt, n, n_blocks)
            if any(back_count.get(u, 0) >= max_back for u in backs):
                continue
            if start_in <= n <= end_out or _is_insertion_node(graph, n):
                neigh.append(n)
        if neigh:
            last_node = nxt
            last_stack = neigh
        else:
            backtracks += 1
            last_stack = []
            nback = 0
            last_node = nxt
            while not last_stack and stacks:
                last_stack = stacks.pop()
                nback += 1
                a = -2 * nback - 1
                edge = path[a:] if a + 2 == 0 else path[a : a + 2]
                if len(edge) == 2:
                    for u in _back_nodes(edge[0], edge[1], n_blocks):
                        back_count[u] = back_count.get(u, 0) - 1
                last_node = path[-2 * nback - 1]
            if not last_stack:
                return paths, "ok"
            path = path[: -2 * nback]
        stacks.append(last_stack)
        dest = last_stack.pop()
        for u in _back_nodes(last_node, dest, n_blocks):
            back_count[u] = back_count.get(u, 0) + 1
        cur = dest
    if backtracks > max_backtrack:
        return paths, "skipped_backtrack"
    return paths, "ok"


def _is_insertion_node(graph: Graph, node: int) -> bool:
    block = node // 2
    return block < len(graph.blocks) and graph.blocks[block].is_insertion


def _reference_blocks(blocks: list[Block]) -> int:
    """Reference blocks are the prefix. Novel insertions are appended after them."""
    n = 0
    for block in blocks:
        if block.is_insertion:
            break
        n += 1
    return n


def decompose(graph: Graph, min_support: int) -> list[tuple[int, int]]:
    """Cut points and minimal subgraphs. orig:sv_inference.py::decompose_graph.

    Returns inclusive ``(start_block, end_block)`` pairs. orig uses
    ``end_out = 2 * end + 1``, so ``end`` is the last block of the locus.
    Loci of only two reference blocks are discarded by the caller.
    """
    n = _reference_blocks(graph.blocks)
    if n <= 1:
        return []
    cut = {0, n - 1}
    spanned = [False] * n
    start_in = 0
    end_out = 2 * (n - 1) + 1
    # Paper: B_j is spanned when any edge joins B_i to B_k with i < j < k.
    # Edges into an appended insertion block are not reference spans.
    for edge in graph.edges.values():
        if edge.support < min_support:
            continue
        b1, b2 = edge.u // 2, edge.v // 2
        if b1 >= n or b2 >= n:
            continue
        if graph.blocks[b1].is_insertion or graph.blocks[b2].is_insertion:
            continue
        left, right = sorted((b1, b2))
        for j in range(left + 1, right):
            spanned[j] = True
    for b in range(1, n - 1):
        b_in = 2 * b
        b_out = 2 * b + 1
        inn = [v for v in graph.neighbors(b_in, min_support) if start_in <= v <= end_out]
        out = [v for v in graph.neighbors(b_out, min_support) if start_in <= v <= end_out]
        if (not inn or max(inn) < b_in) and (not out or min(out) > b_out):
            cut.add(b)
        for node in inn + out:
            if node > b_out + 2:
                block_idx = node // 2
                if block_idx >= n or graph.blocks[block_idx].is_insertion:
                    continue
                for s in range(b + 1, block_idx):
                    spanned[s] = True
            elif node < b_in - 2:
                block_idx = node // 2
                if block_idx >= n or graph.blocks[block_idx].is_insertion:
                    continue
                for s in range(block_idx + 1, b):
                    spanned[s] = True
    ordered = [c for c in sorted(cut) if not spanned[c]]
    minimal: list[int] = []
    cur: list[int] | None = None
    for cp in ordered:
        if cur is None:
            cur = [cp, cp]
        elif cp == cur[1] + 1:
            cur[1] = cp
        elif not minimal:
            minimal.append(cur[1])
            cur = [cp, cp]
        else:
            minimal.extend(cur)
            cur = [cp, cp]
    if cur is not None:
        minimal.append(cur[0])
    return [(minimal[i], minimal[i + 1]) for i in range(0, len(minimal) - 1, 2)]


def expand_locus(
    blocks: list[Block],
    start: int,
    end: int,
    width: int,
) -> tuple[int, int]:
    """Extend a locus so each flank covers about ``width`` bases.

    ``start`` and ``end`` are inclusive block indexes, matching
    orig:sv_inference.py where ``end_out = 2 * end + 1`` and the log prints
    ``blocks[end].end``. ``width`` is the insertion search width
    (1.1 × the 99% insert quantile).
    """
    anchor_l = blocks[start].end
    first = start
    while first > 0 and not blocks[first - 1].is_insertion and blocks[first].start > anchor_l - width:
        if blocks[first].gap_left:
            break
        first -= 1
    anchor_r = blocks[end].start
    last = end
    while (
        last + 1 < len(blocks)
        and not blocks[last + 1].is_insertion
        and blocks[last].end < anchor_r + width
    ):
        if blocks[last + 1].gap_left:
            break
        last += 1
    return first, last
