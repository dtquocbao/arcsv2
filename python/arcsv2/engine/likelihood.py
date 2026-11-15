"""Haplotype likelihood and diploid genotyping (§4.6–§4.7).

The normalizer is orig:compute_normalizing_constant: the gap-free rectangle
sum, plus the gap and novel-insertion corrections. Per-read likelihood sums
``f(raw_offset + block_distance) · pmapped`` over placements whose block path
is a subpath of the haplotype (orig:get_block_distances_between_nodes). ``EdgeRead.offset``
is the raw offset, not the insert already adjusted by one shortcut.
A hanging read with no insertion overlap contributes ``class_probs[1]`` when
the mate is unmapped and ``class_probs[3]`` when the mate is distant. The
insertion-overlap term in orig:compute_hanging_edge_likelihood is applied
when the path contains a novel insertion. The overlap uses the mate query
length when that line was in the file, and only the first placement (``pr[0]``).

Diploid likelihood at allele fraction 1/2 follows the paper:
``Π [2π + (1-π)(p1 + p2)] / (G1 + G2)^n``.
Other fractions use orig:sv_inference.py::diploid_likelihood_frac.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from arcsv2.engine.graph import Block, Edge, EdgeRead, Graph, PathLayout, _flip
from arcsv2.engine.library import LibraryModel
from arcsv2.engine.names import path_nodes_to_tokens
from arcsv2.engine.normalizer import (
    HapBlock,
    compute_normalizing_constant,
    get_insertion_overlap_positions,
)


def haplotype_length(path: tuple[int, ...], graph: Graph) -> int:
    total = 0
    for node in path[0::2]:
        total += len(graph.blocks[node // 2])
    return total


def _hap_blocks(blocks: list[Block]) -> list[HapBlock]:
    converted: list[HapBlock] = []
    for block in blocks:
        if block.is_insertion:
            converted.append(HapBlock(0, block.insertion_len, is_de_novo=True))
        else:
            converted.append(HapBlock(block.start, block.end))
    return converted


def normalizer(
    path: tuple[int, ...],
    graph: Graph,
    lib: LibraryModel,
    hap: list[HapBlock] | None = None,
) -> float:
    """G_θ including gap and insertion corrections. An empty path returns 0."""
    if len(path) < 2:
        return 0.0
    return compute_normalizing_constant(
        path,
        hap if hap is not None else _hap_blocks(graph.blocks),
        lib.cdf_sum,
        lib.class_probs,
        lib.read_len_short,
        lib.read_len_long,
    )


def _cdf(lib: LibraryModel, t: int) -> float:
    if t < 0 or not lib.cdf:
        return 0.0
    if t >= len(lib.cdf):
        return lib.cdf[-1]
    return lib.cdf[t]


def _path_has_insertion(path: tuple[int, ...], blocks: list[HapBlock]) -> bool:
    """True when a novel insertion block appears on ``path`` (not merely in the graph)."""
    for node in path[0::2]:
        if blocks[node // 2].is_insertion():
            return True
    return False


def _hanging_probability(
    edge: Edge,
    read: EdgeRead,
    path: tuple[int, ...],
    blocks: list[HapBlock],
    lib: LibraryModel,
    *,
    path_has_ins: bool | None = None,
    overlap_cache: dict[int, object] | None = None,
) -> float:
    """One hanging read. Unmapped mates use class_probs[1]; distant mates use class_probs[3].

    When the path has a novel insertion, add the overlap term from
    orig:compute_hanging_edge_likelihood. Only the first placement of the
    anchored node is used, matching ``pr[0]`` in that function.
    """
    class_index = 3 if read.distant else 1
    base = lib.class_probs[class_index] if len(lib.class_probs) > class_index else 0.0
    anchored_p = 0.0
    if lib.class_probs:
        anchored_p = lib.class_probs[0]
        if len(lib.class_probs) > 1:
            anchored_p += lib.class_probs[1]
        if len(lib.class_probs) > 3:
            anchored_p += lib.class_probs[3]
    # Only insertions on this path change hanging likelihood. Scanning the full
    # chromosome block list (old behavior) forced the slow path after any novel
    # insertion was placed elsewhere; empty overlap intervals then returned base.
    if path_has_ins is None:
        path_has_ins = _path_has_insertion(path, blocks)
    if not path_has_ins:
        return base
    rlen = read.mate_len if read.mate_len > 0 else lib.read_len_long
    offset = read.offset
    if overlap_cache is not None and rlen in overlap_cache:
        intervals = overlap_cache[rlen]
    else:
        intervals, _, _ = get_insertion_overlap_positions(path, blocks, rlen)
        if overlap_cache is not None:
            overlap_cache[rlen] = intervals
    anchored = edge.v if edge.v % 2 == 1 else edge.u
    which = [i for i, node in enumerate(path) if node == anchored]
    if not which or not list(intervals):
        return base
    idx = which[0]
    block_5 = sum(len(blocks[path[i] // 2]) for i in range(0, idx - 1, 2))
    if idx % 2 != path[idx] % 2:
        block_5 += len(blocks[path[idx] // 2])
    mass = 0.0
    if idx % 2 == 1:
        overall = -block_5 + offset + rlen
        for interval in intervals:
            mass += max(0.0, _cdf(lib, interval.hi + overall - 1) - _cdf(lib, interval.lo + overall))
    else:
        overall = block_5 + offset
        for interval in intervals:
            mass += max(0.0, _cdf(lib, -interval.lo + overall - 1) - _cdf(lib, -interval.hi + overall))
    return mass * anchored_p * 0.5 + (1.0 - mass) * base


@dataclass
class ScoredPath:
    path: tuple[int, ...]
    length: int
    log_g: float
    # Per-read probabilities aligned with the read table. Missing reads are 0.
    read_p: list[float]
    # Cached exp(log_g) for diploid mixes; 0 when log_g is 0 / unset.
    g: float = 0.0


def score_path(
    path: tuple[int, ...],
    graph: Graph,
    libs: list[LibraryModel],
    read_index: list[tuple[Edge, int]],
    hap: list[HapBlock] | None = None,
) -> ScoredPath:
    """``read_index`` lists ``(edge, read_offset_in_edge.reads)`` in a fixed order."""
    length = haplotype_length(path, graph)
    if hap is None:
        hap = _hap_blocks(graph.blocks)
    # One normalizer per library; paper mode uses a single pooled library.
    log_g = 0.0
    if libs:
        g_raw = normalizer(path, graph, libs[0], hap=hap)
        log_g = math.log(g_raw) if g_raw > 0 else 0.0
    # Match prior log_diploid: exp(log_g), not the raw normalizer (float path).
    g_val = math.exp(log_g) if log_g else 0.0
    layout = PathLayout(path, graph.blocks)
    path_has_ins = _path_has_insertion(path, hap)
    # Placement distances depend only on endpoints + adj; reuse across reads.
    dist_cache: dict[
        tuple[int, int, tuple[int, ...] | None, tuple[int, ...] | None], list[int]
    ] = {}
    # Hanging probs for identical (edge ends, offset, mate_len, distant, lib).
    hang_cache: dict[tuple, float] = {}
    # Insertion-overlap intervals depend only on path + read length.
    overlap_cache: dict[int, object] = {}
    read_p: list[float] = []
    for edge, ridx in read_index:
        read = edge.reads[ridx]
        lib = libs[min(read.lib_id, len(libs) - 1)]
        if read.hanging:
            if not path_has_ins:
                class_index = 3 if read.distant else 1
                read_p.append(
                    lib.class_probs[class_index] if len(lib.class_probs) > class_index else 0.0
                )
                continue
            rlen = read.mate_len if read.mate_len > 0 else lib.read_len_long
            key_h = (edge.u, edge.v, read.offset, rlen, read.distant, read.lib_id)
            prob = hang_cache.get(key_h)
            if prob is None:
                prob = _hanging_probability(
                    edge,
                    read,
                    path,
                    hap,
                    lib,
                    path_has_ins=True,
                    overlap_cache=overlap_cache,
                )
                hang_cache[key_h] = prob
            read_p.append(prob)
            continue
        key = (edge.u, edge.v, read.adj1, read.adj2)
        dists = dist_cache.get(key)
        if dists is None:
            dists = layout.placement_distances(edge.u, edge.v, read.adj1, read.adj2)
            dist_cache[key] = dists
        total = 0.0
        pmf = lib.pmf
        n_pmf = len(pmf)
        off = read.offset
        for dist in dists:
            t = off + dist
            if 0 <= t < n_pmf:
                total += pmf[t]
        read_p.append(total * read.pmapped)
    return ScoredPath(path, length, log_g, read_p, g=g_val)


def log_haploid(scored: ScoredPath, pi: float, n_reads: int | None = None) -> float:
    """``Σ log(π + (1-π) p_i) - n log G``."""
    total = 0.0
    for p in scored.read_p:
        total += math.log(pi + (1.0 - pi) * p)
    n = len(scored.read_p) if n_reads is None else n_reads
    return total - n * scored.log_g


def log_diploid(
    a: ScoredPath,
    b: ScoredPath,
    pi: float,
    rho: float = 0.5,
    n_reads: int | None = None,
) -> float:
    """Diploid log-likelihood. ``rho`` is the allele fraction of haplotype ``b``.

    ``n_reads`` is the full locus count used in the normalizer. When ``read_p``
    is only the informative subset, the missing reads are charged ``log 2`` at
    ``rho = 0.5``, matching orig:diploid_likelihood.
    """
    n = len(a.read_p) if n_reads is None else n_reads
    dropped = n - len(a.read_p)
    g1 = a.g if a.g else (math.exp(a.log_g) if a.log_g else 0.0)
    g2 = b.g if b.g else (math.exp(b.log_g) if b.log_g else 0.0)
    if abs(rho - 0.5) < 1e-12:
        total = 0.0
        for p1, p2 in zip(a.read_p, b.read_p):
            total += math.log(2.0 * pi + (1.0 - pi) * (p1 + p2))
        if dropped > 0:
            total += dropped * math.log(2.0)
        # (G1 + G2)^n = exp(n * log(G1 + G2))
        if g1 + g2 <= 0:
            return total
        return total - n * math.log(g1 + g2)
    rho1 = 1.0 - rho
    rho2 = rho
    total = 0.0
    for p1, p2 in zip(a.read_p, b.read_p):
        total += math.log(pi + (1.0 - pi) * (rho1 * p1 + rho2 * p2))
    mix = rho1 * g1 + rho2 * g2
    if mix <= 0:
        return total
    return total - n * math.log(mix)


def reference_path(start: int, end_inclusive: int) -> tuple[int, ...]:
    """Nodes from block ``start`` through block ``end_inclusive``."""
    nodes: list[int] = []
    for b in range(start, end_inclusive + 1):
        nodes.extend((2 * b, 2 * b + 1))
    return tuple(nodes)


@dataclass
class GenotypeCall:
    path1: tuple[int, ...]
    path2: tuple[int, ...]
    rho: float
    loglik: float
    score_vs_ref: float
    score_vs_next: float
    second_paths: tuple[tuple[int, ...], tuple[int, ...]]
    n_paths: int
    n_informative: int


def _is_ref(path: tuple[int, ...], ref: tuple[int, ...]) -> bool:
    return path == ref


# HOM vs HET, and a lower copy of the same collapsed block sequence, are
# preferred only inside this band. Adjacency filtering at chr20:2659058 puts
# ABBC ~50 nats ahead of ABBBC, so a wider copy-family band is not applied.
LIKELIHOOD_BAND = 0.5


def _max_block_copies(path: tuple[int, ...]) -> int:
    """Most times any block is visited. A two-copy tandem duplication is 2."""
    counts: dict[int, int] = {}
    for node in path[0::2]:
        block = node // 2
        counts[block] = counts.get(block, 0) + 1
    return max(counts.values()) if counts else 0


def _collapsed(path: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    """Block sequence with adjacent copies removed. Orientation is kept."""
    folded: list[tuple[int, int]] = []
    for token in path_nodes_to_tokens(path):
        if folded and folded[-1] == token:
            continue
        folded.append(token)
    return tuple(folded)


def _same_copy_family(
    item: tuple[float, tuple[int, ...], tuple[int, ...], float],
    peak: tuple[float, tuple[int, ...], tuple[int, ...], float],
    ref: tuple[int, ...],
) -> bool:
    """True when both genotypes are the same blocks and orientations, apart from copy number.

    The homozygous reference is not treated as a lower copy of a duplication
    unless it is itself the best genotype.
    """
    item_sig = frozenset((_collapsed(item[1]), _collapsed(item[2])))
    peak_sig = frozenset((_collapsed(peak[1]), _collapsed(peak[2])))
    if item_sig != peak_sig:
        return False
    if item[1] == item[2] == ref and not (peak[1] == peak[2] == ref):
        return False
    return True


def select_in_band(
    candidates: list[tuple[float, tuple[int, ...], tuple[int, ...], float]],
    ref: tuple[int, ...],
) -> tuple[tuple[float, tuple[int, ...], tuple[int, ...], float], tuple[float, tuple[int, ...], tuple[int, ...], float] | None]:
    """Pick the called genotype and the best genotype outside the band.

    A lower copy number of the same collapsed sequence wins only inside
    ``LIKELIHOOD_BAND`` of the peak. Homozygosity is preferred only inside
    ``LIKELIHOOD_BAND`` of the best genotype at that copy number.
    Other types keep the higher log-likelihood.
    """
    peak_item = max(candidates, key=lambda item: (item[0], item[1], item[2]))
    peak = peak_item[0]
    family = [
        item
        for item in candidates
        if peak - item[0] < LIKELIHOOD_BAND and _same_copy_family(item, peak_item, ref)
    ]
    if family:
        def _copies(item: tuple[float, tuple[int, ...], tuple[int, ...], float]) -> int:
            return max(_max_block_copies(item[1]), _max_block_copies(item[2]))

        min_copies = min(_copies(item) for item in family)
        at_min = [item for item in family if _copies(item) == min_copies]
        best_ll = max(item[0] for item in at_min)
        near = [item for item in at_min if best_ll - item[0] < LIKELIHOOD_BAND]
        near.sort(
            key=lambda item: (
                0 if item[1] == item[2] else 1,
                -item[0],
                item[1],
                item[2],
                item[3],
            )
        )
        chosen = near[0]
    else:
        band = [item for item in candidates if peak - item[0] < LIKELIHOOD_BAND]
        band.sort(key=lambda item: (-item[0], item[1], item[2], item[3]))
        chosen = band[0]
    outside = [
        item
        for item in candidates
        if peak - item[0] >= LIKELIHOOD_BAND and (item[1], item[2]) != (chosen[1], chosen[2])
    ]
    outside.sort(key=lambda item: (-item[0], item[1], item[2], item[3]))
    return chosen, (outside[0] if outside else None)


def genotype_locus(
    paths: list[tuple[int, ...]],
    ref: tuple[int, ...],
    graph: Graph,
    libs: list[LibraryModel],
    *,
    pi: float,
    top_k: int = 50,
    allele_fractions: tuple[float, ...] = (0.5,),
) -> GenotypeCall | None:
    """Top-50 diploid search (§4.7). Returns None when no read is informative."""
    if ref not in paths:
        paths = [ref, *paths]
    # Stable unique.
    uniq: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for path in paths:
        if path not in seen and path:
            seen.add(path)
            uniq.append(path)
    # Reads whose blocks both sit in the locus. A read outside the locus is
    # omitted; scoring it as 0 would still change the -n log G term.
    locus_blocks = {node // 2 for node in ref}
    read_index: list[tuple[Edge, int]] = []
    for edge in sorted(graph.edges.values(), key=lambda e: e.key()):
        if edge.u // 2 not in locus_blocks or edge.v // 2 not in locus_blocks:
            continue
        for i in range(len(edge.reads)):
            if not edge.reads[i].in_likelihood:
                continue
            read_index.append((edge, i))
    if not read_index:
        return None
    # HapBlock list is path-independent; build once for all score_path calls.
    hap = _hap_blocks(graph.blocks)
    scored = {path: score_path(path, graph, libs, read_index, hap=hap) for path in uniq}
    ref_s = scored[ref]
    n_locus = len(read_index)

    def best_for(path: tuple[int, ...]) -> float:
        alt = scored[path]
        hom = log_diploid(alt, alt, pi, 0.5, n_locus)
        het = log_diploid(ref_s, alt, pi, 0.5, n_locus)
        return max(hom, het)

    ranked = sorted(uniq, key=lambda p: (-best_for(p), p))
    kept = []
    for path in ranked:
        if path == ref or len(kept) < top_k:
            kept.append(path)
        if len(kept) >= top_k and ref in kept:
            break
    if ref not in kept:
        kept.append(ref)

    # Informative reads differ across the kept paths. The product uses every
    # locus read, so the normalizer's n is the full count, not this subset.
    informative = []
    for i in range(len(read_index)):
        values = {scored[p].read_p[i] for p in kept}
        if len(values) > 1:
            informative.append(i)

    def ll(p1: tuple[int, ...], p2: tuple[int, ...], rho: float) -> float:
        return log_diploid(scored[p1], scored[p2], pi, rho, n_locus)

    candidates: list[tuple[float, tuple[int, ...], tuple[int, ...], float]] = []
    for i, p1 in enumerate(kept):
        for p2 in kept[i:]:
            fracs = (1.0,) if p1 == p2 else allele_fractions
            for rho in fracs:
                # Order paths so the comparison does not depend on kept order.
                left, right = (p1, p2) if p1 <= p2 else (p2, p1)
                used_rho = rho if left == p1 else (1.0 - rho if p1 != p2 else rho)
                if p1 == p2:
                    used_rho = 0.5
                    left = right = p1
                candidates.append((ll(left, right, used_rho), left, right, used_rho))
    # Copy-number ties use the lower copy, then homozygosity. Other types
    # keep the higher likelihood. score_vs_next uses the best genotype outside
    # the band.
    (best_ll, b1, b2, brho), runner = select_in_band(candidates, ref)
    ref_ll = ll(ref, ref, 0.5)
    if runner is not None:
        next_ll, n1, n2, _rho = runner
        next_paths = (n1, n2)
    else:
        next_ll = best_ll
        next_paths = (b1, b2)
    return GenotypeCall(
        path1=b1,
        path2=b2,
        rho=brho,
        loglik=best_ll,
        score_vs_ref=best_ll - ref_ll,
        score_vs_next=best_ll - next_ll,
        second_paths=next_paths,
        n_paths=len(uniq),
        n_informative=len(informative),
    )


def block_partner(node: int) -> int:
    return _flip(node)
