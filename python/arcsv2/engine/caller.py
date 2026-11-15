"""Evidence → calls. Ties §4.3–§4.8 together for one chromosome region."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from arcsv2.config import OUTPUT_COLUMNS, CallConfig
from arcsv2.engine.breakpoints import (
    Breakpoint,
    Drp,
    clamp_breakpoint,
    cluster_drps,
    cluster_to_intervals,
    fdr_keep,
    merge_breakpoints,
    relocate_drps,
    score_cluster,
)
from arcsv2.engine.classify import (
    affected_length,
    classify_simple,
    complex_type,
    simplify_blocks_diploid,
)
from arcsv2.engine.evidence import SplitHit
from arcsv2.engine.graph import (
    Block,
    Edge,
    EdgeRead,
    Graph,
    block_distance,
    block_seq_to_path,
    create_blocks,
    decompose,
    enumerate_paths,
    expand_locus,
    find_block,
    in_node,
    out_node,
    overlapping_vertices,
    read_vertex,
)
from arcsv2.engine.library import LibraryModel, fit_library, pi_robust, pool_inserts
from arcsv2.engine.likelihood import (
    genotype_locus,
    log_diploid,
    log_haploid,
    reference_path,
    score_path,
)
from arcsv2.engine.names import (
    block_idx_to_name,
    path_nodes_to_tokens,
    reference_tokens,
    tokens_to_string,
)
from arcsv2.progress import Progress


@dataclass
class PairRecord:
    """One model-eligible pair. Positions are 0-based half-open aligned intervals."""

    start1: int
    end1: int
    start2: int
    end2: int
    reverse1: bool
    reverse2: bool
    mapq1: int
    mapq2: int
    insert: int
    lib_id: int = 0
    hanging: bool = False
    qname_hash: int = 0
    # True when the mate is on another chromosome or farther than max_pair_distance.
    distant: bool = False
    # Query length of the hanging mate. 0 means the mate line was not in the file.
    mate_len: int = 0


@dataclass
class CallRow:
    values: dict[str, str]

    def as_dict(self) -> dict[str, str]:
        return {col: self.values.get(col, "NA") for col in OUTPUT_COLUMNS}


def _fmt_score(value: float) -> str:
    return f"{value:.2f}"


def _drp_from_pair(pair: PairRecord, lib: LibraryModel, index: int) -> Drp | None:
    """Classify a completed same-chrom pair (orig:pecluster.process_discordant_pair).

    Hanging pairs are never DRPs in the original — they stay graph evidence only.
    Opposite-orientation pairs whose insert sits in ``[z_ins, z_del]`` are
    concordant (early return, before Del/Dup/Ins). Else Del, then Dup (outward
    or nested), then Ins (short insert; overlapping mates allowed).
    Same-orientation pairs are InvL / InvR.
    """
    if pair.hanging:
        return None
    # Same orientation → inversion.
    if pair.reverse1 == pair.reverse2:
        if pair.reverse1 and pair.reverse2:
            positions = sorted((pair.start1, pair.start2))
            return Drp("InvR", positions[0], positions[1], pair.insert, index)
        if (not pair.reverse1) and (not pair.reverse2):
            positions = sorted((pair.end1, pair.end2))
            return Drp("InvL", positions[0], positions[1], pair.insert, index)
        return None
    # Opposite orientation. Concordant insert → not a DRP (orig:pecluster:39–41).
    if lib.z_ins_cut <= pair.insert <= lib.z_del_cut:
        return None
    # "First" = expected-left = forward for FR libraries.
    if not pair.reverse1 and pair.reverse2:
        fwd_start, fwd_end = pair.start1, pair.end1
        rev_start, rev_end = pair.start2, pair.end2
    else:
        fwd_start, fwd_end = pair.start2, pair.end2
        rev_start, rev_end = pair.start1, pair.end1
    if pair.insert > lib.z_del_cut and fwd_end < rev_start:
        return Drp("Del", fwd_end, rev_start, pair.insert, index)
    # Dup before Ins (orig:pecluster.py:57–65). Outward or nested FR geometry.
    if (fwd_start > rev_start) or (fwd_end > rev_end):
        return Drp("Dup", rev_start, fwd_end, pair.insert, index)
    if pair.insert < lib.z_ins_cut:
        # Overlapping mates keep (fwd_end, rev_start) even when fwd_end > rev_start.
        return Drp("Ins", fwd_end, rev_start, pair.insert, index)
    return None


def _fit_libs(
    pairs: list[PairRecord] | Sequence,
    cfg: CallConfig,
) -> list[LibraryModel]:
    by_lib: dict[int, list[int]] = {}
    for pair in pairs:
        if pair.hanging:
            continue
        if pair.mapq1 < cfg.min_mapq or pair.mapq2 < cfg.min_mapq:
            continue
        # Opposite-orientation signed inserts feed the KDE (orig:process_insert_len
        # appends only when 0 ≤ ilen ≤ cap). Same-orient → insert 0, skipped.
        if pair.reverse1 != pair.reverse2 and pair.insert >= 0:
            by_lib.setdefault(pair.lib_id, []).append(pair.insert)
    if cfg.sample_insert_sizes:
        cap = max(2, int(cfg.insert_sample_size))
        for lib_id, inserts in list(by_lib.items()):
            if len(inserts) > cap:
                by_lib[lib_id] = inserts[:cap]
    for lib_id, inserts in list(by_lib.items()):
        if len(inserts) < 2:
            continue
        med = sorted(inserts)[len(inserts) // 2]
        by_lib[lib_id] = [x for x in inserts if 0 <= x <= 3 * med] or inserts
    if not by_lib:
        # Fall back to every finite insert so a DRP-only fixture still fits.
        for pair in pairs:
            if not pair.hanging and pair.insert > 0:
                by_lib.setdefault(pair.lib_id, []).append(pair.insert)
    if cfg.pool_libraries:
        inserts = pool_inserts(list(by_lib.values()))
        return [fit_library(inserts, name="pooled", log_c=cfg.discordance_log_c, robustness=cfg.robustness_parameter)]
    libs = []
    for lib_id, inserts in sorted(by_lib.items()):
        libs.append(
            fit_library(
                inserts,
                name=f"lib{lib_id}",
                lib_id=lib_id,
                log_c=cfg.discordance_log_c,
                robustness=cfg.robustness_parameter,
            )
        )
    return libs


def discover_breakpoints(
    pairs: list[PairRecord] | Sequence,
    lib: LibraryModel,
    chrom_len: int,
    cfg: CallConfig,
    rng: random.Random,
    extra: list[Breakpoint] | None = None,
    progress: Progress | None = None,
) -> tuple[list[Breakpoint], list[str]]:
    logs: list[str] = []
    grouped: dict[str, list[Drp]] = {}
    if progress is not None:
        progress.stage("discover", f"scan {len(pairs):,} pairs for DRPs")
    for i, pair in enumerate(pairs):
        if pair.hanging:
            continue
        if min(pair.mapq1, pair.mapq2) < cfg.min_mapq:
            continue
        drp = _drp_from_pair(pair, lib, i)
        if drp is not None:
            grouped.setdefault(drp.type, []).append(drp)
    pieces: list[Breakpoint] = list(extra or [])
    distance = lib.insert_median + cfg.cluster_distance_sd * max(lib.insert_scale, 1.0)
    # One counter for the whole chromosome so Del/Inv/… clusters never share an id.
    cluster_id = 0
    for kind, group in sorted(grouped.items()):
        if progress is not None:
            progress.stage("discover", f"{kind}: {len(group):,} DRPs + {cfg.null_reps} null reps")
        clusters, cap_logs = cluster_drps(
            group,
            distance,
            ins_slop=cfg.drp_insertion_slop,
            ins_cluster_slop=cfg.drp_insertion_cluster_slop,
            max_cluster_size=cfg.max_pecluster_size,
        )
        logs.extend(cap_logs)
        scores = [score_cluster(cl, lib.insert_median, max(lib.insert_scale, 1.0)) for cl in clusters]
        null_scores: list[float] = []
        max_insert = lib.insert_median * cfg.insert_max_mu_multiple if kind == "Del" else None
        for rep in range(cfg.null_reps):
            # Independent stream per replicate, still determined by cfg.seed.
            moved = relocate_drps(group, chrom_len, rng, max_insert=max_insert)
            null_clusters, _ = cluster_drps(
                moved,
                distance,
                ins_slop=cfg.drp_insertion_slop,
                ins_cluster_slop=cfg.drp_insertion_cluster_slop,
                max_cluster_size=cfg.max_pecluster_size,
            )
            null_scores.extend(
                score_cluster(cl, lib.insert_median, max(lib.insert_scale, 1.0)) for cl in null_clusters
            )
            if progress is not None and (rep + 1) % max(1, cfg.null_reps // 2) == 0:
                progress.stage("discover", f"{kind}: null rep {rep + 1}/{cfg.null_reps}")
        keep = set(fdr_keep(scores, null_scores, cfg.null_reps, cfg.target_fdr))
        n_pass = len(keep)
        n_fail = len(clusters) - n_pass
        logs.append(
            f"{kind}: {len(group)} DRPs, {len(clusters)} clusters, "
            f"{n_pass} passing {n_fail} failing"
        )
        for idx, cluster in enumerate(clusters):
            if idx not in keep:
                continue
            left, right = cluster_to_intervals(cluster, cluster_id=cluster_id)
            cluster_id += 1
            pieces.extend((left, right))
    clamped = [item for item in (clamp_breakpoint(bp, chrom_len) for bp in pieces) if item is not None]
    result = merge_breakpoints(clamped, min_split_or_clip=cfg.min_edge_support)
    if progress is not None:
        progress.stage("discover done", f"{len(result):,} breakpoints")
    return result, logs


def _read_offset(start: int, end: int, reverse: bool, block: Block) -> int:
    """Offset of one read contained in ``block``. orig:get_overlapping_blocks."""
    if reverse:
        return end - block.start
    return block.start - start


def _pair_path(v1: int, v2: int) -> list[int]:
    """orig:add_pair combined_full_path for two single-block reads."""
    return [v1 ^ 1, v1, v2, v2 ^ 1]


def _adjacency_insert(pair: PairRecord, blocks: list[Block], b1: int, b2: int) -> int:
    """Insert implied by joining the two reads. orig:add_pair ``dist + offset``."""
    v1 = read_vertex(b1, pair.reverse1)
    v2 = read_vertex(b2, pair.reverse2)
    offset = _read_offset(pair.start1, pair.end1, pair.reverse1, blocks[b1])
    offset += _read_offset(pair.start2, pair.end2, pair.reverse2, blocks[b2])
    path = _pair_path(v1, v2)
    return block_distance(path, blocks, 1, len(path) - 2) + offset


def _is_outward(pair: PairRecord) -> bool:
    """Reverse read sits upstream of the forward read (tandem-duplication orientation)."""
    if pair.reverse1 == pair.reverse2:
        return False
    rev_start = pair.start1 if pair.reverse1 else pair.start2
    fwd_start = pair.start2 if pair.reverse1 else pair.start1
    return rev_start < fwd_start


def _dup_without_insert(pair: PairRecord, blocks: list[Block], block: int) -> int:
    """Insert implied if the same-block duplication edge is not used.

    orig:add_pair alt_combined_path. The duplication edge is kept only when
    this value is below the lower insert quantile.
    """
    v1 = read_vertex(block, pair.reverse1)
    v2 = read_vertex(block, pair.reverse2)
    offset = _read_offset(pair.start1, pair.end1, pair.reverse1, blocks[block])
    offset += _read_offset(pair.start2, pair.end2, pair.reverse2, blocks[block])
    alt = _pair_path(v1, v2)[2:]
    if len(alt) < 2:
        return offset
    return block_distance(alt, blocks, 1, len(alt) - 2) + offset


def _ref_path_insert(pair: PairRecord, blocks: list[Block], b1: int, b2: int) -> int:
    """Insert implied following reference through intervening blocks (no skip).

    Analogous to orig:add_pair ``alt_dist`` for tandem dups, but for a novel
    skip ``B_i.out → B_k.in`` with ``|i−k| > 1``: walk the reference path
    through the middle blocks instead of the direct adjacency.
    """
    left, right = sorted((b1, b2))
    v1 = read_vertex(b1, pair.reverse1)
    v2 = read_vertex(b2, pair.reverse2)
    offset = _read_offset(pair.start1, pair.end1, pair.reverse1, blocks[b1])
    offset += _read_offset(pair.start2, pair.end2, pair.reverse2, blocks[b2])
    # Canonical FR skip: leave left.out, enter right.in via reference middles.
    if v1 == out_node(left) and v2 == in_node(right):
        path = [v1 ^ 1, v1]
        for b in range(left + 1, right):
            path.extend([in_node(b), out_node(b)])
        path.extend([v2, v2 ^ 1])
        return block_distance(path, blocks, 1, len(path) - 2) + offset
    if v1 == out_node(right) and v2 == in_node(left):
        path = [v1 ^ 1, v1]
        for b in range(right - 1, left, -1):
            path.extend([out_node(b), in_node(b)])
        path.extend([v2, v2 ^ 1])
        return block_distance(path, blocks, 1, len(path) - 2) + offset
    # Non-canonical spanning adjacency: fall back to the direct insert.
    return _adjacency_insert(pair, blocks, b1, b2)


def _pair_supported(
    pair: PairRecord,
    blocks: list[Block],
    b1: int,
    b2: int,
    v1: int,
    v2: int,
    lib: LibraryModel,
    *,
    skip_edge_ref_gate: bool = False,
) -> bool:
    """Whether the between-mate adjacency increments edge support.

    Base rule (orig:add_pair): implied insert under the adjacency in
    ``[q_low, q_high]``. Same-block ``B.out↔B.in``: no support unless the
    pair is outward (reverse read upstream of the forward read) and the
    no-dup insert ``alt_dist < q_low`` (one-sided; orig:78–100). The outward
    test is not in literal ``add_pair``; without it, non-outward pairs open
    pair-only tandem edges. Novel skip edges optionally require the
    reference-path insert outside the window when ``skip_edge_ref_gate`` is set
    (default off = literal add_pair).
    """
    bi, bk = v1 // 2, v2 // 2
    same_block_dup = bi == bk and abs(v1 - v2) == 1
    implied = _adjacency_insert(pair, blocks, bi, bk)
    if not (lib.q_low <= implied <= lib.q_high):
        return False
    if same_block_dup:
        # orig:add_pair (sv_parse_reads.py:105–132) has no orientation test.
        # A non-outward pair whose raw offset is below q_low still has a
        # duplication-path insert inside the quantile when the block is about
        # one library insert long, and that opens a pair-only tandem edge.
        # chr20:59.4–64.4 Mb, each v2fix5 change toggled: requiring the reverse
        # read upstream of the forward read drops pair-only DUPs 13 → 2.
        # LIKELIHOOD_BAND and the pe column do not. See docs/PARITY.md.
        if not _is_outward(pair) or _dup_without_insert(pair, blocks, bi) >= lib.q_low:
            return False
        return True
    if skip_edge_ref_gate and abs(bi - bk) > 1:
        ref_implied = _ref_path_insert(pair, blocks, bi, bk)
        if lib.q_low <= ref_implied <= lib.q_high:
            return False
    return True


def _add_hanging(graph: Graph, pair: PairRecord, blocks: list[Block], lib: LibraryModel, cfg: CallConfig) -> None:
    """Store a mate-unmapped read on its block edge. Support stays 0.

    orig:add_hanging_pair. The edge is not a path adjacency.
    """
    if pair.mapq1 < cfg.min_mapq:
        return
    block_idx = find_block(blocks, pair.start1)
    if block_idx is None:
        return
    anchor = in_node(block_idx)
    graph.add_edge(
        anchor,
        anchor + 1,
        support=0,
        read=EdgeRead(
            offset=_read_offset(pair.start1, pair.end1, pair.reverse1, blocks[block_idx]),
            lib_id=0 if cfg.pool_libraries else pair.lib_id,
            pmapped=lib.pmapped,
            hanging=True,
            mapq=pair.mapq1,
            kind="pe",
            distant=pair.distant,
            mate_len=pair.mate_len,
        ),
    )


def _split_agrees(is_first: bool, reverse: bool, start: int, end: int, block: Block, leeway: int) -> bool:
    """orig:get_blocked_alignment junction check. One-sided, within ``leeway``."""
    if is_first != reverse:
        return block.end - end <= leeway
    return start - block.start <= leeway


def _split_support_edges(hit: SplitHit, blocks: list[Block], cfg: CallConfig) -> list[tuple[int, int]]:
    """Within-read edges of a split, from its block path.

    orig:add_pair walks ``block_seq_to_path`` two vertices at a time and calls
    ``add_support`` on each pair. ``find_block(start)`` misses a tandem edge when
    one half begins a few bases into the previous block.
    """
    parsed = _split_path(hit, blocks, cfg.split_read_leeway)
    if parsed is not None:
        verts, _offset = parsed
        path = block_seq_to_path(verts)
        edges = [(path[i], path[i + 1]) for i in range(0, len(path) - 1, 2)]
        if edges:
            return edges
    b1 = find_block(blocks, hit.start1)
    b2 = find_block(blocks, hit.start2)
    if b1 is None or b2 is None:
        return []
    if not _split_agrees(True, hit.reverse1, hit.start1, hit.end1, blocks[b1], cfg.split_read_leeway):
        return []
    if not _split_agrees(False, hit.reverse2, hit.start2, hit.end2, blocks[b2], cfg.split_read_leeway):
        return []
    v1 = read_vertex(b1, hit.reverse1)
    v2 = read_vertex(b2, hit.reverse2) ^ 1
    if v1 == v2:
        return []
    return [(v1, v2)]


def _add_split(
    graph: Graph,
    hit: SplitHit,
    blocks: list[Block],
    lib: LibraryModel,
    cfg: CallConfig,
    *,
    in_likelihood: bool = True,
) -> None:
    """Support each within-read edge of a split. orig:add_pair's path walk.

    Same-block Dup+/Dup− splits become the tandem-dup edge ``B.out ↔ B.in``.
    The PE insert-quantile gate does not apply to splits (Methods p. 21).
    A joined split is not scored (the pair carries its path); an unjoined split
    is scored on one edge so the evidence is not repeated.
    """
    edges = _split_support_edges(hit, blocks, cfg)
    if not edges:
        return
    b1 = find_block(blocks, hit.start1)
    b2 = find_block(blocks, hit.start2)
    if b1 is None or b2 is None:
        return
    offset = _read_offset(hit.start1, hit.end1, hit.reverse1, blocks[b1])
    offset += _read_offset(hit.start2, hit.end2, hit.reverse2, blocks[b2])
    scored = False
    for v1, v2 in edges:
        if v1 == v2:
            continue
        left_b, right_b = sorted((v1 // 2, v2 // 2))
        if left_b != right_b and blocks[right_b].start - blocks[left_b].end > cfg.max_edge_span:
            continue
        graph.add_edge(
            v1,
            v2,
            support=1,
            read=EdgeRead(
                offset=offset,
                lib_id=0,
                pmapped=lib.pmapped,
                hanging=False,
                mapq=hit.mapq,
                kind="split",
                supports_adjacency=True,
                in_likelihood=in_likelihood and not scored,
            ),
        )
        scored = True


def _end_vertices(
    blocks: list[Block], start: int, end: int, reverse: bool
) -> tuple[list[int], int] | None:
    """Anchor vertex list and offset. Reverse reads anchor on the rightmost block."""
    verts = overlapping_vertices(blocks, start, end, reverse)
    if not verts:
        block_idx = find_block(blocks, start)
        if block_idx is None:
            return None
        return [read_vertex(block_idx, reverse)], _read_offset(start, end, reverse, blocks[block_idx])
    return verts, _read_offset(start, end, reverse, blocks[verts[0] // 2])


def _adj_of(verts: list[int]) -> tuple[int, ...] | None:
    if len(verts) <= 1:
        return None
    return block_seq_to_path(verts)


def _segment_agrees(
    verts: list[int],
    *,
    is_first: bool,
    reverse: bool,
    start: int,
    end: int,
    blocks: list[Block],
    leeway: int,
) -> bool:
    """Junction check against the split-facing end of this segment.

    orig:get_blocked_alignment uses the last vertex of the first segment and
    the first vertex of the second. Failure is ``difference > leeway``.

    A 1-bp uncovered CI gap before ``block.start`` (or after ``block.end``)
    can place the split-facing coordinate one base outside the block; treat
    that as within leeway so Dup halves still form ``[B,B,…]`` paths.
    """
    vertex = verts[-1] if is_first else verts[0]
    block = blocks[vertex // 2]
    if is_first != reverse:
        # Allow end one base into the following 1-bp gap (block.end == gap start).
        return block.end - end <= leeway
    # Allow start one base before block.start (in the preceding 1-bp gap).
    return start - block.start <= leeway


def _split_path(
    hit: SplitHit, blocks: list[Block], leeway: int
) -> tuple[list[int], int] | None:
    """Query-order block path of both split segments, plus the first segment's offset."""
    first = overlapping_vertices(blocks, hit.start1, hit.end1, hit.reverse1)
    second = overlapping_vertices(blocks, hit.start2, hit.end2, hit.reverse2)
    if not first or not second:
        return None
    if not _segment_agrees(
        first, is_first=True, reverse=hit.reverse1, start=hit.start1, end=hit.end1, blocks=blocks, leeway=leeway
    ):
        return None
    if not _segment_agrees(
        second, is_first=False, reverse=hit.reverse2, start=hit.start2, end=hit.end2, blocks=blocks, leeway=leeway
    ):
        return None
    offset = _read_offset(hit.start1, hit.end1, hit.reverse1, blocks[first[0] // 2])
    return first + second, offset


def _resolve_end(
    start: int,
    end: int,
    reverse: bool,
    hits: list[SplitHit],
    blocks: list[Block],
    cfg: CallConfig,
    used_splits: set[int],
) -> tuple[list[int], int, tuple[int, ...] | None, bool] | None:
    """Vertices, offset, adjacency, and whether a split replaced the primary.

    The split flag lets ``build_graph`` skip within-read PE ``add_support`` when
    ``_add_split`` will credit the same edges.
    """
    matched = [hit for hit in hits if hit.primary_start == start and hit.primary_reverse == reverse]
    if len(matched) == 1:
        parsed = _split_path(matched[0], blocks, cfg.split_read_leeway)
        if parsed is not None:
            used_splits.add(id(matched[0]))
            verts, offset = parsed
            return verts, offset, _adj_of(verts), True
    parsed_end = _end_vertices(blocks, start, end, reverse)
    if parsed_end is None:
        return None
    verts, offset = parsed_end
    return verts, offset, _adj_of(verts), False


def _within_read_path_edges(verts: list[int]) -> list[tuple[int, int]]:
    """Within-read edges from a multi-block mate (orig:add_pair badj walk).

    Matches ``block_seq_to_path`` pairs before the between-mate junction. The
    original appends ``flip_parity(last)`` only to form that junction endpoint;
    within-read support does not use it.
    """
    if len(verts) <= 1:
        return []
    path = block_seq_to_path(verts)
    return [(path[i], path[i + 1]) for i in range(0, len(path) - 1, 2)]


def _extended_badj(verts: list[int]) -> tuple[int, ...]:
    """``block_seq_to_path`` plus orig ``flip_parity`` append for multi-block mates."""
    badj = block_seq_to_path(verts)
    if len(badj) > 1:
        return badj + (badj[-1] ^ 1,)
    return badj


def _between_junction(verts1: list[int], verts2: list[int]) -> tuple[int, int] | None:
    """Between-mate support endpoints after flip_parity (orig:add_pair).

    For single-block mates this equals ``(verts1[0], verts2[0])``. For
    multi-block mates the junction is the last extended badj vertex of each
    side — often the tandem ``B.out↔B.in`` while the stored PE edge is a skip.
    """
    badj1 = _extended_badj(verts1)
    badj2 = _extended_badj(verts2)
    if not badj1 or not badj2:
        return None
    return badj1[-1], badj2[-1]


def _between_junction_supported(
    pair: PairRecord,
    verts1: list[int],
    verts2: list[int],
    off1: int,
    off2: int,
    blocks: list[Block],
    lib: LibraryModel,
) -> bool:
    """Whether orig:add_pair credits the flip_parity between-mate junction.

    Insert under the combined path must lie in ``[q_low, q_high]``. Same-block
    junctions also need ``alt_dist < q_low``. The ARC-SV2 outward guard applies
    only when both mates are single-block (pair-only flood on insert-sized
    blocks); multi-block mates match literal add_pair.
    """
    badj1 = _extended_badj(verts1)
    badj2 = _extended_badj(verts2)
    if not badj1 or not badj2:
        return False
    combined = badj1 + tuple(reversed(badj2))
    if len(combined) < 2:
        return False
    full = (combined[0] ^ 1,) + combined + (combined[-1] ^ 1,)
    dist = block_distance(full, blocks, 1, len(full) - 2) + off1 + off2
    if not (lib.q_low <= dist <= lib.q_high):
        return False
    j1, j2 = badj1[-1], badj2[-1]
    same_block = j1 // 2 == j2 // 2 and abs(j1 - j2) == 1
    if same_block:
        m = len(badj1) - 1
        alt_full = full[:m] + full[(m + 2) :]
        if len(alt_full) < 2:
            return False
        alt_dist = block_distance(alt_full, blocks, 1, len(alt_full) - 2) + off1 + off2
        if alt_dist >= lib.q_low:
            return False
        if len(verts1) == 1 and len(verts2) == 1 and not _is_outward(pair):
            return False
    return True


def _add_within_read_pe_support(
    graph: Graph,
    verts: list[int],
    *,
    offset: int,
    lib: LibraryModel,
    cfg: CallConfig,
    mapq: int,
) -> None:
    """Credit within-read edges of a multi-block PE mate (orig:add_pair).

    Unconditional on insert concordance — same as the original's
    ``i != len(badj1)-1`` branch. Between-mate support stays in ``build_graph``.
    """
    for u, v in _within_read_path_edges(verts):
        if u == v:
            continue
        left_b, right_b = sorted((u // 2, v // 2))
        if (
            left_b != right_b
            and graph.blocks[right_b].start - graph.blocks[left_b].end > cfg.max_edge_span
        ):
            continue
        graph.add_edge(
            u,
            v,
            support=1,
            read=EdgeRead(
                offset=offset,
                lib_id=0,
                pmapped=lib.pmapped,
                hanging=False,
                mapq=mapq,
                kind="pe",
                supports_adjacency=True,
                in_likelihood=False,
            ),
        )


def build_graph(
    blocks: list[Block],
    pairs: list[PairRecord] | Sequence,
    lib: LibraryModel,
    cfg: CallConfig,
    splits: list[SplitHit] | None = None,
    progress: Progress | None = None,
) -> Graph:
    graph = Graph(blocks)
    graph.add_reference_edges(cfg.min_edge_support)
    splits_by_q: dict[int, list[SplitHit]] = {}
    for hit in splits or []:
        if hit.qname_hash < 0:
            continue
        splits_by_q.setdefault(hit.qname_hash, []).append(hit)
    used_splits: set[int] = set()
    n_pairs = len(pairs)
    if progress is not None:
        progress.stage("build_graph", f"{n_pairs:,} pairs × {len(blocks):,} blocks")
    report_every = max(100_000, n_pairs // 20) if n_pairs else 1
    for i, pair in enumerate(pairs):
        if progress is not None and i > 0 and i % report_every == 0:
            progress.stage("build_graph", f"{i:,}/{n_pairs:,} pairs")
        if pair.hanging:
            _add_hanging(graph, pair, blocks, lib, cfg)
            continue
        if pair.mapq1 < cfg.min_mapq or pair.mapq2 < cfg.min_mapq:
            continue
        hits = splits_by_q.get(pair.qname_hash, [])
        left = _resolve_end(
            pair.start1, pair.end1, pair.reverse1, hits, blocks, cfg, used_splits
        )
        right = _resolve_end(
            pair.start2, pair.end2, pair.reverse2, hits, blocks, cfg, used_splits
        )
        if left is None or right is None:
            continue
        verts1, off1, adj1, from_split1 = left
        verts2, off2, adj2, from_split2 = right
        v1, v2 = verts1[0], verts2[0]
        if v1 == v2:
            continue
        b1, b2 = v1 // 2, v2 // 2
        span_blocks = sorted((b1, b2))
        if span_blocks[0] != span_blocks[1]:
            span = blocks[span_blocks[1]].start - blocks[span_blocks[0]].end
            if span > cfg.max_edge_span:
                continue
        orient_ok = v1 == read_vertex(b1, pair.reverse1) and v2 == read_vertex(b2, pair.reverse2)
        junction = _between_junction(verts1, verts2)
        store_key = (min(v1, v2), max(v1, v2))
        junction_key = (
            (min(junction[0], junction[1]), max(junction[0], junction[1]))
            if junction is not None
            else None
        )
        # Support on the stored (v1,v2) edge only when it is the between junction
        # (single-block mates). Multi-block mates store PE on (v1,v2) for
        # likelihood but credit support on the flip_parity junction.
        if junction_key == store_key:
            supported = (
                _pair_supported(
                    pair, blocks, b1, b2, v1, v2, lib, skip_edge_ref_gate=cfg.skip_edge_ref_gate
                )
                if orient_ok
                else False
            )
        else:
            supported = False
        implied = _adjacency_insert(pair, blocks, b1, b2) if orient_ok else pair.insert
        if v1 < v2:
            adj_lo, adj_hi = adj1, adj2
        else:
            adj_lo, adj_hi = adj2, adj1
        graph.add_edge(
            v1,
            v2,
            support=1 if supported else 0,
            read=EdgeRead(
                offset=int(off1 + off2),
                lib_id=0 if cfg.pool_libraries else pair.lib_id,
                pmapped=lib.pmapped,
                hanging=False,
                mapq=min(pair.mapq1, pair.mapq2),
                kind="pe",
                discordant=not (lib.q_low <= implied <= lib.q_high),
                supports_adjacency=supported,
                adj1=adj_lo,
                adj2=adj_hi,
            ),
        )
        if (
            junction is not None
            and junction_key != store_key
            and junction[0] != junction[1]
            and _between_junction_supported(
                pair, verts1, verts2, int(off1), int(off2), blocks, lib
            )
        ):
            ju, jv = junction
            left_b, right_b = sorted((ju // 2, jv // 2))
            if (
                left_b == right_b
                or blocks[right_b].start - blocks[left_b].end <= cfg.max_edge_span
            ):
                graph.add_edge(
                    ju,
                    jv,
                    support=1,
                    read=EdgeRead(
                        offset=int(off1 + off2),
                        lib_id=0 if cfg.pool_libraries else pair.lib_id,
                        pmapped=lib.pmapped,
                        hanging=False,
                        mapq=min(pair.mapq1, pair.mapq2),
                        kind="pe",
                        supports_adjacency=True,
                        in_likelihood=False,
                    ),
                )
        # orig:add_pair: after the between-mate edge, walk each mate's within-read
        # path and add_support unconditionally. Splits already get that via
        # _add_split; skip those ends to avoid double-counting.
        if not from_split1 and len(verts1) > 1:
            _add_within_read_pe_support(
                graph,
                verts1,
                offset=int(off1),
                lib=lib,
                cfg=cfg,
                mapq=pair.mapq1,
            )
        if not from_split2 and len(verts2) > 1:
            _add_within_read_pe_support(
                graph,
                verts2,
                offset=int(off2),
                lib=lib,
                cfg=cfg,
                mapq=pair.mapq2,
            )
    for hit in splits or []:
        _add_split(
            graph,
            hit,
            blocks,
            lib,
            cfg,
            in_likelihood=id(hit) not in used_splits,
        )
    graph.prune(cfg.min_edge_support, cfg.max_edge_span)
    if progress is not None:
        progress.stage("build_graph done", f"{len(graph.edges):,} edges")
    return graph


def _support_for_junction(graph: Graph, b_left: int, b_right: int) -> tuple[int, int]:
    """Split and pair support on edges for one novel junction.

    Same-block tandem ``(B, B)`` only counts the within-block edge (orig pe is
    typed Dup DRP ∩ at breakpoints — split-only dups report pe=0). Cross-block
    junctions require both endpoints on the edge. PE counts only reads that
    incremented ``edge.support``.
    """
    pe = 0
    sr = 0
    for edge in graph.edges.values():
        eu, ev = edge.u // 2, edge.v // 2
        if b_left == b_right:
            if not (eu == ev == b_left):
                continue
        elif {eu, ev} != {b_left, b_right}:
            continue
        for read in edge.reads:
            if read.kind == "split":
                sr += 1
            elif read.supports_adjacency:
                pe += 1
    return sr, pe


def _insertion_position(blocks: list[Block], block: int, strand: int) -> int:
    """Junction coordinate. A clip inside the uncovered interval replaces the interval's left edge."""
    if strand > 0:
        clipped = blocks[block].clip_pos
        if clipped is not None:
            return clipped
        return blocks[block].end
    return blocks[block].start


def _blocks_column(blocks: list[Block], start: int, end: int) -> str:
    """``start-end`` for each locus block, in order. An insertion block is ``_``."""
    parts: list[str] = []
    for index in range(start, end + 1):
        block = blocks[index]
        if block.is_insertion:
            parts.append("_")
        else:
            parts.append(f"{block.start}-{block.end}")
    return ",".join(parts)


def _breakpoint_point(left: Block, right: Block) -> tuple[int, str]:
    """Clip position inside the uncovered interval, otherwise its midpoint."""
    lo, hi = left.end, right.start
    if hi < lo:
        return lo, "mid"
    clip = left.clip_pos
    if clip is not None and lo <= clip <= hi:
        return clip, "clip"
    return (lo + hi) // 2, "mid"


def _blocks_point_column(blocks: list[Block], start: int, end: int) -> str:
    """Locus blocks that meet at each breakpoint's point estimate.

    ``blocks`` keeps the uncovered interval. This column closes that interval
    at the soft-clip base when the interval has one, otherwise at its midpoint.
    Insertion blocks stay ``_``.
    """
    indexes = list(range(start, end + 1))
    reference = [index for index in indexes if not blocks[index].is_insertion]
    points = {
        (left, right): _breakpoint_point(blocks[left], blocks[right])
        for left, right in zip(reference, reference[1:])
    }
    closed: dict[int, str] = {}
    for place, index in enumerate(reference):
        block = blocks[index]
        left_edge = block.start if place == 0 else points[(reference[place - 1], index)][0]
        right_edge = (
            block.end if place == len(reference) - 1 else points[(index, reference[place + 1])][0]
        )
        if right_edge < left_edge:
            right_edge = left_edge
        closed[index] = f"{left_edge}-{right_edge}"
    parts = ["_" if blocks[index].is_insertion else closed[index] for index in indexes]
    return ",".join(parts)


def _allele_base_coords(
    blocks: list[Block],
    path: tuple[int, ...] | list[int],
    start: int,
    end: int,
    *,
    flank: int | None = None,
) -> tuple[tuple[int, int], ...]:
    """Ordered ``(ref_pos, strand)`` bases along a path under ``blocks_point`` closure.

    Same closure the labeler uses. With ``flank=None`` (default), the full
    closed spans are used — ABCBCD and A(BC)(BC)D are the same bases when the
    fused block spans the closed B+C intervals. With ``flank=N``, outer blocks
    are trimmed to N bp exactly as ``label._prepared_blocks`` does; fusion that
    pulls an internal kept block into an outer flank then changes the query.
    """
    point_col = _blocks_point_column(blocks, start, end)
    spans: dict[int, tuple[int, int]] = {}
    ordered: list[int] = []
    for index, part in zip(range(start, end + 1), point_col.split(",")):
        part = part.strip()
        if not part or part == "_":
            continue
        left, right = part.split("-", 1)
        spans[index] = (int(left), int(right))
        ordered.append(index)
    if flank is not None and flank > 0 and len(ordered) >= 2:
        left_i, right_i = ordered[0], ordered[-1]
        lo, hi = spans[left_i]
        if hi - lo > flank:
            spans[left_i] = (hi - flank, hi)
        lo, hi = spans[right_i]
        if hi - lo > flank:
            spans[right_i] = (lo, lo + flank)
    coords: list[tuple[int, int]] = []
    for block, strand in path_nodes_to_tokens(path):
        if block < len(blocks) and blocks[block].is_insertion:
            continue
        span = spans.get(block)
        if span is None:
            continue
        left, right = span
        if strand > 0:
            coords.extend((pos, 1) for pos in range(left, right))
        else:
            coords.extend((pos, -1) for pos in range(right - 1, left - 1, -1))
    return tuple(coords)


def _bp_point_source_column(blocks: list[Block], start: int, end: int) -> str:
    """``clip`` or ``mid`` for each internal junction in locus order."""
    indexes = list(range(start, end + 1))
    reference = [index for index in indexes if not blocks[index].is_insertion]
    sources = [
        _breakpoint_point(blocks[left], blocks[right])[1]
        for left, right in zip(reference, reference[1:])
    ]
    return ",".join(sources) if sources else "NA"


def _format_inslen(blocks: list[Block], ins_lengths: dict[int, int]) -> str:
    """Point estimate, or ``>N`` when the tested size is only a lower bound."""
    if not ins_lengths:
        return "NA"
    parts: list[str] = []
    for block_id, length in ins_lengths.items():
        text = str(length)
        if block_id < len(blocks) and blocks[block_id].length_is_lower_bound:
            text = ">" + text
        parts.append(text)
    return ",".join(parts)


def _gap_before(blocks: list[Block], b: int) -> int:
    """Uncovered width on the left of block ``b`` (orig:block_gap on the in-node)."""
    if b <= 0 or b >= len(blocks):
        return 0
    return max(0, blocks[b].start - blocks[b - 1].end)


def _gap_after(blocks: list[Block], b: int) -> int:
    """Uncovered width on the right of block ``b`` (orig:block_gap on the out-node)."""
    if b < 0 or b + 1 >= len(blocks):
        return 0
    return max(0, blocks[b + 1].start - blocks[b].end)


def reference_block_lengths_midpoint(
    blocks: list[Block],
    start: int,
    end: int,
) -> dict[int, int]:
    """Reference block lengths with breakpoint gaps split at their median.

    Methods p.25 / orig:sv_output.py: each internal block boundary is
    ``floor(median(prev.end, next.start))``. ``len_affected`` sums those
    midpoint-extended lengths for affected reference blocks. Uncovered
    ``blocks`` / ``bp_uncertainty`` columns are unchanged.
    """
    idxs = [b for b in range(start, end + 1) if not blocks[b].is_insertion]
    if not idxs:
        return {}
    boundaries = [blocks[idxs[0]].start]
    for i in range(1, len(idxs)):
        prev_b = idxs[i - 1]
        cur_b = idxs[i]
        boundaries.append(int(math.floor((blocks[prev_b].end + blocks[cur_b].start) / 2)))
    boundaries.append(blocks[idxs[-1]].end)
    return {
        idxs[i] - start: boundaries[i + 1] - boundaries[i]
        for i in range(len(idxs))
    }


def _support_for_merged_tandem(
    graph: Graph, orig_blocks: list[int]
) -> tuple[int, int]:
    """Split/PE support for a tandem spanning one or more original blocks.

    Single-block: within-block edge. Multi-block (after simplify fused B..C):
    the wrap ``last_out → first_in`` edge, matching orig Dup support on the
    outer breakpoints of the duplicated segment.
    """
    if not orig_blocks:
        return 0, 0
    if len(orig_blocks) == 1:
        return _support_for_junction(graph, orig_blocks[0], orig_blocks[0])
    return _support_for_junction(graph, orig_blocks[-1], orig_blocks[0])


def _row_for_genotype(
    chrom: str,
    blocks: list[Block],
    graph: Graph,
    start: int,
    end: int,
    geno,
    cfg: CallConfig,
) -> CallRow | None:
    ref = reference_path(start, end)
    alt = geno.path1 if geno.path1 != ref else geno.path2
    if alt == ref and geno.path2 == ref:
        return None
    if geno.score_vs_ref <= 0:
        return None
    # orig:sv_inference.py — simplify unused cuts before path strings / simple
    # classification; rearrangement_next stays on the unsimplified runner-up.
    nb, np1, np2, orig_groups = simplify_blocks_diploid(blocks, geno.path1, geno.path2)
    n_ref = sum(1 for b in nb if not b.is_insertion)
    if n_ref < 1:
        return None
    ref_s = reference_path(0, n_ref - 1)
    alt_s = np1 if np1 != ref_s else np2
    if alt_s == ref_s and np2 == ref_s:
        return None
    out_blocks = nb
    out_start, out_end = 0, n_ref - 1
    # Fusion renames paths (ABCBCD → ABBC) but must not change the predicted
    # allele sequence the labeler builds from blocks_point. Emit point blocks
    # from the simplified locus (not the pre-fusion cut list).
    fused = any(
        len(group) > 1 and not blocks[group[0]].is_insertion for group in orig_groups
    )
    if fused:
        for before_path, after_path, which in (
            (geno.path1, np1, "path1"),
            (geno.path2, np2, "path2"),
        ):
            before = _allele_base_coords(blocks, before_path, start, end)
            after = _allele_base_coords(out_blocks, after_path, out_start, out_end)
            if before != after:
                raise RuntimeError(
                    "simplify_blocks_diploid changed predicted allele sequence "
                    f"on {which} ({len(before)} -> {len(after)} bases under blocks_point)"
                )
        # Labeler also trims each outer block to 150 bp. Fusing a kept internal
        # block into an outer flank (DEL ACD→AC) changes that query even when
        # full allele bases match — revert simplify so flanks stay stable.
        # ABCBCD→ABBC keeps A/D as outer blocks, so it still simplifies.
        label_flank = 150
        flank_changed = any(
            _allele_base_coords(blocks, before_path, start, end, flank=label_flank)
            != _allele_base_coords(
                out_blocks, after_path, out_start, out_end, flank=label_flank
            )
            for before_path, after_path in (
                (geno.path1, np1),
                (geno.path2, np2),
            )
        )
        if flank_changed:
            out_blocks = blocks
            out_start, out_end = start, end
            np1, np2 = geno.path1, geno.path2
            orig_groups = [[i] for i in range(len(blocks))]
            ref_s = reference_path(out_start, out_end)
            alt_s = np1 if np1 != ref_s else np2
            if alt_s == ref_s and np2 == ref_s:
                return None
    ref_tokens = reference_tokens(out_start, out_end)
    alt_tokens = path_nodes_to_tokens(alt_s)
    named_ref = [(b - out_start, s) for b, s in ref_tokens]
    ins_lengths: dict[int, int] = {}
    named_alt: list[tuple] = []
    for b, s in alt_tokens:
        if out_blocks[b].is_insertion:
            named_alt.append(("_", b))
            ins_lengths[b] = len(out_blocks[b])
        else:
            named_alt.append((b - out_start, s))
    lengths = reference_block_lengths_midpoint(out_blocks, out_start, out_end)
    length, _ids = affected_length(named_ref, named_alt, lengths, ins_lengths)
    # Do not drop short events here. orig:do_inference emits them; paper ≥50 bp
    # lives in merge-filter (PARITY D11). CI-shrunk middle blocks of true 54–74 bp
    # DELs were 39–48 bp and were lost at emit (CURSOR_TASK_NEXT_V2FIX6 §4).
    pure_ins = [t for t in named_alt if t[0] != "_"] == named_ref and bool(ins_lengths)
    genome_alt = [t for t in named_alt if t[0] != "_"]
    simple = None if ins_lengths else classify_simple(named_ref, genome_alt)
    # Complex lookup on simplified absolute nodes (locus rooted at 0).
    translated = tuple(
        n
        for n in alt_s
        if (n // 2) < len(out_blocks) and not out_blocks[n // 2].is_insertion
    )
    ctype = None if ins_lengths else (complex_type(translated) or complex_type(alt_s))
    if pure_ins:
        svtype = "INS"
        complextype = "NA"
        num_sv = 1
    elif simple is not None and ctype is None:
        svtype = simple
        complextype = "NA"
        num_sv = 1
    else:
        complextype = ctype or "other"
        # One BND per novel junction.
        num_sv = 0
        svtype_parts: list[str] = []
        for (b1, s1), (b2, s2) in zip(genome_alt, genome_alt[1:]):
            if not (s1 > 0 and s2 > 0 and b2 == b1 + 1):
                num_sv += 1
                svtype_parts.append("BND")
        if num_sv == 0:
            num_sv = 1
            svtype_parts = [simple or "BND"]
        svtype = ",".join(svtype_parts)
    boundaries = [out_blocks[out_start].start]
    uncertainty = [0]
    for b in range(out_start, out_end):
        boundaries.append(out_blocks[b].end)
        uncertainty.append(_gap_after(out_blocks, b))
    boundaries.append(out_blocks[out_end].end)
    uncertainty.append(0)
    junctions: list[tuple[int, int]] = []
    junction_ci: list[tuple[int, int]] = []
    sr_list: list[int] = []
    pe_list: list[int] = []
    for (b1, s1), (b2, s2) in zip(alt_tokens, alt_tokens[1:]):
        if out_blocks[b1].is_insertion or out_blocks[b2].is_insertion:
            if out_blocks[b2].is_insertion and not out_blocks[b1].is_insertion:
                pos = _insertion_position(out_blocks, b1, s1)
                gap = _gap_after(out_blocks, b1) if s1 > 0 else _gap_before(out_blocks, b1)
                junctions.append((pos, pos))
                junction_ci.append((gap, gap))
                sr_list.append(0)
                pe_list.append(0)
            continue
        if s1 > 0 and s2 > 0 and b2 == b1 + 1:
            continue
        pos1 = out_blocks[b1].end if s1 > 0 else out_blocks[b1].start
        pos2 = out_blocks[b2].start if s2 > 0 else out_blocks[b2].end
        junctions.append((pos1, pos2))
        u1 = _gap_after(out_blocks, b1) if s1 > 0 else _gap_before(out_blocks, b1)
        u2 = _gap_before(out_blocks, b2) if s2 > 0 else _gap_after(out_blocks, b2)
        junction_ci.append((u1, u2))
        if b1 == b2 and 0 <= b1 < len(orig_groups):
            sr, pe = _support_for_merged_tandem(graph, orig_groups[b1])
        else:
            # Map simplified endpoints back to outer original blocks.
            o1 = orig_groups[b1][-1] if s1 > 0 else orig_groups[b1][0]
            o2 = orig_groups[b2][0] if s2 > 0 else orig_groups[b2][-1]
            sr, pe = _support_for_junction(graph, o1, o2)
        sr_list.append(sr)
        pe_list.append(pe)
    if not junctions and simple == "DEL":
        missing = [b for b, _s in ref_tokens if b not in {t[0] for t in alt_tokens}]
        if missing:
            junctions.append((out_blocks[missing[0]].start, out_blocks[missing[-1]].end))
            junction_ci.append(
                (_gap_before(out_blocks, missing[0]), _gap_after(out_blocks, missing[-1]))
            )
            o_left = orig_groups[missing[0]][0]
            o_right = orig_groups[missing[-1]][-1]
            flank_l = o_left - 1 if o_left > 0 else o_left
            flank_r = o_right + 1 if o_right + 1 < len(blocks) else o_right
            sr, pe = _support_for_junction(graph, flank_l, flank_r)
            sr_list = [sr]
            pe_list = [pe]
    hom = geno.path1 == geno.path2
    gt = "HOM" if hom else "HET"
    af = 1.0 if hom else geno.rho
    # Runner-up strings stay unsimplified (orig next_best_pathstring).
    next_a = _block_string(path_nodes_to_tokens(geno.second_paths[0]), blocks, start)
    next_b = _block_string(path_nodes_to_tokens(geno.second_paths[1]), blocks, start)
    locus_id = f"{chrom}_{blocks[start].start}-{blocks[end].end}"
    minbp = min(p for pair in junctions for p in pair) if junctions else out_blocks[out_start].start
    maxbp = max(p for pair in junctions for p in pair) if junctions else out_blocks[out_end].end
    values = {
        "chrom": chrom,
        "minbp": str(minbp),
        "maxbp": str(maxbp),
        "id": locus_id,
        "svtype": svtype,
        "complextype": complextype,
        "num_sv": str(num_sv),
        "bp": ",".join(str(x) for x in boundaries),
        "bp_uncertainty": ",".join(str(x) for x in uncertainty),
        "reference": tokens_to_string(named_ref),
        "rearrangement": _block_string(alt_tokens, out_blocks, out_start),
        "len_affected": str(length),
        "filter": "PASS",
        "sv_bp": (
            ";".join(str(a) for a, _b in junctions)
            if svtype == "INS" and junctions
            else (";".join(f"{a},{b}" for a, b in junctions) if junctions else "NA")
        ),
        "sv_bp_uncertainty": (
            ";".join(str(a) for a, _b in junction_ci)
            if svtype == "INS" and junction_ci
            else (";".join(f"{a},{b}" for a, b in junction_ci) if junction_ci else "NA")
        ),
        "gt": gt,
        "af": f"{af:.6f}",
        "inslen": _format_inslen(out_blocks, ins_lengths),
        "sr_support": ",".join(str(x) for x in sr_list) if sr_list else "0",
        "pe_support": ",".join(str(x) for x in pe_list) if pe_list else "0",
        "score_vs_ref": _fmt_score(geno.score_vs_ref),
        "score_vs_next": _fmt_score(geno.score_vs_next),
        "rearrangement_next": f"{next_a}/{next_b}",
        "num_paths": str(geno.n_paths),
        # Must match simplified reference/rearrangement names (fusion-safe).
        "blocks": _blocks_column(out_blocks, out_start, out_end),
        "blocks_point": _blocks_point_column(out_blocks, out_start, out_end),
        "bp_point_source": _bp_point_source_column(out_blocks, out_start, out_end),
    }
    return CallRow(values)


# Shared genotype context for ProcessPool workers. With the fork start method
# the parent binds these before creating the pool so children inherit the graph
# via copy-on-write and each task only sends (start, end). Spawn falls back to
# an initializer that pickles the graph once per worker.
_GENO_GRAPH: Graph | None = None
_GENO_BY_BLOCK: dict[int, list] | None = None
_GENO_CHROM: str = ""
_GENO_CFG: CallConfig | None = None
_GENO_LIBS: list[LibraryModel] | None = None
_GENO_PI: float = 0.0


def _bind_geno_context(
    graph: Graph,
    chrom: str,
    cfg: CallConfig,
    libs: list[LibraryModel],
    pi: float,
    by_block: dict[int, list] | None = None,
) -> None:
    global _GENO_GRAPH, _GENO_BY_BLOCK, _GENO_CHROM, _GENO_CFG, _GENO_LIBS, _GENO_PI
    _GENO_GRAPH = graph
    _GENO_BY_BLOCK = by_block if by_block is not None else _edges_by_block(graph)
    _GENO_CHROM = chrom
    _GENO_CFG = cfg
    _GENO_LIBS = libs
    _GENO_PI = pi


def _clear_geno_context() -> None:
    global _GENO_GRAPH, _GENO_BY_BLOCK, _GENO_CHROM, _GENO_CFG, _GENO_LIBS, _GENO_PI
    _GENO_GRAPH = None
    _GENO_BY_BLOCK = None
    _GENO_CHROM = ""
    _GENO_CFG = None
    _GENO_LIBS = None
    _GENO_PI = 0.0


def _geno_worker_init(
    graph: Graph,
    chrom: str,
    cfg: CallConfig,
    libs: list[LibraryModel],
    pi: float,
) -> None:
    _bind_geno_context(graph, chrom, cfg, libs, pi)


def _evaluate_locus_span(start: int, end: int) -> tuple[list[str], CallRow | None]:
    """Worker entry: locus view from the bound shared graph, then genotype."""
    assert _GENO_GRAPH is not None and _GENO_CFG is not None and _GENO_LIBS is not None
    view = _locus_view(_GENO_GRAPH, start, end, _GENO_BY_BLOCK)
    return _evaluate_locus(view, _GENO_CHROM, start, end, _GENO_CFG, _GENO_LIBS, _GENO_PI)


def _evaluate_locus(
    graph: Graph,
    chrom: str,
    start: int,
    end: int,
    cfg: CallConfig,
    libs: list[LibraryModel],
    pi: float,
) -> tuple[list[str], CallRow | None]:
    """Genotype one locus. Safe to run in another process."""
    notes: list[str] = []
    n_blocks = end - start + 1
    if n_blocks <= 2 and not _touches_insertion(graph, start, end):
        notes.append(f"discard locus {start}-{end}: only {n_blocks} reference blocks")
        return notes, None
    paths, status, mes_extra = _enumerate_with_support(graph, cfg, start, end)
    blocks = graph.blocks
    notes.append(f"locus\t{blocks[start].start}\t{blocks[end].end}\t{status}\t{len(paths)}\t{mes_extra}")
    if status != "ok":
        return notes, None
    ref = reference_path(start, end)
    geno = genotype_locus(
        paths,
        ref,
        graph,
        libs,
        pi=pi,
        top_k=cfg.top_paths,
        allele_fractions=tuple(cfg.allele_fractions),
    )
    if geno is None:
        return notes, None
    return notes, _row_for_genotype(chrom, blocks, graph, start, end, geno, cfg)


def _enumerate_with_support(
    graph: Graph,
    cfg: CallConfig,
    start: int,
    end: int,
) -> tuple[list[tuple[int, ...]], str, int]:
    """Path search, raising edge support by 2 when the path cap is hit.

    orig:sv_inference.py retries while ``mes_extra <= max_mes_extra``. The
    paper default of ``max_mes_extra`` is 0, so the first overflow skips.
    """
    extra = 0
    while True:
        paths, status = enumerate_paths(
            graph,
            max_back=cfg.max_back_edges,
            max_paths=cfg.max_paths,
            min_support=cfg.min_edge_support + extra,
            start=2 * start,
            end=2 * end + 1,
        )
        if status != "skipped_npaths" or extra + 2 > cfg.max_mes_extra:
            return paths, status, extra
        extra += 2


def insertion_length_is_bound(chosen: int, sizes: list[int], q_high: int) -> bool:
    """True when pairs cannot separate this size from a longer insertion.

    The search grid stops near the insert-size range. The largest tested size,
    and any size past ``q_high``, is a lower bound rather than a point estimate.
    """
    if chosen <= 0 or not sizes:
        return False
    if chosen >= sizes[-1]:
        return True
    return q_high > 0 and chosen > q_high


def _stamp_clip_positions(blocks: list[Block], breakpoints: list[Breakpoint]) -> None:
    """Remember the clip coordinate in the uncovered interval that contains it."""
    for bp in breakpoints:
        if bp.clip_pos is None:
            continue
        support = bp.n_clip_left + bp.n_clip_right + bp.n_softclip
        for i in range(len(blocks) - 1):
            if blocks[i].is_insertion or blocks[i + 1].is_insertion:
                continue
            if blocks[i].end <= bp.clip_pos <= blocks[i + 1].start and support >= blocks[i].clip_reads:
                blocks[i].clip_pos = bp.clip_pos
                blocks[i].clip_reads = support


def _stamp_clip_positions_from_softclips(blocks: list[Block], clips) -> None:
    """Refine point estimates with any soft clip inside an uncovered interval.

    Breakpoint discovery needs two clipped reads to open an interval. Labeling
    still prefers a single clip in that interval over the midpoint. Clip
    positions are sorted once and looked up with ``bisect``.

    ``clips`` may be Softclip objects or a lean int/ndarray of positions (M4).
    """
    import bisect
    from collections import Counter

    from arcsv2.engine.breakpoints import Softclip

    if clips is None or len(clips) == 0:
        return
    if hasattr(clips, "dtype"):
        positions = sorted(int(x) for x in clips)
    elif isinstance(clips[0], Softclip):
        positions = sorted(clip.pos for clip in clips)
    else:
        positions = sorted(int(x) for x in clips)
    if not positions:
        return
    for i in range(len(blocks) - 1):
        if blocks[i].is_insertion or blocks[i + 1].is_insertion:
            continue
        lo, hi = blocks[i].end, blocks[i + 1].start
        if hi <= lo:
            continue
        left = bisect.bisect_left(positions, lo)
        right = bisect.bisect_right(positions, hi)
        if left >= right:
            continue
        counts = Counter(positions[left:right])
        pos, support = counts.most_common(1)[0]
        if support >= blocks[i].clip_reads:
            blocks[i].clip_pos = pos
            blocks[i].clip_reads = support


def _pair_spans_boundary(graph: Graph, after: int) -> bool:
    """True when a discordant pair has ends on both sides of the insertion.

    A concordant pair can step across the reference breakpoint without
    containing the novel sequence. Hanging and clipped reads do not identify
    the length, so a homozygous insertion longer than the fragment stays a
    lower bound.
    """
    for edge in graph.edges.values():
        left, right = edge.u // 2, edge.v // 2
        if left > right:
            left, right = right, left
        if left <= after < right and any(
            read.kind == "pe" and not read.hanging and read.discordant for read in edge.reads
        ):
            return True
    return False


def _edges_by_block(graph: Graph) -> dict[int, list]:
    """Each edge is listed once under every block it touches."""
    index: dict[int, list] = {}
    for edge in graph.edges.values():
        for block in {edge.u // 2, edge.v // 2}:
            index.setdefault(block, []).append(edge)
    return index


def _edge_in_locus(graph: Graph, edge, start: int, end: int) -> bool:
    """An edge is needed to genotype ``[start, end]``, including an appended insertion."""
    ends = (edge.u // 2, edge.v // 2)
    if any(start <= block <= end for block in ends):
        return True
    for block, other in ((ends[0], ends[1]), (ends[1], ends[0])):
        if block < len(graph.blocks) and graph.blocks[block].is_insertion and start <= other <= end:
            return True
    return False


def _locus_view(graph: Graph, start: int, end: int, by_block: dict[int, list] | None = None) -> Graph:
    """Same block indexes, with only the edges this locus genotypes.

    ``by_block`` is built once per chromosome. Workers copy edge shells so
    concurrent genotype tasks cannot touch shared mutable edge state. CSR
    read columns stay shared read-only (M4); only the Edge shell is private.
    """
    view = Graph(list(graph.blocks), edge_reads=graph.edge_reads)
    index = by_block if by_block is not None else _edges_by_block(graph)
    seen: set[tuple[int, int]] = set()
    chosen = []
    for block in range(start, end + 1):
        chosen.extend(index.get(block, ()))
    for block, edges in index.items():
        if block < len(graph.blocks) and graph.blocks[block].is_insertion:
            chosen.extend(edges)
    for edge in chosen:
        key = edge.key()
        if key in seen or not _edge_in_locus(graph, edge, start, end):
            continue
        seen.add(key)
        view.edges[key] = edge.shell_copy()
    return view


def insertion_test_sizes(width: int) -> list[int]:
    """orig:do_inference ``10 ** arange(1.5, 1 + log10(width), 0.5)``."""
    if width <= 0:
        return []
    limit = 1.0 + math.log10(max(width, 1))
    sizes: list[int] = []
    step = 1.5
    while step < limit - 1e-9:
        sizes.append(int(10**step))
        step += 0.5
    return sizes


def _blocks_near(blocks: list[Block], idx: int, width: int) -> tuple[int, int]:
    """Inclusive start and exclusive end of blocks within ``width`` of block ``idx``'s right end."""
    pos = blocks[idx].end
    first = idx
    while first > 0 and blocks[first - 1].end > pos - width:
        if blocks[first].gap_left:
            break
        first -= 1
    last = idx + 1
    while last < len(blocks) and not blocks[last].is_insertion and blocks[last].start < pos + width:
        if blocks[last].gap_left:
            break
        last += 1
    return first, last


def _boundary_index(blocks: list[Block], bp: Breakpoint) -> int | None:
    for b in range(len(blocks) - 1):
        if blocks[b].is_insertion or blocks[b + 1].is_insertion or blocks[b + 1].gap_left:
            continue
        if bp.start < blocks[b + 1].start and bp.end > blocks[b].end:
            return b
    return None


def _reads_touching(graph: Graph, first: int, last_inclusive: int) -> list:
    sorted_edges = sorted(graph.edges.values(), key=lambda item: item.key())
    return _reads_touching_sorted(graph, first, last_inclusive, sorted_edges)


def _reads_touching_sorted(
    graph: Graph,
    first: int,
    last_inclusive: int,
    sorted_edges,
    by_block: dict | None = None,
) -> list:
    if by_block is not None:
        seen: set[tuple[int, int]] = set()
        local = []
        for block in range(first, last_inclusive + 1):
            for edge in by_block.get(block, []):
                key = edge.key()
                if key in seen:
                    continue
                seen.add(key)
                local.append(edge)
        local.sort(key=lambda item: item.key())
        edges = local
    else:
        edges = sorted_edges
    index = []
    for edge in edges:
        if by_block is None:
            if edge.u // 2 > last_inclusive and edge.v // 2 > last_inclusive:
                continue
            if edge.u // 2 < first and edge.v // 2 < first:
                continue
        if edge.u // 2 >= len(graph.blocks) or edge.v // 2 >= len(graph.blocks):
            continue
        if graph.blocks[edge.u // 2].is_insertion or graph.blocks[edge.v // 2].is_insertion:
            continue
        # Keep only edges that touch the window (same filter as the linear scan).
        if edge.u // 2 > last_inclusive and edge.v // 2 > last_inclusive:
            continue
        if edge.u // 2 < first and edge.v // 2 < first:
            continue
        for i in range(len(edge.reads)):
            if not edge.reads[i].in_likelihood:
                continue
            index.append((edge, i))
    return index


def _graph_plus_insertion(graph: Graph, after: int, size: int, support: int) -> tuple[Graph, int]:
    """Temporary graph with a novel insertion after block ``after``.

    ``score_path`` / ``normalizer`` only read ``blocks``, so edges are not copied.
    """
    del support  # kept for call-site parity with committed insertions
    blocks = list(graph.blocks)
    new = len(blocks)
    blocks.append(Block(0, 0, is_insertion=True, insertion_len=size))
    return Graph(blocks), new


def _path_through(start: int, end_inclusive: int) -> tuple[int, ...]:
    nodes: list[int] = []
    for b in range(start, end_inclusive + 1):
        nodes.extend((in_node(b), out_node(b)))
    return tuple(nodes)


def _remap_node(node: int, lo: int) -> int:
    return in_node(node // 2 - lo) if node % 2 == 0 else out_node(node // 2 - lo)


def _remap_adj(adj: tuple[int, ...] | None, lo: int) -> tuple[int, ...] | None:
    if adj is None:
        return None
    return tuple(_remap_node(n, lo) for n in adj)


def _remap_edge_read(read: EdgeRead, lo: int) -> EdgeRead:
    if read.adj1 is None and read.adj2 is None:
        return read
    return EdgeRead(
        offset=read.offset,
        lib_id=read.lib_id,
        pmapped=read.pmapped,
        hanging=read.hanging,
        mapq=read.mapq,
        kind=read.kind,
        discordant=read.discordant,
        distant=read.distant,
        mate_len=read.mate_len,
        supports_adjacency=read.supports_adjacency,
        adj1=_remap_adj(read.adj1, lo),
        adj2=_remap_adj(read.adj2, lo),
        in_likelihood=read.in_likelihood,
    )


def _remap_edge_for_insertion(edge: Edge, lo: int) -> Edge:
    """Private edge shell with remapped endpoints/adj; CSR stays shared."""
    remapped = Edge(
        u=_remap_node(edge.u, lo),
        v=_remap_node(edge.v, lo),
        support=edge.support,
        reference=edge.reference,
        store=edge.store,
        read_start=edge.read_start,
        read_end=edge.read_end,
        _reads=None if edge.store is not None else list(edge._reads or ()),
    )
    if edge.store is not None:
        overrides: dict[int, tuple[tuple[int, ...] | None, tuple[int, ...] | None]] = {}
        for i in range(len(remapped.reads)):
            read = edge.reads[i]
            if read.adj1 is None and read.adj2 is None:
                continue
            overrides[i] = (_remap_adj(read.adj1, lo), _remap_adj(read.adj2, lo))
        remapped._adj_override = overrides or None
    else:
        remapped._reads = [_remap_edge_read(r, lo) for r in (edge._reads or ())]
    return remapped


def _insertion_test_window(
    graph: Graph,
    first: int,
    last: int,
    after: int,
    reads: list[tuple[Edge, int]],
) -> tuple[list[Block], list[Block], int, tuple[int, ...], tuple[int, ...], list[tuple[Edge, int]]]:
    """Local blocks + remapped paths/reads for insertion scoring.

    Includes one flanking ref block on each side when present so ``block_gap`` /
    ``node_gap`` at the window edges match the full-graph scores. Ref scoring
    uses ``ref_blocks`` (no insertion); test scoring mutates ``test_blocks[-1]``.
    """
    lo = first - 1 if first > 0 else first
    hi = last
    if last < len(graph.blocks) and not graph.blocks[last].is_insertion:
        hi = last + 1
    ref_blocks = list(graph.blocks[lo:hi])
    ins = Block(0, 0, is_insertion=True, insertion_len=0)
    test_blocks = ref_blocks + [ins]
    new_local = len(test_blocks) - 1
    first_l, after_l, last_l = first - lo, after - lo, last - lo
    ref_path = _path_through(first_l, last_l - 1)
    test_path = (
        _path_through(first_l, after_l)
        + (in_node(new_local), out_node(new_local))
        + _path_through(after_l + 1, last_l - 1)
    )
    # Remap each distinct edge once. CSR columns stay shared; adj overrides
    # live only on the private shell (avoids the M3 share bug).
    edge_map: dict[int, Edge] = {}
    local_reads: list[tuple[Edge, int]] = []
    for edge, ridx in reads:
        remapped = edge_map.get(id(edge))
        if remapped is None:
            remapped = _remap_edge_for_insertion(edge, lo)
            edge_map[id(edge)] = remapped
        local_reads.append((remapped, ridx))
    return ref_blocks, test_blocks, new_local, ref_path, test_path, local_reads


def place_novel_insertions(
    graph: Graph,
    breakpoints: list[Breakpoint],
    libs: list[LibraryModel],
    pi: float,
    cfg: CallConfig,
    progress: Progress | None = None,
) -> list[tuple[int, int]]:
    """Test novel insertions at clip or insertion-pair boundaries.

    Sizes and the haploid/heterozygous comparison follow orig:do_inference.
    A winning insertion is appended and linked around the boundary.
    """
    if not libs or not breakpoints:
        return []
    width = int(1.1 * libs[0].q_high) if libs[0].q_high else 1000
    sizes = insertion_test_sizes(width)
    if not sizes:
        return []
    added: list[tuple[int, int]] = []
    # Boundaries are decided on the reference prefix. Appending does not move them.
    n_ref = sum(1 for b in graph.blocks if not b.is_insertion)
    candidates = [
        bp
        for bp in breakpoints
        if (bp.n_clip_left > 0 and bp.n_clip_right > 0) or "Ins" in bp.drp_types
    ]
    if progress is not None:
        progress.stage(
            "insertions",
            f"{len(candidates)} candidate sites, {len(sizes)} sizes, {len(graph.edges)} edges",
        )
    # Sort edges once; _reads_touching used to re-sort per site.
    sorted_edges = sorted(graph.edges.values(), key=lambda item: item.key())
    by_block = _edges_by_block(graph)
    report_every = max(1, len(candidates) // 20) if candidates else 1
    for i, bp in enumerate(candidates):
        if progress is not None and i > 0 and i % report_every == 0:
            progress.stage("insertions", f"{i}/{len(candidates)} sites ({len(added)} kept)")
        after = _boundary_index(graph.blocks[:n_ref], bp)
        if after is None:
            continue
        first, last = _blocks_near(graph.blocks, after, width)
        if last <= first:
            continue
        reads = _reads_touching_sorted(
            graph, first, last - 1, sorted_edges, by_block=by_block
        )
        if not reads:
            continue
        ref_blocks, test_blocks, _new_local, ref_path, test_path, local_reads = (
            _insertion_test_window(graph, first, last, after, reads)
        )
        # Ref path must not see the novel insertion block (hanging-read branch).
        ref_scored = score_path(ref_path, Graph(ref_blocks), libs, local_reads)
        ref_ll = log_haploid(ref_scored, pi)
        best_ll = ref_ll
        best_size = 0
        best_hom = True
        # Cache (size -> hom, het); second pass is a lookup with the same == tie rule.
        score_cache: dict[int, tuple[float, float]] = {}
        temp = Graph(test_blocks)
        ins_block = test_blocks[-1]
        for size in sizes:
            ins_block.insertion_len = size
            test_scored = score_path(test_path, temp, libs, local_reads)
            hom = log_haploid(test_scored, pi)
            het = log_diploid(ref_scored, test_scored, pi, 0.5)
            score_cache[size] = (hom, het)
            if hom >= best_ll and hom >= het:
                best_ll = hom
                best_size = size
                best_hom = True
            elif het > best_ll and het > hom:
                best_ll = het
                best_size = size
                best_hom = False
        if best_size <= 0 or best_ll <= ref_ll:
            continue
        # Prefer the smallest size that matches the winning zygosity, as upstream does.
        winner = best_ll
        chosen = best_size
        for size in sizes:
            hom, het = score_cache[size]
            score = hom if best_hom else het
            if score == winner:
                chosen = size
                break
        new = len(graph.blocks)
        graph.blocks.append(
            Block(
                0,
                0,
                is_insertion=True,
                insertion_len=chosen,
                length_is_lower_bound=(
                    insertion_length_is_bound(chosen, sizes, libs[0].q_high)
                    or not _pair_spans_boundary(graph, after)
                ),
            )
        )
        graph.add_edge(out_node(after), in_node(new), support=cfg.min_edge_support, reference=True)
        graph.add_edge(out_node(new), in_node(after + 1), support=cfg.min_edge_support, reference=True)
        added.append((after, chosen))
    if progress is not None:
        progress.stage("insertions done", f"{len(added)} insertions placed")
    return added


def _touches_insertion(graph: Graph, start: int, end: int) -> bool:
    for edge in graph.edges.values():
        ends = (edge.u // 2, edge.v // 2)
        if not any(i < len(graph.blocks) and graph.blocks[i].is_insertion for i in ends):
            continue
        if any(start <= i <= end and not graph.blocks[i].is_insertion for i in ends):
            return True
    return False


def _block_string(tokens: list[tuple[int, int]], blocks: list[Block], start: int) -> str:
    parts: list[str] = []
    for block, strand in tokens:
        if block < len(blocks) and blocks[block].is_insertion:
            parts.append("_")
            continue
        parts.append(block_idx_to_name(block, start))
        if strand < 0:
            parts.append("'")
    return "".join(parts)


def call_evidence(
    pairs: list[PairRecord] | Sequence,
    *,
    chrom: str,
    chrom_len: int,
    cfg: CallConfig | None = None,
    gaps: list[tuple[int, int]] | None = None,
    libs: list[LibraryModel] | None = None,
    breakpoints: list[Breakpoint] | None = None,
    extra_breakpoints: list[Breakpoint] | None = None,
    splits: list[SplitHit] | None = None,
    softclips=None,
    progress: Progress | None = None,
) -> tuple[list[CallRow], list[str]]:
    """Call SVs from pair records already extracted from one region."""
    cfg = cfg or CallConfig()
    rng = random.Random(cfg.seed)
    if progress is not None:
        progress.stage("fit_libs", f"{len(pairs):,} pairs")
    libs = libs or _fit_libs(pairs, cfg)
    lib = libs[0]
    # One π for the sample: median across libraries (orig:compute_pi_robust).
    pi = pi_robust([item.pmf for item in libs], cfg.robustness_parameter)
    if breakpoints is None:
        breakpoints, logs = discover_breakpoints(
            pairs, lib, chrom_len, cfg, rng, extra=extra_breakpoints, progress=progress
        )
    else:
        logs = []
    intervals = [(bp.start, bp.end) for bp in breakpoints]
    blocks = create_blocks(intervals, chrom_len, gaps)
    if progress is not None:
        progress.stage(
            "blocks",
            f"{len(breakpoints):,} breakpoints -> {len(blocks):,} blocks "
            f"(chrom_len={chrom_len:,})",
        )
    if len(blocks) < 2:
        return [], logs
    # R2: rust block prep + build_graph when engine=rust. Missing native is an
    # error — only --engine python uses the pure-Python graph builder.
    graph = None
    if cfg.engine == "rust":
        from arcsv2 import native as native_mod

        if not native_mod.build_graph_available():
            raise RuntimeError(
                "engine=rust requires arcsv2_native.build_graph "
                "(rebuild arcsv2:dev or pass --engine python)"
            )
        from arcsv2.engine.graph_native import build_graph_rust

        if progress is not None:
            progress.stage(
                "build_graph",
                f"{len(pairs):,} pairs x {len(blocks):,} blocks (rust)",
            )
        graph = build_graph_rust(
            blocks,
            pairs,
            lib,
            cfg,
            splits,
            breakpoints=breakpoints,
            softclips=[] if softclips is None else softclips,
        )
        if progress is not None:
            progress.stage("build_graph done", f"{len(graph.edges):,} edges")
    if graph is None:
        _stamp_clip_positions(blocks, breakpoints)
        _stamp_clip_positions_from_softclips(
            blocks, [] if softclips is None else softclips
        )
        graph = build_graph(blocks, pairs, lib, cfg, splits, progress=progress)
    # M3/M4: pairs/splits/softclips are fully consumed into the graph. Drop them
    # before insertions/genotype so fork workers inherit a leaner address space.
    import gc

    del pairs
    splits = None
    softclips = None
    gc.collect()
    try:
        import ctypes

        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass
    if progress is not None:
        progress.stage("insertions", "test novel insertion sites")
    inserted = place_novel_insertions(graph, breakpoints, libs, pi, cfg, progress=progress)
    gc.collect()
    try:
        import ctypes

        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass
    width = int(1.1 * lib.q_high) if lib.q_high else 1000
    loci = decompose(graph, cfg.min_edge_support)
    for after, _size in inserted:
        loci.append((after, after + 1))
    expanded: list[tuple[int, int]] = []
    for start, end in loci:
        expanded.append(expand_locus(blocks, start, end, width))
    # Open intervals: loci that only share a cut-point block stay separate.
    # orig:sv_inference.py merges pyinter.open(start, end), so a shared endpoint
    # is not an overlap. Chaining on ``<=`` collapses a chromosome into one locus.
    expanded.sort()
    merged: list[tuple[int, int]] = []
    for start, end in expanded:
        if merged and start < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    rows: list[CallRow] = []
    pending = [
        (start, end)
        for start, end in merged
        if not (end - start + 1 <= 2 and not _touches_insertion(graph, start, end))
    ]
    for start, end in merged:
        if (start, end) not in pending:
            logs.append(f"discard locus {start}-{end}: only {end - start + 1} reference blocks")
    if progress is not None:
        progress.stage("genotype", f"{len(pending)} loci (of {len(merged)} merged)")
    workers = cfg.threads if cfg.threads > 1 and len(pending) > 1 else 1
    if workers == 1:
        evaluated = []
        for i, (start, end) in enumerate(pending):
            if progress is not None and (i == 0 or (i + 1) % 25 == 0 or i + 1 == len(pending)):
                progress.stage(
                    "genotype",
                    f"locus {i + 1}/{len(pending)} blocks {start}-{end}",
                )
            evaluated.append(_evaluate_locus(graph, chrom, start, end, cfg, libs, pi))
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Fork inherits the bound graph (no per-locus pickle). Spawn pickles once
        # per worker via initializer. Workers build a locus view in-process.
        # M3: collect after releasing pairs so children CoW a leaner heap.
        import gc

        gc.collect()
        start_methods = mp.get_all_start_methods()
        use_fork = "fork" in start_methods
        ctx = mp.get_context("fork" if use_fork else "spawn")
        by_block = _edges_by_block(graph)
        _bind_geno_context(graph, chrom, cfg, libs, pi, by_block=by_block)
        initargs = () if use_fork else (graph, chrom, cfg, libs, pi)
        initializer = None if use_fork else _geno_worker_init
        in_flight_cap = max(workers * 2, workers + 1)
        try:
            with ProcessPoolExecutor(
                max_workers=workers,
                mp_context=ctx,
                initializer=initializer,
                initargs=initargs,
            ) as pool:
                evaluated = [None] * len(pending)
                future_to_i: dict = {}
                next_i = 0
                done = 0
                while next_i < len(pending) or future_to_i:
                    while next_i < len(pending) and len(future_to_i) < in_flight_cap:
                        start, end = pending[next_i]
                        future_to_i[pool.submit(_evaluate_locus_span, start, end)] = next_i
                        next_i += 1
                    for fut in as_completed(future_to_i):
                        idx = future_to_i.pop(fut)
                        evaluated[idx] = fut.result()
                        done += 1
                        if progress is not None and (
                            done == 1 or done % 25 == 0 or done == len(pending)
                        ):
                            progress.stage(
                                "genotype", f"locus {done}/{len(pending)} done"
                            )
                        break
            evaluated = list(evaluated)
        finally:
            _clear_geno_context()
    # Emit in subgraph (start, end) order so threads=1/4/8 write identical tables.
    for (_start, _end), (notes, row) in sorted(
        zip(pending, evaluated), key=lambda item: item[0]
    ):
        logs.extend(notes)
        if row is not None:
            rows.append(row)
    rows.sort(key=lambda r: (r.values["chrom"], int(r.values["minbp"]), r.values["rearrangement"]))
    if progress is not None:
        progress.stage("genotype done", f"{len(rows)} calls")
    return rows, logs


def write_calls(rows: list[CallRow], path: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            data = row.as_dict()
            handle.write("\t".join(data[col] for col in OUTPUT_COLUMNS) + "\n")
