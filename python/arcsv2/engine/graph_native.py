"""R2/M4: attach native CSR columns as Graph.edge_reads (no EdgeRead objects)."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from arcsv2.config import CallConfig
from arcsv2.engine.graph import (
    KIND_PE,
    Block,
    Edge,
    EdgeReadStore,
    Graph,
)
from arcsv2.engine.library import LibraryModel
from arcsv2.engine.pairs import Pairs


def _blocks_to_arrays(blocks: list[Block]) -> dict[str, np.ndarray]:
    n = len(blocks)
    start = np.empty(n, dtype=np.int64)
    end = np.empty(n, dtype=np.int64)
    is_ins = np.zeros(n, dtype=np.uint8)
    ins_len = np.zeros(n, dtype=np.int32)
    gap_left = np.zeros(n, dtype=np.uint8)
    clip_pos = np.full(n, -1, dtype=np.int64)
    clip_reads = np.zeros(n, dtype=np.int32)
    for i, b in enumerate(blocks):
        start[i] = b.start
        end[i] = b.end
        is_ins[i] = 1 if b.is_insertion else 0
        ins_len[i] = b.insertion_len
        gap_left[i] = 1 if b.gap_left else 0
        clip_pos[i] = -1 if b.clip_pos is None else int(b.clip_pos)
        clip_reads[i] = b.clip_reads
    return {
        "start": start,
        "end": end,
        "is_insertion": is_ins,
        "insertion_len": ins_len,
        "gap_left": gap_left,
        "clip_pos": clip_pos,
        "clip_reads": clip_reads,
    }


def _splits_to_arrays(splits: Sequence | None) -> dict[str, np.ndarray]:
    if not splits:
        empty_i64 = np.empty(0, dtype=np.int64)
        empty_u8 = np.empty(0, dtype=np.uint8)
        return {
            "start1": empty_i64,
            "end1": empty_i64,
            "reverse1": empty_u8,
            "start2": empty_i64,
            "end2": empty_i64,
            "reverse2": empty_u8,
            "mapq": empty_u8,
            "primary_start": empty_i64,
            "primary_reverse": empty_u8,
            "qname_hash": empty_i64,
        }
    n = len(splits)
    start1 = np.empty(n, dtype=np.int64)
    end1 = np.empty(n, dtype=np.int64)
    reverse1 = np.zeros(n, dtype=np.uint8)
    start2 = np.empty(n, dtype=np.int64)
    end2 = np.empty(n, dtype=np.int64)
    reverse2 = np.zeros(n, dtype=np.uint8)
    mapq = np.zeros(n, dtype=np.uint8)
    primary_start = np.empty(n, dtype=np.int64)
    primary_reverse = np.zeros(n, dtype=np.uint8)
    qname_hash = np.empty(n, dtype=np.int64)
    for i, hit in enumerate(splits):
        start1[i] = hit.start1
        end1[i] = hit.end1
        reverse1[i] = 1 if hit.reverse1 else 0
        start2[i] = hit.start2
        end2[i] = hit.end2
        reverse2[i] = 1 if hit.reverse2 else 0
        mapq[i] = hit.mapq & 0xFF
        primary_start[i] = hit.primary_start
        primary_reverse[i] = 1 if hit.primary_reverse else 0
        qname_hash[i] = hit.qname_hash
    return {
        "start1": start1,
        "end1": end1,
        "reverse1": reverse1,
        "start2": start2,
        "end2": end2,
        "reverse2": reverse2,
        "mapq": mapq,
        "primary_start": primary_start,
        "primary_reverse": primary_reverse,
        "qname_hash": qname_hash,
    }


def _bp_clip_arrays(breakpoints: Sequence | None) -> tuple[np.ndarray, np.ndarray]:
    if not breakpoints:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    pos_list: list[int] = []
    sup_list: list[int] = []
    for bp in breakpoints:
        if bp.clip_pos is None:
            continue
        pos_list.append(int(bp.clip_pos))
        sup_list.append(int(bp.n_clip_left + bp.n_clip_right + bp.n_softclip))
    return (
        np.asarray(pos_list, dtype=np.int64),
        np.asarray(sup_list, dtype=np.int32),
    )


def _softclip_pos_array(softclips) -> np.ndarray:
    from arcsv2.engine.breakpoints import Softclip

    if softclips is None or len(softclips) == 0:
        return np.empty(0, dtype=np.int64)
    if isinstance(softclips, np.ndarray):
        return np.asarray(softclips, dtype=np.int64)
    if isinstance(softclips[0], Softclip):
        return np.asarray([clip.pos for clip in softclips], dtype=np.int64)
    return np.asarray(softclips, dtype=np.int64)


def _pairs_as_soa(pairs: Pairs | Sequence) -> Pairs:
    if isinstance(pairs, Pairs):
        return pairs
    return Pairs.from_records(pairs)


def graph_from_columns(blocks: list[Block], cols: dict) -> Graph:
    """Attach native CSR columns; Edge shells index into the shared store."""
    # Apply stamped clip fields returned by rust block prep.
    bcp = cols["block_clip_pos"]
    bcr = cols["block_clip_reads"]
    for i, block in enumerate(blocks):
        cp = int(bcp[i])
        block.clip_pos = None if cp < 0 else cp
        block.clip_reads = int(bcr[i])

    store = EdgeReadStore(
        offset=np.asarray(cols["read_offset"]),
        lib_id=np.asarray(cols["read_lib_id"]),
        pmapped=np.asarray(cols["read_pmapped"]),
        mapq=np.asarray(cols["read_mapq"]),
        kind=np.asarray(cols["read_kind"]),
        flags=np.asarray(cols["read_flags"]),
        mate_len=np.asarray(cols["read_mate_len"]),
        adj1_ptr=np.asarray(cols["read_adj1_ptr"]),
        adj1_len=np.asarray(cols["read_adj1_len"]),
        adj2_ptr=np.asarray(cols["read_adj2_ptr"]),
        adj2_len=np.asarray(cols["read_adj2_len"]),
        adj_flat=np.asarray(cols["adj_flat"]),
    )

    edge_u = cols["edge_u"]
    edge_v = cols["edge_v"]
    edge_support = cols["edge_support"]
    edge_reference = cols["edge_reference"]
    edge_read_ptr = cols["edge_read_ptr"]

    graph = Graph(blocks, edge_reads=store)
    n_edges = len(edge_u)
    for ei in range(n_edges):
        u = int(edge_u[ei])
        v = int(edge_v[ei])
        key = (min(u, v), max(u, v))
        r0 = int(edge_read_ptr[ei])
        r1 = int(edge_read_ptr[ei + 1])
        graph.edges[key] = Edge(
            key[0],
            key[1],
            support=int(edge_support[ei]),
            reference=bool(edge_reference[ei]),
            store=store,
            read_start=r0,
            read_end=r1,
        )
    return graph


def build_graph_rust(
    blocks: list[Block],
    pairs: Pairs | Sequence,
    lib: LibraryModel,
    cfg: CallConfig,
    splits: Sequence | None = None,
    *,
    breakpoints: Sequence | None = None,
    softclips=None,
) -> Graph:
    """Call native build_graph and attach CSR columns (no EdgeRead rebuild)."""
    from arcsv2 import native as native_mod

    soa = _pairs_as_soa(pairs)
    ba = _blocks_to_arrays(blocks)
    sa = _splits_to_arrays(splits)
    bp_pos, bp_sup = _bp_clip_arrays(breakpoints)
    sc_pos = _softclip_pos_array(softclips)

    cols = native_mod.build_graph(
        ba["start"],
        ba["end"],
        ba["is_insertion"],
        ba["insertion_len"],
        ba["gap_left"],
        ba["clip_pos"],
        ba["clip_reads"],
        soa.start1[: len(soa)],
        soa.end1[: len(soa)],
        soa.start2[: len(soa)],
        soa.end2[: len(soa)],
        soa.insert[: len(soa)],
        soa.mate_len[: len(soa)],
        soa.qname_hash[: len(soa)],
        soa.mapq1[: len(soa)],
        soa.mapq2[: len(soa)],
        soa.lib_id[: len(soa)],
        soa.flags[: len(soa)],
        sa["start1"],
        sa["end1"],
        sa["reverse1"],
        sa["start2"],
        sa["end2"],
        sa["reverse2"],
        sa["mapq"],
        sa["primary_start"],
        sa["primary_reverse"],
        sa["qname_hash"],
        bp_pos,
        bp_sup,
        sc_pos,
        float(lib.pmapped),
        int(lib.q_low),
        int(lib.q_high),
        int(cfg.min_mapq),
        int(cfg.min_edge_support),
        int(cfg.max_edge_span),
        int(cfg.split_read_leeway),
        bool(cfg.pool_libraries),
        bool(cfg.skip_edge_ref_gate),
    )
    # Rust copied the pair/split arrays; drop Python SoA buffers now so peak
    # RSS is not (pairs + CSR) while Edge shells are attached.
    import gc

    if isinstance(pairs, Pairs):
        pairs.release()
    del sa, ba, bp_pos, bp_sup, sc_pos
    gc.collect()
    return graph_from_columns(blocks, cols)
