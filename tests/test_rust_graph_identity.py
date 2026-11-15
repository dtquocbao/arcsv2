"""Bit-identity: Rust build_graph vs Python on Edge / EdgeRead fields."""

from __future__ import annotations

import os
import random

import pytest

from arcsv2.config import CallConfig
from arcsv2.engine.caller import (
    _fit_libs,
    _stamp_clip_positions,
    _stamp_clip_positions_from_softclips,
    build_graph as build_graph_py,
    discover_breakpoints,
)
from arcsv2.engine.graph import EdgeRead, Graph, create_blocks
from arcsv2.engine.graph_native import build_graph_rust
from arcsv2.engine.pairs import Pairs
from arcsv2.samio import load_region


def _read_tuple(r: EdgeRead) -> tuple:
    return (
        r.offset,
        r.lib_id,
        r.pmapped,
        r.hanging,
        r.mapq,
        r.kind,
        r.discordant,
        r.distant,
        r.mate_len,
        r.supports_adjacency,
        r.adj1,
        r.adj2,
        r.in_likelihood,
    )


def _compare_graphs(py: Graph, rs: Graph) -> list[str]:
    diffs: list[str] = []
    if len(py.blocks) != len(rs.blocks):
        diffs.append(f"n_blocks {len(py.blocks)} != {len(rs.blocks)}")
        return diffs
    for i, (a, b) in enumerate(zip(py.blocks, rs.blocks)):
        if a.clip_pos != b.clip_pos or a.clip_reads != b.clip_reads:
            diffs.append(
                f"block[{i}] clip py=({a.clip_pos},{a.clip_reads}) "
                f"rs=({b.clip_pos},{b.clip_reads})"
            )
            if len(diffs) > 20:
                return diffs
    py_keys = list(py.edges.keys())
    rs_keys = list(rs.edges.keys())
    if py_keys != rs_keys:
        first = "none"
        for i in range(max(len(py_keys), len(rs_keys))):
            pk = py_keys[i] if i < len(py_keys) else None
            rk = rs_keys[i] if i < len(rs_keys) else None
            if pk != rk:
                first = f"i={i} py={pk} rs={rk}"
                break
        diffs.append(
            f"edge key order/set: py={len(py_keys)} rs={len(rs_keys)} first_diff={first}"
        )
    shared = [k for k in py_keys if k in rs.edges]
    missing = [k for k in py_keys if k not in rs.edges]
    extra = [k for k in rs_keys if k not in py.edges]
    if missing:
        diffs.append(f"missing in rust ({len(missing)}): {missing[:5]}")
    if extra:
        diffs.append(f"extra in rust ({len(extra)}): {extra[:5]}")
    for k in shared:
        pe = py.edges[k]
        re = rs.edges[k]
        if pe.support != re.support or pe.reference != re.reference:
            diffs.append(
                f"edge{k} support/ref py=({pe.support},{pe.reference}) "
                f"rs=({re.support},{re.reference})"
            )
        if len(pe.reads) != len(re.reads):
            diffs.append(f"edge{k} n_reads py={len(pe.reads)} rs={len(re.reads)}")
            continue
        for ri, (a, b) in enumerate(zip(pe.reads, re.reads)):
            ta, tb = _read_tuple(a), _read_tuple(b)
            if ta != tb:
                diffs.append(f"edge{k} read[{ri}] py={ta} rs={tb}")
                if len(diffs) > 40:
                    return diffs
    return diffs


def _build_both(bam: str, chrom: str, start: int | None, end: int | None, threads: int):
    from arcsv2 import native as native_mod

    if not native_mod.rust_available() or not native_mod.build_graph_available():
        pytest.skip("rust build_graph not available")

    cfg = CallConfig(engine="rust", columnar_pairs=True, threads=threads)
    pairs, extras, splits, _chrom_name, chrom_len, softclips = load_region(
        bam, chrom, cfg, start, end
    )
    assert isinstance(pairs, Pairs)
    libs = _fit_libs(pairs, cfg)
    lib = libs[0]
    breakpoints, _logs = discover_breakpoints(
        pairs, lib, chrom_len, cfg, random.Random(cfg.seed), extra=extras
    )
    intervals = [(bp.start, bp.end) for bp in breakpoints]
    blocks_py = create_blocks(intervals, chrom_len, gaps=None)
    blocks_rs = create_blocks(intervals, chrom_len, gaps=None)
    if len(blocks_py) < 2:
        pytest.skip("too few blocks in window")

    _stamp_clip_positions(blocks_py, breakpoints)
    _stamp_clip_positions_from_softclips(
        blocks_py, [] if softclips is None else softclips
    )
    g_py = build_graph_py(blocks_py, pairs, lib, cfg, splits)

    g_rs = build_graph_rust(
        blocks_rs,
        pairs,
        lib,
        cfg,
        splits,
        breakpoints=breakpoints,
        softclips=[] if softclips is None else softclips,
    )
    return g_py, g_rs


@pytest.mark.skipif(
    __import__("importlib.util").util.find_spec("arcsv2_native") is None,
    reason="arcsv2_native not installed",
)
@pytest.mark.parametrize("rust_threads", [1, 8])
def test_rust_python_graph_identity_chr20_window(rust_threads: int):
    bam = "/data/work/HG01884/chr20.bam"
    chrom, start, end = "chr20", 10_000_000, 11_000_000
    g_py, g_rs = _build_both(bam, chrom, start, end, rust_threads)
    diffs = _compare_graphs(g_py, g_rs)
    assert not diffs, (
        f"Graph mismatch (threads={rust_threads}):\n  " + "\n  ".join(diffs[:30])
    )


@pytest.mark.skipif(
    __import__("importlib.util").util.find_spec("arcsv2_native") is None,
    reason="arcsv2_native not installed",
)
def test_rust_python_graph_identity_chr20_full_optional():
    """Full chr20 identity (slow). Enable with ARCSV2_R2_FULL=1."""
    if os.environ.get("ARCSV2_R2_FULL") != "1":
        pytest.skip("set ARCSV2_R2_FULL=1 to run full chr20 identity")
    bam = "/data/work/HG01884/chr20.bam"
    g_py, g_rs = _build_both(bam, "chr20", None, None, 8)
    diffs = _compare_graphs(g_py, g_rs)
    assert not diffs, "Full chr20 Graph mismatch:\n  " + "\n  ".join(diffs[:50])
