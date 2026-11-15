"""Bit-identity: Rust BAM loader vs Python samio.load_region on Pairs SoA."""

from __future__ import annotations

import numpy as np
import pytest

from arcsv2.config import CallConfig
from arcsv2.engine.pairs import Pairs
from arcsv2.samio import load_region


def _arrays_equal(a: Pairs, b: Pairs) -> list[str]:
    diffs: list[str] = []
    if len(a) != len(b):
        diffs.append(f"len {len(a)} != {len(b)}")
        return diffs
    for name in (
        "start1",
        "end1",
        "start2",
        "end2",
        "insert",
        "mate_len",
        "qname_hash",
        "mapq1",
        "mapq2",
        "lib_id",
        "flags",
    ):
        xa = getattr(a, name)
        xb = getattr(b, name)
        if not np.array_equal(xa, xb):
            mism = int(np.sum(xa != xb))
            # First mismatch index for debugging.
            idx = int(np.flatnonzero(xa != xb)[0]) if mism else -1
            diffs.append(
                f"{name}: {mism} mismatches"
                + (f" first@{idx} py={xa[idx]} rust={xb[idx]}" if idx >= 0 else "")
            )
    return diffs


@pytest.mark.skipif(
    __import__("importlib.util").util.find_spec("arcsv2_native") is None,
    reason="arcsv2_native not installed",
)
@pytest.mark.parametrize("rust_threads", [1, 8])
def test_rust_python_pairs_identity_chr20_window(rust_threads: int):
    from arcsv2 import native as native_mod

    if not native_mod.rust_available():
        pytest.skip("rust engine not available")

    bam = "/data/work/HG01884/chr20.bam"
    # 1 Mb window near the start of chr20 (0-based half-open).
    chrom, start, end = "chr20", 10_000_000, 11_000_000
    cfg_py = CallConfig(engine="python", columnar_pairs=True, threads=1)
    cfg_rs = CallConfig(engine="rust", columnar_pairs=True, threads=rust_threads)

    pairs_py, extras_py, splits_py, c1, e1, clips_py = load_region(
        bam, chrom, cfg_py, start, end
    )
    pairs_rs, extras_rs, splits_rs, c2, e2, clips_rs = load_region(
        bam, chrom, cfg_rs, start, end
    )
    assert isinstance(pairs_py, Pairs)
    assert isinstance(pairs_rs, Pairs)
    assert c1 == c2 == chrom
    diffs = _arrays_equal(pairs_py, pairs_rs)
    assert not diffs, f"Pairs SoA mismatch (rust threads={rust_threads}):\n  " + "\n  ".join(
        diffs
    )
    assert len(splits_py) == len(splits_rs)
    assert len(clips_py) == len(clips_rs)
    assert len(extras_py) == len(extras_rs)


@pytest.mark.skipif(
    __import__("importlib.util").util.find_spec("arcsv2_native") is None,
    reason="arcsv2_native not installed",
)
def test_rust_python_pairs_identity_chr20_full_optional():
    """Full chr20 identity (slow). Enable with ARCSV2_R1_FULL=1."""
    import os

    if os.environ.get("ARCSV2_R1_FULL") != "1":
        pytest.skip("set ARCSV2_R1_FULL=1 to run full chr20 identity")

    from arcsv2 import native as native_mod

    if not native_mod.rust_available():
        pytest.skip("rust engine not available")

    bam = "/data/work/HG01884/chr20.bam"
    cfg_py = CallConfig(engine="python", columnar_pairs=True, threads=1)
    cfg_rs = CallConfig(engine="rust", columnar_pairs=True, threads=8)
    pairs_py, _, _, _, _, _ = load_region(bam, "chr20", cfg_py)
    pairs_rs, _, _, _, _, _ = load_region(bam, "chr20", cfg_rs)
    diffs = _arrays_equal(pairs_py, pairs_rs)
    assert not diffs, "Full chr20 Pairs mismatch:\n  " + "\n  ".join(diffs)
