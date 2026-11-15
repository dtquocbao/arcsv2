"""Rust extension ``arcsv2_native`` (required when ``engine=rust``).

Import may fail if the wheel is missing; callers must use ``--engine python``
in that case. There is no silent fallback from rust → python.
"""

from __future__ import annotations

from typing import Any

_NATIVE: Any | None
try:
    import arcsv2_native as _NATIVE
except ImportError:  # pragma: no cover - wheel missing until image rebuild
    _NATIVE = None


def rust_available() -> bool:
    return _NATIVE is not None and bool(getattr(_NATIVE, "rust_engine_available", lambda: False)())


def native_version() -> str | None:
    if _NATIVE is None:
        return None
    return str(_NATIVE.native_version())


def load_bam_pairs(
    path: str,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    reference: str | None = None,
    max_pair_distance: int = 2_000_000,
    threads: int = 1,
    filter_high_nm: bool = False,
    max_nm: int = 5,
    filter_both_mates_clipped: bool = False,
) -> dict[str, Any]:
    """Call the Rust BAM loader. Raises RuntimeError if the extension is missing."""
    if _NATIVE is None or not hasattr(_NATIVE, "load_bam_pairs"):
        raise RuntimeError("arcsv2_native.load_bam_pairs is not available")
    return _NATIVE.load_bam_pairs(
        path,
        chrom,
        start,
        end,
        reference,
        int(max_pair_distance),
        int(threads),
        bool(filter_high_nm),
        int(max_nm),
        bool(filter_both_mates_clipped),
    )


def build_graph_available() -> bool:
    return _NATIVE is not None and hasattr(_NATIVE, "build_graph")


def build_graph(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Call the Rust build_graph kernel. Raises RuntimeError if missing."""
    if _NATIVE is None or not hasattr(_NATIVE, "build_graph"):
        raise RuntimeError("arcsv2_native.build_graph is not available")
    return _NATIVE.build_graph(*args, **kwargs)
