"""Annotation features (§7.2 Tier C). Tracks are BED intervals, not hard filters."""

from __future__ import annotations

from arcsv2.merge_filter import _frac, _overlap_bp


def annotate_span(
    start: int,
    end: int,
    *,
    segdups: list[tuple[int, int]] | None = None,
    simple_repeats: list[tuple[int, int]] | None = None,
    repeats: dict[str, list[tuple[int, int]]] | None = None,
    low_mappability: list[tuple[int, int]] | None = None,
) -> dict[str, float]:
    repeats = repeats or {}
    features = {
        "frac_segdup": _frac(start, end, segdups or []),
        "frac_simple_repeat": _frac(start, end, simple_repeats or []),
        "frac_low_mappability": _frac(start, end, low_mappability or []),
    }
    for name in ("SINE", "LINE", "LTR", "satellite"):
        features[f"frac_{name}"] = _frac(start, end, repeats.get(name, []))
    # Overlap of a ±0 window is the breakpoint point itself; callers pass the
    # uncertainty-expanded window when they have one.
    features["bp_overlap_segdup"] = float(_overlap_bp(start, end, segdups or []) > 0)
    return features
