"""Paper merge-filter (§4.9). Interval overlap is computed on sorted arrays."""

from __future__ import annotations


def _overlap_bp(start: int, end: int, intervals: list[tuple[int, int]]) -> int:
    """Overlap of half-open ``[start, end)`` with a sorted interval list."""
    if end <= start or not intervals:
        return 0
    total = 0
    for a, b in intervals:
        if b <= start:
            continue
        if a >= end:
            break
        total += max(0, min(end, b) - max(start, a))
    return total


def _frac(start: int, end: int, intervals: list[tuple[int, int]]) -> float:
    span = end - start
    if span <= 0:
        return 0.0
    return _overlap_bp(start, end, intervals) / span


def inslen_is_lower_bound(text: str | None) -> bool:
    """True when ``inslen`` is a minimum (``>N``) rather than a point estimate."""
    return ">" in (text or "")


def passes_filter(
    start: int,
    end: int,
    affected: int,
    *,
    simple_repeats: list[tuple[int, int]] | None = None,
    segdups: list[tuple[int, int]] | None = None,
    blacklist: list[tuple[int, int]] | None = None,
    gaps: list[tuple[int, int]] | None = None,
    min_affected: int = 50,
    length_is_lower_bound: bool = False,
) -> tuple[bool, str]:
    """Return whether one call survives §4.9, and the first reason it failed.

    A censored insertion length is a minimum, so it is not dropped for being short.
    """
    if affected < min_affected and not length_is_lower_bound:
        return False, "short"
    if _frac(start, end, simple_repeats or []) > 0.85:
        return False, "simple_repeat"
    if _frac(start, end, segdups or []) > 0.50:
        return False, "segdup"
    if _overlap_bp(start, end, blacklist or []) > 0:
        return False, "blacklist"
    if _overlap_bp(start, end, gaps or []) > 0:
        return False, "gap"
    return True, "PASS"


def reciprocal_overlap(a0: int, a1: int, b0: int, b1: int) -> float:
    inter = max(0, min(a1, b1) - max(a0, b0))
    if inter == 0:
        return 0.0
    return min(inter / (a1 - a0), inter / (b1 - b0))


def merge_calls(
    records: list[dict],
    *,
    n_samples: int,
) -> list[dict]:
    """Merge calls with ≥ 50% reciprocal overlap and the same class.

    Allele frequency is allele count / (2 × samples). ``records`` need
    ``chrom, start, end, sv_class, gt``.
    """
    grouped: dict[tuple, list[dict]] = {}
    for rec in records:
        grouped.setdefault((rec["chrom"], rec.get("sv_class") or rec.get("complextype") or rec.get("svtype")), []).append(rec)
    merged: list[dict] = []
    for (_chrom, _cls), group in grouped.items():
        group.sort(key=lambda r: (r["start"], r["end"]))
        clusters: list[list[dict]] = []
        for rec in group:
            placed = False
            for cluster in clusters:
                anchor = cluster[0]
                if reciprocal_overlap(rec["start"], rec["end"], anchor["start"], anchor["end"]) >= 0.5:
                    cluster.append(rec)
                    placed = True
                    break
            if not placed:
                clusters.append([rec])
        for cluster in clusters:
            # One dose per sample. A repeated call from the same sample keeps
            # the larger genotype (homozygous = 2, heterozygous = 1).
            per_sample: dict[str, int] = {}
            for i, rec in enumerate(cluster):
                sample = rec.get("sample")
                key = str(sample) if sample not in (None, "") else f"#{i}"
                dose = 2 if rec.get("gt") == "HOM" else 1
                per_sample[key] = max(per_sample.get(key, 0), dose)
            ac = sum(per_sample.values())
            out = dict(cluster[0])
            out["allele_count"] = ac
            out["allele_frequency"] = ac / (2 * n_samples) if n_samples else 0.0
            out["n_samples_called"] = len(per_sample)
            merged.append(out)
    merged.sort(key=lambda r: (r["chrom"], r["start"], r["end"]))
    return merged
