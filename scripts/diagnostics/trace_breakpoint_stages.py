#!/usr/bin/env python3
"""Trace ARC-SV2 breakpoint stages on a small region (current caller path)."""

from __future__ import annotations

import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path

from arcsv2.config import CallConfig
from arcsv2.engine.breakpoints import (
    Breakpoint,
    clamp_breakpoint,
    cluster_drps,
    cluster_to_intervals,
    fdr_keep,
    merge_breakpoints,
)
from arcsv2.engine.caller import (
    _drp_from_pair,
    _evaluate_locus,
    _fit_libs,
    _touches_insertion,
    build_graph,
    place_novel_insertions,
)
from arcsv2.engine.graph import create_blocks, decompose, expand_locus
from arcsv2.engine.library import pi_robust
from arcsv2.pipeline import fasta_length
from arcsv2.samio import load_region


def fmt_bp(bp: Breakpoint) -> str:
    return (
        f"[{bp.start},{bp.end}) clip={bp.n_softclip} split={bp.n_split} "
        f"drp={bp.n_drp} types={bp.drp_types} clusters={bp.cluster_ids}"
    )


def dump(label: str, bps: list[Breakpoint], focus: tuple[int, int] | None, out: list[str]) -> None:
    out.append(f"## {label} (n={len(bps)})")
    if focus is None:
        show = bps
    else:
        lo, hi = focus
        show = [bp for bp in bps if not (bp.end < lo or bp.start > hi)]
        out.append(f"   (showing {len(show)} overlapping focus {lo}-{hi})")
    for bp in show[:80]:
        out.append(f"   {fmt_bp(bp)}")
    if len(show) > 80:
        out.append(f"   ... {len(show) - 80} more")
    out.append("")


def discover_raw(pairs, lib, chrom_len, cfg, rng, extras):
    grouped: dict[str, list] = defaultdict(list)
    for pair in pairs:
        drp = _drp_from_pair(pair, lib)
        if drp is not None:
            grouped[drp.kind].append(drp)
    pieces: list[Breakpoint] = list(extras)
    logs: list[str] = []
    cluster_id = 0
    for _kind, items in grouped.items():
        clusters, clog = cluster_drps(items, lib, cfg, rng)
        logs.extend(clog)
        keep = fdr_keep(clusters, lib, cfg, rng)
        for idx, cluster in enumerate(clusters):
            if idx not in keep:
                continue
            left, right = cluster_to_intervals(cluster, cluster_id=cluster_id)
            cluster_id += 1
            pieces.extend((left, right))
    return pieces, logs, grouped


def run_pipeline(
    sam: str,
    reference: str,
    gaps_bed: str | None,
    chrom: str,
    start0: int,
    end0: int,
    focus: tuple[int, int],
) -> str:
    cfg = CallConfig()
    chrom_len = fasta_length(reference, chrom) or end0
    pairs, extras, splits, _name, _inf, softclips = load_region(
        sam, chrom, cfg, start0, end0, reference
    )
    lines: list[str] = [
        f"# region {chrom}:{start0 + 1}-{end0}",
        f"pairs={len(pairs)} extras={len(extras)} splits={len(splits)} softclips={len(softclips)}",
        "",
    ]
    dump("clip/split extras (pre-DRP)", extras, focus, lines)
    libs = _fit_libs(pairs, cfg)
    lib = libs[0]
    pi = pi_robust([item.pmf for item in libs], cfg.robustness_parameter)
    rng = random.Random(cfg.seed)
    pieces, logs, grouped = discover_raw(pairs, lib, chrom_len, cfg, rng, extras)
    lines.append(
        "DRP counts: " + ", ".join(f"{k}={len(v)}" for k, v in sorted(grouped.items()))
    )
    for log in logs[:20]:
        lines.append(f"  cluster_log: {log}")
    dump("after FDR DRP + extras (pre-merge)", pieces, focus, lines)

    clamped = [item for item in (clamp_breakpoint(bp, chrom_len) for bp in pieces) if item is not None]
    breakpoints = merge_breakpoints(clamped, min_split_or_clip=cfg.min_edge_support)
    dump("after merge_breakpoints (precedence)", breakpoints, focus, lines)

    gaps = []
    if gaps_bed and Path(gaps_bed).is_file():
        from arcsv2.beds import read_bed

        gaps = read_bed(gaps_bed, chrom)
    intervals = [(bp.start, bp.end) for bp in breakpoints]
    blocks = create_blocks(intervals, chrom_len, gaps)
    lines.append(f"## blocks n={len(blocks)}")
    for i, b in enumerate(blocks):
        if b.end < focus[0] or b.start > focus[1]:
            continue
        lines.append(f"   block{i} [{b.start},{b.end}) len={b.end - b.start}")
    lines.append("")

    if len(blocks) < 2:
        lines.append("CALLS: none (fewer than 2 blocks)")
        return "\n".join(lines)

    graph = build_graph(blocks, pairs, lib, cfg, splits)
    place_novel_insertions(graph, breakpoints, libs, pi, cfg)
    width = int(1.1 * lib.q_high) if lib.q_high else 1000
    loci = decompose(graph, cfg.min_edge_support)
    expanded = sorted(expand_locus(blocks, s, e, width) for s, e in loci)
    merged_loci: list[tuple[int, int]] = []
    for s, e in expanded:
        if merged_loci and s < merged_loci[-1][1]:
            merged_loci[-1] = (merged_loci[-1][0], max(merged_loci[-1][1], e))
        else:
            merged_loci.append((s, e))
    lines.append(f"## loci n={len(merged_loci)}")
    for s, e in merged_loci:
        if blocks[e].end < focus[0] or blocks[s].start > focus[1]:
            continue
        lines.append(f"   locus blocks {s}-{e} span {blocks[s].start}-{blocks[e].end}")
    pending = [
        (s, e)
        for s, e in merged_loci
        if not (e - s + 1 <= 2 and not _touches_insertion(graph, s, e))
    ]
    rows = []
    for s, e in pending:
        _notes, row = _evaluate_locus(graph, chrom, s, e, cfg, libs, pi)
        if row is not None:
            rows.append(row)
            lines.append(
                f"CALL {row.values.get('svtype')} {row.values.get('minbp')}-"
                f"{row.values.get('maxbp')} gt={row.values.get('gt')} "
                f"sr={row.values.get('sr_support')} pe={row.values.get('pe_support')} "
                f"rearr={row.values.get('rearrangement')}"
            )
    lines.append(f"CALLS total={len(rows)}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--gaps", default="")
    parser.add_argument("--chrom", default="chr20")
    parser.add_argument(
        "--center",
        type=int,
        required=True,
        help="focus midpoint (1-based ok; treated as coord)",
    )
    parser.add_argument("--pad", type=int, default=20_000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    center = args.center
    start0 = max(0, center - args.pad)
    end0 = center + args.pad
    focus = (center - 2_000, center + 2_000)
    text = run_pipeline(
        args.bam,
        args.reference,
        args.gaps or None,
        args.chrom,
        start0,
        end0,
        focus,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text + "\n", encoding="utf-8")
    sys.stdout.write(text + "\n")


if __name__ == "__main__":
    main()
