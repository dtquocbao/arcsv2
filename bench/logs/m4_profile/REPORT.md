# M4.1 baseline profile (R2 engine, no code change)

git: d8ab20b (caller R2)
date: 2026-10-05

## Full chr1 HG01884 @ threads=8 (engine=rust)

Died during rust build_graph (process became zombie; likely OOM while
rebuilding EdgeRead objects). Stage marks before death:

| Stage | Wall (s) | Peak RSS (MB) |
|:--|--:|--:|
| after_load | 934 | **12986** (12.7 GB) |
| after_fit_libs | 1034 | 12986 |
| after_blocks | 1105 | 12986 |

| Metric | Value |
|:--|--:|
| pairs | 29,449,918 |
| bytes/pair (SoA payload) | **53.0** |
| pairs payload | 1.56 GB |
| splits / extras | 10,086 / 85,664 |
| breakpoints / blocks | 67,231 / 67,232 |

Load RSS is ~8x the pair payload: take-record Python dicts (cigar/qual) and
Softclip objects dominate before graph build. EdgeRead rebuild on top of that
is what pushed the process over the edge.

## Window chr1:30-35 Mb (threads=1, tracemalloc)

| Stage | Wall (s) | Peak RSS (MB) |
|:--|--:|--:|
| after_load | 16.4 | 319 |
| after_fit_libs | 85 | 323 |
| after_blocks | 143 | 323 |
| after_build_graph | 231 | **631** |
| after_drop_pairs | 264 | 805 (tracemalloc overhead) |
| after_insertions | 359 | 805 |

| Metric | Value |
|:--|--:|
| pairs | 612,791 |
| bytes/pair | **53.0** |
| edges | 4,200 |
| EdgeRead count | **665,942** |
| bytes/EdgeRead | **448** |
| EdgeRead total est | **284.5 MB** |

Scaled EdgeRead count for full chr1: ~32.0M x 448 B ≈ **13.3 GB** object
footprint if R2 conversion kept.

### Top allocation sites after build_graph (tracemalloc, window)

| # | Size (MB) | Count | Site |
|--:|--:|--:|:--|
| 1 | 116.9 | 1,331,884 | graph_native.py:195 EdgeRead(...) |
| 2 | 19.0 | 623,246 | EdgeRead field boxes (offset) |
| 3 | 15.2 | 665,942 | EdgeRead pmapped boxes |
| 4 | 5.5 | 4,193 | edge.reads.append |
| 5 | 5.2 | 134,543 | adj tuple materialization |

## Implication

1. Keep R2 CSR as storage; drop EdgeRead objects (saves ~13 GB on chr1).
2. Free takes immediately after softclip/split extract; avoid retaining Softclip
   objects when only positions are needed (cuts the 12.7 GB load floor).
