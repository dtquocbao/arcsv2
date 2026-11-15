# M1 memory profile

git: 9475836 (caller-v1.4-fast lineage / C5 docs HEAD)
date: 2026-10-05

## Summary

Python `PairRecord` objects dominate memory. Full chr20 HG01884 load holds
**8,087,727** pairs at **~451 bytes/pair** (getsizeof object + attrs), ~3.5 GB
just for the pair list estimate, **RSS 3.93 GB** after load alone. Target for
M2 is **< 100 bytes/pair**.

## Full chr20 load (no tracemalloc)

| Metric | Value |
|:--|--:|
| alignments | 16,625,580 |
| pairs | 8,087,727 |
| PairRecord sizeof | 56 B |
| attrs (getsizeof) | 395 B |
| est bytes/pair | 451 |
| pairs total est | 3546 MB |
| RSS after load | 3932 MB |
| wall | 266.6 s |
| splits | 4,413 |
| clip/split bps | 33,288 |

Artifacts: `scratch/m1_profile/full_load.json`

## Stage RSS (5 Mb window chr20:30-35 Mb, threads=1, no tracemalloc)

Normal `arcsv2 call` with 1 Hz VmRSS sampler. Window has ~750k pairs.

| Stage | Peak RSS (MB) | Notes |
|:--|--:|:--|
| load | 715 | BAM evidence |
| fit_libs | 728 | |
| block_prep (discover+blocks) | 699 | |
| build_graph | 774 | |
| insertions | 795 | |
| genotype | 829 | peak for window |

Wall for window call: **145.5 s** (load 41 s, build_graph ~14 s, insertions ~58 s, genotype ~27 s).

Artifacts: `scratch/m1_profile/rss_samples.tsv`, `rss_call.log`

## Top 10 allocation sites (tracemalloc, window, after load)

From `snap_loaded.json` (749,918 pairs; tracemalloc adds overhead so RSS here is
higher than the RSS-only window run).

| # | Size (MB) | Count | Site |
|--:|--:|--:|:--|
| 1 | 142.6 | 1,495,672 | `samio.emit_pair` -> `PairRecord(...)` |
| 2 | 47.1 | 1,542,888 | `samio._from_pysam` / BAM decode |
| 3 | 22.8 | 747,836 | PairRecord field `end2` / int boxes |
| 4 | 22.8 | 747,836 | PairRecord field `end1` |
| 5 | 22.8 | 747,836 | `qname_hash` / sha256 path |
| 6 | 21.7 | 712,219 | `insert_length` temps |
| 7 | 17.1 | 558,927 | `_lib_id` hashing |
| 8 | 7.6 | 124,532 | softclip / evidence.py |
| 9 | 6.8 | 111,858 | softclip take(a) |
| 10 | 6.8 | 111,032 | softclip take(b) |

After `blocks` the top-10 list is unchanged: pair objects still dominate before
`build_graph` grows the edge/read structures.

## Bytes per pair (window, after load)

- n_pairs: 749,918
- est bytes/pair: **451**
- pairs total est: 329 MB
- RSS (with tracemalloc): 1763 MB after load / 2046 MB at blocks

## Implication for M2

Replace per-pair Python objects with numpy SoA (~52 B payload/row). Free BAM
`pending` mate buffer after load (already cleared in the M2 path). Keep thin
`PairView` only where attribute access is required.
