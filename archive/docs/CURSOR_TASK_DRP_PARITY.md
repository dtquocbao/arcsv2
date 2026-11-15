# Cursor task: DRP classification parity (inversion, insertion, duplication)

## Where we are (chr20, HG01884 + NA19700, labeled with the same haplotype labeler)
| | before fix | after precedence + Del mean (`results_v2fix`) | original |
|:--|--:|--:|--:|
| Validated calls | 106 | **157** | 243 |
| Deletion precision | 24% | **35%** | 67% |
| Duplication validated / labeled | 5/36 | 3/18 | 66/164 |
| Original true calls recovered | 85/243 | **131/243** | — |
| ARC-SV2-only calls (precision) | 671 (8%) | 597 (9%) | — |
| HG01884 call time | 44 min | **106 min** | 24 min |

Still missed: 63 original true duplications (52 with split reads), 48 deletions.
ARC-SV2-only false calls: 261, almost all pair-only (sr=0) ~100 bp deletions, low score_vs_ref (median 2.5).

## The divergence is now in DRP classification (stage 1–3), HG01884 chr20
| Type | original DRPs → clusters passing | ARC-SV2 DRPs → clusters passing |
|:--|:--|:--|
| Del | 24,373 → 298 pass / 2,747 fail | 22,026 → 279 / 2,288 (similar) |
| Dup | 11,047 → 295 / 179 | **8,968** → 250 / 133 (fewer DRPs) |
| Ins | 9,695 → **19** / 0 | 4,808 → **475** / 0 |
| InvL | **705** → 47 / 0 | **7,622** → 1,047 / 0 |
| InvR | **578** → 45 / 0 | **7,397** → 1,015 / 0 |

ARC-SV2 finds 11–13× more inversion DRPs and 25× more insertion clusters. That is ~2,500 extra
breakpoints (26,249 vs 15,919 before), which plausibly explains the extra pair-only calls, the over-merged
locus at 5962348 (it chained through InvL clusters 1,2), and the 2.4× slowdown. Duplication DRPs are
under-counted, matching the missed duplications.

Sources (original, read only): `bamparser_streaming.py` (parse_bam, process_inverted, process_discordant /
how pairs are typed, mapq/flags/read-through filters), `invertedreads.py`, `pecluster.py` (type assignment,
compatibility rules per type incl. insertion microhomology h=20, min cluster size, `cluster_pairs`, null).
Cell Methods p. 17–19 (DRP types, likelihood-ratio cutoffs, clustering, FDR).

## Tasks
1. For each type, list the exact rule that makes a pair a DRP of that type in the original (orientation,
   insert cutoffs, mapq, duplicates/secondary/supplementary, same-chrom distance, "read-through" /
   overlapping mates, inverted-pair handling) vs ARC-SV2 (`samio.py`, `engine/evidence.py`,
   `caller.discover_breakpoints`). Cite file:line.
2. On chr20:5.9–6.1 Mb (HG01884) dump DRP read names per type from both tools and diff the sets:
   which pairs ARC-SV2 calls InvL/InvR/Ins that the original does not, and why (one rule per bucket).
   Same for Dup pairs the original has and ARC-SV2 lacks.
3. Same diff for clustering: why 4,808 Ins DRPs give 475 clusters vs the original's 19 (min size?
   compatibility? singleton removal? largest-clique only?).
4. Fix ARC-SV2 to match the original for each rule found (flags only if the Methods say otherwise),
   with unit tests. Then rerun chr20 both samples (`-j 1`), relabel, export to `export/v2fix2/`, and
   fill the table above (DRPs/clusters per type, calls, validated, precision, recovered TPs, time).
5. Commit before the run; record the row in `bench/RESULTS.md` with the revision.

Do not touch the original repo.

## Review notes after c6b9dce (do after the v2fix2 run finishes; don't interrupt it)
Two typing rules still differ from `pecluster.process_discordant_pair` / `bamparser_streaming.process_insert_len`:

1. **Concordant early return missing.** Original returns None first for any opposite-orientation pair with
   `min_insert <= ilen <= max_insert` (pecluster.py:39–41), *before* the Dup geometry test. ARC-SV2
   `_drp_from_pair` (caller.py:113–128) tests Dup geometry regardless of insert, so nested / read-through
   pairs with a normal insert become Dup DRPs.
2. **Insert length definition.** Original `ilen` is signed: `rev.reference_end − fwd.reference_start`,
   plus read-trim and 5′ soft-clip adjustments (bamparser_streaming.py:350–373). Outward Dup pairs get a
   small or negative `ilen`. ARC-SV2 uses `abs(TLEN)` (samio.py:151). This changes Ins/Dup typing and the
   insert KDE used for all cutoffs.

Port both exactly (unit tests: outward pair → Dup with negative ilen; nested pair with normal ilen → None;
5′-clipped pair ilen matches original), then compare per-type DRP counts with the original's log
(Del 24,373 / Dup 11,047 / Ins 9,695 / InvL 705 / InvR 578 on HG01884) before another full rerun.
