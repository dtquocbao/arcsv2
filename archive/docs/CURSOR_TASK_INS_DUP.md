# Cursor task: insertion-pair explosion + missing split-read tandem duplications

## Result of c6b9dce (export/v2fix2, chr20 HG01884 + NA19700)
| | v2fix | **v2fix2** | original |
|:--|--:|--:|--:|
| Validated calls | 157 | 151 | 243 |
| Deletion precision | 35% | 36% | 67% |
| Duplications validated/labeled | 3/18 | 2/18 | 66/164 |
| Original TPs recovered | 131 | 132 | — |
| HG01884 call time | 106 min | 81 min | 24 min |

DRP parity (HG01884), DRPs → clusters passing:
| Type | ARC-SV2 v2fix2 | original | status |
|:--|:--|:--|:--|
| Del | 22,026 → 279 | 24,373 → 298 | close |
| Dup | 11,376 → 295 | 11,047 → 295 | **matches** |
| InvL / InvR | 639 → 43 / 507 → 40 | 705 → 47 / 578 → 45 | **matches** |
| Ins | **266,834 → 596** | 9,695 → 19 | **27× too many DRPs** |

## Status (693fade → results_v2fix3)
Implemented A + B; unit tests green; commit `693fade`. Discover probe (HG01884 chr20):
| Type | probe DRPs | original DRPs |
|:--|--:|--:|
| Ins | **7,106** | 9,695 |
| Dup | 10,797 | 11,047 |
| Del | 22,217 | 24,373 |

Ins is in the right order of magnitude (was 266k). Full `label.v2fix3.yaml` batch running (`arcsv2-v2fix3-batch`). Fill the results table + passing-cluster counts from `graph_log.txt` when done.

## A. Insertion DRPs (done)
Removing the non-overlap gate made every overlapping / read-through pair an Ins DRP. Ported:
1. `pecluster.process_discordant_pair:39–41` — early return for opposite-orientation with
   `z_ins ≤ ilen ≤ z_del`, before Del/Dup/Ins (`caller._drp_from_pair`).
2. `bamparser_streaming.process_insert_len` — signed outer insert + soft-clip adjust
   (`samio.insert_length`). KDE keeps `ilen ≥ 0`.
3. `filter_read_through` left **False** (original default).

## B. Tandem duplications supported by split reads (done in code)
`_add_split` no longer drops `b1 == b2`. Same-block Dup+/Dup− → `B.out ↔ B.in` without the PE 95% insert gate. Recovery of the five small HOM/HET dups is checked on the v2fix3 label run.

## C. Then
Unit tests for A and B ✓, commit ✓, rerun chr20 both samples (`-j 1`) → `results_v2fix3/` **in progress**.
Remaining after that: 231 ARC-SV2-only false calls (212 pair-only deletions, median score_vs_ref ~2.5) —
analyse after A/B since Ins clusters may be feeding them.
