# Cursor task: duplications called with 3 copies (last big gap)

## Status after v2fix4 (chr20, HG01884 + NA19700)
| | v2fix3 | **v2fix4** | original |
|:--|--:|--:|--:|
| Validated calls | 160 | 151 | 243 |
| Original TPs recovered | 162 | 163 | — |
| Deletion precision | 33% | **50%** | 67% |
| Duplication precision | 6/50 | **7/68 (10%)** | 66/164 (40%) |
| Tool-only false calls | 310 | **152** | — |
| HG01884 call time | 60 min | **44 min** | 24 min |

Deletions are close; duplications are now the gap (7 vs 66 validated). 38 ARC-SV2 DUPs are ≥3-copy.
The traced loci are still 3-copy **and still carry pair support the original does not have**:
| locus (HG01884) | ARC-SV2 v2fix4 | original |
|:--|:--|:--|
| 3067783 | ABCCCD, sr 27, **pe 31** | ABBC HOM, sr 27, pe 0 |
| 2659058 | ABBBC, sr 11, **pe 28** | ABBC HET, sr 11, pe 0 |
| 306278 | ABBBCD, sr 2, **pe 27** | ABBC HET, sr 2, pe 0 |

## Status: v2fix5 (batch finished)
The chr20 batch finished (`results_v2fix5`, exit 0). The live comparison is the window section at the bottom. Implemented:
1. Same-block gate literal `add_pair` (no `_is_outward`); pe junction match only counts the within-block edge → pe=0 for split-only dups.
2. `COPY_FAMILY_BAND=5.0` prefers lower copy number; HOM vs HET still uses the narrow 0.5-nat band.
3. `CallConfig.skip_edge_ref_gate` (default True).
4. Unit tests (`test_split_only_tandem_dup_reports_pe_zero_and_two_copies`, copy-band).
5. chr20 batch finished. Export is `export/v2fix5/`.
6. Open: original `-v 2` at weak DEL loci 23857951 / 13358681 / 63458085.

## Review of the v2fix5 changes (after the run; don't interrupt it)
`COPY_FAMILY_BAND = 5.0` is a new tie-break, not in the original. The original has no band: it takes the
argmax, and its likelihood clearly prefers 2 copies. Original `arcsv_out.tab` (HG01884):
| locus | original call | score_vs_ref | score_vs_next | next |
|:--|:--|--:|--:|:--|
| 2659058 | ABBC HET | 60.55 | **8.0** (52.54) | ABBBC/ABC |
| 3067783 | ABBC HOM | 110.84 | 3.69 | ABCCD/ABCD |
| 306278 | ABBC HET | 8.77 | 8.77 | ABCD/ABCD |
ARC-SV2 dump at 306278: best 3-copy beats 2-copy by **1.24 nats** and beats ref by only 2.63, with sr=2.
At 2659058 ARC-SV2 vs_ref 19.85 vs original 60.55. So the **likelihood values differ**; the band only hides it.

Next step — likelihood oracle test (decisive):
1. For one locus (2659058), fix the blocks and the read set, and score the same candidate haplotypes
   (ABC, ABBC, ABBBC; HET/HOM combos) with the original's functions (`sv_inference` likelihood /
   `compute_edge_likelihood`, normalizer, hanging-read terms) and with ARC-SV2 `score_path`/`log_diploid`.
2. Compare per-path log-likelihoods and per-read contributions; find the first term that differs
   (within-block read placements / `which_dup`, normalizer for duplicated blocks, n in G^n, robust π —
   see PARITY.md item 6).
3. Fix that term, then set `COPY_FAMILY_BAND` back to `LIKELIHOOD_BAND` (or remove it) and confirm the
   locus picks ABBC with score_vs_ref close to 60.
Keep `skip_edge_ref_gate` documented as an ARC-SV2 deviation in PARITY.md.

## v2fix5 result (export/v2fix5)
Dup calls 68 → 159 labeled, still **10% precision** (16/159). 117 of the 143 false DUPs have **no split
reads** (pair-only, median score_vs_ref 13, ~160 bp, mostly HET). The band fixed 2659058 and 306278;
3067783 is still ABCCCD. Complex calls rose 30 → 52 (4 validated).

## Which change made the pair-only DUPs
HG01884 chr20:59,444,167–64,444,167 (the 5 Mb with the most pair-only DUPs). One discovery, then
each v2fix5 change toggled. Script: `scratch/gap_analysis/ablate_v2fix5_window.py`.

| setting | calls | DUP | pair-only DUP (sr=0) |
|:--|--:|--:|--:|
| v2fix5 | 88 | 18 | 13 |
| outward check restored | 78 | **7** | **2** |
| `COPY_FAMILY_BAND` = 0.5 | 88 | 18 | 13 |
| `skip_edge_ref_gate` off | 105 | 15 | 10 (DELs 27→52) |
| old pe reporter (any touching edge) | 88 | 18 | 12 |

The pair-only calls come from dropping the outward test on `B.out↔B.in`, not from the pe
column and not from the copy-number band. The old reporter only inflates printed `pe`
(median 3→22) on the same 18 DUPs. `skip_edge_ref_gate` off adds weak deletions.
The outward test is restored in `_pair_supported`. It is not in literal `add_pair`; recorded
in `PARITY.md` item 11 and `EDGE_SUPPORT.md`. `skip_edge_ref_gate` is item 10.

## Likelihood oracle (2659058, same blocks and reads)
Window chr20:2.4–2.9 Mb, outward check on. The locus is blocks
2652733–2659058 / 2659059–2659109 / 2659110–2663187 (ABC). ARC-SV2 calls
**ABBBC HET, score_vs_ref 24.93**. Original full-chr call is ABBC HET, 60.55.

`compute_normalizing_constant` and `compute_edge_likelihood` from
`/data/orig-src`, fed this graph's blocks, insert pmf, and non-hanging reads,
match ARC-SV2 to the printed digits:

| path | G | haploid (n=1295, no hanging) | HET vs ref (all 1297 reads) |
|:--|--:|--:|--:|
| ABC | 10008.9 | −19866.160 | 0 |
| ABBC | 10059.9 | −19847.453 | 13.90 |
| ABBBC | 10110.9 | −19843.585 | **24.93** |

3-copy leads 2-copy by **11.0 nats**, so `COPY_FAMILY_BAND` (5.0) does not apply
and is not what this window calls. The diploid wrapper is not the gap either:
haploid HOM and ρ=0.5 HET reduce to `log_diploid` when p and G match, and they do.
The 60-vs-25 gap is upstream of `score_path` (breakpoints, which reads sit on
which edges, or the chromosome-wide library). `COPY_FAMILY_BAND` stays until
that upstream difference is fixed; removing it here would not produce ABBC or a
score near 60. Script: `scratch/gap_analysis/oracle_dup_likelihood.py`.

## Review of the ablation + oracle (done; result below)
Good result: the likelihood **code** matches the original given identical inputs. So the 25-vs-60 gap at
2659058 is in the **inputs**. The library-sampling guess below was tested and is not the cause;
the same-window section is the current result. The original does sample ~1M pairs genome-wide
(`extract_approximate_library_stats`) and then refits the pmf on the region.

Do next, at 2659058:
1. Run the original on the same window with `-v 2`; take its `logging/lib1_insert_pmf.txt`, its blocks,
   and its read→edge assignments (`add_pair` / `parse_pair` extracts already in `_orig_extract/`).
2. Diff against ARC-SV2 in this order: (a) insert pmf (and q_low/q_high, π), (b) block boundaries,
   (c) reads per edge — names present in one tool only, with the filter that dropped them
   (MAPQ, read-through, distant, hanging, `min_mapq_reads`).
3. Feed the original's inputs into ARC-SV2 `score_path`: if ABBC wins with ~60, the input that differs is the
   cause. Fix it (likely: genome-wide sampled library like the original), then retest the band at 0.5.
Keep the outward check (PARITY item 11) — the ablation shows it removes the pair-only DUPs.

## Same-window comparison (2659058)
Original `-v 2` on chr20:2400000-2900000 (HG01884), the oracle window.
Output: `scratch/gap_analysis/loci/orig_2659058_win/`. Script: `run_orig_2659058_win.sh`.

The original calls **ABBC HET, score_vs_ref 57.74**, next ABBBC/ABC at a margin of 50.86
(full chromosome was 60.55 / 8.0). sr 11, pe 0. Blocks:

| | A | B | C |
|:--|:--|:--|:--|
| original | 2652732–2659058 | 2659058–2659109 | 2659109–2663187 |
| ARC-SV2 | 2652733–2659058 | 2659059–2659109 | 2659110–2663187 |

Same three blocks (the start differs by the 1-based breakpoint coordinate).
7 edges, **1287 reads**, 78 informative. Reads with likelihood 0:
ABC 13, ABBC 17, ABBBC 23. ARC-SV2 on this window still calls ABBBC HET at 24.93
(1297 reads, 5 zeros on every path).

(a) Insert pmf is not the gap. Window fit q=[268, 691], π=1.40e-6.
The chromosome-wide `lib1_insert_pmf.txt` has q=[269, 689], π=1.22e-6, KL=0.0003.
Rescoring the ARC-SV2 graph with that pmf (`swap_pmf_score.py`) leaves ABBBC ahead:
24.96 vs ABBC 14.52. A genome-wide library will not produce ABBC or a score near 60.

(b) Block boundaries match, above.

(c) Read placement does not. Two original rules are missing from `score_path`:

- A read that overlaps several blocks stores `adj = block_seq_to_path` of that
  overlap (`add_pair`). A reverse read is anchored on its rightmost block
  (`get_overlapping_blocks` reverses the list). `placement_distances` keeps every
  opposite-parity placement and ignores `is_adj_satisfied`.
- A split is one read's block path (first segment + second segment), an adjacency
  constraint on the pair edge with its mate. ARC-SV2 stores the two halves as their
  own edge. The 11 Dup splits at 2659058 pass the junction check, but several have
  a start a few bases into block A (`2659051`), so `find_block(start)` connects A–B
  instead of requiring the B tandem path. Rebuilding endpoints from overlap intervals
  (`adj_rescore.py`) recovers **1287 reads and 78 adjacency reads**, and still calls
  ABBBC (with the filter applied, both alts fall below the reference). The split is
  not yet joined to its mate, so this is not the original's edge list.

Do next (done): each read stores its block path, a split is both halves on the
pair with its mate, and a placement is kept only when `is_adj_satisfied`.
Retest of this window (`scratch/gap_analysis/retest_2659058.py`), band still 5.0:

| | rearrangement | gt | score_vs_ref | vs next | zeros ABC / ABBC / ABBBC | sr |
|:--|:--|:--|--:|--:|:--|--:|
| original window | ABBC | HET | 57.74 | 50.86 | 13 / 17 / 23 | 11 |
| ARC-SV2 retest | ABBC | HET | 60.51 | 50.05 | 13 / 15 / 21 | 3 |

1276 likelihood reads, 70 with an adjacency, 0 splits scored a second time
(11 split edges are support-only). `COPY_FAMILY_BAND` is back to 0.5. The
50-nat margin means the band does not choose this call. Outward check stays.

The printed `sr` follows that path. Window retest (`retest_dup_windows.py`), HG01884:

| locus | original | this window |
|:--|:--|:--|
| 2659058 | ABBC HET, sr 11, pe 0, 57.74 (full chr 60.55) | ABBC HET, sr 11, pe 0, 60.51, next ABBBC/ABC by 50.05 |
| 3067783 | ABBC HOM, sr 27, pe 0, 110.84, next ABCCD/ABCD by 3.69 | ABCCD HOM, sr 27, pe 0, 104.84, next ABCCD/ABCD by 3.93 |
| 306278 | ABBC HET, sr 2, pe 0, 8.77, next ABCD/ABCD | ABCBCD HET, sr 4, pe 0, 10.44, next ABCD/ABCD |

2659058 matches, including sr. 3067783 has the original split count and a score near 111, but the winning collapsed family is still the extra C (ABCCD), not ABBC. 306278 is a different rearrangement (ABCBCD, sr 4). `COPY_FAMILY_BAND` is 0.5. Outward check stays. Full chr20 is `label.v2fix6.yaml`.
