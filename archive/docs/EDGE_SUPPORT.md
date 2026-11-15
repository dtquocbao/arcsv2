# Pair → edge support (CURSOR_TASK_EDGE_SUPPORT)

Original: `sv_parse_reads.py` `GenomeGraph.add_pair`.
ARC-SV2: `caller.build_graph` / `_pair_supported`.

## Rules

| Case | Original | ARC-SV2 |
|:--|:--|:--|
| Record pair on edge | Always | Always |
| Between-mate support | Implied insert under that adjacency ∈ `[q_low, q_high]` | Same |
| Same-block `B.out↔B.in` | Additionally: no-dup insert `alt_dist < q_low` (one-sided). No orientation test | Same alt_dist test, **and** the reverse read must sit upstream of the forward read |
| Novel skip `\|i−k\| > 1` | Literal `add_pair` has **no** extra gate | **Also** require insert under the reference path through middle blocks to fall **outside** `[q_low, q_high]` |
| Within-read (split) | Always support | Always support |
| `pe_support` column | Typed DRP intersection at BPs (`sv_classify.supporting_read_count`) | Counts PE reads with `supports_adjacency` (not every stored within-block pair) |

## Why the skip gate

Literal `add_pair` credits a concordant FR pair that jumps an ~80 bp middle block to the skip edge (implied insert under the skip is still in the 95% window). That floods `edge.support`, path search emits weak pair-only DELs (score_vs_ref ~2.5, pe tens, sr=0) that the original never calls. Extending the Methods tandem-dup idea (“insert under *no* novel adjacency outside 95%”) to skip edges matches original *behavior*.

## Same-block dup pe inflation

Tandem-dup edges are mostly split-supported. Concordant within-block pairs are stored with `support=0`.
`_support_for_junction(B,B)` must match only the within-block edge — counting any edge that *touches* B
(ref flanks) inflated pe to tens while original reported 0. That reporter does not create calls: on
chr20:59.4–64.4 Mb the old and new reporters both emit 18 DUPs; only the printed `pe` column changes
(median 3 vs 22 on the pair-only subset).

## Outward orientation

Literal `add_pair` will credit a non-outward pair when the duplication-path insert is inside
`[q_low, q_high]` and the no-dup offset is below `q_low`. On a block about one insert long, a
concordant FR pair does that and opens a tandem edge. Those edges are the v2fix5 pair-only DUP
calls (same window, each change toggled):

| setting | DUP calls | pair-only (sr = 0) |
|:--|--:|--:|
| v2fix5 (no outward check) | 18 | 13 |
| outward check restored | 7 | 2 |
| `COPY_FAMILY_BAND` = 0.5 | 18 | 13 |
| `skip_edge_ref_gate` off | 15 | 10 |
| old pe reporter | 18 | 12 |

The outward check is what removes them. It is not in `add_pair`.

## Config

`CallConfig.skip_edge_ref_gate` (default `False`) toggles the novel-skip ref-concordance gate.
Default is literal `add_pair` (gate off). Set `True` to restore the ARC-SV2-only veto.

## flip_parity between junction (kept) + point cuts (ablated off)

Original merges are zero-width, so blocks are contiguous. ARC-SV2 half-open softclip
CIs leave 1-bp gaps. v2fix9 forced point cuts (`breakpoint_block_cuts` → contiguous
blocks) and recovered same-block DUP support, but collapsed INS and fragmented
tandem DUPs into ABCBC complex calls. Default path is gapped again
(`use_point_block_cuts=False`); point cuts remain as an ablation toggle.

Original `add_pair` appends `flip_parity` before the between-mate junction. For
multi-block mates that junction is often `B.out↔B.in` while the stored PE edge is a
skip. ARC-SV2 credits that junction (`_between_junction_supported`) with the
`alt_dist < q_low` rule. The outward orientation guard applies only when **both**
mates are single-block (keeps the pair-only DUP flood at 2 on chr20:59.4–64.4 Mb).
Same-block support under gapped geometry needs a 1-bp gap fix in vertex/path
checks — not contiguous blocks.

## Why original skips the weak pair-only DELs (23857951 / 13358681 / 63458085)

Original `-v 2` on ±20 kb (HG01884) **does** credit skip-edge PE literally and **does**
open a locus with a DEL path — it still emits nothing:

| site | skip edge support | locus / paths | best genotype | next (DEL) | `arcsv_out.tab` |
|:--|--:|:--|:--|:--|:--|
| 23857951–079 | `5.out→7.in` **27** | blocks 5–8; ABCD vs ACD | ABCD HOM | ACD HET (−7 nats) | empty |
| 13358681–745 | `1.out→3.in` **23** | blocks 1–4; ABCD vs ACD | ABCD HOM | ACD HET (−11 nats) | empty |
| 63458085–141 | `66.out→68.in` **25** | blocks 66–68; ABC vs AC | ABC HOM | AC HET (−22 nats) | INS elsewhere only |

Softclip BPs pass `min_bp_support`; FDR/path-search are not the killer. Filter is
**genotype: reference preferred** (`score_vs_ref ≤ 0`). ARC-SV2 has the same emit
guard; with `skip_edge_ref_gate=False`, window `call_evidence` also emits 0 near
these three. Prefer gate off (PARITY item 10) over a narrowed clip-aware veto.
