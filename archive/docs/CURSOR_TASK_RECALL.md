# Cursor task: recall — affected length, then the never-called SVs

## v2fix7 status (chr20, HG01884 + NA19700)
Labels are **identical to v2fix6** (182 validated, 65% precision, 169/243 original TPs, 48 tool-only FPs),
so the build_graph speed fix preserved output. HG01884 call 43 min (build_graph 2.9 min). Merge-filter BEDs
are still empty (min-size 50 only), same as the original's labels — comparison is like for like.

## Finding: short filtered original TPs (len_affected gap)
Originally noted as ~30 of 74 missed TPs; recount on v2fix7 labels: **35 of 75** missed orig TPs
are present in `results_v2fix7/*/chr20/arcsv_out.tab` with `len_affected` < 50, so
`merge-filter --min-size 50` removes them. ARC-SV2 measured them smaller than the original:

| sample | call | original len_affected | ARC-SV2 (raw) | mid-gap fix | ARC-SV2 rearr |
|:--|:--|--:|--:|--:|:--|
| HG01884 | DEL 613786–613840 | 54 | 47 | **54** | ACD |
| HG01884 | DEL 6606020–6606074 | 54 | 27 | **55** | ACD |
| HG01884 | DEL 9168924–9168977 | 53 | 48 | **53** | AC |
| HG01884 | DUP 11021371–11021437 | 66 | 35 | **66** | ABBC |
| HG01884 | DEL 14072417–14072480 | 63 | 23 | **63** | AC |
| HG01884 | DUP 30602119–30602182 | 63 | 49 | **63** | ABBC |

### Root cause
Methods p. 25 / `orig:sv_output.py`: before `sv_affected_len`, reference blocks are rebuilt with
internal boundaries at `floor(median(prev.end, next.start))` (`blocks_midpoints`). Affected length
sums those midpoint-extended lengths.

ARC-SV2 `caller._row_for_genotype` used `len(block) = end - start` with breakpoint intervals left
uncovered, so the deleted/duplicated middle block was short by ~half of each flanking gap
(e.g. DEL 613786: 47 + 3 + 4 from gaps of 7 → 54).

Likelihood / normalizer already use `genome_blocks_gaps` medians; this was an **output/measurement**
bug only. Merge-filter `--min-size 50` stays.

### Fix (implemented)
- `reference_block_lengths_midpoint` in `python/arcsv2/engine/caller.py` — used for `len_affected` only.
- Unit test `test_reference_block_lengths_include_median_gaps` (synthetic 613786 + the 6 geometries).
- Recompute script: `scratch/gap_analysis/recompute_len_affected_midpoint.py`
  - **35** short-filtered missed TPs; **34/35** now ≥50 (one fragmented DUP stays at 49).
  - Of all 181 raw calls with len < 50, 48 recompute to ≥50.

### Rerun note
Existing `arcsv_out.tab` still has old lengths. Either rewrite `len_affected` from the `blocks`
column (no re-call) then merge+label → `export/v2fix8/`, or full chr20 recall. Do not start
full chr20 yet until asked.

## Then: never called (not in raw output) — DONE (diagnostic)

Earlier estimate "~44" was wrong. Against `results_v2fix7/*/chr20/arcsv_out.tab` (raw):

| bucket | n | notes |
|:--|--:|:--|
| missed orig TPs (vs v2fix7 label=1) | **75** | 243 − 169 − (sample split) = 75 |
| called but short (`len_affected` < 50) | **35** | merge-filter drop; midpoint fix recovers 34 |
| called ≥50 but label=0 / not TP | **20** | present in raw; fail labeler (`low_mapq`, `clipped`, `large_gap`, `low_identity`, …) |
| **never in raw** | **20** | this section |

Never-called mix: 10 DEL + 10 DUP; HG01884 9 / NA19700 11. Of the 20: **5 in a graph_log `ok` locus**, **15 no locus**.

Script: `scratch/gap_analysis/diag_never_called_v2fix7.py` (Docker `arcsv2:dev`).
Outputs: `scratch/gap_analysis/never_called_v2fix7.{out,json}`.

### Ten dumps (5 in-locus + 5 no-locus)

| sample | event | type | orig len | sr/pe | locus | local forced gt | verdict |
|:--|--:|:--|--:|:--|:--|:--|:--|
| HG01884 | 709820–709945 | DUP | 125 | 0/0 | ok p=2 | DEL ABCDEH vs_ref=27.7 | **locus_ok_silent_full** (wrong geom; no ABBC) |
| NA19700 | 55751957–55752029 | DEL | 72 | 1/0 | ok p=2 | DEL ABE vs_ref=25.7 | **locus_ok_silent_full** |
| HG01884 | 45267979–45268083 | DEL | 104 | 1/0 | ok p=2 | DEL AE vs_ref=22.1 | **locus_ok_silent_full** |
| NA19700 | 62993597–62993683 | DEL | 86 | 0/0 | ok p=2 | DEL AG vs_ref=48.7 | **locus_ok_silent_full** |
| HG01884 | 63571651–63571976 | DEL | 325 | 0/5 | ok p=3 | DEL ABL vs_ref=79.7 | **locus_ok_silent_full** |
| NA19700 | 52764531–52764583 | DUP | 52 | 1/0 | none | — | **bp_present_weak_novel** (same-block support=1) |
| NA19700 | 18675742–18675855 | DUP | 113 | 1/5 | none | — | **bp_present_weak_novel** (softclips=51; nov support≤1) |
| HG01884 | 40558523–40558662 | DUP | 139 | 1/0 | none | — | **bp_present_weak_novel** (same-block support=1) |
| NA19700 | 19154873–19154926 | DEL | 53 | 0/0 | none | — | **bp_present_weak_novel** (skip support=1; softclips=27) |
| NA19700 | 22505500–22505590 | DEL | 90 | 0/0 | none | — | **bp_present_weak_novel** (skip support=0) |

### Root causes

**A. No-locus, breakpoints present, novel edges too weak (13/20 cheap-class `bp_novel_edge_but_no_locus`; 5/5 no-locus dumps)**  
Discover finds softclip (often split) BPs at the site. Local rebuild shows skip and/or same-block novel edges, but `edge.support` is usually 0–1 while `min_edge_support=2`, so `decompose()` never opens a locus. Classic example: NA19700 DUP 52764531 has the tandem same-block edge with split support=1 and skip support=1 — both ignored for decompose. Softclip evidence is abundant (9–51 clips at site).

**B. In-locus `ok` but full chr20 emits nothing (5/20)**  
`graph_log` has an `ok` locus (2–3 paths) covering the site, yet no raw call shares that locus id or lies within 5 kb. Full-run genotype therefore returned no emit (ref preferred / no informative row). Local force-genotype of the published span can invent a winning **larger DEL** (often wrong geometry for DUP sites — extra BPs fragment the event; ABBC never appears). `call_evidence` on the same window still emits 0 rows, so the forced-gt "would emit" is a window artifact, not a full-run miss of an already-built call.

**C. Missing breakpoint (2/20 cheap-class)**  
HG01884 DEL 62961107–62961290 (orig pe=8) and DEL 64090911–64091153 (orig pe=6): no covering BP in local discover.

### Cheap classify all 20
`locus_with_novel_edge` 5 | `bp_novel_edge_but_no_locus` 13 | `missing_breakpoint` 2

### Recommended next fixes (by leverage)
1. **DONE (v2fix8):** len_affected midpoint rewrite + merge/label — **+24** orig TPs (24/35 short TPs label=1); see Rerun.
2. **Novel-edge support for clip/split-driven small events** — largest lever on the 20: credit softclip/split so skip and same-block edges reach `min_edge_support` (or lower the decompose gate for clip-backed junctions). Targets ~13 no-locus sites (and may help silent in-locus loci).
3. **DUP block fragmentation** — merge adjacent clip intervals that split a short tandem-dup block so ABBC remains a path (overlaps the short-DUP / genotype-ref issues from v2fix6 dumps).
4. **In-locus silent genotype (5 sites)** — after (2)–(3), re-trace; likely ref-vs-wrong-DEL competition on large loci rather than a separate emit bug.
5. **Labeler/filter on the 20 called-≥50 label=0** — and on the **10** short TPs that now pass merge but label≠1; `low_mapq` / `clipped` / `large_gap` / `low_identity` / `ambiguous`.

Do **not** implement engine fixes yet unless trivial. Do **not** start a full chr20 re-call until asked.

## Rerun — DONE (v2fix8, no re-call)

Rewrote `len_affected` from the `blocks` column on v2fix7 raw tabs into
`/data/results_v2fix8/{sample}/chr20/` (INS/`_` rows left unchanged; tab rewrite
cannot recover insertion lengths). Then snakemake `label.v2fix8.yaml` `-j 1`
(merge + label + dataset only; call skipped). Wall **84 s**. Export:
`export/v2fix8/`.

| metric | v2fix7 | v2fix8 |
|:--|--:|--:|
| Validated | 182 | **207** (+25) |
| Precision | 182/282 (65%) | 207/321 (65%) |
| Original TPs recovered | 168/243 (docs 169) | **192/243** (+24) |
| Short filtered TPs (of 35) now label=1 | 0 | **24** |
| Merged calls (HG01884 / NA19700) | 248 / 218 | 270 / 244 (+48 raw≥50) |

Of the 35 short missed TPs: **34** recompute to ≥50; **24** become labeled TPs;
**1** stays at 49 (fragmented DUP); **10** pass merge-filter but label ≠ 1
(`ambiguous` / `clipped` / `low_identity` / …). Validated +25 vs TP recovery +24:
one extra tool-side validated call that is not an original TP match.

Scripts: `scratch/gap_analysis/rewrite_len_affected_tabs.py`,
`run_v2fix8_snakemake.sh`, `export_and_summarize_v2fix8.sh`, `summarize_vs_orig.py`.

## v2fix8 result (len_affected midpoint rewrite, merge+label only)
Validated 207 (was 182), precision 64.5% (unchanged), del 151/214 (71%), dup 39/67 (58%), complex 17/40 (42%).
Original TPs recovered 198/243 (±100 bp, same class); 25 validated calls the original does not make;
56 ARC-SV2-only false calls. Original: 243 validated, 54% precision.

## Edge-support parity — IMPLEMENTED (v2fix9 blocked on DUP same-block)

### Weak DELs (why original never calls) — DONE

Original `-v 2` on HG01884 ±20 kb around 23857951 / 13358681 / 63458085
(`scratch/gap_analysis/run_orig_weak_del_trace.sh`). **Literal skip PE is credited** and
DEL paths are scored; genotype still prefers reference:

| site | skip support | paths | best | next (DEL) | call |
|:--|--:|:--|:--|:--|:--|
| 23857951 | 27 | ABCD vs ACD | ABCD HOM | ACD HET (−7 nats) | none |
| 13358681 | 23 | ABCD vs ACD | ABCD HOM | ACD HET (−11 nats) | none |
| 63458085 | 25 | ABC vs AC | ABC HOM | AC HET (−22 nats) | none (INS elsewhere) |

Not missing BPs / min_bp_support / FDR / path skip. Filter = **ref preferred**
(`score_vs_ref ≤ 0`). ARC-SV2 already has that emit guard. See `docs/EDGE_SUPPORT.md`.

### Gate choice — A (literal `add_pair`, gate default False)

`CallConfig.skip_edge_ref_gate` default **`False`** (PARITY item 10). Evidence: with gate
off, window `call_evidence` emits **0** near the three weak DELs; pair-only DUP flood on
chr20:59.4–64.4 Mb stays at **2** (not the old 13). Option B (clip-aware narrowed gate)
not needed.

### Within-read PE `add_support` — DONE

`build_graph` credits within-read edges of multi-block PE mates via
`_add_within_read_pe_support` / `_within_read_path_edges` (orig `badj` walk before the
between-mate junction). Skips ends resolved from splits (those use `_add_split`). Unit
test: `test_multiblock_pe_mate_adds_within_read_novel_support`. Does **not** raise
same-block Dup support (still split-only).

### Validation before v2fix9 — PASS

Script: `scratch/gap_analysis/validate_v2fix9_windows.py` + `quick_dup_flood_check.py`.

| site | loci open | emit | expected | result |
|:--|:--|:--|:--|:--|
| NA19700 DUP 52764531 | yes | ABBC HET vs_ref=0.11 sr=1 pe=1 | ABBC | **PASS** |
| HG01884 DUP 40558523 | yes | ABBC HET vs_ref=12.72 sr=1 pe=6 | ABBC | **PASS** |
| NA19700 DEL 19154873 | yes | AC HOM vs_ref=6.49 pe=22 | AC | **PASS** |
| NA19700 DUP 18675742 | yes | BND ABCBCBCDE vs_ref=53.14 | any | **PASS** (fragmented vs orig ABBBC) |
| weak DEL ×3 | — | 0 | uncalled | **PASS** |
| 59.4–64.4 Mb pair-only DUP | — | 2 | ≤5 | **PASS** |

### Diagnosis (NA19700 DUP 52764531 ±20 kb)

1. **Dup splits on same-block:** original credits **1** Dup+ within-read (`[5,5,7]` → `(5,4)`) plus **1** multi-block PE between-mate on the flip_parity junction (`r1=[6,4] r2=[3,5]` → between `(4,5)`, `alt_dist=254 < q_low`) → same-block **support=2**. ARC-SV2 had only the Dup+ split (**support=1**); the PE landed on a skip store edge and never credited tandem.
2. **Blocks:** original point BPs → contiguous `B2=[52764531,52764583)`. ARC-SV2 half-open CIs left 1-bp gaps (`B2` started at 52764532).
3. **PE within-read:** original same-block support is **not** only splits — multi-block PE between-mate (after `flip_parity`) is the second credit. Soft clips are not edge support.

### Fix (v2fix9 engine)

- `breakpoint_block_cuts` + `create_blocks` point cuts → contiguous blocks (Option A).
- `_between_junction` / `_between_junction_supported`: credit flip_parity between-mate when it differs from `(v1,v2)`; same-block uses `alt_dist < q_low`; outward flood guard only when both mates are single-block (Option C from orig trace).
- Keep `skip_edge_ref_gate=False`, within-read PE `add_support`, outward check for single-block tandem. No lowering `min_edge_support`; no clip-as-support.

**v2fix9:** commit `ef23352`. Configs ready (`label.v2fix9.yaml`, `v2fix9_batch_inner.sh`,
`scratch/start_v2fix9_batch.ps1`). Docker `arcsv2-v2fix9-batch` **not started** from this
session (auto-review blocked the launch); run `scratch/start_v2fix9_batch.ps1` locally.

## v2fix9 result — regression, v2fix8 stays the best run
| chr20, both samples | v2fix8 | v2fix9 | original |
|:--|--:|--:|--:|
| Validated | 207 | 198 | 243 |
| Precision | 64.5% | 56.6% | 54.1% |
| Deletion | 151/214 (71%) | 128/196 (65%) | 172/255 |
| Duplication | 39/67 (58%) | 46/67 (69%) | 66/164 |
| Complex | 17/40 (42%) | 24/87 (28%) | 5/30 |
| Insertion calls | 181 | **29** | 256 |
| Original INS matched (±200 bp) | 128/256 | **19/256** | — |
| Orig TPs recovered | 198 | 174 | — |
| Tool-only FPs | 56 | 81 | — |

What broke (from the labeled tables):
- 45 of 87 complex calls sit on original **duplications**; typical rearrangements `ABCBCD`, `ABCDCDE`,
  `ABCBCBCD` — a tandem dup fragmented into 2–3 blocks and reported as BND/complex.
- 24 v2fix8 true DELs have **no call** in v2fix9 (e.g. NA19700 DEL 61783826–61785044 `AE` → `ACDE`).
- Insertions collapsed 181 → 29 (insertions placed 196 → 31 per sample).
These match the **point-cut / contiguous-block change** (Option A), not the gate or the PE support change.
Methods p. 20–21 say gaps between blocks exist when `b_j − a_j > 1`; the original keeps breakpoint intervals.

## v2fix10 gate — PASS on validated (full chr20 started)

Point cuts disabled by default (`CallConfig.use_point_block_cuts=False`). `call_evidence`
passes half-open CI intervals to `create_blocks` again. `breakpoint_block_cuts` kept for
ablation/tests only. Kept: `skip_edge_ref_gate=False`, within-read PE `add_support`,
flip_parity between-junction, outward check for single-block tandem.

### Step 1 — Ablation (script `scratch/gap_analysis/ablate_point_cuts_windows.py`)

Same discovery; A = gapped (current minus point cuts); B = point cuts ON; C = existing tabs.

| window | setting | DEL | DUP | complex | INS | total |
|:--|:--|--:|--:|--:|--:|--:|
| NA19700 61.5–64.5 Mb | A_gapped | 22 | 11 | 16 | 20 | 69 |
| (lost DEL 61783826) | B_point_cuts | 22 | 9 | 10 | 7 | 48 |
| | C_v2fix8_tab | 19 | 11 | 9 | 24 | 63 |
| | C_v2fix9_tab | 20 | 9 | 14 | 7 | 50 |
| HG01884 60–65 Mb | A_gapped | 28 | 16 | 19 | 22 | 85 |
| (ABCBC complex) | B_point_cuts | 27 | 15 | 21 | 3 | 66 |
| | C_v2fix8_tab | 22 | 16 | 10 | 27 | 75 |
| | C_v2fix9_tab | 27 | 15 | 23 | 3 | 68 |

Gapped restores INS (20/22 vs point-cut 7/3) and the lost DEL geometry
(`61783826 AE len=810` vs point-cut `ACDE len=162`). Complex stays **above** v2fix8
(16 vs 9; 19 vs 10) — between-junction / within-read PE still fragment some tandem
DUPs into ABCBC… even without point cuts.

### Step 2 — 52764531 without point cuts — PASS

With gapped blocks: same-block support=2 (1 Dup+ split + 1 between-junction PE),
locus opens, emits `ABBC HET`. Gap handling: `overlapping_vertices` / `_in_one_bp_gap`
for coordinates in a 1-bp CI gap; `_segment_agrees` documents the leeway for gap
abutment. Unit test: `test_gapped_ci_blocks_same_block_support_from_split_and_pe`.

Validation (`scratch/gap_analysis/validate_gapped_v2fix10.py`):
| site | result |
|:--|:--|
| NA19700 DUP 52764531 | PASS ABBC sameblock_sup=2 |
| HG01884 DUP 40558523 | PASS ABBC sameblock_sup=7 |
| NA19700 DEL 19154873 | PASS AC |
| NA19700 DUP 18675742 | PASS (fragmented ABCBCBCDE) |
| weak DEL ×3 | PASS uncalled |
| 59.4–64.4 Mb pair-only DUP | PASS (=1 ≤5) |

### Step 3 — ml-label windows vs v2fix8 (decision metric = validated + precision)

Call-count ablation above is **not** the gate. Gate script:
`scratch/gap_analysis/window_ml_label_v2fix10.sh` →
`scratch/gap_analysis/v2fix10_windows/`. Commit with defaults: `d4c2265`.

| window | arm | validated | n01 | prec | DEL | DUP | complex | orig TP |
|:--|:--|--:|--:|--:|:--|:--|:--|:--|
| NA19700 61.5–64.5 Mb | A_current (gapped) | **15** | 40 | 37.5% | 9/20 | 4/8 | 2/12 | 12/18 |
| | B_v2fix8 | 14 | 30 | 46.7% | 8/16 | 4/7 | 2/7 | 11/18 |
| HG01884 60–65 Mb | A_current | **23** | 50 | 46.0% | 16/28 | 5/6 | 2/16 | 21/30 |
| | B_v2fix8 | 20 | 36 | 55.6% | 13/21 | 5/7 | 2/8 | 19/30 |
| **aggregate** | A_current | **38** | 90 | 42.2% | — | — | — | 33/48 |
| | B_v2fix8 | 34 | 66 | 51.5% | — | — | — | 30/48 |

**DECISION (historical):** passed the then-stated gate on validated only (+4). Precision was
**worse** (−9.3 pp). **Retrospective: treat as EXPERIMENT** — full chr20 confirmed
precision loss; **v2fix8 remains best**. Complex *call counts* remain inflated vs v2fix8
(ABCBC trace below); next lever is `simplify_blocks_diploid` (v2fix11), not point cuts.

## v2fix10 result — EXPERIMENTAL (window precision drop; do not prefer over v2fix8)

Full chr20 batch finished (`arcsv2-v2fix10-batch` exit 0, wall ~2 h 12 min).

| chr20, both samples | v2fix8 | v2fix10 (experiment) | original |
|:--|--:|--:|--:|
| Validated | **207** | 222 | 243 |
| Precision (0/1) | **64.5%** (207/321) | 59.0% (222/376) | 54.1% |
| Deletion | 151/214 (71%) | 156/218 (72%) | — |
| Duplication | 39/67 (58%) | 43/71 (61%) | — |
| Complex | 17/40 (42%) | 23/87 (26%) | — |

Window ml-label (decision metric) had already shown **+4 validated** but **−9.3 pp precision**
(38/90 = 42.2% vs v2fix8 34/66 = 51.5%). Full chr20 confirms validated up, precision down
(~5.5 pp). Treat v2fix10 as an **experiment only** — **v2fix8 stays the best production run**.
Do not start further full-chr20 recalls from this tip unless a later window gate recovers
precision.

### Step 4 — full chr20 v2fix10 **DONE (experiment)**

`label.v2fix10.yaml` + `v2fix10_batch_inner.sh` → docker `arcsv2-v2fix10-batch`
wrote `/data/results_v2fix10` / `/data/labels_v2fix10`. Finished exit 0.

## ABCBC complex vs original (±20 kb) — HG01884 60–65 Mb

Script: `scratch/gap_analysis/trace_abcbc_orig.sh` → `loci/orig_abcbc_*` + `v2_abcbc_*`.

| center | original near call | ARC-SV2 (pre-v2fix11) | takeaway |
|:--|:--|:--|:--|
| 60792115 | **DUP ABBC** len=68; path next ABCBCD | **BND ABCBCD** len=68; extra cut @60792144 (clip×2) | Same BLOCKS; orig `simplify_blocks_diploid` → ABBC |
| 64177566 | **DUP ABBC** len=118 (nearby ABBBC + complex) | **BND ABCDCDE** len=118 | Same; simplify fuses unused internal cuts |
| 63028037 | complex `ABCDC_DCDEFH` / `AB_CDBCDBCDFH` | `ABCDBCDCDE_F` + DEL AC | Both complex; not the ABBC case |

**Summary:** original **keeps the same extra softclip cuts in BLOCKS** (e.g. 60792144
clip-only). Genotype path is still `ABCBCD` / `ABCDCDE`. Before writing the tab,
`simplify_blocks_diploid` (`orig:sv_validate.py`) fuses adjacent blocks whose mutual
adjacency is only used as reference (degree-1) in **both** haplotypes → output
`ABC`/`ABBC` DUP. ARC-SV2 lacked that step → BND. Not a point-cuts issue.

## v2fix11 — port `simplify_blocks_diploid` — window gate **PASS** (flank-preserving)

### Rule citation
- **File:** `orig:arcsv/sv_validate.py::simplify_blocks_diploid` (also `simplify_blocks`)
- **Call site:** `orig:arcsv/sv_inference.py` after genotype, before `sv_output` path strings
- **Logic:** build neighbor sets from both haplotype paths; if `B_out`’s only neighbor is
  `C_in` and vice versa, drop the cut and fuse `[B.start, C.end)` into one block; remap paths.
  Turns `ABCD`/`ABCBCD` → `ABC`/`ABBC`.
- **Classify note:** original `classify_paths` still runs on the *unsimplified* path (multi-block
  `DUP:TANDEM` in `sv_classify.py`); tab `reference`/`rearrangement`/`bp` use the simplified
  blocks. ARC-SV2 applies simplify then `classify_simple` (single-block DUP) — equivalent for
  these sites.

### Port
- `python/arcsv2/engine/classify.py::simplify_blocks_diploid`
- Wired in `caller._row_for_genotype` (strings, `bp`, `len_affected`, simple class); support for
  fused tandem maps via `orig_groups`; `rearrangement_next` stays unsimplified.
- **Sequence guard:** full `blocks_point` allele bases must match pre/post fuse (RuntimeError).
- **Flank guard:** if 150 bp outer-block trim (labeler) would change allele coords, **revert**
  simplify for that row. Keeps ABCBCD→ABBC; blocks DEL ACD→AC that broke labeling.
- Tests: `test_simplify_blocks_diploid_collapses_unused_cut_to_abbc`,
  `test_simplify_blocks_diploid_keeps_novel_cut`,
  `test_simplify_blocks_diploid_preserves_allele_base_coords`,
  `test_simplify_del_acd_changes_flank_trimmed_allele`. Point cuts remain **OFF**.

### Lost-TP diagnosis (seqfix 35 vs v2fix10 38)
Report: `scratch/gap_analysis/lost_tp_seqfix_report.md`. Net −3 = 5 lost − 2 gained.
HG01884 losses (the aggregate drop): DEL 62270773, BND 63442864, DEL 63553544 — full
allele identical but flank-trimmed queries differed after end-adjacent fuse.

### Window ml-label — **PASS**
Script: `scratch/gap_analysis/window_ml_label_v2fix11.sh` → `v2fix11_seqfix_windows/`.

| aggregate | validated | n01 | precision | orig TP |
|:--|--:|--:|--:|--:|
| A_current (v2fix11 flank-safe) | **38** | 90 | **42.2%** | 33/48 |
| B_v2fix8 | 34 | 66 | 51.5% | 28/48 |

Gate `validated ≥ 38` and `precision ≥ 42%`: **PASS**. 60792115 still **DUP ABBC**.
Full chr20: `label.v2fix11.yaml` + `scratch/start_v2fix11_batch.ps1`.
