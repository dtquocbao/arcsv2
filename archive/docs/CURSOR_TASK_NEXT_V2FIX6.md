# Cursor task: after v2fix6 — speed, labeler guard, filters, recall

## v2fix6 result (chr20, HG01884 + NA19700, same labeler)
| | before | v2fix5 | **v2fix6** | original |
|:--|--:|--:|--:|--:|
| Validated calls | 106 | 156 | **182** | 243 |
| Overall precision (labeled) | 24% | 33% | **65%** | 54% |
| Deletion precision | 24% | 52% | **70%** (138/197) | 67% |
| Duplication precision | 14% | 10% | **56%** (31/55) | 40% |
| Complex precision | 33% | 8% | **43%** (13/30) | 17% |
| Original TPs recovered | 85 | 164 | **169**/243 | — |
| TPs the original misses | — | — | **20** | — |
| ARC-SV2-only false calls | — | 233 | **48** | — |
| HG01884 call time | 44 min | 49 min | **7.4 h** | 24 min |
Precision is now at or above the original in every class. Remaining: recall (182 vs 243) and speed.

## 1. Speed regression (blocker) — DONE
`build_graph` was O(pairs × blocks) because `overlapping_vertices` scanned every
block per read. It now uses `find_block` + a short walk. 5 Mb window
(`scratch/gap_analysis/bench_build_graph_5mb.py`): **build_graph 7.1 s**
(pairs 629k, blocks 2849). Target ≤ 3 min on full chr20 should hold.

## 2. Labeler guard — DONE
Length is estimated from block spans **before** `haplotype_sequences` /
mappy (`MAX_LABEL_QUERY_LEN = 1_000_000`, reason `query_too_long`). Each
`aligner.map` is capped at 60 s (`LABEL_MAP_TIMEOUT_S`, reason `timeout`).
Unit tests: `test_megabase_ref_query_is_skipped`,
`test_two_megabase_call_is_skipped_before_sequence_build`.

## 3. Paper merge-filter in label.smk — DONE (wiring); numbers still needed
Methods p. 27: drop calls with >85% simple-repeat/low-complexity/satellite/STR overlap, >50% segdup overlap,
<50 bp affected, blacklist overlap, or reference-gap overlap.

**Done**
- `scripts/fetch_filter_beds.py` downloads UCSC simpleRepeat/rmsk(Low_complexity|Satellite|Simple_repeat)/microsat,
  genomicSuperDups, and 10x `svblacklist.bed`; writes `/data/ref/{simple_repeats,segdups,blacklist}.bed`
  (gaps remain from `make_ref_beds.py` N-runs). Also `fetch_data.py --filter-beds`.
- Config keys `simple_repeats` / `segdups` / `blacklist` in `label.yaml` (paths filled),
  `label.v2fix*.yaml` (empty = keep published min-size-only numbers), and `write_label_auto_yaml`.
- `label.smk` `merge`: when any annotation key is non-empty, passes those flags plus `--gaps`.
  Empty annotation keys => min-size only (prior v2fix* behavior).
- Same CLI flags for the original in `docs/ORIGINAL_COMPARISON.md` §5.
- CLI still ignores missing BED paths (`is_file()`).

**Still needed for with-vs-without numbers (no full chr20 re-call)**
1. Inside container: `python scripts/fetch_filter_beds.py --dest /data/ref` (rmsk is large; once).
2. Without filter: current v2fix6 labeled tables (min-size only) — already recorded above.
3. With filter — re-merge + relabel only (reuse `arcsv_out.tab` for both tools):
   ```bash
   BEDS="--simple-repeats /data/ref/simple_repeats.bed --segdups /data/ref/segdups.bed --blacklist /data/ref/blacklist.bed --gaps /data/ref/gaps.bed"
   for s in HG01884 NA19700; do
     arcsv2 merge-filter /data/results_v2fix6/$s/chr20/arcsv_out.tab $BEDS -o /data/results_v2fix6/$s/merged.paper.tab
     # same $BEDS on orig v2cols.tab -> /data/orig/$s/merged.paper.tab
     # ml-label + ml-dataset into labels_v2fix6_paper / orig/labels_paper
   done
   ```
   Or set the three keys in `label.v2fix6.yaml` to the `/data/ref/*.bed` paths and re-run
   `merge`+`label`+`dataset` only (`--forcerun merge`).
4. Record both sides in `bench/RESULTS.md` (call counts + labeled precision per class).

## 4. Recall: 74 original true calls ARC-SV2 misses — DONE (diagnostic)
Inventory (export labeled tables + `results_v2fix6/*/chr20/graph_log.txt`):
orig TP 243; missed vs v2 TP **75** (≈74); no ARC-SV2 call within 2 kb **55**; of those,
**40 in a locus** (39 `ok` with 2–3 paths) and **15 with no locus**.

Script: `scratch/gap_analysis/diag_recall_misses.py` (Docker `arcsv2:dev`).
Outputs: `scratch/gap_analysis/recall_misses_v2fix6.{out,json}`.

### 4a. Ten in-locus dumps (ok, 2–3 paths)
| sample | event | type | orig len | called | vs_ref | len_aff | verdict |
|:--|--:|:--|--:|:--|--:|--:|:--|
| HG01884 | 613786–613840 | DEL | 54 | ACD/ACD | 6.56 | **47** | **min_affected_length_drop** |
| HG01884 | 5644313–5644364 | DUP | 51 | ABCD/ABCD | 0 | 0 | genotype_prefers_ref |
| HG01884 | 9168924–9168977 | DEL | 53 | AC/AC | 5.32 | **48** | **min_affected_length_drop** |
| HG01884 | 50452086–50452147 | DEL | 61 | AD/AD | 39.34 | 259 | would_emit (window; not on full chr20) |
| NA19700 | 9168924–9168977 | DEL | 53 | AC/AC | 5.44 | **48** | **min_affected_length_drop** |
| NA19700 | 9716117–9716171 | DEL | 54 | AC/AC | 5.30 | **45** | **min_affected_length_drop** |
| NA19700 | 53500266–53500342 | DEL | 76 | AC/AC | 8.27 | **39** | **min_affected_length_drop** |
| HG01884 | 30602119–30602182 | DUP | 63 | ABC/ABC | 0 | 0 | genotype_prefers_ref |
| NA19700 | 46420064–46420127 | DEL | 63 | AC/AC | 6.94 | **41** | **min_affected_length_drop** |
| HG01884 | 11021371–11021437 | DUP | 66 | ABC/ABC | 0 | 0 | genotype_prefers_ref |

**Root-cause fractions (10 dumps):**
- **6/10 min_affected_length drop** — alt wins (vs_ref 5–8) but middle-block `len_affected` is 39–48 bp after BP CI shrinks the event; emit filter at 50 kills the call (PARITY D11; original `do_inference` does not). Skip edges + splits are present (e.g. DEL 613786 skip support=22, sr≈21).
- **3/10 genotype prefers ref** — all DUPs. Extra BPs split the short duplicated interval (e.g. 5644313 → blocks 43 bp + 6 bp); paths are multi-copy of the fragments (ABCBCD…), not ABBC; ref is best by ~3.6–4 nats.
- **1/10 window/emit mismatch** — local genotype would emit a larger DEL; full-chr20 table has no nearby call (secondary).

### 4b. Fifteen no-locus sites
- **14/15**: breakpoints **present** at the site (usually softclip-only ends; some with DRP). Discover stage is fine; locus decompose / novel-edge survival fails (no `ok` locus in graph_log).
- **1/15**: truly missing breakpoint (HG01884 DEL 64090911–64091153).

### Recommended next engine fix
1. **DONE for v2fix7:** emit no longer applies `min_affected_length=50` in
   `_row_for_genotype` (match original D11). Paper ≥50 bp stays in merge-filter.
2. **DUP recall (secondary):** merge adjacent clip intervals that fragment a short tandem-dup block so ABBC remains a path; otherwise genotype correctly prefers ref.
3. **No-locus (later):** softclip BPs without surviving novel edges — separate diagnose of split wiring / min_edge_support for clip-only small events.

## 5. Then — in progress (v2fix7)
Commit speed + labeler + emit-filter + merge-filter wiring; rerun chr20
(`label.v2fix7.yaml`, `-j 1`); relabel; export `export/v2fix7/`; record
`bench/RESULTS.md`. Paper with-filter numbers still need BED fetch + re-merge.
