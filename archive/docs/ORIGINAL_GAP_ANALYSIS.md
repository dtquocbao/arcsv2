# Why ARC-SV2 underperforms the original ARC-SV (gap analysis)

**Scope:** chr20 HG01884 + NA19700 (same BAMs, same merge-filter, same haplotype labeler).  
**Sources:** original `../../src/arcsv/arcsv/` (READ ONLY), ARC-SV2 `python/arcsv2/`, exports under `export/`, thesis `Chapter 1.pdf`, papers `nihms-2024296.pdf` / `nihms-605372.pdf`, roadmap `../ARC-SV2 - From-Scratch Implementation Roadmap.md`.  
**Rule:** facts vs hypotheses are marked. Engine code was not modified.

---

## Summary (top causes, ranked by estimated impact)

1. **Breakpoint merge lacks original precedence (proven causal for missing true DELs).**  
   Original (`breakpoint_merge.py:63–71`, `__add__:32–45`, `bp_mergefun_precedence:261–294`): softclip/split BPs have precedence **2**, PE-only BPs precedence **1**. When a PE interval touches a clip/split interval, support is absorbed but **coordinates stay on the clip/split interval**.  
   ARC-SV2 (`breakpoints.merge_breakpoints:593–615`): **no precedence** — overlapping/touching intervals are always unioned.  
   **Single-locus proof (HG01884 ±20 kb, see below):** at DEL 5962348–5962677, ARC-SV2 has separate left/right split endpoints before merge, then **one** union interval `[5954931,5972086)` after `merge_breakpoints` → no deleted middle block → **0 calls**. Scratch precedence merge recovers `DEL … sr=41 pe=50 rearr=AC`. Same pattern for all five original-only true DELs. Removing `merge_intervals_for_min_block` or `absorb_nearby_drp` alone does **not** recover them.  
   **Cell 2024 Methods p20–21:** “merge all overlapping” breakpoint intervals; keep if DRP or ≥2 split/clip — **no precedence rule**. Not in the chapter provided (`Chapter 1.pdf` is background only). Roadmap has no precedence.

2. **`merge_intervals_for_min_block` (ARC-SV2-only; secondary).**  
   Joins breakpoints when the intervening reference would be `< min_affected_length` (50 bp) (`breakpoints.py:618`; `caller.py`). Original has **no** such step. On the five true DELs it did not change the outcome once union-merge had already collapsed both ends; it can still over-merge dense PE/clip fields elsewhere.

3. **Pair-only false calls are not fixed by precedence alone.**  
   Flip experiments on ARC-SV2-only false loci (28549495, 61390333, 64235350): high-PE `sr=0` calls **remain** under precedence. So cause 1 explains missing trues; false PE calls need additional FDR/scoring or post-filters.

4. **`max_paths=1000` (paper) vs original code 2000 — downstream of over-merge.**  
   Documented in Cell p25. Four `skipped_npaths` loci on full chr20. Dense union-merged graphs make path explosion more likely (**hypothesis** linking to cause 1).

5. **Runtime/memory:** in-memory full-chr pairs (~5.1 GB / ~44 min) vs original streaming (~1.8 GB / ~24.5 min). See Performance.

---

## Stage map

| Stage | Original (file:line) | ARC-SV2 (file:line) |
|:--|:--|:--|
| 1 BAM / read filters | `call_sv.call_sv:103`; `bamparser_streaming.parse_bam:98`; `softclip.process_softclip:39`; `splitreads.parse_splits:68` | `pipeline.call_sam`; `samio.load_region:235`; `evidence.softclips_from_alignment:95`; `evidence.parse_split:254` |
| 2 Library / insert model | `bamparser` KDE `pmf_kernel_smooth:459`; LR cutoffs `:134`; approx stats `:25` | `library.fit_library:262`; `caller._fit_libs:120`; `discordance_cutoffs` log C=4.5 |
| 3 DRP cluster + FDR | `pecluster.apply_discordant_clustering:82`; `cluster_pairs:175`; `fdr_discordant_clusters:412` | `caller.discover_breakpoints:162`; `breakpoints.cluster_drps:327`; `fdr_keep:420` |
| 4 Softclip / split | softclip merge in `breakpoint_merge:107–118`; splits collected in bamparser | `evidence.clip_breakpoints:142`; `cluster_softclips:130`; `split_breakpoints:359` |
| 5 Breakpoint merge | `breakpoint_merge.merge_breakpoints:75`; precedence `:63–71`, `:261–294` | `breakpoints.merge_breakpoints:593`; `merge_intervals_for_min_block:618`; `absorb_nearby_drp_breakpoints:650` (no precedence) |
| 6 Blocks / graph | `sv_parse_reads.parse_reads_with_blocks:246`; `GenomeGraph.add_pair:26` | `graph.create_blocks:60`; `caller.build_graph:343` |
| 7 Paths / genotype | `sv_inference.do_inference:26`; `get_paths_iterative:861`; `max_paths=2000` | `graph.decompose`; `enumerate_paths`; `likelihood.genotype_locus`; `max_paths=1000` |
| 8 Insertions | `sv_inference:70–181`; `sv_inference_insertions.py` | `caller.place_novel_insertions:969` |
| 9 Output / filter | `sv_output.py`; `sv_filter.py`; `filter_output.py` | `caller.write_calls`; `merge_filter.py` |

---

## Difference table

| # | Stage | Original behavior | ARC-SV2 behavior | Impact (evidence) | Thesis/paper status | Question for Dr. Bo? |
|--:|:--|:--|:--|:--|:--|:--|
| D1 | 5 | Clip/split precedence 2 over PE 1; PE absorbs without moving interval (`breakpoint_merge.py:32–45,63–71,261–294`) | Interval union only (`breakpoints.py:593–615`) | Proven: 5/5 missing true DELs recovered by scratch precedence; lost at `merge_breakpoints` | **Partially documented** — Cell p4 “merged with” softclip/split; p20 “merge all overlapping”; no precedence. Not in the chapter provided. | Yes |
| D2 | 5 | No 100 bp DRP absorb pass | `absorb_nearby_drp_breakpoints(..., 100)` (`breakpoints.py:650`) | Not sufficient alone to explain missing trues (flip: still 0 calls) | **Not documented** in Cell 2024 Methods or the chapter provided. Searched: absorb, nearby, 100 bp. | Yes |
| D12 | 5 | No min-block breakpoint join | `merge_intervals_for_min_block` when gap `< min_affected_length` (50) (`breakpoints.py:618`) | Not sufficient alone for the 5 true DELs; can over-merge elsewhere | **Not documented** in Cell 2024 Methods or the chapter provided. Original has no equivalent. | Yes |
| D3 | 3 | `max_pecluster_size = phys_cov × 1.5` (`bamparser:298–299`; HG01884 → ≈84) | Fixed `80` (`config.py:42`) | 32 oversized components discarded on HG01884 (`graph_log`) | **Partially documented** — coverage-adaptive FDR (Cell p18–19); formula 1.5×phys not found. | Yes |
| D4 | 7 | `max_paths=2000` (`arcsv_call_options.py:69`) | `max_paths=1000` paper mode (`config.py:55`) | 4 skipped_npaths loci on HG01884 | **Documented** Cell p25 (1,000). Roadmap notes paper 1k vs code 2k. | No (ask only if paper mode should match code) |
| D5 | 5 | Separate `min_bp_support=2` vs `min_edge_support=2` | Reuses `min_edge_support` for BP keep + edges | Same numeric default; naming conflation | **Not documented** as two distinct parameters in Cell 2024 Methods or the chapter provided. | Yes |
| D6 | 3 | Del LR MLE: `min(max_span, mean(z)−μ)` (`pecluster.lr_del:303–309`) | `D=min(min(z)−m, min b−max a)` (`score_cluster:410–413`) | May change which Del clusters pass FDR (hypothesis) | **Partially documented** — Cell p18–19 heuristic score + FDR; exact MLE not found. | Yes |
| D7 | 1–2 | Approx chunk sampling then full-region re-KDE | Pools all region pairs (`pool_libraries=True`) | Runtime/RSS; slight cutoff shifts (hypothesis) | **Partially documented** — library-specific insert distribution (Cell p4); chunk approx not found. | Yes |
| D8 | 6–7 | Second-pass streaming graph build | One in-memory pass over all pairs × all blocks | RSS 1.8→5.1 GB; BAM load 8.5 min | Implementation detail; not a method claim | No |
| D9 | 4–5 | Softclip merge ≤5 bp, min support 2; MH ≤25 | Same defaults (`config.py:24–37`) | Aligned | Softclip/split described (Cell p4, p20) | No |
| D10 | 3 | FDR target 0.1, null reps 4 | Same | Aligned | **Documented** Cell p18–19 | No |
| D11 | 9 | Live call path does not apply `min_simplesv_size` | `min_affected_length=50` at emit + merge-filter | Mild size filter difference | Evaluation contexts use ≥50 bp; live code path differs | Low priority |

---

## Per-call traces (HG01884; from `scripts/diagnostics/`)

Candidates: `scratch/gap_analysis/HG01884_call_trace_candidates.tsv`.  
Full text: `scratch/gap_analysis/HG01884_call_traces.txt`.

### Original-only true (ABSENT in ARC-SV2)

| Locus | Orig call | ARC-SV2 | Notes |
|:--|:--|:--|:--|
| DEL 5962348–5962677 | HOM, sr=41, pe=29, score=478, rearr=AC | no nearby call | High SR+PE; first divergence cannot be “genotype only” |
| DEL 2883362–2883684 | HOM, sr=33, pe=9, score=363, rearr=AC | absent | |
| DEL 61995575–61995895 | HOM, sr=28, pe=12, score=320, rearr=AC | absent | |
| DEL 32311430–32311735 | HOM, sr=26, pe=10, score=285, rearr=AC | absent | |
| DEL 57219054–57219328 | HOM, sr=26, pe=11, score=278, rearr=AC | absent | |

**Fact:** absence means no emitted call within ±200 bp.  
**Hypothesis:** breakpoints never formed as a genotype-able deleted block (merge precedence / over-merge / path skip), not a late filter.

### ARC-SV2-only false (no original nearby)

| Locus | ARC-SV2 call | Orig | Notes |
|:--|:--|:--|:--|
| DEL 28549495–28549577 | HET, sr=0, pe=94, rearr=ABCDEFH | none | Pair-only; complex multi-block |
| DUP 61390333–61390777 | HET, sr=0, pe=72, rearr=ABCCD | none | |
| DUP 64235350–64235471 | HET, sr=0, pe=68, rearr=ABCDDEFGHIJKL | none | Many blocks |
| DEL 38728711–38728795 | HET, sr=0, pe=56, rearr=ABD, score=0.22 | none | Weak margin |
| DUP 18022682–18022889 | HET, sr=0, pe=50, rearr=ABCDEEFGH | none | |

**Fact:** these are pair-only, high-PE, often multi-block.  
**Hypothesis:** PE clusters that original FDR rejected or never merged into a callable simple locus.

### Window chosen for stage counts

`chr20:63000000–66000000` (3 Mb): 11 orig-only true + 19 v2-only false (HG01884).  
From full-chr exports restricted to the window (`scratch/gap_analysis/window_63_66_counts.txt`):

| Tool | Calls | sr=0 | By type (abbrev.) |
|:--|--:|--:|:--|
| Original | 58 | 36 | DEL 23, DUP 14, INS 8, BND-complex 11, mixed 2 |
| ARC-SV2 | 42 | 38 | DEL 19, INS 16, DUP 4, BND 3 |

Full per-stage clip/DRP counts for the window were not re-run end-to-end; per-locus stage dumps under `scratch/gap_analysis/loci/` cover the five true + three false foci (±20 kb). Scripts: `scripts/diagnostics/pick_gap_window.py`, `window_stage_counts.py`, `trace_breakpoint_stages.py`.

**First stage that diverges materially (proven):** stage **5 (breakpoint merge)** — same softclip/split/DRP ingredients, different merge algebra → different blocks/loci/calls. Stages 1–4 defaults largely match; stage 7 path skip is a downstream consequence.

---

## Performance

| Phase | Original (HG01884 chr20) | ARC-SV2 (HG01884 chr20) | Design cause |
|:--|:--|:--|:--|
| Wall | 24:29 (`*.time.txt`) | ~44.4 min (`call.log` done ~2660 s on probe; Snakemake sample later ~32–44 min) | See rows below |
| Peak RSS | 1.79 GB | ~5.1 GB | Hold all pairs + full block list in RAM |
| BAM / evidence | ~1.6 min (task brief) | ~8.5 min (`call.log` load_bam) | Python pair materialization; no approx-stats shortcut |
| Discover / FDR | included in parse+cluster | ~1 min after load | Similar algorithm; more DRPs listed |
| Graph build | second pass ~6.5 min (brief) | ~1 min after binary-search `find_block` fix | Was pathological O(pairs×blocks) before fix |
| Insertions | part of inference | ~18 min (3632 sites) | Per-site likelihood grid; many candidates after dense BPs |
| Genotype | SV inference ~15 min (brief) | ~14 min (2070 loci) | More loci from denser breakpoints; some `skipped_npaths` |

---

## Documentation search log (for “not documented” claims)

Texts extracted with `pypdf` to `scratch/pdf_text/` (`Chapter 1.txt`, `nihms-2024296.txt`, `nihms-605372.txt`, plus `*_paged.txt` via `pages[i].extract_text()`).  
**Note:** `Chapter 1.pdf` is thesis **background** only (no ARC-SV method text). Method claims are checked against Cell 2024 Methods (nihms-2024296) and the chapter provided.  
Terms searched (case-insensitive): precedence, absorb, nearby, soft-clip/softclip, split read, discordant, DRP, FDR / false discovery, merge, breakpoint, insert-size, clustering, 1,000 / 1000 / 2,000, path, subgraph, min_bp_support, min_edge_support, physical coverage, null, quantile, microhomology, min_block / min_affected.  
**Real Cell Methods pages used** (author manuscript): softclip/split thresholds **p17**; Del score / FDR / clique clustering **p18–19**; breakpoint merge **p20–21**; path limit 1,000 **p25**; ML labeler **p25–26**; merge-filter **p27**.

---

## Single-locus traces

Scratch diagnostic only (`scripts/diagnostics/trace_breakpoint_stages.py`); engine tree unchanged. Modes: `default` (union), `precedence`, `no_minblock`, `no_absorb`, `precedence_no_minblock_no_absorb`. Original ±20 kb with `-v 2`. Raw dumps: `scratch/gap_analysis/loci/`.

### Worked example — DEL 5962348–5962677 (HG01884)

**Pre-merge (both tools):** left/right split endpoints present separately, e.g. `[5962341,5962355)` and `[5962670,5962684)` plus Del PE clusters overlapping each end.

| Step | ARC-SV2 default (union) | Original (`-v 2`) / ARC-SV2 precedence |
|:--|:--|:--|
| After merge | **One** interval `[5954931,5972086)` clip=11 split=82 drp=22 | Two ends kept: `(5962341,5962355)` PE=29 absorbed; `(5962670,5962684)` PE=29; precedence: `[5962341,5962356)` + `[5962670,5962685)` |
| After min-block / absorb | Still one fat interval (n_bp 6–7) | Ends remain separate |
| Blocks / call | No deleted middle block → **0 calls** | Middle block ~314 bp → `DEL 5962341–5962685 sr=41 pe=50 rearr=AC` |

**Exact loss step (fact):** `merge_breakpoints` union — PE CI spans both clip/split ends and swallows them into one interval before genotyping.

### Flip experiments — five original-only true DELs

| Locus | default | no_minblock | no_absorb | precedence | precedence + no_minblock + no_absorb |
|:--|:--|:--|:--|:--|:--|
| 2883362 | 0 | 0 | 0 | **DEL sr=33 pe=36 AC** | **same** |
| 5962348 | 0 | 0 | 0 | **DEL sr=41 pe=50 AC** | **same** |
| 32311430 | 0 | 0 | 0 | **DEL sr=26 pe=29 ABD** | **same** |
| 57219054 | 0 | 0 | 0 | **DEL sr=26 pe=30 ABD** | **DEL sr=26 pe=29 ABCDEG** |
| 61995575 | nearby PE-only only; **true DEL absent** | same | same | **DEL 61995568–61995905 sr=28 pe=30 AC** | **same** |

**Conclusion:** precedence alone is **necessary and sufficient** to recover all five missing true DELs on these windows. Removing min-block or absorb alone never recovers them.

### Flip experiments — three ARC-SV2-only false calls

| Locus | default | precedence | precedence + no_minblock + no_absorb |
|:--|:--|:--|:--|
| DEL 28549495 (sr=0 pe=94) | present | **still present** (+ extra DUP) | **still present** |
| DUP 61390333 (sr=0 pe~72) | present | **still present** | **still present** |
| DUP 64235350 (sr=0 pe=68) | present | becomes BND-complex (sr=0 pe high) | still PE-driven complex |

**Conclusion:** precedence does **not** remove these pair-only false calls. Separate FDR/scoring/post-filter work remains.

### InvL note at DEL 5962348 (report only)

ARC-SV2 union dump shows the fat interval `[5954931,5972086)` with `types=('Del','InvL') clusters=(0,1,2)`.  
Original `-v 2` on the same ±20 kb window (`scratch/gap_analysis/loci/orig_5962348/verbose.log`) clusters **Del / Ins / Dup only**, then goes straight to merge — **no InvL/InvR section** (0 Inv discordants in that window under the original classifier). Full-chr original log still has InvL 47 pass / InvR 45 pass, so Inv clustering works globally; those InvL DRPs that expand the ARC-SV2 union at this locus are **local to ARC-SV2’s InvL calls**, not present as InvL clusters in the original window. Precedence still recovers the true DEL by keeping split ends; InvL is an extra over-merge ingredient under union, not the sole cause.

### Revised ranked causes (after single-locus proof)

1. **Missing true DELs:** no clip/split precedence in `merge_breakpoints` (proven; 5/5).  
2. **`merge_intervals_for_min_block` / 100 bp absorb:** ARC-SV2-only; not causal alone for the five trues; keep as secondary over-merge risk (D2, D12).  
3. **False PE-only multi-block calls:** survive precedence; need FDR / Del-score (mean vs min-z) / filters.  
4. **`max_paths=1000` skips:** documented; may worsen under over-merged graphs (**hypothesis**).  
5. **Runtime/memory:** in-memory design vs streaming.

---

## Scripts

- `scripts/diagnostics/pick_gap_window.py` — choose window; write orig-only true / v2-only false candidates  
- `scripts/diagnostics/trace_gap_calls.py` — compare candidates to both `arcsv_out.tab` files  
- `scripts/diagnostics/trace_breakpoint_stages.py` — ±20 kb stage dumps + merge flips (scratch precedence)  
- `scripts/diagnostics/window_stage_counts.py` — restrict full-chr call tables to a genomic window  
- `scripts/diagnostics/run_*_loci.sh` — batch original / ARC-SV2 locus runs  

Outputs under `scratch/gap_analysis/` (including `loci/`).
