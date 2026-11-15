# Questions for Dr. Bo (undocumented / partially documented ARC-SV behaviors)

Derived from `docs/ORIGINAL_GAP_ANALYSIS.md`. Only items that are **not** or only **partially** described in Cell 2024 Methods (nihms-2024296) or the chapter provided (`Chapter 1.pdf` is background only). Grouped by stage; most important first. Max 15.

---

## Breakpoint merge (highest priority)

1. **Precedence when merging PE vs softclip/split breakpoints.**  
   The original code assigns softclip/split breakpoints precedence 2 and PE-only breakpoints precedence 1 (`arcsv/breakpoint_merge.py` `Breakpoint.precedence` ~63–71, `__add__` ~32–45, `bp_mergefun_precedence` ~261–294). When a PE interval touches a clip/split interval, support is combined but the **coordinates stay on the clip/split interval**. Cell Methods (p20–21) say merge all overlapping intervals / keep if DRP or ≥2 split/clip — **no precedence**.  
   **Single-locus proof (HG01884):** on five original-only true ~300 bp DELs, ARC-SV2 union-merge collapses both ends into one interval → 0 calls; scratch precedence recovers `rearr=AC` with high SR on all five. Was precedence intentional method, and should ARC-SV2 reproduce it (e.g. as a flag, default on for parity)?

2. **Distance for that PE→clip/split absorb.**  
   Final merge uses `max_distance=0` (intervals must touch or overlap). Is “touch only” the intended rule for combining PE CIs with clip/split points, or should a small positive slop (bp) be allowed?

3. **100 bp DRP absorb after merge.**  
   ARC-SV2 adds `absorb_nearby_drp_breakpoints(..., max_distance=100)` after the main merge. The original tree has no equivalent pass. Cell Methods do not describe a 100 bp absorb. Flip experiments: removing absorb alone does **not** recover the five missing true DELs. Should this exist, and if so what distance (suggest: flag, default off)?

4. **`merge_intervals_for_min_block` (ARC-SV2-only).**  
   ARC-SV2 joins breakpoints when the intervening reference would be `< min_affected_length` (50 bp). Original has no such step; not in Cell Methods. Alone it does not explain the five missing trues (union already collapsed them), but it further reduces BP counts on dense loci. Keep, drop, or flag (default off)?

5. **`min_bp_support` vs `min_edge_support`.**  
   Original defaults both to 2 but as **separate** options (`arcsv_call_options.py` 40–41, 60). Non-PE breakpoints need `(supp_split + supp_clip) ≥ min_bp_support`. ARC-SV2 reuses one `min_edge_support` for both breakpoint keep and graph edges. Should these remain two knobs, and are both “2” the values to keep for paper mode?

---

## Discordant clustering / FDR

6. **Max PE-cluster size formula.**  
   Original sets `max_pecluster_size = physical_coverage × 1.5` (`bamparser_streaming.py` ~298–299; HG01884 phys cov 56× → ~84). Components larger than this are kept only if fully mutually compatible, else dropped (`pecluster.py` ~202–210). ARC-SV2 hard-codes 80. Cell text discusses coverage-adaptive FDR/clustering (p18–19) but not the 1.5× physical-coverage cap. Is `1.5 × phys_cov` the method, and should paper mode use it?

7. **Deletion cluster score (LR) formula.**  
   Original `lr_del` uses MLE size `min(max_span, mean(insert) − μ)` (`pecluster.py` ~303–309). ARC-SV2 uses `D = min(min(z) − median, min(b) − max(a))` (`engine/breakpoints.py` `score_cluster`). Cell (p18–19) uses **mean** z − μ. Which formula is intended (ARC-SV2 `min z` looks like a bug vs Methods)?

8. **Null simulation details.**  
   Original relocates discordant pairs along gap-free chromosome length, `pecluster_null_reps=4`, `pecluster_target_fdr=0.1`. Del nulls cap insert at `3μ`. Are these values and the relocation model part of the published method, or tunable implementation choices?

---

## Paths / genotyping

9. **Path cap: paper 1,000 vs code 2,000.**  
   Cell Methods skip subgraphs with >1,000 paths (p25). Original code default is `max_paths=2000` (`arcsv_call_options.py:69`). Which should a “paper mode” reimplementation use when comparing accuracy to the published caller binary?

10. **`max_mes_extra` default 0.**  
    Original can raise `min_edge_support` and retry when paths explode, but default `max_mes_extra=0` skips immediately. Was disabling escalation intentional for published runs?

---

## Library model

11. **Approximate insert sampling vs full-region KDE.**  
    Original uses `approx_stats_nreads=1e6` / `nchunks=1000` then re-estimates the KDE on all pairs in the region. ARC-SV2 fits one pooled KDE on all loaded pairs. Is the two-step approx procedure required for parity, or only an optimization?

12. **Pooling libraries.**  
    ARC-SV2 defaults `pool_libraries=True`. Original is set up for multi-library but these BAMs behave as one. For multi-RG samples, should paper mode pool or keep per-library models?

---

## Soft clips / splits (confirmation)

13. **Softclip filters and merge.**  
    Defaults: mapq≥10, ≥5 clipped bases, median Q≥15, merge within 5 bp, min support 2, microhomology ≤25 bp (Cell p17). Confirm these numeric thresholds are the ones used for the published HGSVC/1KGP calls.

14. **Single split endpoint as a breakpoint.**  
    A lone split contributes `n_split=1` and is dropped unless combined with another split/clip or PE (`min_bp_support=2`). Is that intentional (require two pieces of local evidence), or should a single high-MAPQ split suffice?

---

## Output filters

15. **`min_simplesv_size=50` on the live call path / merge-filter BEDs.**  
    Original `do_inference` does not apply `min_simplesv_size`; ARC-SV2 applies `min_affected_length=50` when emitting. Cell p27 merge-filter also uses repeats/segdup/blacklist/gaps — should paper-mode labeling always apply those BEDs?
