# Questions for Dr. Bo — ARC-SV paper vs code

Source checked: Zhou et al., Cell 2024 (author manuscript, `nihms-2024296.pdf`), STAR Methods
"Probabilistic modeling in ARC-SV" (manuscript pp. 17–25) and "Machine learning models" (p. 26).
Code: SUwonglab/arcsv @ dd6f795.

Context: ARC-SV2 was re-implemented from the Methods. On chr20 (HG01884, NA19700) the original code is
much more accurate (deletion precision 67% vs 24%). Several differences are places where the **code does
something the Methods do not describe**. I would like to confirm which behavior is intended.

1. **Breakpoint merging (most important).** Methods p. 20–21: "We merge all overlapping breakpoint
   intervals and exact breakpoints ... and keep merged intervals with a supporting DRP cluster or at least
   two supporting split/clipped reads." The code instead ranks soft-clip/split breakpoints above
   DRP-only ones: when they overlap, it keeps the **soft-clip/split coordinates** and only adds the DRP
   support (`breakpoint_merge.py`: `Breakpoint.precedence`, `bp_mergefun_precedence`). ARC-SV2 followed the
   Methods (plain union). On five HG01884 true ~300 bp DELs, union collapses both ends into one interval
   (0 calls); scratch precedence recovers high-SR `AC` deletions on all five. Is the precedence rule the
   intended method?

2. **ARC-SV2-only post-merge steps (not in Methods or original).** (a) `absorb_nearby_drp` at 100 bp;
   (b) `merge_intervals_for_min_block` when gap <50 bp. Alone they do not recover those five trues, but
   they further over-merge. Keep behind flags (default off)?

3. **Deletion size D.** Methods p. 18–19 use min(**mean** z − μ, …). Original code uses mean; ARC-SV2
   uses **min z**. Should we fix ARC-SV2 to mean?

4. **DRP cluster size cap.** Methods p. 18: largest cliques, no size limit. Code caps at 1.5 × physical
   coverage (`bamparser_streaming.py`, `pecluster.py`). Is the cap intended, and is 1.5 the published value?

5. **Path limit.** Methods p. 25: >1,000 paths skipped; code default `max_paths = 2000`. Which was used
   for published 1000GP/HPRC calls?

6. **Filters for the ML training data.** Methods p. 27: >85% simple-repeat, >50% segdup, <50 bp,
   blacklist, gaps. Were ML labels (p. 26) built **after** this merge-filter? Which BED files?

7. **ML features (Table S1) / Arthur thesis.** Could you share Table S1, the 31/11 HPRC split, and the
   thesis chapters on breakpoint discovery and the likelihood model?
