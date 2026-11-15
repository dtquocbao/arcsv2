# Cursor follow-up: prove the breakpoint-merge mechanism on single loci

Review of `docs/ORIGINAL_GAP_ANALYSIS.md` found three gaps. Fix them; do not change engine code except in a
scratch copy for debug output.

## A. New evidence to explain
For the five original-only true deletions (HG01884: 2883362, 5962348, 32311430, 57219054, 61995575),
ARC-SV2 has **no locus at all** covering the deletion in `export/results/HG01884/chr20/graph_log.txt`
(e.g. loci end at 32310932 and the next starts at 32323153). So the call is lost before genotyping:
breakpoints or blocks never form there.

Candidate mechanism to test:
1. `merge_breakpoints` (`engine/breakpoints.py:593`) unions overlapping intervals with no precedence, so a wide
   PE confidence interval can swallow the clip/split breakpoints at both ends of a ~300 bp deletion.
2. `merge_intervals_for_min_block` (`engine/breakpoints.py:618`, called in `caller.py:221` and `:1136` with
   `min_affected_length=50`) then joins breakpoints closer than 50 bp. **The original has no such step.**
   Together these can collapse a deletion's two breakpoints into one interval → no deleted block → no locus.
3. `absorb_nearby_drp_breakpoints` (100 bp, ARC-SV2 only) may add to it.

## B. Tasks
1. For each of the 5 loci, run ARC-SV2 on a ±20 kb region (`--region chr20:start-end`) inside the container
   with a debug dump of breakpoints after each step: clip clusters, splits, DRP clusters (with intervals),
   after `merge_breakpoints`, after `merge_intervals_for_min_block`, after `absorb_nearby_drp_breakpoints`,
   final blocks. Run the original on the same region with `-v 2` (it prints `all_bp`, `just merged bp`,
   `merged and filtered bp`). Put both side by side per locus.
2. State the exact step where the deletion's breakpoints disappear or merge. Then flip one thing at a time in
   the scratch copy (a) precedence merge like `bp_mergefun_precedence`, (b) remove min-block merge,
   (c) remove 100 bp absorb, and report whether the locus and call appear.
3. Do the same for 3 ARC-SV2-only false calls (28549495, 61390333, 64235350).
4. Complete the planned stage-count table on chr20:63–66 Mb (both tools, same window).

## C. Documentation corrections
- `Chapter 1.pdf` is the thesis background chapter (no ARC-SV method text). Say "not in Cell 2024 Methods or
  the chapter provided", not "not in the thesis".
- Replace approximate page numbers with real ones (extract per page with pypdf `pages[i].extract_text()`).
- Add D12: `merge_intervals_for_min_block` exists only in ARC-SV2.

## Deliverable
Append a section "Single-locus traces" to `docs/ORIGINAL_GAP_ANALYSIS.md` with the step-by-step breakpoint
dumps, the flip experiments, and a revised ranked cause list.

## D. Paper check results (Cell 2024 Methods, nihms-2024296 pp. 17–27) — use these, not the approximate pages
| Item | Paper | Original code | ARC-SV2 | Action |
|:--|:--|:--|:--|:--|
| Breakpoint merge | union of overlapping intervals; keep if DRP or ≥2 split/clip (p. 20–21) | precedence: clip/split coords win over DRP | union (matches paper) | paper ≠ code; test precedence as a flag; ask Dr. Bo |
| 100 bp DRP absorb | not described | absent | present | ARC-SV2 invention; put behind a flag, default off |
| Min-block merge (<50 bp) | not described | absent | present | ARC-SV2 invention; flag, default off |
| Deletion size D (eq. 2) | min(**mean** z − μ, min b − max a) (p. 18–19) | mean (`pecluster.lr_del`) | **min z** (`breakpoints.score_cluster`) | **bug: fix to mean** |
| DRP cluster cap | none (largest cliques) (p. 18) | 1.5 × phys coverage | fixed 80 | add `coverage × 1.5` mode |
| Path limit | 1,000 (p. 25) | 2,000 | 1,000 | keep paper; add code value as option |
| Soft-clip/split thresholds | ≥5 bp, Q≥15, MAPQ≥10, merge 5 bp, MH 25 bp; splits MAPQ≥20 (p. 17) | same | same | none |
| Validation labeler | 150 bp flanks, full query coverage, no clipping, MAPQ 60, >90% identity, gap <49 bp (p. 25–26) | — | `label.py` matches | none |
| Merge-filter | repeats >85%, segdup >50%, <50 bp, blacklist, gaps (p. 27) | `filter-merge` | `merge-filter` supports BEDs, but `label.smk` passes none | add BEDs to label.smk |

Run the flip experiments in section B with: (a) precedence on, (b) min-block off, (c) absorb off, (d) Del mean fix,
alone and combined, and re-score chr20 with the labeler.
