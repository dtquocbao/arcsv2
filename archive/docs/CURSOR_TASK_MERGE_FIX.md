# Cursor task: adopt the original's breakpoint merge, re-measure chr20

Proven in `docs/ORIGINAL_GAP_ANALYSIS.md` (single-locus traces): union merge collapses both ends of
~100–500 bp deletions. Why: DRP breakpoint CIs are ~200 bp wide and extend toward the deleted segment
(left CI rightward, right CI leftward), so for a 322 bp DEL at chr20:2883362 the left CI
[2883365,2883569) and right CI [2883496,2883682) overlap and union into one interval. The original's
clip/split precedence keeps the precise ends.

## 1. Engine changes (with tests)
- `merge_mode: "precedence" | "union"` in `CallConfig`; default **precedence** (port
  `bp_mergefun_precedence` semantics exactly: precedence 2 = any split or clip, 1 = DRP only; merge
  same-level intervals first with distance 0, then fold lower-level intervals that touch into existing
  higher-level ones without moving coordinates; then drop intervals with no DRP and split+clip < 2).
- Flags, default **off**: `absorb_nearby_drp` (100 bp), `min_block_merge` (<50 bp).
- Fix `score_cluster` deletion size to Methods eq. 2 / `pecluster.lr_del`: `min(mean(z) − μ, min b − max a)`.
- Unit tests: the 2883362 geometry (two overlapping DRP CIs + split ends) yields two breakpoints under
  precedence and one under union; Del score matches `lr_del` on a fixed cluster.
- Investigate (report only): at chr20:5962348 the union interval grew to [5954931,5972086) through InvL
  clusters (clusters 1,2). Check whether those InvL clusters also exist in the original (`-v 2` logs in
  `scratch/gap_analysis/loci/orig_5962348/`).

## 2. Re-run and re-label chr20 (both samples, one at a time)
Same pipeline as before (`label.smk`, `-j 1`), output to `/data/results_v2fix/`. Then compare against
`export/orig/` with the same script used for the original-vs-ARC-SV2 table: calls, validated calls,
precision per class, shared vs tool-only, sr=0 share, wall time, RSS.
Record the row in `bench/RESULTS.md` and commit first so the row has a git revision.

## 3. False positives (analysis only, next step)
Precedence does not remove the pair-only false calls (28549495, 61390333, 64235350). For HG01884 compare
DRP clusters passing FDR per type: original log says Del 298 pass / 2747 fail, Dup 295/179, InvL 47,
InvR 45, Ins 19. Print ARC-SV2's equivalent counts and list the first stage where they differ.
