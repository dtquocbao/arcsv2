# Cursor task: cleanup + performance before the ML phase

Baseline: **v2fix11** (chr20, HG01884 + NA19700): validated 222, precision 59.0%,
221/243 original TPs, ~50 min/sample call, `export/v2fix11/`.

**Golden rule for this whole task: no behavior change.** Every part ends with a
re-run that must produce byte-identical `merged.tab` (both samples) and an
identical `dataset.tsv` label table vs `export/v2fix11/`. If anything differs,
stop and report the diff; do not "fix forward".

Work in this order: Part A (freeze) → Part B (cleanup) → Part C (perf). One
commit per numbered step so any step can be reverted alone.

---

## Part A — Freeze the baseline (do first)

1. Commit everything that is currently uncommitted as `v2fix11: simplify_blocks_diploid (flank-preserving)`.
2. Tag it `caller-v1` (annotated; message = the v2fix11 numbers above).
3. Add `scripts/check_golden.py`: given a results dir + labels dir, compare
   against `export/v2fix11/` (sorted `merged.tab` rows, `dataset.tsv` rows,
   validated/precision summary). Exit non-zero on any difference, print the
   first 20 differing rows. Every later step runs this.
4. Add a fast golden check too: a small window set (the loci from
   CURSOR_TASK_RECALL: 2659058, 5962348, 52764531, 60792115 + the 3 weak-DEL
   loci) whose v2fix11 calls are saved under `tests/golden/`. Runs in minutes;
   use it after every step, the full chr20 check at the end of each Part.

## Part B — Cleanup

**Archive, don't delete.** Everything removed from the working tree stays in
git history at `caller-v1`. Move to `archive/` only what is worth browsing later;
`git rm` the rest. Do not run `rm -rf` on data dirs outside git.

### B1. Root clutter
- `git rm`: `_tmp_find_labeled.ps1`, empty `arcsv2/` dir.
- `ref.fa`, `*.cram.crai`: confirm nothing in `workflows/` or `scripts/` uses
  the repo-root copies (data lives in `/data`); then remove and add to `.gitignore`.
- `pilot/`, `pilot-large/`, `work/`, `results/`, `labels/`, `haps/`, `target/`,
  `.snakemake/`: already gitignored. Leave on disk; list them in README
  "local scratch you can delete" (~175 MB).

### B2. scratch/ (15 MB)
- `scratch/gap_analysis/`: keep only the final reports (`*report*.md`,
  oracle/ablation summary tables) → `archive/gap_analysis/`. Dump scripts,
  batch logs, `loci/` raw dumps → `git rm`.
- `scratch/pdf_text/`, `scratch/*.sh`, `scratch/*.ps1` (do_commit_drp,
  do_drp_tests, do_smoke_tests, run_drp_tests, run_in_docker, run_merge_tests,
  start_v2fix10_batch) → `git rm`. `scratch/` should end empty and be gitignored.

### B3. scripts/
- Keep: `arcsv2.ps1/.sh`, `fetch_data.py`, `fetch_filter_beds.py`,
  `make_ref_beds.py`, `simulate.py`, `simulate_haplotypes.py`,
  `orig_to_arcsv2.py`, new `check_golden.py`.
- Move to `archive/scripts/`: `check_label_pilot.py`, `pilot_large.py`,
  `make_dashboard_demo.py` (unless the dashboard is still used — say which).
- `scripts/diagnostics/`: keep `trace_breakpoint_stages.py` and
  `window_stage_counts.py` (useful for future debugging, make them run against
  current code); archive the rest (`compare_labeled.py`, `pick_gap_window.py`,
  `run_*_loci.sh`, `trace_gap_calls.py`).

### B4. workflows/config/
- Keep: `data.yaml`, `samples.tsv.example`, `label.yaml` (= v2fix11 settings,
  make it the default), `label.smoke.yaml`.
- `label.v2fix.yaml` … `label.v2fix11.yaml`, `label.pilot.yaml`,
  `label.large.yaml`: `git rm` (they live in history at `caller-v1`). Keep one
  `label.caller-v1.yaml` that reproduces the golden run exactly.

### B5. docs/ — consolidate to:
| Keep / new | Content |
|:--|:--|
| `README.md` | install, run, one-command chr20 benchmark, golden check |
| `docs/METHOD.md` | pipeline as implemented, stage by stage |
| `docs/PARITY.md` | single list: every deviation from paper and from original code, with flag + reason. Merge in DRP_PARITY, EDGE_SUPPORT, ORIGINAL_GAP_ANALYSIS conclusions |
| `docs/ORIGINAL_COMPARISON.md` | orig vs caller-v1 benchmark table |
| `docs/CROSS_PLATFORM.md`, `docs/LABELS_PLAN.md` | as is |
| `docs/HISTORY.md` (new) | one table: v2fix → change → validated/precision → commit. Taken from bench/RESULTS.md prose |
| `bench/RESULTS.md` | timing table only |

Archive to `archive/docs/`: all `CURSOR_TASK_*.md` (except this one),
`reviews/`, `DASHBOARD.md`, `PRODUCTION_PLAN.md`, `DRP_PARITY.md`,
`EDGE_SUPPORT.md`, `ORIGINAL_GAP_ANALYSIS.md`. `QUESTIONS_FOR_DR_BO*.md`:
keep only `_SHORT`.

### B6. Dead flags and code paths (`config.py`, engine)
For each flag, either delete the flag **and** the non-default branch, or keep
it and document it in PARITY.md. Proposed:
- Delete (failed experiments): `absorb_nearby_drp`, `min_block_merge`,
  `use_point_block_cuts` (+ point-cut code), any `COPY_FAMILY_BAND` remnants,
  `merge_mode="union"` branch (keep precedence only; note the paper/code
  difference in PARITY.md).
- Delete if never read outside config (grep first, report): `adaptive_support`,
  `path_search`, `beam_width`, `refine_breakpoints`, `gc_correct`.
- Keep and document: `skip_edge_ref_gate` (default False = original),
  `max_paths`, `max_pecluster_size`, `min_edge_support`, `threads`,
  `min_affected_length`.
- Run `ruff` + `vulture --min-confidence 80 python/` and remove unreferenced
  helpers it finds (report the list in the PR).
- `core/arcsv-core` (Rust) is not imported by Python. Leave it, but add one
  line to README: "not wired in yet; see Part C step 5".

### B7. Tests
- Keep tests that guard parity fixes; delete tests for removed flags.
- Add `tests/test_parity_loci.py` using `tests/golden/` windows:
  2659058 → DUP ABBC; 5962348 → DEL; 60792115 → ABBC; the 3 weak-DEL loci
  → not called; 52764531 → same-block support == 2.
- Mark slow ones `@pytest.mark.slow`; default `pytest` < 2 min.

**Part B gate:** full chr20 re-run (`label.caller-v1.yaml`) →
`check_golden.py` passes. Commit, tag `caller-v1.1-clean`.

---

## Part C — Performance (each step: golden check must pass)

### C0. Profile first (no code change)
Run HG01884 chr20 under `py-spy record --native -o prof.svg` (or cProfile)
plus the existing `Progress` stage timings. Record a table in
`bench/RESULTS.md`: stage → wall time → peak RSS (read BAM / evidence,
breakpoints, DRP clustering, build_graph, decomposition+paths, genotype,
insertions, output). Do the steps below in order of what this shows; report the
table before starting C1.

### C1. Insertion testing (`caller.py::place_novel_insertions`)
Currently every candidate site scores every size **twice** (the second loop
re-scores to find the smallest matching size) and `_graph_plus_insertion` copies
the full block list per size.
- Cache `(size → hom, het)` from the first loop; the second loop becomes a
  lookup. Keep the exact `==` tie rule and size order.
- Build the test graph on a local window (blocks `first..last` + the insertion)
  instead of copying `graph.blocks` per size.
- Keep the sequential `graph.blocks.append` / `add_edge` at the end exactly as
  now. `_pair_spans_boundary` must see the same graph state as before (sites
  are processed in order and later sites see earlier insertions) — so only
  parallelise the scoring, never the commit step, and only if C0 shows it matters.

### C2. Threads end to end
- `workflows/label.smk` hard-codes `--threads 1` (line ~152). Use
  `threads: config.get("threads", 4)` and `--threads {threads}`.
- Confirm genotyping's ProcessPool gives identical output for threads=1/4/8
  (sort results by subgraph id before writing).
- Record wall time at 1/4/8 threads.

### C3. Faster evidence collection (`samio.py`, `evidence.py`)
- `pysam.AlignmentFile(..., threads=N)` for BGZF/CRAM decompression.
- Cheap pre-filter per read before building a `SamRead`: skip unmapped,
  secondary, duplicate, QC-fail, and concordant proper pairs with no soft clip,
  no SA tag and insert within library range — but **still count them** for
  coverage / insert-size stats that need them. Read `query_qualities` only when
  the read is clipped.
- Insert-size distribution: if the code currently uses all pairs, test
  sampling (e.g. first 1M proper pairs). This **may** change the PMF, so it
  goes behind a flag and is only kept if golden passes; otherwise leave it.
- Store per-read evidence as numpy columnar arrays instead of Python objects
  where it is the memory hot spot (C0 decides).
- Evidence cache: write the extracted evidence to `npz`/Parquet keyed on
  (BAM path + size + mtime, region, evidence-relevant config, code version).
  Re-runs that only change downstream params skip BAM reading entirely.

### C4. Per-chromosome parallelism (whole-genome prep)
Snakemake rule per (sample, chrom), then a merge rule. Not needed for chr20
benchmarks but must be in place before the ML-phase scale-up.

### C5. Later (do not start without asking)
Optional noise filters (smoove-style) behind flags; Rust kernels via PyO3 for
whatever C0 shows is still the top stage.

**Part C gate:** golden check passes; `bench/RESULTS.md` has before/after rows
(wall time per stage, peak RSS, threads). Target: chr20 call < 15 min/sample at
8 threads, peak RSS < 2.5 GB. Commit, tag `caller-v1.2-fast`.

---

## Report back
1. Part B: list of removed / archived files and deleted flags; golden PASS.
2. C0 profile table.
3. Per C step: wall time, RSS, golden PASS/FAIL.
