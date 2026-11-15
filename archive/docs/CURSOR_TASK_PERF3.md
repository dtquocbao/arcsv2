# Cursor task: PERF3 — Docker prune → R2 Rust graph build → chr1 memory

Start: tag `caller-v1.5-lean` + R1 (Rust BAM loader, `engine=rust` default).

| Run (8 threads) | Wall | Peak RSS | Dominant stages |
|:--|--:|--:|:--|
| chr20 HG01884, rust | 10.2 min | not recorded | BAM 183 s, block prep ~200 s, build_graph ~130 s |
| chr1 HG01884, python | 48.5 min | **9.84 GB** | BAM 1047 s, build_graph **1265 s** |

This is the last performance task before the ML phase. Do the steps in order:
**D1 → R2 → M4 → gate**. Stop and report after each step.

## Rules (unchanged)
- `check_golden.py` vs `export/v2fix11/`:
  - `--fast` after each step
  - full chr20 (both samples) at the gate
  - never loosen it
- Rust must be **bit-identical** to Python:
  - keep the same operation order for floats; no unordered parallel float sums
  - edges and reads must come out in the same order as the Python code
- New paths go behind a flag until they pass golden. Then flip the default and keep `--engine python` as the fallback.
- **Every timing row records wall time AND peak RSS** (R1's gate had no RSS figure).
- Keep `scratch/` clean. When the task ends, move useful logs to `bench/logs/` and delete the rest. Remove the stray `HG01884.final.cram.crai` from the repo root.

---

## D1. Docker prune (host only — never the data volume)
`docker system df` shows ~19 GB of reclaimable images and 15.6 GB of build cache.

1. Confirm no container is running a job (`docker ps`).
2. Run these, in this order:
   - `docker container prune -f`
   - `docker image prune -a -f --filter "until=24h"`, but keep `arcsv2:dev` (the current image) and its base `mambaorg/micromamba`
   - `docker builder prune -f --filter "until=24h"`
3. **Forbidden:** `docker volume prune`, `docker system prune --volumes`, `docker compose down -v`, and any `rm` inside `arcsv2-data`.
4. Report `docker system df` before and after, and confirm `arcsv2-data` mounts with `/data/work`, `/data/ref` and `/data/haps` intact.

The next `docker compose build` will be slower because the cache is gone. That's fine.

## R2. Rust kernel: block prep + build_graph
Port the code from `blocks` → `build_graph done` to `arcsv2_native`. It reads the M2 Pairs arrays and the R1 take-records directly.

**Scope**
- Block prep: clip/split stamping onto blocks (the stage between the `blocks` and `build_graph` progress lines).
- `graph.build_graph`, covering:
  - `add_reference_edges`
  - the `splits_by_q` / `used_splits` logic
  - `_resolve_end` and `overlapping_vertices`, including the 1-bp-gap handling
  - `_between_junction`, the `max_edge_span` gate and `_add_hanging`
  - adjacency paths (`adj1`/`adj2` via `block_seq_to_path`)
  - `supports_adjacency` and the flip_parity / between-junction support credit
  - the discordant / distant flags
  - `pmapped` (float: compute it exactly as Python does)

**Output**
- Rust returns edges in columnar form, CSR-style:
  - edge arrays: `u`, `v`, `support`, `reference`, `read_ptr`
  - flat per-read arrays: offset, lib_id, pmapped (f64), mapq, kind, the flag bits, mate_len
  - `adj` encoded as a flat int array with offsets (−1 or an empty span means None)
- Edge order and per-edge read order must match Python's insertion order.

**Python side**
- For R2, rebuild the existing `Graph` / `Edge` / `EdgeRead` objects from those arrays, so everything downstream is unchanged and golden can be checked.
- Do this conversion only for now; M4 removes it.

**Parallelism**
- Assigning pairs to blocks (binary search) and resolving each end are per-pair work. Use rayon with indexed `par_iter` that keeps row order.
- Count edge support in a deterministic pass: either a serial pass over the ordered results, or per-thread integer tallies merged in a fixed order.
- `used_splits` is a sequential dependency (the first pair to use a split wins). Keep that part serial, in pair order.

**Identity test**
- New `tests/test_rust_graph_identity.py`.
- On a 1 Mb window, compare the Python vs Rust graphs at threads 1 and 8: the same edges, supports and reference flags, and every EdgeRead field equal (floats exactly equal).
- Also run it on full chr20 once and save the log.

**Gate**
- Golden `--fast` + full chr20 PASS with `engine=rust`.
- Record wall time and RSS, plus the block prep / build_graph stage times, before vs after.

## M4. chr1 memory reduction
Goal: **chr1 peak RSS ≤ 5 GB** at threads=8, so that 2 chromosome jobs fit on a 16 GB host and 5 fit on 32 GB.

1. **Profile chr1 first** (no code change). Take RSS per stage (and memray if it's affordable) for chr1 with the R2 engine. Report:
   - bytes per pair
   - bytes per EdgeRead and the EdgeRead count
   - the top allocation sites

   Expected: per-read `EdgeRead` dataclass objects and `Edge.reads` lists dominate after build_graph, now that pairs are already columnar.
2. **Columnar edge reads.**
   - Keep the R2 CSR arrays as the real storage and drop the `EdgeRead` objects.
   - Give `Edge` a lightweight read view (indices into the CSR arrays).
   - Update the consumers to work on arrays/slices: `likelihood.score_path`, placement distances, hanging-read handling, insertion scoring, sr/pe columns in output, `_locus_view`.
   - Mutable per-locus state must live in per-locus arrays, never in the shared CSR. This avoids the M3 edge-sharing bug.
   - Golden must PASS.
3. **Free memory early.**
   - Drop the pair arrays and take-records once build_graph and insertion prep no longer need them. Check which later stages still read pairs (hanging, insertions).
   - Call `gc.collect()` at stage boundaries.
   - Report RSS after each change.
4. **Genotype workers.** With fork, workers should inherit the CSR arrays read-only. Pass `(start, end)` only, and copy only the per-locus mutable slices. RSS at t8 should stay at about the t1 level.
5. **Snakemake memory.** Set `call_chrom`'s `mem_mb` from the chromosome length, using the measured MB per Mb plus 20%, instead of a flat 8000 that understates chr1–5. Document the formula.

**Gate (end of PERF3)**
- Full chr20, both samples, golden PASS (engine=rust).
- chr1 HG01884 at 8 threads: report wall time + peak RSS.
- Update the whole-genome estimates:
  - serial run
  - packed run on 16 GB and 32 GB hosts, using the new per-chromosome memory figures
- Tag `caller-v1.6-ready`.

Then stop. The next task is the ML phase.
