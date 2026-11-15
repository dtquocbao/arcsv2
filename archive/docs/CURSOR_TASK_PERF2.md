# Cursor task: PERF2 — Docker cleanup, memory, multi-threading, Rust kernels

Start: tag `caller-v1.4-fast` (C5b). HG01884 chr20 @ 8 threads: **12.4 min, 5.57 GB**, golden PASS.

Stage profile (c5final, 745 s): BAM load **260 s** · fit_libs 15 s · discover 6 s ·
block prep (blocks → build_graph) **198 s** · build_graph **130 s** · insertions 90 s ·
genotype 40 s.

## Rules (same as before, plus Rust)
- `check_golden.py` must PASS vs `export/v2fix11/`: `--fast` after each step,
  full chr20 for both samples at each gate. Never loosen tolerances.
- **Stop and report after each numbered step.** Don't run ahead.
- One commit per step. Every new path sits behind a config flag until it passes
  golden; then flip the default and keep the old path as a fallback.
- Rust kernels must give **bit-identical** results to Python:
  - no `-ffast-math` style flags
  - no parallel float reductions with an unordered sum (rayon `sum()` over floats);
    reduce in the same order Python uses
  - integers are fine to parallelise freely
- Targets at the gate (chr20, threads=8):
  - peak RSS **< 2.5 GB** (main goal; needed for the whole genome)
  - wall time **< 8 min/sample**

---

## D0. Docker cleanup (do first)

**Never delete the `arcsv2-data` volume.** It holds the BAMs, reference,
assemblies and every result. No `docker system prune --volumes`, no
`docker volume prune`, no `docker compose down -v`.

1. Report first: `docker system df -v`, `docker ps -a`, `docker images -a`.
2. Remove stopped containers: `docker container prune -f`. This includes the
   old `arcsv2-*-batch` runs. Check that nothing is still running a job first.
3. Remove dangling and old image layers: `docker image prune -f`. Remove old
   tagged images except `arcsv2:dev` and its base. Then run `docker builder prune -f`.
4. Inside the volume (`/data`): list the sizes of `results_*` and `labels_*`. Delete
   only intermediate run dirs that are already exported to `export/` or are
   superseded:
   - **Keep:** `work/`, the references, the assemblies, `results_c5final`,
     `labels_caller_v1`, the original-ARC-SV results, and anything `check_golden` reads.
   - Post the list, with sizes, of what you plan to delete, then **wait for my OK**
     before removing anything in `/data`.
5. Report the space freed. Note: the Rust build in R0 will grow the image, so check
   the new size there.

## M1. Memory profile per stage (no code change)
Use `memray` or tracemalloc snapshots at the end of each stage (load, block prep,
build_graph, insertions, genotype). Report:
- the top 10 allocation sites
- bytes per pair
- RSS at each stage

Expected: ~8.1M Python pair/read objects dominate.

## M2. Columnar pair storage (Python, numpy)
Replace the per-pair Python objects with columnar numpy arrays. This is the
data layout the Rust kernels will read and write without copying.
- Use a `Pairs` struct-of-arrays:
  - positions int64
  - insert / lengths int32
  - strand / flags uint8
  - lib id uint8
  - clip fields
- Keep thin per-read views only where code really needs them.
- Free BAM-side temporaries after load.
- Target: < 100 bytes/pair. Golden PASS.

## M3. Lean genotype workers
The workers fork-share the graph and pair arrays read-only (no pickling of
reads). RSS at threads=8 must be about the same as at threads=1.

## T1. Multi-threading audit (Python side)
Measure what each stage does with threads=1/4/8. Then:
- **BAM load:**
  - Use pysam `threads=` for decompression (already in).
  - Try region-split loading: N processes, each `fetch`es a slice of the
    chromosome and returns arrays. Mates that cross slice borders are
    joined by qname in a final merge.
  - Golden PASS required. If the merge is fragile, skip this and let R1 cover it.
- **Block prep / build_graph:** vectorise with numpy (`searchsorted`, `bincount`).
  numpy releases the GIL, so a thread pool over chromosome chunks works without
  copying memory.
- **Insertions (90 s):** parallel scoring is only worth it if it's cheap. Score sites
  in parallel against the pre-insertion graph, then commit in order, and re-score
  serially any site whose window contains an already-committed insertion.
- **Per-chromosome parallelism (C4):** check that the Snakemake per-chrom
  rules run concurrently. Use `resources: mem_mb` so jobs don't exceed RAM.
- Report the per-stage speed-up for 1/4/8 threads.

## R0. Rust toolchain in Docker
- Add `rust` + `maturin` to the image. conda-forge in `environment.yml` is the
  simplest; a multi-stage Dockerfile also works.
- Create a PyO3 + `numpy` crate (`arcsv2-native` or extend `core/arcsv-core`)
  that builds a Python wheel installed in the image.
- Python imports it optionally: `try: import arcsv2_native except ImportError:
  use Python`. Config: `engine: python | rust`, default `python`.
- Keep `cargo test` in CI. Update `docs/CROSS_PLATFORM.md`: Windows works
  through Docker; native Windows uses the Python fallback.
- Audit the existing `core/arcsv-core` modules (library, drp, paths, affected,
  names). Are they tested against the current Python (caller-v1.4)? Reuse what
  still matches; delete or archive what is stale. Report which.

## R1. Rust kernel: BAM → columnar evidence (biggest win: time + memory)
- Use `rust-htslib` (CRAM + multithreaded BGZF; links htslib) or `noodles`.
  Pick the one that reads the project's CRAM/BAM in the image; justify the choice.
- Apply the same filters as `samio.py` / `evidence.py`:
  - MAPQ thresholds, secondary/dup/QC/supplementary skips
  - soft-clip quality rules and SA parsing
  - insert stats
- Output the M2 arrays directly into numpy (no per-read Python objects).
- Use rayon over chromosome regions with a deterministic merge order.
- Gate: arrays identical to the Python loader on chr20 (add a test that
  compares them); then golden PASS with `engine: rust`.

## R2. Rust kernel: block prep + build_graph
- Do pair → block assignment, clip stamping and edge-support counting in Rust
  over the M2 arrays.
- Parallelise with rayon over pairs. Edge counts are integers, so merge the
  per-thread tallies.
- Return the edges and supports in the same sorted order as Python.
- Golden PASS.

## R3. (Only if R1–R2 leave it in the top 2) Rust `score_path` / likelihood
- Floats: keep the exact Python operation order and the `log` / `exp` calls.
- Compare the per-locus log-likelihoods for **exact** equality on the golden
  windows before running the full chr20.
- If exact equality isn't reachable, stop and report. Don't change genotype
  outputs.

## Gate
- Full chr20, both samples, threads=8, `engine: rust`: golden PASS.
- Run once with `engine: python` too, to prove the fallback still passes.
- Add rows to `bench/RESULTS.md` for python vs rust × threads 1/8: per-stage wall
  time + peak RSS.
- Tag `caller-v1.5-lean`.
- After the gate, estimate whole-genome time and RAM per sample.

## Housekeeping
- `scratch/` has ~80 run scripts and logs from Part C/C5:
  - Move the useful ones (c0 profile svg, gate logs, c5final logs) to
    `bench/logs/`.
  - Empty the rest, including the `*.bak` and `kill_*` scripts.
- Clean `scratch/` at the end of every task from now on.
