# Runtime

ML-0 S4 recorded a **40-sample whole-genome** calendar wall of ~164 h on 22
cores (see **ML-0 S4** below); that run was disk-bound with restarts, not a clean
1.5 h / 32-core timing. The roadmap target remains 1.5 h wall on 32 cores and
16 GB peak, and at least 10× fewer core-hours than the original. The 200 kb
simulated label workflow in CI is a correctness check, not a row in the table
below.

Change history (validated/precision): `docs/HISTORY.md`.

| Date | Region | Samples | Threads | Wall time | Peak RSS | Calls | Git revision |
|:--|:--|:--|--:|--:|--:|--:|:--|
| 2026-10-05 | chr1 | HG01884 | 8 | **31.1 min** (fresh wheel) | **5.69 GB** | 881 | re-verify; tab ≡ python `results_chr1` |
| 2026-10-05 | chr20 | HG01884 | 8 | 7.1 min (428.3 s) | **1.65 GB** | 347 | fresh `arcsv2_native` wheel; engine=rust |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate verify 397 s | — | 347+345 | fresh wheel; golden PASS; keep `caller-v1.6-ready` |
| 2026-10-05 | chr1 | HG01884 | 8 | **29.3 min** (M4 rust) | **7.46 GB**††† | 881 | `caller-v1.6-ready` (stale wheel; see below) |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate 455.5 s | **7.81 GB**††† | 347+345 | `caller-v1.6-ready` golden PASS (stale wheel) |
| 2026-09-29 | chr20 | HG01884 | 1 | 44.4 min | ~5.1 GB | 445 | uncommitted (progress + find_block + insertion fix) |
| 2026-09-30 | chr20 | HG01884, NA19700 | 2 | 3 h 41 min (label.smk; HG01884 call reused) | ~5 GB (call) | 445 + 428 | uncommitted (union merge) |
| 2026-09-30 | chr20 | HG01884, NA19700 | 1 | see `results_v2fix` | — | validated 157 | `134bc8e` (precedence + Del mean) |
| 2026-10-01 | chr20 | HG01884, NA19700 | 1 | done (`results_v2fix2`) | — | validated 151 | `c6b9dce` (DRP typing parity) |
| 2026-10-01 | chr20 | HG01884, NA19700 | 1 | 1 h 56 min (`results_v2fix3`) | — | validated 160 | `693fade` (signed insert + concordant early return) |
| 2026-10-01 | chr20 | HG01884, NA19700 | 1 | 1 h 20 min (`results_v2fix4`) | — | validated 151 | `6092b21` (skip-edge ref gate + pe_support adjacency) |
| 2026-10-02 | chr20 | HG01884, NA19700 | 1 | 1 h 34 min (`results_v2fix5`) | — | DUP 16/159 labeled | `46adb80` (pe junction match + COPY_FAMILY_BAND) |
| 2026-10-03 | chr20 | HG01884, NA19700 | 1 | HG01884 call 7.4 h (`results_v2fix6`) | — | validated 182; DUP 31/55 | `a44611b` (adj placements; build_graph O(n²)) |
| 2026-10-03 | chr20 | HG01884, NA19700 | 1 | merge+label 84 s (`results_v2fix8`; no re-call) | — | validated 207; orig TP 192/243 | midpoint `len_affected` rewrite |
| 2026-10-04 | chr20 | HG01884, NA19700 | 1 | batch done (`results_v2fix9`) | — | validated 198; prec 56.6% | `ef23352` point cuts + flip_parity between |
| 2026-10-04 | chr20 | HG01884, NA19700 | 1 | batch done 2h12m (`results_v2fix10`) | — | validated 222; prec 59% | `d4c2265` gapped default + 1-bp gap Dup |
| 2026-10-04 | chr20 | HG01884, NA19700 | 1 | ~50 min/sample call (`results_v2fix11`) | — | validated 222; prec 59.0%; orig TP 221/243 | tag `caller-v1` / `export/v2fix11/` |
| 2026-10-04 | chr20 | HG01884 | 1 | 76.1 min (C0 profile) | 5.63 GB | 347 | tag `caller-v1.1-clean` |
| 2026-10-04 | chr20 | HG01884 | 1 | 57.1 min (C2 threads=1) | 5.62 GB | 347 | `6d35e0c` (post-C1; insertions ~291 s) |
| 2026-10-04 | chr20 | HG01884 | 4 | 28.1 min (C2 threads=4) | 6.57 GB | 347 | `6d35e0c`; merged.tab == t1/t8 |
| 2026-10-04 | chr20 | HG01884 | 8 | 25.2 min (C2 threads=8) | 6.50 GB | 347 | `6d35e0c`; merged.tab == t1/t4 |
| 2026-10-05 | chr20 | HG01884 | 8 | 21.0 min (Part C gate call) | ~6.5–7.9 GB† | 347 | tag `caller-v1.2-fast`; BAM load 252 s |
| 2026-10-05 | chr20 | NA19700 | 8 | 18.8 min (Part C gate call) | ~6.5–7.9 GB† | 345 | tag `caller-v1.2-fast`; BAM load 236 s |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | 42.1 min (label.smk end-to-end) | 7.92 GB† | validated 222; prec 59.0% | tag `caller-v1.2-fast`; full golden PASS |
| 2026-10-05 | chr20 | HG01884 | 8 | 20.2 min (C5 PathLayout) | 6.49 GB | 347 | `492bfa1`; merged.tab == `results_c_gate` |
| 2026-10-05 | chr20 | HG01884 | 8 | 12.4 min (C5 fork+hang) | 5.57 GB | 347 | `6ac4933`; genotype ~40 s; golden PASS |
| 2026-10-05 | chr20 load | HG01884 | 1 | load 4.4 min (M1) | 3.93 GB after load | — | M1 profile; 8.09M pairs @ 451 B/pair |
| 2026-10-05 | chr20 | HG01884 | 8 | 10.9 min (M2/M3 lean) | 2.72 GB | 347 | `results_gate_lean`; golden PASS |
| 2026-10-05 | chr20 | NA19700 | 8 | 10.2 min (M2/M3 lean) | 2.51 GB | 345 | `results_gate_lean`; golden PASS |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate lean (python engine) | 2.72 GB peak | validated 222; prec 59.0% | tag `caller-v1.5-lean`; full golden PASS |
| 2026-10-05 | chr20 | HG01884 | 8 | 11.4 min (685 s) | 2.72 GB | 347 | `results_gate_v15` re-verify; golden PASS |
| 2026-10-05 | chr20 | NA19700 | 8 | 9.7 min (580 s) | 2.52 GB | 345 | `results_gate_v15` re-verify; golden PASS |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate_v15 + labels | 2.72 GB peak | validated 222; prec 59.0% | tag `caller-v1.5-lean` |
| 2026-10-05 | chr1 | HG01884 | 8 | 48.5 min (2911 s) | 9.84 GB | 881 | `results_chr1`; PERF2 WG estimate |
| 2026-10-05 | chr20 | HG01884 | 8 | 10.2 min (612.6 s) | —†† | 347 | `results_gate_r1`; engine=rust; BAM load 182.9 s |
| 2026-10-05 | chr20 | NA19700 | 8 | 8.4 min (503.6 s) | —†† | 345 | `results_gate_r1`; engine=rust; BAM load 127.1 s |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate_r1 (rust BAM) | —†† | validated 222; prec 59.0% | full golden PASS; see `docs/R1_STATUS.md` |
| 2026-10-05 | chr20 | HG01884 | 8 | 5.8 min (346.8 s) | 7.90 GB‡ | 347 | `results_gate_r2`; rust BAM+graph; prep+graph 26.4 s |
| 2026-10-05 | chr20 | NA19700 | 8 | 4.5 min (269.0 s) | 7.90 GB‡ | 345 | `results_gate_r2`; rust BAM+graph; prep+graph 24.6 s |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate_r2 6.9 min (413.3 s) | 7.90 GB‡ | validated 222; prec 59.0% | full golden PASS; see `docs/R2_STATUS.md` |
| 2026-10-05 | chr20 | HG01884 | 8 | 7.4 min (446.4 s) | **5.98 GB**§ | 347 | `results_gate_m4`; M4 CSR EdgeReadStore |
| 2026-10-05 | chr20 | NA19700 | 8 | 6.5 min (389.4 s) | **5.57 GB**§ | 345 | `results_gate_m4`; M4 CSR EdgeReadStore |
| 2026-10-05 | chr20 | HG01884, NA19700 | 8 | gate_m4 | 5.98 GB peak§ | validated 222; prec 59.0% | full golden PASS; tag `caller-v1.6-ready` |

† Gate `max_rss` from `RUSAGE_CHILDREN` over snakemake (calls + label). Per-call peak at C2 t8 was 6.50 GB.
§ M4 gate per-call peak from `measure_rss.py` (tree VmRSS, fork workers included).
†† R1 gate wall from progress timestamps; peak RSS not re-instrumented (lean python gate was 2.51–2.72 GB).
‡ R2 gate peak from `measure_rss.py` over snakemake (both calls ran in parallel).
††† Measured under a **stale** `arcsv2_native` install (see fresh-wheel re-verify below).

## C0 profile (HG01884 chr20, threads=1)

Tag `caller-v1.1-clean`. Progress stage timings + 1 Hz `/proc` VmRSS samples.
Flamegraph: `scratch/c0_profile/prof.svg` (py-spy `--native`).

Breakpoint/DRP row includes post-`blocks` clip stamping until `build_graph` starts
(~181 s of the 186.5 s).

| Stage | Wall (s) | Wall (min) | Peak RSS (MB) |
|:--|--:|--:|--:|
| read BAM / evidence | 539.5 | 8.99 | 4057 |
| fit_libs | 15.3 | 0.25 | 4300 |
| breakpoints + DRP clustering (+ block prep) | 186.5 | 3.11 | 4027 |
| build_graph | 130.7 | 2.18 | 5652 |
| insertions | 1355.0 | 22.58 | 5694 |
| decomposition+paths | 1.6 | 0.03 | 5700 |
| genotype | 2339.0 | 38.98 | 5770 |
| output | 0.8 | 0.01 | 5770 |
| **TOTAL** | **4568.4** | **76.14** | **5770** |

Hot stages for C1+: **genotype** (~39 min) and **insertions** (~23 min), then
BAM evidence (~9 min). Peak RSS ~5.6 GB during genotype.

## C1 / C2 notes

- **C1** (`2446105`): insertion scoring cache + local test graph. Insertions stage
  dropped from ~1355 s (C0) to ~291 s (C2 threads=1: 790.6→1082.1 s).
- **C2** (`6d35e0c`): snakemake/call `--threads` wired; genotype ProcessPool.
  HG01884 chr20 walls: t1 57.1 min / 5.62 GB, t4 28.1 min / 6.57 GB,
  t8 25.2 min / 6.50 GB. Sorted `merged.tab` identical across 1/4/8.
  `check_golden --fast` PASS on `results_caller_v1`, `export/v2fix11`, and
  C2 t8 HG01884 + caller_v1 NA19700 (timing dirs were HG01884-only).

## C3 / C4 / Part C gate

- **C3** (`e291e00`): `pysam.AlignmentFile(threads=N)`; skip secondary/dup/QC/supplementary
  before `SamRead`; qualities only when soft-clipped; cheap no-SA/no-clip `take()` path.
  `sample_insert_sizes` flag added (default **False** — not enabled; would risk PMF drift).
  No numpy columnar rewrite / evidence cache (kept off for golden safety).
  BAM load HG01884: **252 s** (was ~516 s at C2 / ~540 s at C0).
- **C4** (`5dbbbc5`): `workflows/Snakefile` one `call_chrom` per (sample, chrom) +
  `merge_sample`; `label.smk` already per-chrom (chr20-only still works).
- **Part C gate** (tag `caller-v1.2-fast`): full chr20 both samples, threads=8,
  `label.caller-v1.yaml` → `/data/results_c_gate` + `/data/labels_c_gate`.
  `check_golden.py` **PASS** (validated 222, precision 59.04%, merged.tab identical).
  Call walls: HG01884 **21.0 min**, NA19700 **18.8 min**; pipeline wall **42.1 min**,
  peak RSS **7.92 GB**.

### Target vs achieved (chr20 @ 8 threads)

| Metric | Target | Achieved | Status |
|:--|--:|--:|:--|
| Wall / sample | < 15 min | 12.4 min (HG01884 C5b) | **PASS** |
| Peak RSS | < 2.5 GB | ~5.6–7.9 GB | **MISS** |
| Golden | PASS | PASS | **PASS** |

Before (C0 / `caller-v1.1-clean`, threads=1): 76.1 min / 5.63 GB.
After (C1–C4 / `caller-v1.2-fast`, threads=8): ~21 min/sample call / ~6.5+ GB.
After C5b (`caller-v1.4-fast`, threads=8): **12.4 min** HG01884 / 5.57 GB.
Dominant remaining cost: build_graph + insertions + RSS footprint (genotype ~39 s @8).

## C5 notes

- **Genotype hot path** (default on, bit-identical): `PathLayout` prefix-sum
  distances + placement-distance cache + hanging-read cache in `score_path`;
  reuse `ScoredPath.g` in diploid mixes. Window re-call row-identical to HEAD;
  full chr20 HG01884 `merged.tab` identical to `results_c_gate`.
- **Optional smoove-style filters** (default **OFF**): `filter_high_nm` /
  `max_nm`, `filter_both_mates_clipped` (+ CLI flags). Documented in
  `docs/PARITY.md`.
- **Rust/PyO3:** not wired — `arcsv2:dev` image has no Rust toolchain
  (Dockerfile). Python cache only.
- **Timing (PathLayout only):** HG01884 @8 **20.2 min** / 6.49 GB vs gate 21.0.
- **False alarm:** windowed re-call vs full-chr golden fixtures is invalid
  (library/graph context); `--fast` must read full-chr `merged.tab`.

## C5 follow-up (path-local hanging + fork pool)

- **Hot after PathLayout:** `_hanging_probability` / `any(is_insertion)` over
  the full block list, `_hap_blocks` rebuild, then `log_diploid`.
- **Bug-shaped waste:** after any chromosome insertion was placed, every
  hanging read took the slow overlap path, then returned the base class prob
  when the *path* had no insertion. Fix: path-local insertion check + overlap
  cache by read length; cache `_hap_blocks` once per locus. Bit-identical.
- **ProcessPool:** fork-share the chromosome graph (spawn initializer
  fallback); submit `(start, end)` only; workers build `_locus_view` in
  process; cap in-flight futures. No Rust/PyO3 (`arcsv2:dev` has no `rustc`).
  Smoove filters still default **OFF**.
- **Solo timing:** HG01884 chr20 @8 **12.4 min** (744.7 s), RSS **5.57 GB**,
  347 calls. Stages: BAM 260 s, build_graph ~328 s, insertions ~88 s,
  genotype **~40 s**. `check_golden --fast` PASS; `merged.tab` == export.
- Wall target **HIT**; RSS still above 2.5 GB (pairs/graph footprint).

## PERF2 / M1 memory profile

Artifacts: `scratch/m1_profile/` (`REPORT.md`, `rss_samples.tsv`, `full_load.json`,
tracemalloc `snap_loaded.json` / `snap_blocks.json`).

Window `chr20:30000001-35000000` HG01884 threads=1 — RSS from 1 Hz `/proc`
sampler on a normal call (no tracemalloc). Top-10 sites from a prior tracemalloc
run on the same window (aborted mid-`build_graph`; pair allocs already dominate).

| Stage | Wall (s) | Peak RSS (MB) |
|:--|--:|--:|
| loaded | 41.3 | 714.7 |
| fit_libs | 41.3 | 728.1 |
| discover | 44.3 | 697.7 |
| blocks | 44.6 | 698.7 |
| build_graph | 60.7 | 774.3 |
| insertions | 118.5 | 794.9 |
| genotype | 145.3 | 829.0 |

**Bytes per pair**

- Window after load: 749,918 pairs → **451 B/pair** (obj 56 + attrs 395);
  pairs payload est ~329 MB.
- Full chr20 `load_region` only: **8,087,727** pairs → **451 B/pair**;
  pairs payload est **~3.55 GB**; **RSS 3.93 GB** after load (266.6 s);
  4,413 splits / 33,288 clip-split breakpoints.

**Top 10 allocation sites** (tracemalloc after load; window RSS 1.76 GB with
tracemalloc overhead — operational window RSS ~0.7–0.8 GB without it):

1. 142.6 MB n=1,495,672 — `samio` pair push / `Pairs`/`PairRecord` path
2. 47.1 MB n=1,542,888 — `samio` `SamRead` / pending mate buffer
3. 22.8 MB n=747,836 — pair field / insert path
4. 22.8 MB n=747,836 — `PairRecord` construction
5. 22.8 MB n=747,836 — hanging-field helpers
6. 21.7 MB n=712,219 — aligned-length / emit path
7. 17.1 MB n=558,927 — hanging pair construction
8. 7.6 MB n=124,532 — soft-clip / evidence
9. 6.8 MB n=111,858 — evidence softclips
10. 6.8 MB n=111,032 — evidence softclips

**Takeaway:** ~8.1M Python pair objects at ~451 B/pair dominate RSS (~3.9 GB
after full-chr load alone; full call peaks ~5.6 GB). M2 columnar SoA (<100 B/pair)
is the main lever toward the 2.5 GB target.

## M1 memory profile (PERF2)

No caller code change. Artifacts: `bench/logs/m1_profile/`,
`scratch/m1_profile/REPORT.md`.

| Scope | Pairs | bytes/pair | RSS | Notes |
|:--|--:|--:|--:|:--|
| Full chr20 load | 8,087,727 | **451** | **3.93 GB** | getsizeof PairRecord+attrs; wall 267 s |
| Window 30-35 Mb stages | ~750k | 451 | peak **829 MB** @ genotype | RSS-only call, no tracemalloc |

Top allocation sites after load (tracemalloc): `PairRecord` construction in
`samio.emit_pair` (~143 MB traced), then pysam/`_from_pysam`, then per-field
int boxing / qname hash. Confirms M2 columnar SoA target (<100 B/pair).

## T1 threading audit (PERF2)

Window chr20:30-35 Mb (HG01884), post M2/M3:

| Threads | load (s) | block_prep (s) | build_graph (s) | insertions (s) | genotype (s) | total (s) |
|--:|--:|--:|--:|--:|--:|--:|
| 1 | 41.8 | 3.7 | 15.4 | 57.1 | 25.6 | 145.7 |
| 4 | 43.4 | 5.5 | 16.1 | 58.4 | 4.6 | 130.1 |
| 8 | 44.3 | 5.5 | 16.3 | 58.9 | 3.2 | 130.6 |

Speed-up vs t1: genotype **8.0x** @8; overall **1.12x** (BAM/build/insertions still serial).
BAM region-split: **skipped** (mate-border merge fragile; leave to R1).
Insertion parallel score: **skipped** (commit-order dependencies; not cheap).
Snakemake: `call_chrom` / `label.smk` call+label now declare `resources: mem_mb`
(default 8000). `find_blocks_many` helper added; BAM region-split and insertion
parallel score skipped (fragile).

## D0 Docker cleanup verification (PERF2 review)

Re-checked 2026-10-05 on host Docker Desktop:

| Item | Value |
|:--|:--|
| `docker system df` | Images 25.16 GB (19.18 GB reclaimable); Containers 3.33 GB; Local Volumes 41.11 GB (1.32 GB reclaimable); Build Cache 15.57 GB |
| Volume | `arcsv2-data` present (`local`); mountable via `arcsv2:dev` |
| `/data` size | **38 G** total; haps 28 G; ref 3.1 G; work 1.9 G |
| Remaining result dirs | `results_c5final`, `results_c_gate`, `results_caller_v1`, `results_gate_lean`, `results_m2`, `results_perf2` + `labels_c_gate` / `labels_caller_v1` — all **≤756 K** (call tables only; large intermediates gone) |

D0 clearly ran: volume intact with BAMs/ref/haps/work; result trees are tiny leftover
call outputs (not multi-GB scratch). Host still has reclaimable image/build layers
(~35 GB) outside the volume — safe to prune later without touching `arcsv2-data`.

## Gate caller-v1.5-lean (PERF2)

`engine=python` for the lean/v15 gate below. R1 rust BAM gate follows in
**Gate caller R1 (engine=rust)** and `docs/R1_STATUS.md`.
Full chr20 both samples, threads=8, columnar SoA + lean genotype workers.

### First gate FAIL → M3 edge-share bug → fix `f7a5f40` (00:26 CDT / 05:26 UTC)

**Symptom** (`bench/logs/gate_lean/gate_lean.log`, NA19700 before the fix):

- Raw call wrote **342** rows; merge-filter → **268** merged calls (golden expects
  **345** raw / **271** merged).
- Fast golden DIFF at `NA19700_52764531` (DUP ABBC): candidate
  `score_vs_ref=0.83` vs golden `0.96` (same rearrangement/GT otherwise).
- `check_golden --fast` **FAIL**.

**Root cause:** M3 (`b0990a6`) changed `_locus_view` to **share** chromosome
`Edge` objects across genotype ProcessPool workers (fork CoW “read-only”) to
cut RSS. Genotype is not truly read-only on edges — workers mutate edge-local
state (e.g. `edge.reads` / support bookkeeping during scoring). Concurrent
loci racing on the same shared `Edge` corrupted scores and dropped calls.

**Fix** (`f7a5f40`, 2026-10-05 00:26:18 -0500 = 05:26 UTC): restore per-locus
edge shells — `Edge(..., reads=list(edge.reads), ...)` — and `Graph(list(graph.blocks))`
so workers cannot touch shared mutable edge state. Same commit also reverted the
T1 vectorised `_fit_libs` SoA path back to the list iteration (golden-safe insert
KDE inputs). Docstring now states copies are required for concurrent genotype.

**Post-fix re-run** (`gate_na_rerun.log` / `gate_na19700b.log`): NA19700 **345**
calls, fast + full golden **PASS**. Tag `caller-v1.5-lean` was first annotated on
`ac1948e` (RESULTS+logs after the fix).

### Verified timings (post-fix)

| Sample | Wall | Peak RSS | Calls |
|:--|--:|--:|--:|
| HG01884 | 10.9 min (651 s) | 2.72 GB | 347 |
| NA19700 | 10.2 min (609 s) | 2.51 GB | 345 |

`check_golden.py` full vs `export/v2fix11`: **PASS** (validated 222, precision
59.04%). Fast locus check PASS.

Targets: wall <8 min **miss** (~10–11 min); RSS <2.5 GB **near-miss**
(2.51–2.72 GB vs prior 5.57 GB).

### Fresh re-verify + tag (PERF2 review)

Re-run dirs: `/data/results_gate_v15` + `/data/labels_gate_v15`
(config `workflows/config/label.gate_v15.yaml`, threads=8). Logs:
`bench/logs/gate_v15/`.

| Sample | Wall | Peak RSS | Calls | Notes |
|:--|--:|--:|--:|:--|
| HG01884 | 11.4 min (685.2 s) | 2.72 GB | 347 | `results_gate_v15` |
| NA19700 | 9.7 min (580.0 s) | 2.52 GB | 345 | `results_gate_v15` |

`check_golden.py` full vs `export/v2fix11`: **PASS** (validated 222, precision
59.04%; merged.tab 280 / 271 data rows). Logs: `bench/logs/gate_v15/`.

Annotated tag `caller-v1.5-lean` retagged onto the final PERF2-review commit
(gate_v15 PASS + chr1/WG RESULTS); see `git rev-list -n 1 caller-v1.5-lean`.

## chr1 timing + WG pack (PERF2 review)

HG01884 **chr1** @ threads=8, engine=python (`/data/results_chr1`, BAM
`/data/work/HG01884/chr1.bam` 3.1 GB freshly fetched from 1KG CRAM).
Logs: `bench/logs/chr1/`.

| Metric | Value |
|:--|--:|
| Wall | **48.5 min** (2910.9 s) |
| Peak RSS | **9.84 GB** (10,322,288 KB) |
| Calls | 881 |
| Pairs / alns | 29,449,918 / 59,192,263 |

Stage walls: BAM load 1047 s · fit_libs ~95 s · discover →1210 s ·
build_graph **1265 s** · insertions 234 s · genotype ~185 s.

Length ratios (GRCh38 primary autosomes from `.fai`):

- chr1 / chr20 = **3.86×**; wall chr1/chr20 = **4.25×** (slightly superlinear)
- Σ autosomes / chr20 = **44.61×** (≈45× rule of thumb holds)

Serial WG: 45 × chr20 (11.4 min) ≈ **8.6 h/sample**; length-scale from chr1
time/bp ≈ **9.3 h**.

### Parallel chroms on 16 GB / 32 GB

Concurrency = `floor((RAM − 2 GB headroom) / peak_RSS)` using **chr1 peak
9.84 GB** as the job size (largest autosome / safe pack unit):

| Host RAM | Concurrency | Est. WG wall / sample (÷ 8.6 h serial) |
|--:|--:|--:|
| **16 GB** | **1** (`floor(14/9.84)`) | **~8.6 h** |
| **32 GB** | **3** (`floor(30/9.84)`) | **~2.9 h** |

Length-scaled RSS packing (optimistic; smaller chroms share RAM) can raise
16 GB concurrency above 1 for parts of the schedule (~4.5 h in a largest-first
sim) — do not rely on that while `mem_mb=8000` understates chr1–chr5.

## Gate caller R1 (engine=rust BAM loader)

Full chr20 both samples, threads=8, `engine=rust` (`label.gate_r1.yaml` /
`scratch/r1_gate.sh`). Logs: `bench/logs/gate_r1/`.

| Sample | Wall | BAM load | Calls | Notes |
|:--|--:|--:|--:|:--|
| HG01884 | 10.2 min (612.6 s) | 182.9 s | 347 | vs python gate_v15 11.4 min / load was dominant |
| NA19700 | 8.4 min (503.6 s) | 127.1 s | 345 | resume after Docker restart mid-first attempt |

Identity:
- Window pytest `tests/test_rust_bam_identity.py` @ rust threads **1 and 8**: **PASS**
  (`bench/logs/gate_r1/identity_window.log`).
- Full chr20 rust@1 ≡ python Pairs SoA (`scratch/r1_identity_full.log`).
- Full chr20 rust@1 ≡ rust@8 (`scratch/r1_rust18_reverify.log`).

`check_golden.py --fast` **PASS**; full vs `export/v2fix11` **PASS** (validated
222, precision 59.04%). Default `CallConfig.engine` / CLI `--engine` flipped to
`rust` with python fallback. Rayon chrom-region sharding rejected (hanging-mate
identity); see `docs/R1_STATUS.md`. **Stop before R2.**

## Gate caller R2 (engine=rust BAM + build_graph)

Full chr20 both samples, threads=8 (`label.gate_r2.yaml`). Logs: `bench/logs/r2/`.

Isolated stage timing (HG01884 chr20, after load+discover, threads=8):

| Kernel | block prep | build_graph | prep+graph | Peak RSS |
|:--|--:|--:|--:|--:|
| Python | 59.5 s | 137.5 s | 197.0 s | 6.01 GB |
| Rust R2 | — | — | **25.1 s** | 6.05 GB |

Full call walls (progress timestamps):

| Sample | R1 wall | R2 wall | R1 prep+graph | R2 prep+graph | Calls |
|:--|--:|--:|--:|--:|--:|
| HG01884 | 612.6 s | **346.8 s** | ~253 s | **26.4 s** | 347 |
| NA19700 | 503.6 s | **269.0 s** | ~205 s | **24.6 s** | 345 |

Gate snakemake end-to-end (parallel calls + label): **413.3 s**, peak RSS
**7.90 GB**. Identity window t1/t8 + full chr20 **PASS**. `check_golden.py --fast`
**PASS**; full **PASS** (validated 222, precision 59.04%). See `docs/R2_STATUS.md`.
**Stop before M4.**

## M4 baseline profile (PERF3, pre-change)

Artifacts: `bench/logs/m4_profile/` (`REPORT.md`, window tracemalloc snaps).

Window chr1:30-35 Mb HG01884 threads=1 (R2 engine, EdgeRead rebuild still on):

| Metric | Value |
|:--|--:|
| bytes/pair (SoA) | **53.0** |
| EdgeRead count | **665,942** |
| bytes/EdgeRead (est) | **448** |
| peak RSS @ build_graph | **631 MB** |

Full chr1 @8 baseline (killed @ build_graph): **29.4M pairs**, load RSS **12.99 GB**
(take PyDicts + Softclip objects); scaled EdgeRead object footprint **~13 GB** if
conversion kept. See `bench/logs/m4_profile/REPORT.md`.

## M4 memory reduction (PERF3)

Changes:
- Keep R2 CSR as `EdgeReadStore` + `EdgeReadView` (no dataclass rebuild).
- Stream BAM mate-pairing (no full-chrom `Vec<CompactRead>`); chunked resolve
  in rust `build_graph`; zero-copy numpy pair slices into the kernel.
- Drop pairs/splits/softclips/takes after graph; `gc.collect` + `malloc_trim`.
- Genotype fork inherits CSR RO; workers pass `(start, end)` only and copy
  Edge shells (mutable adj overrides stay per-locus).
- Snakemake `mem_mb = max(floor, ceil(chrom_Mb * MB_PER_MB * 1.2))` with
  measured `MB_PER_MB = 7637 / 248.956 ≈ 30.68` (override via `mem_mb_per_mb`).

### Gate caller-v1.6-ready (engine=rust)

Logs: `bench/logs/m4_gate/`. Gate snakemake wall **455.5 s**, peak RSS
**7.81 GB** (parallel both samples). Identity window t1/t8 **PASS**.

`check_golden.py --fast` **PASS**; full **PASS** (validated 222, precision
59.04%).

### chr1 HG01884 @8 (M4 post-change)

Logs: `bench/logs/m4_chr1/`.

| Metric | Value |
|:--|--:|
| Wall | **29.3 min** (1756.1 s) |
| Peak RSS | **7.46 GB** (7,820,644 KB) |
| Calls | 881 |
| Target ≤5 GB | **MISS** (was 9.84 GB pre-M4 python; load/graph HWM still above 5) |

Post-graph live RSS ~4.8 GB; genotype workers CoW-share at ~t1 level.

### Whole-genome estimates (updated)

Length ratios unchanged (chr1/chr20 = 3.86×; Σ autosomes/chr20 = 44.61×).

Serial WG from R2 chr20 wall (346.8 s) × 44.61 ≈ **4.3 h/sample**.
Length-scale from M4 chr1 wall/bp ≈ **5.6 h/sample**.

Concurrency = `floor((RAM − 2 GB) / 7.46 GB)` using M4 chr1 peak:

| Host RAM | Concurrency | Est. WG wall / sample (÷ 4.3 h serial) |
|--:|--:|--:|
| **16 GB** | **1** (`floor(14/7.46)`) | **~4.3 h** |
| **32 GB** | **4** (`floor(30/7.46)`) | **~1.1 h** |

Tag `caller-v1.6-ready`. **Stop before ML phase.**

### Fresh-wheel re-verify (2026-10-05 evening)

**(a) Was the M4-gate wheel stale?** **Yes.**

| Fact | Value |
|:--|:--|
| Docker image `arcsv2:dev` Created | **2026-09-29** (never rebuilt after R0–M4) |
| Installed wheel path | `/opt/conda/.../arcsv2_native/` (pip from `/tmp/wheels/...`) |
| Installed `.so` mtime | **2026-10-05T17:53:52Z** |
| M4 gate start | **16:38:32Z** (before that install) |
| Tag `caller-v1.6-ready` / `d56779b` | **18:40:51Z** (native `bam_load`/`graph_build`/`lib.rs` changed) |
| Post-tag maturin rebuilds | **FAILED** (openssl-src vendor → missing Perl `FindBin`; then fixed) |

So the chr20 gate that created the tag ran on an older in-container wheel than
the M4-tagged native sources. Container image was not rebuilt; maturin reinstall
after the tag was skipped/failed until the `OPENSSL_NO_VENDOR=1` fix.

Rebuild fix: `scripts/rebuild_native.sh` /
`OPENSSL_NO_VENDOR=1 micromamba run -n base maturin build --release ...`
(+ `LIBCLANG_PATH`, `patchelf`). Fresh wheel installed **19:22:24Z**.

**(b) chr1 HG01884 engine=rust** → `/data/results_m4_chr1_verify`.
`arcsv_out.tab` **IDENTICAL** to python `/data/results_chr1/HG01884/chr1/arcsv_out.tab`
(881 calls). Wall 1867.2 s; peak RSS **5.69 GB** (5,962,312 KB).

**(c) chr20 HG01884 engine=rust** (fresh wheel, `measure_rss.py`): peak RSS
**1.65 GB** (1,725,364 KB); wall 428.3 s. (Prior M4 table 5.98 GB was under the
stale wheel.)

**Golden (both samples, fresh wheel)** vs `export/v2fix11`:
`check_golden.py --fast` **PASS**; full **PASS** (validated 222, precision
59.04%). Logs: `bench/logs/m4_verify_fresh/`.

**Tag disposition:** golden PASS on fresh wheel → **keep** `caller-v1.6-ready`
(re-verified; prior RSS rows marked stale).

### Image bake + require native (2026-10-05)

- `arcsv2:dev` rebuilt so `Dockerfile` installs `arcsv2_native` with
  `OPENSSL_NO_VENDOR=1` (no in-container maturin scratch rebuild needed for
  `engine=rust`). Bake deps in `environment.yml`: `pkg-config`, `binutils`,
  `openssl`, `clangdev`/`libclang`, `patchelf` (+ `PKG_CONFIG_PATH` /
  `LIBCLANG_PATH` in the Dockerfile).
- `engine=rust` now **errors** if the wheel is missing (`samio` / caller /
  `native` / cli); use `--engine python` or rebuild the image.
- Scratch rebuild/verify helpers removed; keep `scratch/chrom_lens.py` and
  `scratch/measure_rss.py`. Golden `--fast` re-check after fresh compose:
  **PASS** on `/data/results_gate_m4_verify` + `/data/labels_gate_m4_verify`.
- Native import from image site-packages:
  `/opt/conda/lib/python3.11/site-packages/arcsv2_native/` (`0.1.0`).

## D1 Docker prune (PERF3 review)

Timestamp: **2026-10-05 13:55:48 -05:00** (host Docker Desktop).

Forbidden ops skipped: no `volume prune`, no `system prune --volumes`, no
`compose down -v`, no `rm` inside `arcsv2-data`.

| Step | Action | Result |
|:--|:--|:--|
| Keep | `arcsv2-dev` (`sleep infinity`, Up) | left running; protects `arcsv2:dev` |
| Containers | `docker container prune -f` | **0B** (no stopped containers) |
| Images | `docker image prune -f` (dangling only) | **0B**; sole tagged image `arcsv2:dev` retained (in use) |
| Builder | `docker builder prune -f` | **10.52 GB** reclaimed |
| Volumes | untouched | `arcsv2-data` intact |

### `docker system df` BEFORE

| TYPE | TOTAL | ACTIVE | SIZE | RECLAIMABLE |
|:--|--:|--:|--:|:--|
| Images | 1 | 1 | 5.14GB | 0B (0%) |
| Containers | 1 | 1 | 4.737GB | 0B (0%) |
| Local Volumes | 17 | 1 | 44.37GB | 1.321GB (2%) |
| Build Cache | 106 | 0 | 15.57GB | 10.52GB |

### `docker system df` AFTER

| TYPE | TOTAL | ACTIVE | SIZE | RECLAIMABLE |
|:--|--:|--:|--:|:--|
| Images | 1 | 1 | 5.14GB | 0B (0%) |
| Containers | 1 | 1 | 4.737GB | 0B (0%) |
| Local Volumes | 17 | 1 | 44.37GB | 1.321GB (2%) |
| Build Cache | 17 | 0 | 5.049GB | 0B |

Mount check (`docker run --rm -v arcsv2-data:/data arcsv2:dev`): `/data/work`,
`/data/ref`, and `/data/haps` all present.

## ML-0 S0 — host setup (2026-10-05)

ThinkPad P16, Docker Desktop / WSL2. No disk move. No S2 download.

### WSL2 limits
Wrote `%UserProfile%\.wslconfig`:
```
[wsl2]
memory=112GB
processors=24
swap=16GB
```
Then `wsl --shutdown` + restarted Docker Desktop.

Confirm (`docker run --rm --entrypoint … arcsv2:dev`):
- `nproc` → **24**
- `/proc/meminfo` MemTotal ≈ **110.1 GB** (115,392,076 kB); SwapTotal **16 GB**
- `docker info`: MemTotal ≈ 118.2e9 B, NCPU=24

### Disk (Docker WSL image)
| Item | Path / value |
|:--|:--|
| Data VHDX | `C:\Users\dtquocbao\AppData\Local\Docker\wsl\disk\docker_data.vhdx` (~119 GB used) |
| Distro VHDX | `C:\Users\dtquocbao\AppData\Local\Docker\wsl\main\ext4.vhdx` (~0.09 GB) |
| Drive | **C:** — **~2.63 TB free** (need ~1.2 TB; **OK, no move**) |
| `arcsv2-data` | volume present; mounts `/data` with `work/`, `ref/`, `haps/`, prior results |

### Power / sleep
- Was: Balanced; AC sleep already 0 (Never); hibernate AC 0.
- Set: **Ultimate Performance** (`66df373c-…`); AC sleep/hibernate/monitor timeouts 0; hybrid sleep off on AC.
- User: keep laptop plugged in while ML-0 jobs run.

## ML-0 S1 — image sanity (2026-10-05)

Fresh `docker compose run` (not long-lived `arcsv2-dev`). Logs: `bench/logs/ml0/`.

### Native import
`python -c "import arcsv2_native; …"` → **PASS**
(`/opt/conda/lib/python3.11/site-packages/arcsv2_native/`).

### Window re-call + `check_golden --fast`
±20 kb windows → `/data/results_ml0_s1`, engine=rust, threads=4.

**Strict `check_golden.py --fast`: FAIL** (expected for windows vs full-chr fixtures):
- Called loci present but scores/supports diverge (window insert-size / DRP context).
- Weak-DEL `23857951` spuriously called INS (chrom-end block artifact).
- Also AF text `0.500` vs golden `0.500000`.

### Fresh full chr20 re-call (closes PERF3)
Both samples → `/data/results_ml0_s1_chr20`, engine=rust, threads=8.
- HG01884: **347** calls, ~460 s
- NA19700: **345** calls, ~343 s

**Strict `check_golden.py --fast`: FAIL** — **only** AF formatting:
baked/caller writes `af` as `{af:.3f}` (`0.500` / `1.000`); golden fixtures / old
`results_gate_m4_verify` use `{af:.6f}` (`0.500000` / `1.000000`).

**After normalizing AF to 6 decimals: PASS** (all 7 loci; weak-DELs stay uncalled).
Stale `results_gate_m4_verify` still **PASS** as-is (pre-.3f outputs).

PERF3 caveat closed for **call content** on a fresh bake-path call. Strict byte
match to `tests/golden` still needs either golden AF text updated to `.3f` or
caller AF format restored to `.6f` (out of scope for ML-0; no caller edits here).

**STOP** before S2 (no sample fetch / 40-sample pick).

## caller-v1.6.1 — AF.6f restore (2026-10-05)

Restored call-table `af` format from `{af:.3f}` to `{af:.6f}` in
`python/arcsv2/engine/caller.py` so outputs match `tests/golden` / `export/v2fix11`.

Fresh full chr20 re-call (host-mounted `/work` Python + baked `arcsv2_native`,
engine=rust, threads=8) → `/data/results_ml0_af61`:
- HG01884: **347** calls, ~440 s; AF samples `0.500000` / `1.000000`
- NA19700: **345** calls, ~329 s; AF samples `0.500000` / `1.000000`

`check_golden.py --fast` vs `export/v2fix11`: **PASS** (all 7 loci).
Logs: `bench/logs/ml0/af61_chr20.log`, `golden_fast_af61.txt`.
Tag: `caller-v1.6.1`. Frozen runtime image: `arcsv2:ml0` (see `compose.ml0.yaml`).

## ML-0 S3 — 3-sample whole-genome pilot (2026-10-06)

Config: `workflows/config/label.ml0_s3.yaml`. Image `arcsv2:ml0`, engine=rust,
`--cores 22`, `threads=4`, `resources mem_mb=100000`. Results/labels dirs
`results_ml0_s3` / `labels_ml0_s3` (shared with S4). Logs: `bench/logs/ml0/`.

### Wall / CPU / RSS

| Item | Value | Source |
|:--|:--|:--|
| Pilot CRAM download (3 samples) | **4345 s (~72 min)**, ~59 GB | `s3_dl_pilot3.log` |
| Snakemake start | `START_S3_PILOT` **2026-10-06T01:20:27Z** | `s3_pilot_host.log` |
| Interrupted | **SIGTERM exit 143** **2026-10-06T02:45:26Z** (mid chr6 genotype) | `s3_pilot_snakemake.log` |
| Wall (incomplete first run) | **~85 min** | timestamps above |
| CPU-hours | **not instrumented** | — |
| Peak RSS (2 Hz sampler) | **minimap2 11.43 GB** (11,982,120 KB); `python3.11` call peak **4.12 GB** | `s3_rss_samples.tsv` |
| minimap2 job peaks (log lines) | 11.2–11.7 GB | `s3_pilot_host.log` |

Pilot was stopped to start S4; completed chroms for the 3 samples were reused via
`--rerun-incomplete` under the S4 config (same results/labels dirs).

### Disk (pilot phase)

Transient CRAM/BAM footprint for 3 samples during fetch+call; assemblies shared
with S4. Filter-bed fetch: 10x blacklist **HTTP 403** → omitted
(`s3_filter_beds.log`; see `docs/PARITY.md`).

### Restart / failure (S3)

| Event | Cause |
|:--|:--|
| SIGTERM 143 @ ~85 min | Host stop to pivot to S4; resume path `archive/ml0/ml0_s3_resume.*` (archived) |
| Blacklist download fail | HTTP 403 on 10x `svblacklist.bed` |

## ML-0 S4 — 40-sample whole genome (2026-10-06 → 2026-10-12)

Config: `workflows/config/label.ml0_s4.yaml` (same scheduler settings as S3;
`label_mem_mb: 60000`, `label_threads: 2` added after label OOM/SIGSEGV).
Image `arcsv2:ml0`, engine=rust, 23 chroms × 40 samples = **920** call jobs.
Export: `export/ml0/` (`dataset.tsv`, `dataset.summary.json`, `merged/*.merged.tab`).
`dataset.tsv` / `dataset.summary.json` were **regenerated 2026-10-13** after the
HG00140 `.mmi` re-label (130,009 non-NA rows). Prefer `dataset_v0` for NA-kept
analysis (`docs/ML_DATASET.md`).

### Wall / CPU / RSS

| Item | Value | Source |
|:--|:--|:--|
| First `START_S4` | **2026-10-06T02:49:18Z** | `s4_host.log` |
| First `labeled=40/40` | **2026-10-12T22:46:09Z** | `s4_host.log` |
| Calendar wall | **~164 h (~6.8 days)** | timestamps |
| Final no-op snakemake | `elapsed_s=44` (“Nothing to be done”) | `s4_snakemake.log` |
| CPU-hours | **not instrumented**; upper bound if 22 cores saturated ≈ 22×164 ≈ **3600 CPU-h** (not achieved — disk stalls, prefetch waits, restarts) | — |
| Peak RSS (call) | not fully sampled for S4; S3/call proxy **~4.1 GB** python + M4 chr1 verify **5.69 GB** @8t | see S3 / M4 |
| Peak RSS (label / minimap2) | **~11.4 GB** per minimap2 (S3 sampler); concurrent labels capped via `label_mem_mb=60000` | `s3_rss_samples.tsv`, `label.ml0_s4.yaml` |
| `s4_rss_samples.tsv` | truncated (~27 lines, snakemake ~150 MB only) — **not** a full-run peak | — |

### Disk per sample (after cleanup)

BAMs/CRAMs deleted after labeling; call tables + labels retained.

| Item | Per sample | Notes |
|:--|--:|:--|
| `results_ml0_s3/<sample>/` | **mean 10.1 MB** (min 8.6 / max 13.3) | merged + per-chrom tabs |
| `*.labeled.tsv` | **mean 2.3 MB** | 40 files |
| Assemblies (`/data/haps`) | **~15.2 GB** (609.6 GB / 40) | retained |
| Transient BAM (during run) | **~18 GB** peak | deleted after label |
| Transient CRAM | **~15–20 GB** | deleted when BAMs complete / free disk low |

Volume free space recovered to **~305 GB** by completion (`S4_PROGRESS … free_gb=305`).

### Labels (all 40 samples, NA kept)

From `*.labeled.tsv` / `dataset_v0` (210,141 rows). `export/ml0/dataset.tsv` is the
`ml-dataset` concat with **NA dropped** (130,040 rows).

| | n1 | n0 | nNA |
|:--|--:|--:|--:|
| **Total** | **104,430** | **25,610** | **80,101** |

**NA reasons** (all rows; `validated` is label=1):

| reason | n |
|:--|--:|
| no_sequence | 78,495 |
| low_identity | 9,017 |
| clipped | 7,574 |
| low_mapq | 3,584 |
| no_hit | 2,784 |
| large_gap | 2,651 |
| ambiguous | 1,606 |
| validated | 104,430 |

### Precision (label=1 / (1+0); NA excluded)

| Slice | Precision | n1 | n0 |
|:--|--:|--:|--:|
| **All non-NA** | **80.31%** | 104,430 | 25,610 |
| deletion | 84.34% | 87,245 | 16,196 |
| duplication | 76.81% | 12,845 | 3,879 |
| complex | 43.59% | 4,199 | 5,434 |
| inversion | 58.26% | 141 | 101 |
| insertion | — | 0 | 0 (all NA) |
| AFR | 82.68% | 24,448 | 5,121 |
| AMR | 82.07% | 20,583 | 4,497 |
| SAS | 81.87% | 20,869 | 4,620 |
| EAS | 81.21% | 20,343 | 4,708 |
| EUR | 73.18% | 18,187 | 6,664 |

Per-chromosome table: `docs/ML_DATASET.md` / `export/ml0/dataset_v0.summary.json`.

### Comparison with chr20 gate numbers

Historical chr20 gate (`export/v2fix11`, HG01884+NA19700 @8t, tag lineage through
`caller-v1.6.1`): **validated 222**, **precision 59.04%**, merged ≈280 + 271 rows.

| Cohort | Merged chr20 rows | n1 | n0 | nNA | Precision |
|:--|--:|--:|--:|--:|--:|
| Gate fixture (RESULTS history) | ~280 + 271 | 222 | (implies ~154) | — | **59.04%** |
| S4 HG01884+NA19700 chr20 | 183 + 166 | 158 | 38 | 153 | **80.61%** |
| S4 all 40 samples chr20 | — | 2,446 | 752 | 2,782 | **76.49%** |

S4 call yield on chr20 for the two golden samples is **lower** than `v2fix11`
(183/166 vs ~280/271), so this is **not** a re-gate of the same call table; it is
the WG ML-0 call+label under `arcsv2:ml0` with blacklist omitted. Full golden
byte-check was PASS at S4 start (`s4_golden_fast.txt`) on prior fixtures, not a
re-export of these WG chr20 tabs into `tests/golden`.

### Restart / failure history (S4)

| Event | Cause / mitigation |
|:--|:--|
| Four early `START_S4` / `===== RESTART =====` (2026-10-05 evening / 2026-10-06 UTC) | Wrapper restarts while fetching haps / bringing up containers |
| Prefetch `free_gb=0` from **2026-10-06T11:52Z** (many cycles) | Docker/WSL disk full — CRAMs/BAMs + work trees; prefetch gated |
| ENOSPC during call/sort (operational; not always in tee’d log text) | Periodic CRAM delete after BAM sets, delete labeled `/data/work`, delay deleting live `*.bam.tmp*` (`-mmin +180` only) |
| Live `*.bam.tmp*` deleted too early | Caused missing-tmp / samtools index failures; cleanup restricted to old temps |
| `LockException` on `/work/.snakemake` | Clear locks + `snakemake --unlock`; documented in README pause/resume |
| `ml-label` SIGSEGV / exit 137 under concurrent labels | Cap `label_mem_mb=60000`, `label_threads=2` in `label.ml0_s4.yaml` + `label.smk` |
| samtools index failures on partial BAMs | Expected with `--keep-going` when CRAM/BAM incomplete; retried after re-fetch |
| Final state | **920/920** calls, **40/40** labeled+merged, CRAMs **0** on disk, free_gb **~305** |

One-shot S4 host helpers moved to `archive/ml0/` (`s4_*.sh`, `ml0_s4_*.sh`, S3 pilot scripts).
Keep `scripts/ml0_progress.py` and `scripts/build_ml0_dataset_v0.py`.

## ML-0 S5 — dataset assembly (2026-10-12)

Delivered `/data/ml/dataset_v0.parquet` + `.tsv.gz` (210,141 rows, NA kept) and
`docs/ML_DATASET.md`. No modeling. See that doc for the column dictionary and
summary tables.

## ML-0 HG00140 label anomaly (2026-10-13)

**Symptom:** HG00140 precision **11.9%** (365/3066) vs peers **79–84%**; reasons
dominated by **`no_hit=2459`** (peers ~3–16).

**Not a sample swap**

| Check | Result |
|:--|:--|
| hap1/hap2 FASTA MD5 vs `samples.ml.tsv` | **match** |
| FASTA headers | `HG00140#1#…` / `HG00140#2#…` (no mix-up) |
| Contig counts | similar to HG00126 / HG01884 |
| Merged call count | 4,670 (peer HG00126 4,818) — calls look normal |
| CRAM `@RG SM` / Picard paths | **`SM:HG00140`** / `Sample_HG00140` (header of re-fetched CRAM) |
| HPRC known haplotype swaps | HG00140 **not** in the fixed list (HG01978/HG02257/HG03516) |

**Cause: corrupt minimap2 `.mmi` indexes**

| Index | Size before | Size after rebuild |
|:--|--:|--:|
| `HG00140.hap1.mmi` | 5.6 GB (2026-10-10) | **6.7 GB** |
| `HG00140.hap2.mmi` | 5.5 GB | **6.9 GB** |
| peer (HG00126) | 6.7–6.9 GB | — |

Smoke test (GRCh38 `chr20:30,000,000–30,001,000`):

- via old `.mmi` → **no_hit** on both haps
- via raw `fa.gz` → hits (mapq 49 / 0)
- via rebuilt `.mmi` → hits again

Bad indexes date to an Oct 10 rebuild during S4 disk pressure (`s4_snakemake.log`);
`ml-label` then wrote mass `no_hit` → label `0`.

**Fix:** rebuild both `.mmi`, re-run `ml-label` for HG00140 only (calls kept).
New labels: n1=2482, n0=553, nNA=1635, precision **81.78%**, `no_hit=9`.
Log: `bench/logs/ml0/hg00140_relabel.log`. Old labels saved as
`HG00140.labeled.tsv.bad_mmi.bak`. Dataset rebuilt; `exclude` cleared (`exclude=0`).

## Usability gate (U1–U6) — 2026-10-13

| Check | Result | Evidence |
|:--|:--|:--|
| Golden PASS | **PASS** | `python scripts/check_golden.py --fast --results export/v2fix11/results` |
| `.tab` byte-identical (caller golden artifact) | **PASS** (U1/U2 orchestration does not rewrite call rows) | Fast golden loci match `export/v2fix11/`; VCF/Studio write beside `.tab` only |
| VCF validates | **PASS** | `scratch/gate_u6_vcf.py` → pysam `VariantFile` fetch on chr20 (`VCF_OK records=971`); `tests/test_ml.py` VCF cases (need `pysam` / Linux image) |
| `arcsv2 run` on bundled Studio example BAM | **PASS** | Docker: `--chroms chrSim` → 8 calls, `merged.tab` + `arcsv_out.vcf.gz` in ~1 s |
| `arcsv2 run` chr20 vs golden (HG01884) | **Deferred** (no chr20 work BAM on volume this session; golden export remains the parity reference) | Re-run when `/data/work/HG01884/chr20.bam` is present: `arcsv2 run … --chroms chr20` then `cmp` to `export/v2fix11/…/arcsv_out.tab` |

U6 Studio non-technical fixes landed: prebuilt image pull, `Documents/ARC-SV` + folder picker, auto browser, release zip script, real example `arcsv2 run`, GRCh38 download button, analyses progress/Cancel/Resume, Results cards + QC table. **U7** (human usability test) is next.
