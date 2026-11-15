# R2 notes (PERF3)

## What shipped

- ``arcsv2_native.build_graph``: block prep (clip stamp from breakpoints +
  softclips) and ``build_graph`` → CSR columnar edges + flat EdgeRead arrays.
- Python rebuilds ``Graph`` / ``Edge`` / ``EdgeRead`` from those arrays
  (``engine/graph_native.py``) so downstream stages are unchanged until M4.
- Wired when ``CallConfig.engine == "rust"``; missing native falls back to the
  Python stamp + ``build_graph`` path. ``--engine python`` keeps full Python.
- Identity: ``tests/test_rust_graph_identity.py`` (1 Mb window @ threads 1 and 8;
  full chr20 via ``ARCSV2_R2_FULL=1``).

## Parallelism

- Indexed rayon ``par_iter`` resolves pair ends (row order preserved).
- ``used_splits`` claims stay serial in pair order.
- Edge insertion + support counts merge serially in pair order (deterministic).

## Gate (engine=rust BAM + graph) — COMPLETE

- Config: ``workflows/config/label.gate_r2.yaml``
- Outputs: ``/data/results_gate_r2``, ``/data/labels_gate_r2``
- Logs: ``bench/logs/r2/``

### Identity

| Check | Result |
|:--|:--|
| Window t1 / t8 | **PASS** (2 passed, 7.4 s) |
| Full chr20 | **PASS** (1 passed, 429.5 s) |

### Stage timing (HG01884 chr20, isolated load→graph, threads=8)

| Path | block prep (s) | build_graph (s) | prep+graph (s) | Peak RSS |
|:--|--:|--:|--:|--:|
| Python (R1 baseline kernel) | 59.5 | 137.5 | 197.0 | 6.01 GB |
| Rust R2 | (inside kernel) | (inside kernel) | **25.1** | 6.05 GB |

Same edge count (92,846). ~7.9× on the prep+graph stage.

### Full call wall (progress timestamps)

| Sample | R1 wall (s) | R2 wall (s) | R1 prep+graph (s) | R2 prep+graph (s) | Calls |
|:--|--:|--:|--:|--:|--:|
| HG01884 | 612.6 | **346.8** | ~253 | **26.4** | 347 |
| NA19700 | 503.6 | **269.0** | ~205 | **24.6** | 345 |

Gate snakemake (both samples + label, parallel calls): wall **413.3 s**, peak RSS
**7.90 GB** (``measure_rss.py`` tree peak).

### Golden

- ``check_golden.py --fast`` **PASS**
- full vs ``export/v2fix11`` **PASS** (validated 222, precision 59.04%)

Stop after R2; do not start M4.
