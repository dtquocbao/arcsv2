# R1 notes (PERF2)

## Crate choice: rust-htslib

Chose ``rust-htslib`` over noodles:
- Links the same htslib stack pysam uses (CRAM + multithreaded BGZF).
- Project BAMs in ``arcsv2:dev`` already validated through pysam; conda
  ``libhts`` / ``htslib.pc`` are present in the image.
- noodles is pure Rust but would be a second CRAM/BAM stack to validate.

Build notes: ``OPENSSL_NO_VENDOR=1``, system ``libssl``, ``LIBCLANG_PATH`` for
bindgen, ``patchelf`` for the maturin wheel. Do **not** use the openssl-src
vendor path (needs Perl ``FindBin``). Rebuild helper:

```bash
OPENSSL_NO_VENDOR=1 micromamba run -n base \
  env OPENSSL_NO_VENDOR=1 LIBCLANG_PATH=/opt/conda/lib \
  maturin build --release -m native/arcsv2-native/Cargo.toml -o target/wheels
micromamba run -n base pip install --force-reinstall --no-deps target/wheels/arcsv2_native-*.whl
```

Or: ``bash scripts/rebuild_native.sh`` (installs clang/libclang if missing).

## What shipped

- ``arcsv2_native.load_bam_pairs`` → M2 Pairs SoA numpy arrays + take-records.
- ``samio.load_region`` with ``engine=rust`` calls the native loader; softclip /
  split discovery reuses ``evidence.py`` on the same take stream (order-matched).
- ``engine=rust`` requires ``arcsv2_native`` (no silent fallback). Use
  ``--engine python`` for the pysam path, or rebuild the image so the wheel is
  baked in.
- Config: ``CallConfig.engine`` / ``arcsv2 call --engine python|rust``.
- **Default flipped to ``rust`` after gate_r1 golden PASS** (python remains
  explicit ``--engine python``).
- Identity test: ``tests/test_rust_bam_identity.py`` (chr20 1 Mb window @ rust
  threads 1 and 8). Re-verified 2026-10-05: **2 passed** (`identity_window.log`).
- Full chr20: rust@1 arrays identical to python (8,087,727 pairs;
  splits/clips/extras counts match — `scratch/r1_identity_full.log`).
- Full chr20: rust@1 ≡ rust@8 bit-identical after serial+BGZF fix
  (`scratch/r1_rust18_reverify.log`).

## Parallelism / rayon

- ``threads=N`` sets htslib BGZF decompression threads (same role as pysam
  ``threads=``). Mate pairing is a single ordered pass so arrays stay
  bit-identical at threads 1 and 8.
- Rayon **region-sharding** was tried and **rejected for the default path**:
  chrom fetch returns ~20k unmapped records interleaved with mates; start-bucket
  shards pull those into the first shard and change hanging-pair fields
  (same ``n_pairs``, different multiset). Documented in ``bam_load.rs``.
  (Literal CURSOR_TASK_PERF2 "rayon over chromosome regions" acceptance is
  therefore a justified deviation: identity > region parallelism.)
- Rayon is used for ordered column packing when building the numpy return
  (indexed ``par_iter``; row order preserved).

## Gate (engine=rust) — COMPLETE

- Config: ``workflows/config/label.gate_r1.yaml`` (``engine: rust``, threads=8).
- Outputs: ``/data/results_gate_r1``, ``/data/labels_gate_r1``.
- Logs: ``bench/logs/gate_r1/`` (``gate_r1.log``, ``golden_fast.txt``,
  ``golden_full.txt``, ``identity_window.log``).
- HG01884: **347** calls / 280 merged; wall **612.6 s** (10.2 min); BAM load
  **182.9 s**.
- NA19700: **345** calls / 271 merged; wall **503.6 s** (8.4 min); BAM load
  **127.1 s** (complete resume after Docker Desktop restart mid-first-attempt).
- ``check_golden.py --fast`` **PASS**; full **PASS** (validated 222, precision
  59.04%).
- Stop after R1; do not start R2.
