# Head-to-head: original ARC-SV vs ARC-SV2 (chr20, HG01884 + NA19700)

Goal: run the original ARC-SV on the exact chr20 BAMs ARC-SV2 used, push its calls through the
same merge-filter and the same haplotype labeler, and compare precision, call counts and runtime.
Merge-filter: `--min-size 50` plus optional paper BEDs from `label.smk` config
(`simple_repeats`, `segdups`, `blacklist`, `gaps` — fetch with `scripts/fetch_filter_beds.py`).

Original source: `..\..\src\arcsv` (fork of SUwonglab/arcsv). Only `arcsv/cli.py` and `setup.py`
differ from upstream `dd6f795`; the algorithm is unchanged. Working-tree changes are line endings only.

## 1. Start a long-lived container with the original mounted (PowerShell, arcsv2 repo root)

```powershell
# Set $ORIG to your local checkout of the original ARC-SV (SUwonglab/arcsv fork).
$ORIG = "$env:ORIG_ARCSV"   # e.g. path\to\Project-08-Psychiatric\src\arcsv
docker compose run -d --name arcsv-orig -v "${ORIG}:/orig:ro" --entrypoint sleep dev infinity
docker exec -it arcsv-orig bash
```

Everything below runs inside that shell. Results live in the `arcsv2-data` volume under `/data/orig/`.

## 2. Install the original in its own env (once)

```bash
micromamba create -y -p /data/envs/arcsv-orig -c conda-forge -c bioconda \
  python=3.10 pysam "numpy<2" scipy "scikit-learn<1.6" matplotlib python-igraph pyinter time
rm -rf /data/orig-src && cp -r /orig /data/orig-src
/data/envs/arcsv-orig/bin/pip install --no-deps /data/orig-src
export OA=/data/envs/arcsv-orig/bin
$OA/arcsv call -h | head -5
```

## 3. Smoke test (a few minutes)

```bash
cd /data/orig-src
$OA/arcsv call -i example/input.bam -r 20:0-250000 -R example/reference.fa -G example/gaps.bed -o /data/orig/example --overwrite
diff <(cut -f1-24 /data/orig/example/arcsv_out.tab) example/expected_output.tab && echo EXAMPLE OK
```

## 4. Run the original on chr20, one sample at a time (timed)

```bash
mkdir -p /data/orig
for s in HG01884 NA19700; do
  $OA/time -v -o /data/orig/$s.time.txt \
    $OA/arcsv call -i /data/work/$s/chr20.bam -r chr20 \
      -R /data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -G /data/ref/gaps.bed \
      -o /data/orig/$s/chr20 --overwrite > /data/orig/$s.call.log 2>&1
done
```

This can take hours. To leave it running, start it with `nohup bash -c '...' &` or just keep the
container up (it survives closing the terminal; re-enter with `docker exec -it arcsv-orig bash`).
Check progress: `tail -f /data/orig/HG01884.call.log`.
Wall time and peak RSS: `grep -E "Elapsed|Maximum resident" /data/orig/*.time.txt`.

## 5. Same merge-filter and labeler as ARC-SV2 (one sample at a time)

Fetch paper BEDs once (if not already under `/data/ref/`):

```bash
python scripts/fetch_filter_beds.py --dest /data/ref
```

Use the same BED flags as `label.smk` (empty config keys = min-size only).

```bash
cd /work
BEDS="--simple-repeats /data/ref/simple_repeats.bed --segdups /data/ref/segdups.bed --blacklist /data/ref/blacklist.bed --gaps /data/ref/gaps.bed"
for s in HG01884 NA19700; do
  python scripts/orig_to_arcsv2.py /data/orig/$s/chr20/arcsv_out.tab /data/orig/$s/chr20/arcsv_out.v2cols.tab
  arcsv2 merge-filter /data/orig/$s/chr20/arcsv_out.v2cols.tab $BEDS -o /data/orig/$s/merged.tab
  arcsv2 ml-label --calls /data/orig/$s/merged.tab \
    --reference /data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --hap1 /data/haps/$s.hap1.fa.gz --hap2 /data/haps/$s.hap2.fa.gz \
    --hap1-index /data/haps/$s.hap1.mmi --hap2-index /data/haps/$s.hap2.mmi \
    --unreliable /data/ref/unreliable.bed --threads 1 --sequential \
    --sample $s -o /data/orig/labels/$s.labeled.tsv
done
arcsv2 ml-dataset /data/orig/labels/HG01884.labeled.tsv /data/orig/labels/NA19700.labeled.tsv \
  -o /data/orig/labels/dataset.tsv
```

With-vs-without paper filter (no full re-call): omit `$BEDS` for min-size-only, or set
`simple_repeats`/`segdups`/`blacklist` to `""` in the label config and re-run `merge`+`label`
for both ARC-SV2 (`results_*`) and original (`/data/orig/`).

`scripts/orig_to_arcsv2.py` adds `blocks`, `blocks_point` and `bp_point_source` from the original
`bp` / `bp_uncertainty` columns so the labeler can build allele sequences. Uncertain junctions get
a +/- u/2 window, so the labeler's breakpoint-shift search treats both callers alike.

## 6. Export for review

```bash
mkdir -p /work/export/orig && cp -r /data/orig/labels /data/orig/*.time.txt /data/orig/*.call.log /work/export/orig/
for s in HG01884 NA19700; do mkdir -p /work/export/orig/$s && cp /data/orig/$s/merged.tab /data/orig/$s/chr20/arcsv_out.tab /work/export/orig/$s/; done
exit
docker rm -f arcsv-orig
```

## What gets compared

- Calls per sample and class (raw and after merge-filter)
- Precision per class, by size and by score (same labeler, same haplotypes)
- Shared vs caller-specific calls, and precision of each group
- Wall time and peak RSS per sample (ARC-SV2 reference: HG01884 44.4 min / ~5.1 GB; NA19700 32 min; 1 thread)

## 7. Benchmark results (2026-11, chr20 HG01884 + NA19700)

Full write-up: [orig vs ARC-SV2 report](https://claude.ai/artifact/9tCLAfkKAsnH4ShswQpuZW).
Reproduce tables with `scripts/analysis/compare_orig_vs_new_chr20.py` and
`scripts/analysis/eval_drbo_xgb.py` (see `scripts/analysis/README.md`).

### Before ML filter (same merge-filter + haplotype labeler)

| Caller | Validated / labeled | Precision | Notes |
|:--|--:|--:|:--|
| Original ARC-SV (`export/orig`) | 243 / 449 | **54.1%** | |
| ARC-SV2 (`export/v2fix11`) | 222 / 376 | **59.0%** | golden freeze |

- Overall Fisher exact *p* ≈ **0.16** (not significant).
- Duplication class *p* ≈ **0.011** (ARC-SV2 higher precision on DUPs).
- **221 / 243** original true positives matched by ARC-SV2 (±100 bp, same sample/class).

### After Dr. Bo's XGBoost filter

| Caller | Kept true / kept calls | Precision | Fisher *p* |
|:--|--:|--:|--:|
| Original | 193 / 215 | **89.8%** | |
| ARC-SV2 | 165 / 201 | **82.1%** | ≈ **0.033** |

XGBoost evaluation needs Dr. Bo's unpublished `model_replication` package under
`third_party/drbo_model/` (gitignored; see `docs/ML_DATASET.md`).
