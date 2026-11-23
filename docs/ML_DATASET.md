# ML-0 dataset v0

Frozen caller: image `arcsv2:ml0` / tag `caller-v1.6.1` (engine=rust, AF six decimals).
Labels: `arcsv2-label-1`, feature schema `arcsv2-f1`.
Source: 40 HPRC×1KGP samples, whole genome (chr1–22 + chrX), after merge-filter + `ml-label`.

**Artifacts**

| Path | Description |
|:--|:--|
| `/data/ml/dataset_v0.parquet` | Primary table (210,141 rows) |
| `/data/ml/dataset_v0.tsv.gz` | Same rows, gzip TSV |
| `/data/ml/dataset_v0.summary.json` | Counts + precision slices |
| `export/ml0/dataset_v0.*` | Host export copy (primary ML table; NA kept) |
| `export/ml0/dataset.tsv` | `ml-dataset` concat (**NA dropped**; regenerated after HG00140 fix) |
| `export/ml0/dataset.summary.json` | Summary for `dataset.tsv` (same regen) |
| `export/ml0/merged/*.merged.tab` | 40 per-sample merged call tables |

Prefer **`dataset_v0`** for analysis that needs NA rows / context columns. Use
`dataset.tsv` only when you want the classic `ml-dataset` training concat.

Large files (`dataset_v0.parquet`, `dataset_v0.tsv.gz`, `dataset.tsv`,
`merged/`) are **not** kept in git; download the `benchmark-2026-11` release
(see `docs/DATA.md`). `export/ml0/dataset_v0.summary.json` stays in git.

Rebuild: `python scripts/build_ml0_dataset_v0.py` inside `arcsv2:ml0` with `/data` mounted.
Unlike `arcsv2 ml-dataset`, this assembly **keeps NA labels**.

## Dr. Bo's XGBoost / `model_replication` (not in git)

Unpublished lab material (pickled models, HG002 demo, preprocessed CSVs). Do
**not** commit it. Place a local copy at:

```
third_party/drbo_model/   # gitignored
```

Then run `scripts/analysis/eval_drbo_xgb.py` with explicit `--model`, `--meta`,
`--features`, `--raw-input`, and `--predictions` paths. Obtain the package from
Dr. Bo; it is not redistributed with ARC-SV2.

## Column dictionary

### Identity

| Column | Type | Source |
|:--|:--|:--|
| `sample` | str | labeled table / sample id |
| `superpopulation` | str | `workflows/config/samples.ml.tsv` |
| `population` | str | samples.ml.tsv |
| `sex` | str | samples.ml.tsv |
| `family_id` | str | samples.ml.tsv |
| `chrom` | str | call `chrom` |
| `start` | int | call `minbp` |
| `end` | int | call `maxbp` |
| `call_id` | str | call `id` |

### Label

| Column | Type | Source |
|:--|:--|:--|
| `label` | `1` / `0` / `NA` | `ml-label` |
| `reason` | str | NA / failure / `validated` reason |
| `gt_concordant` | str/bool | labeler |
| `ref_validated` | str/bool | labeler |
| `identity` | float | assembly alignment identity |
| `mapq` | float | labeler MAPQ |
| `max_indel` | float | labeler max gap / indel |
| `strand` | str | hap hit strand |
| `label_bp_shift` | float | breakpoint shift |
| `best_hap` | str | which haplotype hit |

### Call (raw caller columns)

Kept as emitted in `merged.tab` / labeled table, including:
`svtype`, `complextype`, `sv_class`, `rearrangement`, `gt`, `af`, `score_vs_ref`,
`score_vs_next`, `sr_support`, `pe_support`, `len_affected`, `inslen`,
`bp_uncertainty`, `bp_point_source`, `num_sv`, `num_paths`, `blocks`,
`blocks_point`, `filter`, and the other caller fields present in the labeled TSV.

Derived:

| Column | Type | Source |
|:--|:--|:--|
| `length_is_lower_bound` | bool | `>` in `inslen` (`inslen_is_lower_bound`) |

### Context

| Column | Type | Source |
|:--|:--|:--|
| `frac_simple_repeat` | float | overlap / span vs `/data/ref/simple_repeats.bed` |
| `frac_segdup` | float | `/data/ref/segdups.bed` |
| `frac_blacklist` | float | **always 0** — 10x blacklist omitted (HTTP 403); see `docs/PARITY.md` |
| `frac_gap` | float | `/data/ref/gaps.bed` |
| `gc_flank_500` | float | GC fraction of reference `[start−500, end+500]` |
| `coverage_ratio` | float | **NA** — BAMs deleted after S4 labeling |

## Summary tables

Precision = `n(label=1) / (n(1)+n(0))` (NA excluded). Insertions are all NA (`no_sequence`).

### Label totals

| label | n |
|:--|--:|
| 1 (validated) | 106,547 |
| 0 | 23,462 |
| NA | 80,132 |
| **total** | **210,141** |

Overall non-NA precision: **81.95%** (106,547 / 130,009).  
(After HG00140 `.mmi` rebuild + re-label; see `bench/RESULTS.md`.)

### NA / reason counts

| reason | n |
|:--|--:|
| validated | 106,547 |
| no_sequence | 78,495 |
| low_identity | 9,157 |
| clipped | 7,605 |
| low_mapq | 3,662 |
| large_gap | 2,704 |
| ambiguous | 1,637 |
| no_hit | 334 |

(`no_hit` was 2,784 before the HG00140 `.mmi` fix.)

### By SV class

| sv_class | n1 | n0 | nNA | precision |
|:--|--:|--:|--:|--:|
| deletion | 87,245 | 16,196 | 901 | 84.34% |
| duplication | 12,845 | 3,879 | 195 | 76.81% |
| complex | 4,199 | 5,434 | 407 | 43.59% |
| inversion | 141 | 101 | 103 | 58.26% |
| insertion | 0 | 0 | 78,495 | — |

### By superpopulation

| superpopulation | n1 | n0 | nNA | precision |
|:--|--:|--:|--:|--:|
| AFR | 24,448 | 5,121 | 17,317 | 82.68% |
| AMR | 20,583 | 4,497 | 16,867 | 82.07% |
| SAS | 20,869 | 4,620 | 16,514 | 81.87% |
| EAS | 20,343 | 4,708 | 16,063 | 81.21% |
| EUR | 18,187 | 6,664 | 13,340 | 73.18% |

### By chromosome

| chrom | n1 | n0 | nNA | precision |
|:--|--:|--:|--:|--:|
| chr1 | 7,946 | 1,751 | 7,937 | 81.94% |
| chr2 | 8,802 | 1,384 | 10,629 | 86.41% |
| chr3 | 6,858 | 1,407 | 4,582 | 82.98% |
| chr4 | 8,772 | 1,431 | 4,032 | 85.97% |
| chr5 | 6,784 | 1,166 | 3,693 | 85.33% |
| chr6 | 7,174 | 1,801 | 4,911 | 79.93% |
| chr7 | 6,668 | 1,891 | 4,608 | 77.91% |
| chr8 | 5,707 | 1,062 | 3,317 | 84.31% |
| chr9 | 4,800 | 1,081 | 3,752 | 81.62% |
| chr10 | 4,815 | 1,002 | 3,349 | 82.77% |
| chr11 | 4,531 | 1,429 | 4,084 | 76.02% |
| chr12 | 4,546 | 1,462 | 3,436 | 75.67% |
| chr13 | 3,844 | 810 | 2,552 | 82.60% |
| chr14 | 3,364 | 847 | 1,806 | 79.89% |
| chr15 | 2,514 | 574 | 1,949 | 81.41% |
| chr16 | 2,246 | 907 | 2,164 | 71.23% |
| chr17 | 2,327 | 1,205 | 2,515 | 65.88% |
| chr18 | 3,302 | 598 | 1,593 | 84.67% |
| chr19 | 1,718 | 1,008 | 1,654 | 63.02% |
| chr20 | 2,446 | 752 | 2,782 | 76.49% |
| chr21 | 1,394 | 491 | 1,224 | 73.95% |
| chr22 | 828 | 428 | 948 | 65.92% |
| chrX | 3,044 | 1,123 | 2,584 | 73.05% |

## ml-train wall time

`arcsv2 ml-train --labels export/ml0/dataset.tsv` (130,009 labeled rows, 40 samples):

| Run | Code | Wall | Peak RSS |
|:--|:--|--:|--:|
| Grace job 20653724 (1 core) | previous trainer | 38 min 02 s | 2.37 GB |
| Laptop, this tree | featurize-once trainer | 150.5 s | 2.25 GiB |

Job 20653417 hit the 30 min limit before this path existed. The laptop number is not a Grace measurement; this tree has not been rerun in the SIF.

## Notes

- Blacklist omitted (HTTP 403) — see `docs/PARITY.md`.
- Insertions remain unlabeled (`no_sequence`) until ML-1 flank labeling.
- `coverage_ratio` left NA; BAMs were deleted after labeling to free disk.
- Runtime / restart history: `bench/RESULTS.md` (ML-0 S3 / S4).
- **`exclude` column:** `1` = drop from ML training (currently all `0`). HG00140 had
  been `exclude=1` briefly: precision 11.9% from **corrupt `.mmi` indexes** (assemblies
  MD5/headers OK; raw `fa.gz` mapped, `.mmi` returned mass `no_hit`). Indexes rebuilt
  and the sample re-labeled (precision **81.8%**, `no_hit` 9). See `bench/RESULTS.md`.
