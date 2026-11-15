# Data assets (benchmark-2026-11)

Small evidence stays in git. Large tables and full logs ship as a GitHub Release
named **`benchmark-2026-11`**, with `SHA256SUMS` beside the archives.

Tag (code tip): https://github.com/dtquocbao/arcsv2/releases/tag/benchmark-2026-11  
Maintainer upload (after `gh auth login`):

```powershell
.\scripts\upload_benchmark_release.ps1
# or: bash scripts/upload_benchmark_release.sh
```

Local staging (gitignored): `dist/release-staging/`.

## In git (clone these)

| Path | Role |
|:--|:--|
| `export/orig/` | Original ARC-SV chr20 calls, labels, time files |
| `export/v2fix11/` | ARC-SV2 golden freeze (needed by `check_golden.py`) |
| `export/ml0/dataset_v0.summary.json` | ML-0 dataset summary |
| `tests/golden/` | Golden fixtures |
| `bench/logs/ml0/*.head_tail.log` | First/last 500 lines of ML-0 host logs |

`export/ml0/merged/` and large ML tables are **not** in git; use the release
tarball below (or regenerate from `/data`).

## Release assets

```bash
gh release download benchmark-2026-11 -R dtquocbao/arcsv2 -D /tmp/arcsv2-data
cd /tmp/arcsv2-data && sha256sum -c SHA256SUMS
```

| Asset (filename in release) | Contents | ~size |
|:--|:--|--:|
| `dataset_v0.parquet` | Primary ML table | 45 MB |
| `dataset_v0.tsv.gz` | Same rows, gzip TSV | 24 MB |
| `dataset.tsv` | `ml-dataset` concat (NA dropped) | 56 MB |
| `ml0_merged.tar.gz` | `export/ml0/merged/*.merged.tab` | 19 MB |
| `s4_host.log.gz` | Full ML-0 S4 host log | 4 MB |
| `s4_snakemake.log.gz` | Full ML-0 S4 Snakemake log | 4 MB |
| `export_history.tar.gz` | Superseded `v2fix`…`v2fix10`, `results`, `labels` | 4 MB |
| `SHA256SUMS` | Checksums for the above | |

```
c40c1175073bc52d94ddc1bf6d61bc3bc2604df1b18f5dc3e3fdc4e3d12342ca  dataset.tsv
764425b9dda0389a01c0fe571ead50e9b1507eebd92c48dd59fa1b45d2530479  dataset_v0.parquet
01d355d59c9600b97d6f128ade0059c0beccb6586f7c8992b8a4b8f5d4dfd5a7  dataset_v0.tsv.gz
8d1b18e58b2ac572fbfe0dd4e18723a3e13653732d20e195a82661f688038145  export_history.tar.gz
af54ee7c46627a0e6222c78d1b46ff146fc31f218026e434ca8e6f38cea92246  ml0_merged.tar.gz
ef59a6276f8518231f8196d24a6624398768da1f469f8cf979439ec06e2a61d3  s4_host.log.gz
7d2df19fee00e978c4f20d8e5ba3f7843784a7523b70493b539570283eb7be6e  s4_snakemake.log.gz
```

## Not public

| Path | Why |
|:--|:--|
| `third_party/drbo_model/` | Dr. Bo's unpublished XGBoost / `model_replication` |
| `dist/*.tar`, `*.sif` | Container images (use GHCR) |
| BAM/CRAM/FASTA under `/data` | Too large; fetch with `scripts/fetch_data.py` |

## Reproduce comparisons

```bash
python scripts/analysis/compare_orig_vs_new_chr20.py \
  --orig-labels export/orig/labels \
  --new-labels export/v2fix11/labels
```
