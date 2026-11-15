# Fast golden windows (caller-v1 / v2fix11)

Frozen `merged.tab` slices from `export/v2fix11/` for the recall loci plus three
weak-DEL sites that must stay uncalled:

| file | sample | center | expected |
|:--|:--|--:|:--|
| `HG01884_2659058.merged.tab` | HG01884 | 2659058 | DUP ABBC |
| `HG01884_5962348.merged.tab` | HG01884 | 5962348 | DEL |
| `NA19700_52764531.merged.tab` | NA19700 | 52764531 | DUP ABBC |
| `HG01884_60792115.merged.tab` | HG01884 | 60792115 | DUP ABBC |
| `HG01884_13358681.merged.tab` | HG01884 | 13358681 | no call |
| `HG01884_23857951.merged.tab` | HG01884 | 23857951 | no call |
| `HG01884_63458085.merged.tab` | HG01884 | 63458085 | no call |

See `manifest.json` for slop (500 bp).

## How to run

Full chr20 parity (end of each Part; compares sorted `merged.tab`, `dataset.tsv`,
validated/precision summary to `export/v2fix11/`):

```bash
# inside arcsv2:dev, repo at /work, data volume at /data
python /work/scripts/check_golden.py \
  --results /data/results_v2fix11 \
  --labels /data/labels_v2fix11

# self-test against the freeze
python /work/scripts/check_golden.py \
  --results /work/export/v2fix11/results \
  --labels /work/export/v2fix11/labels
```

Fast locus parity (after every cleanup/perf step; seconds if results already
exist, or re-call ±20 kb windows first — minutes):

```bash
python /work/scripts/check_golden.py --fast --results /data/results_v2fix11
# or after a window re-call that wrote {sample}/merged.tab under OUT:
python /work/scripts/check_golden.py --fast --results "$OUT"
```

Example window re-call (one locus; repeat for each center/sample):

```bash
arcsv2 call \
  --bam /data/work/HG01884/chr20.bam \
  --ref /data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --region chr20:2639058-2679058 \
  --outdir /tmp/golden_2659058 \
  --threads 1
# then copy/symlink merged.tab into a results-shaped tree and run --fast
```

Both modes exit non-zero on any difference and print the first 20 differing rows.
