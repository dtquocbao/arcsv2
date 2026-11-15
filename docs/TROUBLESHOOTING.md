# Troubleshooting (U3 / Studio)

Plain fixes for the most common errors. For cluster runs see `CLUSTER.md`.

## Docker / Studio

| Symptom | Fix |
|:--|:--|
| “Docker is not running” | Start Docker Desktop; wait until it shows Running. [Install link](https://docs.docker.com/desktop/) |
| Permission denied on `/data` | Run `docker compose --profile init run --rm init-data` or use `scripts/arcsv2.ps1 dashboard` (runs init automatically) |
| Blank Studio / no runs | Set `--results` to the folder that contains `arcsv_out.tab` or `merged.tab` (default `/data` in Docker) |

## Input errors (`arcsv2 run` / `arcsv2 call`)

| Message | Fix |
|:--|:--|
| Input file does not exist | Check `--input` path; use absolute paths in batch manifests |
| Index missing | `samtools index your.cram` or `your.bam` |
| CRAM / reference mismatch | Use the same GRCh38 FASTA the reads were aligned to (`-R` / `--reference`) |
| Not coordinate-sorted | `samtools sort -o sorted.bam in.bam && samtools index sorted.bam` |
| Coverage below 10× | Use deeper data or `--chroms chr20` for a pilot |
| Not enough disk / RAM | Free space, smaller `--chroms`, or a larger machine |

## Labeling / `.mmi`

| Message | Fix |
|:--|:--|
| Index too small | Rebuild: `minimap2 -d sample.hap1.mmi sample.hap1.fa.gz` (both haps) |
| Smoke alignment failed | Same — stale or corrupt `.mmi` after a FASTA update |
| `no_hit` above 2% | Rebuild indexes and re-run `ml-label`; see ML-0 HG00140 notes in `bench/RESULTS.md` |

## Batch mode

| Symptom | Fix |
|:--|:--|
| One sample failed | Read `<out>/<sample>/FAILED.txt`; other samples continue |
| Status stuck on running | Check `logs/slurm/` or `studio_batch.log`; re-run `arcsv2 batch` (resumes) |

Add `--debug` to CLI commands to see a full Python traceback.
