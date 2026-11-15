# Dashboard

A local web dashboard for ARC-SV2 results. It reads output files only, so it
runs natively on Windows, macOS and Linux; no `pysam` or `mappy` needed.

## Start

Use a dedicated environment. Don't install into Anaconda's `base`: pip can
upgrade NumPy there to 2.x while conda's older `pyarrow` stays built for NumPy
1.x. The dashboard then fails with `ImportError: numpy.core.multiarray failed to
import`.

```powershell
conda create -n arcsv2 python=3.11 -y
conda activate arcsv2
pip install -e ".[dashboard]"
```

```bash
pip install -e ".[dashboard]"
arcsv2 dashboard                      # bundled simulation demo
arcsv2 dashboard --results results/   # your own runs
```

The browser opens at <http://localhost:8501>. Inside the container, publish the
port and listen on all interfaces:

```powershell
./scripts/arcsv2.ps1 dashboard
# defaults to /data (Compose volume with Snakemake outputs)
# or: docker compose --profile dashboard up dashboard
```

## What it reads

A **results folder** holds one or more runs. A run is any folder with
`arcsv_out.tab` or `merged.tab`. Optional files next to it:

| File | Used for |
|:--|:--|
| `truth.tsv` (`chrom start end type gt [size]`) | Accuracy: correct / wrong genotype / wrong type / missed, and false calls |
| `run.json` | Run details in the sidebar (coverage, runtime, ...) |

Anywhere under the results folder:

| File | Used for |
|:--|:--|
| `dataset.tsv`, `dataset.summary.json` (from `ml-dataset`) | ML tab: label counts per class |
| `manifest.json` (from `ml-train`) | ML tab: model per class, thresholds, CV grouping, chromosome hold-out |

Truth `type` uses the caller's names: `DEL`, `DUP`, `INV`, `INS`, or the complex
type (e.g. `inv.del`). A call matches an event when it overlaps it on the same
chromosome, with ±100 bp slop.

## Tabs

- **Guide.** Step-by-step with a **Windows / macOS·Linux** switch: build → fetch →
  Snakemake via the wrapper (`arcsv2.ps1` or `arcsv2.sh`) → open results →
  record `bench/RESULTS.md`. Same commands as the README quick start.
- **Overview.** Headline numbers (true events, detected, genotype concordance,
  false calls, call time), outcome by SV type and by size, the events not called
  correctly, and a comparison across runs.
- **How it works.** The ten pipeline steps with the module that implements each,
  and a worked example. It draws the reference and the two called haplotypes block
  by block for any call, then explains the call in five steps from its own columns.
  The drawing uses `blocks`, the uncovered intervals. Labels use `blocks_point`.
- **Call explorer.** Filter by type, genotype, score and truth match. Select a row
  to see the same block drawing and explanation. Download the filtered table.
- **ML & labels.** Label counts and the trained model, once `ml-label`, `ml-dataset`
  and `ml-train` have produced files. Until then it points at the Guide tab.
- **Project status.** Milestones from `LABELS_PLAN.md` and `CROSS_PLATFORM.md`,
  detected from the repository, and the review history in `docs/reviews/`.

## Demo data

`python/arcsv2/dashboard/demo/` holds two runs of the same simulated 2 Mb
chromosome: 40 events, DEL / DUP / INV / INS / inv.del, 300 bp-3 kb, mixed
genotypes, at 30x and 15x. Reads are simulated from the two haplotypes and go
through the real SAM input path. Rebuild the demo after caller changes with:

```bash
python scripts/make_dashboard_demo.py
```

`scripts/simulate_haplotypes.py` is the simulator. Use it for any end-to-end test.

## Larger local pilot

`make pilot-large` (or `python scripts/pilot_large.py`) simulates 40 events on 2 Mb at 20x and runs `label.smk`. It is local only, inside the container, and writes `labels/dataset.tsv`. The ML tab reads that file when the dashboard is pointed at the repo: `arcsv2 dashboard --results .`. This is a dress rehearsal, not an M2 timing record.
