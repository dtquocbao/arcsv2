# ARC-SV2: from pilot to a fully operating project

Status: plan, 2026-09-29. Phase 0 and round-13 Compose / `fetch --auto` are in the tree. First real fetch is underway: HG01884 (ACB) + NA19700 (ASW), chr20 only, `data_root: /data`, CRAMs streamed (~6.4 GB assemblies + reference). After download: Snakemake with `label.auto.yaml`, then record wall time and peak RAM in `bench/RESULTS.md`.

"Fully operating" means:
- any list of BAM/CRAM files goes in;
- filtered, genotyped SV calls (table + VCF) come out, with a per-sample QC report;
- a versioned ML model trained on real labels, with a model card;
- reproducible runs: pinned container, recorded versions;
- scheduled retraining when labels or code change.

## 1. Where to run: laptop only or not?

Resource estimates are marked **est.** and must be replaced by measured numbers from M2.

| Item | Size or cost | Note |
|:--|:--|:--|
| GRCh38 reference (NYGC analysis set) + index | ~3.5 GB | once |
| One 1KGP 30x CRAM, whole genome | ~15-20 GB | chr20-22 only, streamed by region: ~1 GB |
| One sample's two HPRC haplotype FASTAs | ~2 GB gz, ~6 GB plain | whole genome |
| minimap2 index per haplotype (default parameters) | ~7-8 GB on disk and in RAM (est.) | two loaded together today: ~16 GB RAM |
| Calling one 30x chromosome (e.g. chr20) | ~15-60 min, 2-8 GB RAM (est.) | 2 Mb simulated at 30x takes 19 s; real repeats add breakpoints |
| Calling one whole genome | ~10-40 core-hours (est.) | `--jobs` per chromosome in parallel |
| Training | minutes | small tables |

| Stage | Laptop only | HPC / cloud |
|:--|:--|:--|
| Development, tests, dashboard | ✅ best place | not needed |
| **M2 pilot** (HG002 + 5 samples, chr20-22) | ✅ feasible: ~40-120 GB disk, ≥24 GB RAM given to Docker/WSL, 1-2 days | ✅ hours |
| **M4 label set** (30-50 samples, whole genome) | ⚠️ ~1-2 TB storage, 500-2,000 core-hours: weeks of a laptop running flat out; downloads ~1 TB over home internet | ✅ days; data is already in AWS (`s3://1000genomes`, `s3://human-pangenomics`), so compute in `us-east-1` avoids downloading at all |
| **Laptop-only compromise** for labels (15-20 samples, a few chromosomes) | ✅ possible | – |
| Effect of that compromise on the model | ❌ inversion and complex classes get only tens of labels; weak per-class models; not comparable to the paper | ✅ enough labels for per-class models |
| **Production** on a cohort (tens to thousands of WGS) | ❌ not practical | ✅ required |
| Dr. Bo's clinical / psychiatric cohort data | ❌ laptops are a poor place for protected data | ✅ institutional HPC under the cohort's data-use agreement and IRB. Do not move it to a personal cloud account without approval |

**Recommendation:**
- Laptop for development and the M2 pilot.
- A cloud VM in AWS `us-east-1` (e.g. 16-32 vCPU, 64-128 GB RAM, spot instances) or the university HPC for M4 and public-data production.
- The university HPC only for Dr. Bo's cohort.
- The same container image everywhere (Docker locally; Apptainer on HPC).
- Rough cloud cost for M4: ~2,000 core-hours ≈ 60-125 h on a 16-32 vCPU VM, in the order of **$50-150 on spot** plus temporary storage (est.). Data stays in-region.

**Before deciding, check your laptop:** RAM, free disk on the drive Docker uses, CPU cores. With under 32 GB RAM, labeling needs the sequential-index change in §3.2.

## 2. Phase 0: finish round 11 / 12 (small)

1. `_stamp_clip_positions_from_softclips`: sorted clips + `bisect` — **done**.
2. Rebuild the dashboard demo after B2 gates pass (`scripts/make_dashboard_demo.py`) — **done** when `make_dashboard_demo.py` is re-run after absorb.
3. Fold discordant-pair-only breakpoints into a nearby clip- or split-supported breakpoint — **done** (clip-aware, mate-aware, ≤100 bp, re-merge; large pilot 30/32).
4. **Timed chr20 run** on one real 30x sample → `bench/RESULTS.md` (wall time, peak RSS, calls, git revision). This replaces every **est.** above.
5. **Parity run** against the original ARC-SV on its example BAM (`tests/oracle/`) → `bench/PARITY.md`.

Sequential labeling (`arcsv2 ml-label --sequential`) uses a two-pass map (hap1 then hap2, each loaded once). `scripts/fetch_data.py --auto` builds `samples.tsv` and `label.auto.yaml` from the pinned 1KGP + HPRC catalogs (both files are gitignored). Manual mode still accepts a filled `samples.tsv`. Checksums: `md5:` / `sha256:` / `none` for index sidecars.

## 3. Phase 1: data acquisition

### 3.1 `scripts/fetch_data.py` (new) and `workflows/config/samples.tsv`

- **Build the sample manifest** by intersecting the HPRC Release 2 sample index with the 1000 Genomes high-coverage (NYGC) index. One row per sample:
  - `sample`, `population`, `sex`;
  - `cram_url`, `crai_url`, `hap1_url`, `hap2_url`;
  - checksums.
  - Pin exact file URLs; never "latest".
- **Reference:** the GRCh38 analysis set used by NYGC, plus `.fai`, gaps BED, and an unreliable-region BED (segdups, centromeres, and assembly QC regions where available).
- **Downloads:**
  - Resumable (`--continue` / range requests).
  - Verify checksums.
  - Write `/data/manifest.lock.json` with the URL, size, checksum and download time of every file.
- **Region streaming:** fetch reads for the requested chromosomes straight from the CRAM URL (`samtools view -T ref URL chrN`), as `label.smk` already does. Don't download whole CRAMs for pilots.
- **Storage layout:** `/data/ref/`, `/data/samples/<id>/{reads,haps}/`, `/data/indexes/<id>/`. All paths come from config; nothing is hard-coded.
- **Cloud mode:** the same script accepts `s3://` and rewrites known public buckets to HTTPS (`1000genomes`, `human-pangenomics`), or you can pass `https://` URLs directly.

### 3.2 Lower labeling memory (needed for a laptop)

- Today `HapAligner` loads both haplotype indexes at once (~16 GB).
- Add a sequential mode: map every call's alt and ref sequences against hap1, free that index, then map against hap2, then combine per call. That's about half the RAM, with identical labels. Add a test that sequential and joint modes give the same table.
- Build indexes once and cache them under `/data/indexes/`. Delete them after labeling if disk is short (config flag).

## 4. Phase 2: labeling at scale

1. **M2 pilot:** HG002 (Q100 assembly) + 5 HPRC samples, chr20-22, inside the container.
   - Record wall time, peak RAM and disk per step in `bench/RESULTS.md`.
   - Open the labels in the dashboard (ML tab).
   - **Gates before scaling:**
     - per-class NA rate < 30%;
     - HG002 deletions ≥ 50 bp broadly consistent with the GIAB HG002 benchmark;
     - no class with fewer than ~50 labels on the pilot's chromosomes;
     - every label 0 spot-checked (10 per class) in the dashboard.
2. **M4 label set:** 30-50 samples across populations, whole genome, on cloud or HPC. Keep 5 samples fully held out (different populations) for final testing.
3. **Dataset versioning:**
   - `labels/dataset.tsv` + `dataset.summary.json` + a content hash.
   - Name each release `labels-vN`, together with `label_version` and `FEATURE_SCHEMA`.
   - Don't commit data to git. Store it in object storage or on HPC, and commit only the manifest and hash.

## 5. Phase 3: training and model selection

`arcsv2 ml-train --dataset labels/dataset.tsv --out models/<version>/` becomes a full model comparison.

- **Per class:** deletion, duplication, inversion, complex.
- **Candidates:**
  - logistic regression (current baseline, scaled);
  - random forest;
  - XGBoost;
  - LightGBM, with small grids and depth ≤ 12.
- **Validation:**
  - Nested CV, `GroupKFold` by sample: the inner loop tunes, the outer loop scores.
  - A chromosome hold-out as a leakage check.
  - Final numbers only on the 5 held-out samples.
- **Selection metric:** precision at recall ≥ 0.90 (paper rule, `choose_threshold`), PR-AUC as tie-break. Report bootstrap 95% CIs. Prefer the simpler model when CIs overlap.
- **Baseline to beat:** a `score_vs_ref` threshold alone. If no model beats it for a class, ship the threshold for that class.
- **Calibration:** isotonic or Platt on out-of-fold scores if `ml_prob` is shown as a probability; report Brier.
- **Bundle `models/<version>/`:**
  - `manifest.json`: feature schema, label version, dataset hash, git revision, container digest, per-class model and threshold, CV and held-out metrics;
  - model files (XGBoost JSON / LightGBM text / skops, never pickle);
  - scaler and PCA;
  - `MODEL_CARD.md`, generated automatically: data, intended use, metrics per class with CIs, known limits (insertions unlabeled, short HOM duplication genotypes, etc.).
- **Tests:** a fixed synthetic dataset trains to a known ranking; `ml-predict` refuses a bundle with the wrong schema or label version.

## 6. Phase 4: production

### 6.1 One command per cohort

`workflows/production.smk` (or `arcsv2 run --samples samples.tsv --model models/vN --out runs/<run_id>`), per sample:

1. Stream or read the BAM/CRAM.
2. Call per chromosome (`--jobs`).
3. Merge-filter.
4. `ml-predict` with the pinned bundle.
5. VCF (`FILTER=PASS` / `MLFAIL`).
6. Per-sample QC report.

Then per cohort: merge, allele frequencies, a cohort QC summary, and the dashboard pointed at `runs/<run_id>`.

### 6.2 Operational requirements

- **Provenance:** every run writes `run.json` with:
  - input file checksums;
  - reference;
  - config;
  - arcsv2 version and git revision;
  - container digest;
  - model bundle version and hash;
  - start and end time, and resources used.
- **Resumable and isolated:** Snakemake restarts where it stopped; one failed chromosome or sample never blocks the rest (`failures.tsv`, as `call_genome` does).
- **Execution profiles:** `profiles/local`, `profiles/slurm` (Apptainer), `profiles/aws` (Batch or a single large VM). The same workflow runs under all of them.
- **QC gates per sample:**
  - insert-size model sane;
  - coverage in range;
  - call count within the expected band;
  - fraction of `MLFAIL` within band.

  Outliers are flagged in the report, not silently kept.
- **Releases:**
  - semantic versions;
  - `CHANGELOG.md`;
  - the container image published with a version tag (e.g. GHCR) and converted for Apptainer;
  - model bundles released separately with their own versions.
- **Regression tests in CI:** the simulated pilots and a golden call table for a fixed simulation. Any change to calls or labels must be deliberate: update the golden files in the same commit and explain it in the changelog.
- **Retraining:** re-run Phase 3 when labels, features or the caller change. Compare with the current bundle on the same held-out samples. Promote only if not worse for any class.
- **Data protection** (Dr. Bo's cohort):
  - Run only where the data-use agreement allows.
  - No sample data in logs, git or the dashboard demo.
  - Pseudonymous sample IDs in shared outputs.

### 6.3 Definition of done ("production-ready")

- [ ] chr20 timing and whole-genome timing recorded; within the budget in `bench/RESULTS.md`.
- [ ] Parity report against the original ARC-SV on its example data.
- [ ] Labels v1 released (M4) with summary and hash.
- [ ] Model bundle v1 with model card, beating the `score_vs_ref` baseline on held-out samples for each class it covers.
- [ ] `production.smk` runs end to end on 3 public samples under the local and one remote profile, with identical outputs.
- [ ] CI green: core tests (Windows, macOS, Linux), container tests, simulated pilots, golden regression.
- [ ] Versioned container + Apptainer image; `CHANGELOG.md`; `docs/OPERATIONS.md` (how to run, resume, add samples, retrain, roll back a model).

## 7. Order of work for Cursor

1. Phase 0 + round-13 Compose and `fetch --auto` are in the tree.
2. **In progress:** `fetch --n-samples 2 --chroms chr20 --yes` (HG01884, NA19700). When it finishes, run Snakemake on `label.auto.yaml`, then fill `bench/RESULTS.md` (wall time, peak RAM, git revision).
3. Model comparison and bundle (§5), tested on the M2 labels.
4. `production.smk`, profiles, provenance, QC report (§6).
5. M4 on cloud or HPC → labels v1 → model v1 → production-ready checklist.

## Sources

- 1000 Genomes on AWS: https://registry.opendata.aws/1000-genomes/ ; NYGC high coverage: https://www.internationalgenome.org/announcements/3202-samples-at-high-coverage-from-NYGC/
- HPRC data (AWS `s3://human-pangenomics`): https://registry.opendata.aws/hpgp-data/ ; Release 2: https://humanpangenome.org/hprc-data-release-2/
