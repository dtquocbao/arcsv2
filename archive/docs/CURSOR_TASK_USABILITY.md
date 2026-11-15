# Cursor task: USABILITY — easy for everyone, solid for 10k-sample cohorts

## Requirements from Dr. Bo
- Speed and user-friendliness come first.
- It must run for **non-technical users**.
- It must be **"rock solid"** before the Mt. Sinai Alzheimer's cohort (~10,000 samples).
- **Outputs:** write **VCF and .tab** by default, as the original ARC-SV does.
- **Size scope:** SVs >50 bp.
- **Compute:** the lab cluster becomes available this week.

## Three kinds of users, one engine
| User | Interface |
|:--|:--|
| Biologist, no command line | Double-click launcher → browser app (no terminal) |
| Bioinformatician | One command: `arcsv2 run ...` |
| Cohort / cluster | Batch manifest + Snakemake on SLURM via Apptainer |

## Rules
- **Do not disturb the frozen ML-0 image or labeled data.**
  - Use the dev image `arcsv2:dev`, never `arcsv2:ml0`.
  - Do not rewrite `/data/labels_ml0_s3/` or `/data/ml/dataset_v0.*` except when
    an explicit ML-0 fix is requested.
  - Use a separate compose project name and low CPU priority (≤4 threads) for
    long local runs.
- **Caller `.tab` output must not change.** `check_golden.py` must PASS after
  every step; `.tab` stays byte-identical. VCF is an additional default output.
- **No speed work** in this task.
- Stop and report after each step.

---

## U1. Outputs: VCF + .tab by default
- `arcsv2 call` / `run` writes `arcsv_out.tab` (unchanged) **and** `arcsv_out.vcf.gz` + `.tbi`, using the existing `vcf.py`.
- VCF 4.2 conventions:
  - Simple SVs: `SVTYPE=DEL|DUP|INV|INS`, `END`, `SVLEN`, `CIPOS`/`CIEND` from breakpoint uncertainty. DRP-only breakpoints get an `IMPRECISE` flag, since those carry insert-size uncertainty.
  - Complex SVs: `SVTYPE=CPX`, following gnomAD-SV style. Add INFO fields `CPX_TYPE` (ARC-SV class, e.g. `ddup.del`), `REARR` (e.g. `ABCBE`), `BLOCKS`, `AFFECTED_LEN`.
  - FORMAT: `GT`, support (`SR`, `PE`) and scores (`SVSCORE`, `SVSCORE_NEXT`).
- Same calls as the `.tab` (>50 bp after merge-filter).
- Validate with `bcftools view` plus a VCF validator, and add a test.
- Document the VCF↔tab field mapping in `docs/OUTPUT.md`.

## U2. One command: `arcsv2 run`
```
arcsv2 run --input sample.cram --reference GRCh38 --out results/sample1
```
**Auto-detection**
- sample name (from `@RG SM`)
- CRAM reference, BAM/CRAM index (build one if missing)
- `chr`-prefix style, reference build (warn on mismatch)
- read length, insert size

**Defaults**
- all primary chromosomes; `--chroms chr20` for a quick test
- `--threads auto`: total threads and memory-aware per-chromosome parallelism, internal so no Snakemake knowledge is needed
- resumable: rerunning the same command continues

**`--reference GRCh38`** downloads and caches the pinned reference with a checksum and a progress bar.

**End of run**, print a plain summary:
- calls by type and size bin
- output file paths
- runtime

## U3. Input checks and plain-language errors
Check before any heavy work:
- file exists and is readable
- index present
- reference matches the CRAM
- enough disk (estimate per sample) and memory
- coverage and read length within the supported range
- no unaligned or unsorted BAMs

Every error message says **what went wrong and how to fix it**, with no stack trace by default (`--debug` shows it). Add tests for each check.

Implemented in `python/arcsv2/input_checks.py`, wired into `arcsv2 run` and `arcsv2 call` (`--debug` for tracebacks).

## U4. Robustness for big batches
- **Batch mode:** `arcsv2 batch --manifest samples.tsv --out results/`
  - Each sample is isolated: one failure never stops the others.
  - Retries with backoff for transient errors.
  - A clear `FAILED.txt` per failed sample with the reason.
  - A batch status table: done / running / failed / pending.
- **Haplotype index (`.mmi`) validation before `ml-label`** (lesson from ML-0 HG00140):
  - For each hap1/hap2 index, check **size vs its FASTA** (reject indexes that are
    far smaller than a healthy rebuild for that assembly — e.g. the corrupt
    HG00140 indexes were ~5.5 GB while peers/rebuilds were ~6.8 GB).
  - Run a **smoke query**: map a short unique GRCh38 window (or a built-in probe
    sequence) through mappy against the `.mmi`; require at least one primary hit
    on each haplotype.
  - After labeling, if **`no_hit` exceeds 2% of labeled calls** (excluding
    `no_sequence` insertions), **fail that sample** with a plain-language
    `FAILED.txt` pointing at corrupt/stale `.mmi` (rebuild with `minimap2 -d`).
  - Wire these checks into `label.smk` / `arcsv2 ml-label` so a bad index cannot
    silently poison a cohort dataset again.
- **Per-sample QC report** (HTML + JSON):
  - coverage
  - insert-size distribution
  - % discordant / clipped
  - calls by type and size
  - runtime and peak memory
  - **flags** when something looks abnormal (e.g. coverage < 15×, odd insert size)
- Atomic outputs (write to a temp file, then rename), so a crash never leaves half files.
- Log the version, git tag, parameters and input checksums in every output folder.

Implemented: `arcsv2 batch`, `python/arcsv2/batch.py`, `mmi_validate.py` (wired into `ml-label`),
`atomic_io.py`, per-sample `run_info.json` + `qc/report.{json,html}`, `batch_status.tsv`.

## U5. Cluster: Apptainer + SLURM
- `apptainer.def` (or `apptainer build` from the Docker image), with docs on running on a cluster without Docker.
- Snakemake SLURM profile (`workflows/profiles/slurm/`):
  - per-rule threads, memory from the chromosome-length formula, time limits
  - retries on failure
- `docs/CLUSTER.md`: copy-paste steps for a new lab member.
- Test with `--dry-run`, plus a 2-sample local run through the same profile logic. The real cluster test happens once access exists.

Implemented: `apptainer.def`, `scripts/build_apptainer.sh`, `workflows/profiles/slurm/`,
`docs/CLUSTER.md`, `scripts/cluster_dryrun.sh`, `tests/test_cluster_profile.py`.

## U6. Non-technical app: "ARC-SV Studio"
Extend the existing Streamlit dashboard (`python/arcsv2/dashboard/`) into a full app.

**Launchers** (double-click)
- `Start ARC-SV.bat` (Windows), `Start ARC-SV.command` (macOS), `start-arcsv.sh` (Linux).
- Each one:
  - checks that Docker Desktop is running, and shows a friendly message with an install link if not
  - pulls or updates the image
  - starts the app container with the user's chosen data folder mounted
  - opens the browser at the app

**Pages**
1. **Setup check:** Docker, CPU, RAM, free disk, reference downloaded. Each line is green/red with a fix hint.
2. **Try the example:** one button runs a bundled small dataset (a chr20 slice, a few minutes) and opens its results. It is the first thing a new user sees.
3. **New analysis:**
   - pick BAM/CRAM files from the mounted folder (file browser, multi-select)
   - choose a mode: Quick test (one chromosome) or Whole genome
   - press Start
   - before starting, show estimated time and disk
4. **My analyses:**
   - progress per sample (stage, %, ETA), cancel/resume, and the plain-language error from U3 on failure
   - read status from the U4 batch status, not by parsing logs
5. **Results:**
   - summary cards: calls by type and size
   - a filterable table (type, size, chromosome, score)
   - per-SV view, reusing the existing haplotype/block chart and showing the rearrangement in plain words (e.g. "segment B is duplicated and inserted after C")
   - QC report from U4
   - downloads: VCF, .tab, Excel/CSV
6. **Help:** short guide, glossary (SV types, genotype, score), FAQ.

**Design**
- Plain language; no jargon without a tooltip.
- Works at laptop screen sizes.
- Long jobs run in the background, so closing the browser does not stop them.

Implemented: `dashboard/studio_app.py`, `studio_ui.py`, launchers (`Start ARC-SV.bat`, `.command`, `start-arcsv.sh`), default `arcsv2 dashboard` (use `--legacy` for the old UI).

U6 non-technical pass (review fixes): pull `ghcr.io/arcsv2/arcsv2:studio` (optional `--build`), default
`Documents/ARC-SV` + folder picker, auto-open browser, release zip via `scripts/package_studio_release.sh`,
**Try example** runs `arcsv2 run` on bundled BAM, Setup GRCh38 download with progress, My analyses stage/%/ETA + Cancel/Resume,
Results summary cards + all types + QC table.

## U7. Docs and a real usability test
**README quick starts**
- (a) App: install Docker Desktop → double-click → Try example, with screenshots.
- (b) Command line: `arcsv2 run`.
- (c) Cluster: link to `CLUSTER.md`.

**`docs/TROUBLESHOOTING.md`:** the top errors from U3 and their fixes.

**Usability test script** (`docs/USABILITY_TEST.md` — template added): a person who has never used a terminal must, without help and in under 30 minutes on Windows:
1. install
2. run the example
3. run one real chr20 BAM
4. download the VCF

Record every place they got stuck. I will find the tester.

## Gate
- Golden PASS; `.tab` byte-identical.
- VCF validates.
- `arcsv2 run` on chr20 (HG01884) gives the same calls as golden.
- On Windows: launcher → example → results works from a fresh clone.
- Batch mode survives an injected failure (one bad input file) and finishes the rest.
- Tag `v2.0.0-beta`.

## Not in scope
- Speed changes (after ML-0 S4)
- Galaxy / Terra wrappers (later)
- ML model integration in the app (after ML-2)
