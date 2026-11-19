# Cursor task: HPC1-FIX — fixes from the 2026-11-18 review

Same rules as HPC1: no caller logic, numeric defaults or golden outputs change.
`check_golden.py --fast` must still PASS. Stop and report after each item.

## Must fix before H1 on Grace
1. **Retries exceed the `short` limit.** `merge` (120 min) and `dataset` (60 min) run on
   `short` (2 h) and runtime is multiplied by `attempt`. Attempt 2 of merge asks for 240 min,
   which sbatch rejects on `short`, so the retry never runs.
   Fix: lower `MERGE_RUNTIME_MIN` to 30 and `DATASET_RUNTIME_MIN` to 30, and cap every
   runtime at the partition limit (helper `_cap(minutes, partition)` with short=120,
   medium=1440, long=10080). Add a test: for attempt 1..3, every rule's runtime is at most
   its partition's limit in the slurm-grace profile.
2. **Generic profile is broken.** `workflows/profiles/slurm/config.yaml` has `attempts: 3`,
   which is not a Snakemake 8 option: `snakemake --profile workflows/profiles/slurm -n`
   exits with the usage message. Remove it (`restart-times: 2` already covers retries). Add a
   test that runs a dry-run with each of the three profiles when snakemake is installed.
3. **GHCR pull fails (403).** Either (Bao decides) make the package Public, or support a token in
   `scripts/hpc/pull_sif.sh`: if `GHCR_TOKEN` is set, export
   `SINGULARITY_DOCKER_USERNAME=dtquocbao` and `SINGULARITY_DOCKER_PASSWORD=$GHCR_TOKEN`.
   Never write the token to a file or log. After the image is reachable, pin by digest (F5.3).
4. **Commit the hand-off.** `docs/CLUSTER_TAMU.md` hand-off section is uncommitted. Fix its
   claim that tag `caller-v1.6.1` lacks `run_label.sh` (the tag is on 203ee14, which adds it).
   Replace "rsync the working tree" with `git clone https://github.com/dtquocbao/arcsv2.git
   $SCRATCH/arcsv2/repo` at the pushed commit, and record that commit in the run manifest.

## H4 must measure calling cost separately
5. H4 as written also runs `index_haps` and `label`, which are not part of a production run.
   In the H4 hand-off, run the call-only target first:
   `bash scripts/hpc/run_label.sh workflows/config/label.grace.h4.yaml /scratch/user/dtquocbao/arcsv2/data/results_h4/HG01884/merged.tab`
   Then run labeling as a second invocation. Add a small script
   `scripts/hpc/su_report.sh <snakemake log>` that collects the job ids from the log and prints
   per-rule totals from `sacct --format=JobID,JobName,Elapsed,AllocCPUS,ReqMem,MaxRSS` with
   SU = elapsed_h × max(AllocCPUS, ReqMem_GB / 7.5). Report fetch + call + merge as
   "SU per genome".

## Missing from HPC1
6. `scripts/hpc/run_label.sh` does not write the run manifest (F5.4): `git describe --always
   --dirty`, SIF sha256, snakemake version, config path, to
   `$SCRATCH/arcsv2/logs/run_manifest.<timestamp>.txt`.
7. `docs/CLUSTER_VISION.md` was not created (F6.3).
8. README line 10 still says the Rust core is "not wired in yet" (F6.4).

## Hardening (small)
9. `mmi_validate.mmi_smoke_aligns` passes when mappy/pysam are missing, when chr20 is absent, or
   on a read error. The validator must fail closed: raise `InputCheckError` in those cases.
   Also move the probe out of the chr20 centromere (GRCh38 chr20:26.44–30.04 Mb); use
   chr20:40,000,000–40,001,000 and add a second probe on chr1.
10. `tests/test_cluster_profile.py::_chrom_mb` falls back to 1.0 Mb when the fixture `.fai` is
    missing, so the chr1 runtime test always hits the floor. Use real GRCh38 lengths
    (`scripts/bench/chrom_lens.py`) in the test.
