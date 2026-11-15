# Producing ML training labels — implementation plan

Status: plan, 2026-09-29. M0, M1, and the simulated label workflow are in the tree. M2 (HG002 + 5 HPRC samples) has not been run. Owner: Bao. Implements roadmap §4.10 and §7 using public data.

A HET call has the alt allele on one haplotype and the ref allele on the other. That is label 1, not ambiguous. Hits on either strand count, since assembly contigs are not oriented to the reference. Sequences come from `blocks_point`, with 150 bp of each outer block next to the event.

## 1. Goal

Get a labeled table of ARC-SV2 calls (1 = real, 0 = false, NA = can't be assessed) for DEL, DUP, INV and complex calls. It trains and evaluates the per-class filter.

- Labels come from **independent truth**: phased long-read assemblies of the same individuals.
- They do **not** come from our simulator. Simulated labels are only for testing the pipeline.

## 2. Data

| Need | Source | Notes |
|:--|:--|:--|
| Short reads (what ARC-SV2 calls from) | 1000 Genomes high-coverage Illumina (NYGC, ~30x, GRCh38 CRAM) | Public CRAMs; stream by region, don't download whole files |
| Truth haplotypes | HPRC Release 2: phased assemblies for 200+ individuals, many of them 1000 Genomes samples | Maternal/paternal FASTA per sample |
| Gold-standard sanity sample | HG002 Q100 diploid assembly (GIAB/T2T) | Near-complete; use as the first end-to-end sample |
| Reference | GRCh38 analysis set, the same build as the CRAMs | CRAM decoding needs the exact reference |

**Sample selection.** Build `workflows/config/label_samples.tsv` by intersecting the HPRC R2 sample list with the 1000 Genomes high-coverage index. Record for each sample the CRAM URL, index URL, hap1/hap2 FASTA URLs, sex and population. Pin exact file URLs and checksums in the manifest rather than "latest".

**Scale.** Two phases:

- **Pilot:** HG002 plus 5 HPRC samples, chr20-22.
- **Training set:** about 30-50 samples on chr1-22.

Keep 5 samples, from different populations, fully held out for the final test.

## 3. Label definition (paper §4.10, already coded in `label.py`)

For each call:

1. **Build the predicted haplotype.** Read `blocks_point`: each breakpoint is the soft-clip base when the interval has one, otherwise the midpoint, and neighboring blocks meet there. Cut the first and last blocks to 150 bp next to the event. A shorter outer block is kept and padded from the reference until that side is 150 bp. Rearrange those sequences with `sequence_from_call`. `blocks` stays the uncovered intervals, for uncertainty and the dashboard.
2. **Build the reference haplotype** for the same window, in block order. This is the negative control.
3. **Align both to hap1 and to hap2** with minimap2 through `mappy`, using default parameters and no preset, as the paper did.
4. **An alignment passes** (`is_validated`), on either strand, when it:
   - covers the full query,
   - has MAPQ ≥ 60,
   - has identity (matches / block length) > 0.90,
   - has no CIGAR insertion or deletion of 49 bp or more.
5. **Label:**

| Condition | label | reason |
|:--|:--|:--|
| Call contains `_` (novel insertion) or `=` | NA | `no_sequence` |
| Alt passes on ≥1 haplotype, and alt and ref do not both pass on the same haplotype | **1** | `validated` |
| Alt passes nowhere | **0** | first failing criterion (`clipped`, `low_mapq`, …) |
| Alt and ref both pass **on the same haplotype** (repeat or ambiguous window) | NA | `ambiguous` |
| Window overlaps an assembly gap or an unreliable-region BED | NA | `unreliable_assembly` |

6. **Extra columns, not the training label:**
   - `gt_concordant`: HOM needs alt on both haplotypes, HET on exactly one.
   - `best_hap`, `identity`, `mapq`, `max_indel`, `strand`, `ref_validated`, `label_version`.

The label uses only the paper's rule. Genotype concordance is recorded separately so genotype accuracy can be measured without changing the label.

## 4. Code changes

### 4.0 Prerequisites (must land first)

- **BAM/CRAM input through `pysam`** in `samio.py`:
  - `load_region(path, region, reference)` streams from a local or remote CRAM with its index.
  - Same `PairRecord` / clip / split output as the SAM path.
  - Test: the same reads as SAM and as BAM give identical calls.
- **An explicit block list in the call table.** `blocks` is `start-end` per block, in locus order, 0-based half-open, with breakpoint intervals left uncovered. `blocks_point` closes each interval at the clip or the midpoint. `ml-label` reads `blocks_point` and does not reconstruct blocks from `bp`. `FEATURE_SCHEMA` stays `arcsv2-f1`.

### 4.1 `python/arcsv2/label.py`

- `blocks_from_row(row) -> dict[str, tuple[int, int]]` maps block names (`A`, `B`, …, via `names.py`) to reference coordinates.
- `haplotype_sequences(row, fasta, flank=150) -> tuple[str | None, str]` returns `(alt_seq, ref_seq)` from `blocks_point`, with outer blocks cut to `flank` bp. `alt_seq` is `None` for `_` or `=`. It reuses `sequence_from_call`.
- `class HapAligner`:
  - Wraps two `mappy.Aligner` objects built from `hap1.fa` / `hap2.fa`.
  - Takes an optional prebuilt `.mmi` path for each (build once per sample and cache; each index is several GB).
  - `best_hit(seq) -> Hit | None`: the highest-MAPQ primary hit, with `q_start`, `q_end`, `q_len`, `mapq`, `mlen`, `blen` and `cigar`.
- `label_call(row, fasta, aligner, unreliable=None) -> LabelResult`: implements the table in section 3.
  - `LabelResult` gains `ref_validated`, `gt_concordant`, `best_hap`, `identity` and `max_indel`.

### 4.2 CLI: `arcsv2 ml-label` (replace the current hits-table version)

```
arcsv2 ml-label --calls calls.tab --reference GRCh38.fa \
  --hap1 SAMPLE.hap1.fa --hap2 SAMPLE.hap2.fa [--hap1-index h1.mmi --hap2-index h2.mmi] \
  [--unreliable unreliable.bed] --sample SAMPLE -o SAMPLE.labeled.tsv
```

- Output: every input column plus `sample`, `label`, `reason`, `ref_validated`, `gt_concordant`, `best_hap`, `identity`, `mapq`, `max_indel` and `label_version`.
- Keep the old hits-table mode as `ml-label-hits` for tests, or drop it.

### 4.3 CLI: `arcsv2 ml-dataset`

```
arcsv2 ml-dataset SAMPLE1.labeled.tsv SAMPLE2.labeled.tsv ... -o dataset.tsv
```

- Concatenates the tables, adds `sv_class` (`class_of_row`), `group` (= sample) and `chrom`, and drops NA labels.
- Also writes `dataset.summary.json`: per class and sample, the counts of 1/0/NA and the NA reasons.
- Refuses to mix tables with different `FEATURE_SCHEMA` or `label_version`.

### 4.4 Workflow: `workflows/label.smk`

Per sample:

1. `fetch_region`: `samtools view -u -T ref` piped to `samtools sort`, then `samtools index`. The `.bai` is an output of this rule and an input of `call`.
2. `call`: `arcsv2 call` on each chromosome.
3. `merge`: `arcsv2 merge-filter` for the single sample. This keeps the paper's size and repeat filters and the §4.9 behavior.
4. `index_haps`: build `minimap2 -d` indexes once per haplotype (cached).
5. `label`: `arcsv2 ml-label`.
6. `dataset`: `arcsv2 ml-dataset` over all samples.

Resources: indexing is roughly 30 min and a few GB of RAM per haplotype; labeling is minutes per chromosome. Run the pilot before scaling.

### 4.5 Tests: `tests/test_label.py`

Skip the test module when `mappy` isn't installed.

- **Synthetic haplotypes.** Generate a random 50 kb reference. Build hap1 with a deletion, hap2 with an inversion (`scripts/simulate.py::apply_template`), and index both.
  - The call row for the true deletion gives label 1, `gt_concordant=1` for HET.
  - The true inversion gives 1.
  - A deletion at a wrong position (shifted 500 bp) gives 0 with reason `large_gap` or `clipped`.
  - A call with `_` gives NA `no_sequence`.
  - A window duplicated in the reference (a repeat) gives NA `ambiguous`.
- **`blocks_from_row`** round-trips with the caller's `blocks` column.
- **Determinism:** the same inputs give a byte-identical labeled table.

## 5. Quality checks before training

1. **HG002 sanity check.** Label rate for deletions ≥ 50 bp should broadly match GIAB HG002 benchmark recall/precision for the same calls. Investigate large differences before going on.
2. **Simulated truth.** CI runs `label.smk` on a 200 kb chromosome with 5 events at 20x (`workflows/config/label.pilot.yaml`). `scripts/check_label_pilot.py` requires every called non-insertion event to be label 1, prints events that were not called, and requires at least 3 of the 4 non-insertion events to be called. Insertions are NA. `make pilot-large` is the local 2 Mb, 40-event run and is not a training set. The later bar, once M2 exists, is still that true calls are label 1 and injected false calls are 0.
3. **Per-class NA rate.** If more than 30% of a class is NA (usually `ambiguous` in repeats), report it; don't train on the biased remainder silently.
4. **Class balance** per SV type, recorded in `dataset.summary.json`.

## 6. Training and evaluation (after labels exist)

Fix `train.py` first; see round 5 finding 1:
- a scaler inside the model, stored in the manifest;
- `GroupKFold` by sample;
- the same PCA procedure for out-of-fold and final models;
- no silent zeros or silent fallback.

Then compare per class:

- **Candidates:** logistic regression (baseline), random forest, XGBoost and LightGBM, with small grids and depth capped at about 12.
- **CV:** nested, grouped by sample. The inner loop tunes; the outer loop measures.
- **Primary metric:** precision at recall ≥ 0.90 on out-of-fold scores (the paper rule; `choose_threshold` already implements it).
- **Secondary metrics:** PR-AUC, and Brier score if `ml_prob` is shown as a probability.
- **Baseline to beat:** a threshold on `score_vs_ref` alone.
- **Report:** bootstrap 95% CIs per class. When two candidates' CIs overlap, pick the simpler one.
- **Final check** on the 5 held-out samples, and a chromosome hold-out as an extra leakage check.
- **Bundle:** saved models (XGBoost JSON / skops / coefficients), scaler, PCA, thresholds, metrics, dataset checksum and the sample list, all in the manifest.

## 7. Milestones

| # | Deliverable | Done when |
|:--|:--|:--|
| M0 | `pysam` BAM/CRAM input + `blocks` and `blocks_point` | In the tree. SAM and BAM give identical calls |
| M1 | `label.py` + `ml-label` + `tests/test_label.py` | In the tree. Synthetic tests pass in the container; the 200 kb workflow labels called non-insertion events 1 |
| M2 | Pilot: HPRC∩1KGP samples, chr20–22 (HG002 later via `--benchmark`) | **Fetch started** 2026-09-29: HG01884 + NA19700, chr20. Call/label + SUMMARY still required |
| M3 | `train.py` fixes | Grouped CV, scaler, consistent PCA; tests pass |
| M4 | Full label set (30-50 samples, chr1-22) | Dataset + summary, NA rates reviewed |
| M5 | Model comparison + final bundle | Per-class report with CIs; held-out test recorded in `bench/` |

## 8. Risks

- **Assembly errors** make false negatives in the labels, especially in segdups and centromeres. Mitigation: the unreliable-region BED and the NA `ambiguous` rule.
- **Reference build mismatch** between CRAM, reference FASTA and assemblies. Mitigation: pin GRCh38 and check contig names and lengths at start-up.
- **Leakage across samples:** common SVs appear in many samples. Mitigation: grouped CV plus a chromosome hold-out.
- **Compute and storage:** use region streaming and cache the minimap2 indexes. Budget per sample before scaling past the pilot.

## Sources

- HPRC Release 2: https://humanpangenome.org/hprc-data-release-2/ and https://www.internationalgenome.org/data-portal/data-collections/hprc2.html
- HG002 Q100: https://genomeinformatics.github.io/Q100/ and https://www.nist.gov/programs-projects/genome-bottle
