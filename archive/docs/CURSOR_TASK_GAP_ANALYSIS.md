# Cursor task: why ARC-SV2 underperforms the original ARC-SV (gap analysis)

## Context (read first)

On chr20 for HG01884 and NA19700 (same BAMs, same merge-filter, same haplotype labeler), the original
ARC-SV beat ARC-SV2 on every axis:

| | Original | ARC-SV2 |
|:--|--:|--:|
| Wall time / sample | 24 / 23 min | 44 / 32 min |
| Peak RSS | 1.8 GB | ~5.1 GB |
| True calls (validated) | 243 | 106 |
| Deletion precision | 67% (172/255) | 24% (97/403) |
| Duplication precision | 40% (66/164) | 14% (5/36) |

Signals from the labeled calls:
- Calls both tools make: ~70% precision for either. ARC-SV2-only calls: 8% precision (671 calls).
- Precision by split-read support: original sr>0 85% / sr=0 19%; ARC-SV2 sr>0 67% / sr=0 11%.
- ARC-SV2: 78% of labeled calls have sr=0 (original 47%). Only 101 split-supported calls vs 239.
- Original true calls that ARC-SV2 never proposes: median PE support 0, SR support 7.5 (split-read driven).
- ARC-SV2 false calls with no original call nearby: median SR 0, PE 19, ~105 bp deletions (pair-only).
- Original log (HG01884): Del clusters 298 passing / 2747 failing; Dup 295/179; Ins 19/0; InvL 47/0; InvR 45/0.
  ARC-SV2 log: 33,288 clip/split bps, 4,413 splits, 15,919 breakpoints -> 15,859 blocks.

Data for this task (all in this repo):
- `export/results/{HG01884,NA19700}/chr20/arcsv_out.tab`, `call.log`, `graph_log.txt` (ARC-SV2)
- `export/orig/{HG01884,NA19700}/arcsv_out.tab`, `export/orig/*.call.log`, `*.time.txt` (original)
- `export/labels/*.labeled.tsv` (ARC-SV2 labels), `export/orig/labels/*.labeled.tsv` (original labels)

Sources to compare against:
- Original code: `../../src/arcsv/arcsv/` (upstream SUwonglab/arcsv @ dd6f795; only cli.py/setup.py differ). READ ONLY.
- Original defaults: `../../src/arcsv/arcsv/arcsv_call_options.py`
- Thesis chapter: `../../src/Chapter 1.pdf`
- Papers: `../../src/nihms-2024296.pdf`, `../../src/nihms-605372.pdf`
- Spec used to build ARC-SV2: `../ARC-SV2 - From-Scratch Implementation Roadmap.md`
- Known differences already listed: `bench/PARITY.md`, `docs/reviews/`

## Goal

Explain, with evidence, **why ARC-SV2 is not better than the original**, and for every difference
record **whether the thesis/papers describe it**. The output is used to ask Dr. Bo which behaviors
are undocumented implementation details of the original code.

Do NOT change engine code in this task. Analysis and small diagnostic scripts only.

## Steps

### 1. Map the original pipeline (code reading)
Trace `arcsv call` from `call_sv.py:run` through each stage and write down, per stage, the
functions, files and line numbers, every threshold/constant used (with its value from
`arcsv_call_options.py` or `constants.py`), and every filter that drops reads, clips, splits,
clusters, breakpoints, loci or calls. Stages at minimum:
1. BAM parsing and read filters (`bamparser_streaming.py`, `softclip.py`, `splitreads.py`)
2. Library / insert-size model (KDE, cutoffs)
3. Discordant pair clustering and null-distribution FDR (`pecluster.py`)
4. Soft-clip clustering and split-read parsing
5. Breakpoint merge — how clips, splits and DRP clusters combine, precedence, distances,
   minimum support (`breakpoint_merge.py`)
6. Second pass / read parsing into blocks and the graph (`sv_parse_reads.py`)
7. Locus decomposition, path enumeration, likelihood, genotyping (`sv_inference.py`, `convex_diploid.py`)
8. Insertions (`sv_inference_insertions.py`)
9. Output, classification, filters (`sv_output.py`, `sv_classify.py`, `sv_filter.py`, `filter_output.py`)

### 2. Map ARC-SV2 to the same stages
Same table for `python/arcsv2/` (samio.py, engine/evidence.py, engine/breakpoints.py,
engine/caller.py, engine/graph.py, engine/likelihood.py, merge_filter.py, config.py).

### 3. Diff stage by stage
For each stage, list every behavioral difference: missing step, different threshold, different
order, different merge/precedence rule, different filter, different default. Prioritize stages
3–5 (breakpoint discovery): the evidence says ARC-SV2 misses split-read-driven SVs and keeps too
many pair-only calls.

### 4. Measure where the counts diverge (small region)
Pick a 2–5 Mb window of chr20 that contains several original-only true calls and several
ARC-SV2-only false calls (use the labeled TSVs to choose). Run both tools on that window and
record counts at each stage: reads kept, soft clips, split reads, DRPs per type, clusters
passing/failing FDR, merged breakpoints, blocks, loci, calls. Use the original's `logging/`
output and ARC-SV2's `call.log`/`graph_log.txt`; add debug counters to a scratch copy only if needed.
Identify the **first stage** where counts diverge materially.

Then take 5 original-only true calls and 5 ARC-SV2-only false calls and trace each through both
pipelines: which breakpoints each tool proposed there, and at which step the call was created
or dropped.

### 5. Check the thesis/papers for each difference
If you cannot read PDFs directly, convert them first (e.g. `pdftotext -layout` in the container, or
`python -c "import pypdf"` extraction) into `scratch/pdf_text/`, and cite page numbers from the PDF.
For every difference from step 3, search the thesis chapter and both papers. Record:
- `documented` — quote the sentence, with document, page and section
- `partially documented` — the idea is described, but not the value or rule (quote it)
- `not documented` — not found (list the terms searched)
Also note when the ARC-SV2 roadmap says something different from both the thesis and the code.

### 6. Performance differences
Explain the runtime/memory gap by stage (original: parse BAM ~1.6 min, second pass ~6.5 min,
SV inference ~15 min, 1.8 GB; ARC-SV2: BAM load ~8.5 min, insertions ~18 min, genotype ~14 min,
~5.1 GB) with the specific design choices responsible.

## Deliverables

1. `docs/ORIGINAL_GAP_ANALYSIS.md` with:
   - Summary (half a page): the top 3–5 causes of the accuracy gap, ranked by estimated impact,
     each tied to evidence from step 4.
   - Stage map table: stage | original (file:line) | ARC-SV2 (file:line)
   - Difference table: # | stage | original behavior | ARC-SV2 behavior | impact (evidence) |
     thesis/paper status (documented / partially / not documented + quote or search terms) |
     question for Dr. Bo
   - Per-call traces for the 10 calls from step 4
   - Performance section
2. `docs/QUESTIONS_FOR_DR_BO.md`: only the `not documented` and `partially documented` items,
   each as one plain question ("The original code does X with value Y in file:line. The thesis
   does not describe this. Was this intended as part of the method, and is Y the value to use?").
   Group by stage, most important first. Max ~15 questions.
3. Any diagnostic scripts under `scripts/diagnostics/`.

## Rules
- Cite file:line for every code claim and page/section for every thesis/paper claim.
  Do not claim something is undocumented without listing what was searched.
- Separate facts (measured, read in code) from hypotheses.
- Do not modify `../../src/arcsv`. Do not modify ARC-SV2 engine code.
