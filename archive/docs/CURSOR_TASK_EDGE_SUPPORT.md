# Cursor task: pair edge support + genotype parity (after 693fade)

## Where we are (export/v2fix3, chr20 HG01884 + NA19700)
| | v2fix2 | **v2fix3** | original |
|:--|--:|--:|--:|
| Validated calls | 151 | **160** | 243 |
| Original TPs recovered | 132/243 | **162/243** | — |
| Deletion precision | 36% | 33% | 67% |
| Duplications validated/labeled | 2/18 | 6/50 | 66/164 |
| ARC-SV2-only false calls | 231 | 310 | — |
| HG01884 call time | 81 min | **60 min** | 24 min |

**Discovery is now at parity** (HG01884 DRPs → passing clusters):
Del 22,217→286 (orig 24,373→298) · Dup 10,797→292 (11,047→295) · Ins 7,106→15 (9,695→19) ·
InvL 639→43 (705→47) · InvR 507→40 (578→45).
The remaining gap is in **graph edge support / genotyping**.

## Evidence
1. **Duplications are called with too many copies.** 33 of ARC-SV2's 50 labeled DUPs are 3-copy (`ABBBC`,
   `ABCCCD`); only 6 are validated. The original's are mostly 2-copy (`ABBC`). The ARC-SV2 calls carry large
   PE support the original does not have, e.g. (HG01884):
   | locus | ARC-SV2 | original |
   |:--|:--|:--|
   | 2659058–2659109 | ABBBC HET, sr 11, **pe 31** | ABBC HET, sr 11, **pe 0** |
   | 3067783–3067865 | ABCCCD HOM, sr 27, **pe 49** | ABBC HOM, sr 27, **pe 0** |
   | 306278–306466 | ABBBCD HET, sr 2, **pe 55** | ABBC HET, sr 2, pe 0 |
2. **Weak pair-only deletions.** 307 false DEL calls: median score_vs_ref **2.6**, sr = 0, **pe median 14**,
   ~92 bp, HET. True DELs: median score 50, ~300 bp. Examples with many pairs but tiny likelihood gain:
   23857951–23858079 (pe 29, score 2.18), 13358681–13358745 (pe 31, 2.25), 63458085–63458141 (pe 38, 2.28).
   Full list: `scratch/gap_analysis/v2fix3_genotype_loci.tsv`.

**Hypothesis (test it):** ARC-SV2 credits a read pair to a novel adjacency edge whenever the implied insert
under that adjacency is in the 95% range, even when the pair is **also concordant with the reference**. For
SVs ≲ 100–150 bp most concordant pairs satisfy both, so they become "PE support" for small DELs and for the
same-block tandem-dup edge. Methods p. 21 adds, for `B_i.out → B_i.in`, that the insert implied under *no*
duplication must fall **outside** the 95% range; check whether the original applies the analogous condition
to all adjacency edges, and how `sv_parse_reads.parse_reads_with_blocks` / `GenomeGraph.add_pair` choose
which edge(s) a pair supports.

## Status (2026-10-01)
- **Done:** skip-edge ref-concordance gate (`_ref_path_insert` / `_pair_supported`); `EdgeRead.supports_adjacency`;
  `pe_support` counts only adjacency-supporting PE; docs in `EDGE_SUPPORT.md`; unit tests in `test_core.py`.
- **In progress:** chr20 re-run `label.v2fix4.yaml` → `/data/results_v2fix4/` (batch `arcsv2-v2fix4-batch`).
- **Next:** after batch, export `export/v2fix4/`, refresh the table above; revisit copy-number band only if
  3-copy DUPs persist.

## Tasks
1. Write down the original's pair→edge rule (file:line) vs ARC-SV2 `caller.build_graph` / pair handling. ✅
2. For the 3 DUP loci and 3 DEL loci above, dump per locus in both tools: edges with support and which read
   names support them; paths enumerated; top genotypes with log-likelihoods (original `-v 2` or a debug hook
   in a scratch copy). Show the first place the two differ. (scratch dumps; gate confirmed)
3. Fix to match the original (unit test: concordant pair spanning a 80 bp block adds no support to the
   skip edge / same-block dup edge). ✅
4. Check the copy-number band rule (PARITY.md item 6) once edge support matches; only change it if 3-copy
   calls persist.
5. Commit, rerun chr20 (`-j 1`), relabel, export to `export/v2fix4/`, record `bench/RESULTS.md`. (batch running)

## Review of the v2fix4 changes (do after the run; don't interrupt it)
- The **pe_support column fix** (count only pairs that support the adjacency) matches the original's intent. Fine.
- The **skip-edge reference gate is not in the original.** `sv_parse_reads.GenomeGraph.add_pair:61–109`
  credits any pair whose implied insert under the between-read adjacency is in range; only the same-block
  dup edge gets the extra `alt_dist` test. Keep the gate behind a flag (`skip_edge_ref_gate`, default off
  until proven) and record it in PARITY.md as an ARC-SV2 deviation.
- The original was **not** run on the weak-DEL loci, so we still don't know why it doesn't call them.
  Run the original with `-v 2` on ±20 kb around 23857951, 13358681 and 63458085 and compare
  `merged and filtered bp` and blocks with ARC-SV2 at the same window. If the original has no breakpoint pair
  that creates the ~60–130 bp middle block, the real cause is upstream (breakpoints), and the gate only hides it.
- Report v2fix4 twice if possible: gate on vs gate off (same commit), to show what the gate costs in recall.
