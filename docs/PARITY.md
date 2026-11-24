# Parity: paper, original code, and ARC-SV2

Single list of deviations from the Cell 2024 Methods / original `arcsv` code,
with the controlling flag (if any) and reason. Merged from former
`DRP_PARITY.md`, `EDGE_SUPPORT.md`, and `ORIGINAL_GAP_ANALYSIS.md` conclusions.
Oracle setup notes that lived in `bench/PARITY.md` are at the end.

**Baseline:** tag `caller-v1` / `export/v2fix11/` (validated 222, precision 59.0%,
221/243 original TPs on chr20 HG01884+NA19700).

## Intentional / kept flags

| Flag / knob | Default | Paper / original | Why kept |
|:--|:--|:--|:--|
| `skip_edge_ref_gate` | `False` | Literal `add_pair` has **no** ref-concordance veto on novel skip edges | Default matches original. `True` is an ARC-SV2-only veto that reduces weak pair-only DELs |
| `max_paths` | 1000 | Paper 1000; original code 2000 | Paper mode |
| `max_pecluster_size` | 80 | Original ≈ `phys_cov × 1.5` | Fixed paper-scale default |
| `min_edge_support` | 2 | Same numeric default | Also reused as BP keep threshold (original has separate names) |
| `threads` | 1 | N/A | Genotyping pool size |
| `min_affected_length` | 50 | Evaluation ≥50 bp; live original path differs | Merge-filter / size floor |
| `filter_high_nm` | `False` | N/A (smoove-style) | C5 optional; when True drop NM > `max_nm` |
| `filter_both_mates_clipped` | `False` | N/A (smoove-style) | C5 optional; when True drop both-soft-clipped pairs |
| breakpoint merge | always precedence | Paper says “merge all overlapping”; original code uses clip/split precedence over PE | Precedence recovers true DELs lost under naive union; paper/code difference (no `union` toggle) |

## DRP typing (aligned after v2fix2–v3)

| Rule | Original | ARC-SV2 |
|:--|:--|:--|
| Secondary / duplicate / supplementary | Drop | Drop |
| MAPQ ≥ 20 both mates | Always for DRPs | Always for DRPs |
| Hanging / unmapped mate | Never enters DRP dicts | Not classified as DRP |
| Concordant opposite + insert in library range | Return None first | Early return before Del/Dup/Ins |
| Insert length | Signed `process_insert_len` | Same |
| Dup before Ins; nested FR Dup geometry | Yes | Yes |
| Ins cluster post-filter slop | 25 | `drp_insertion_cluster_slop=25` |

## Pair → edge support

| Case | Original | ARC-SV2 |
|:--|:--|:--|
| Between-mate support | Implied insert under adjacency ∈ `[q_low, q_high]` | Same |
| Same-block `B.out↔B.in` | `alt_dist < q_low` only | Same alt_dist **and** reverse upstream of forward (outward). Without outward, pair-only DUP flood |
| Novel skip `\|i−k\| > 1` | No extra gate in `add_pair` | Optional `skip_edge_ref_gate` (default off = original) |
| Within-read (split) | Always support | Always support |
| `pe_support` column | Typed DRP ∩ BPs | Counts PE with `supports_adjacency` |
| Multi-block between junction | `flip_parity` before junction | Same |

Blocks: half-open CI intervals (gapped). Contiguous point cuts were an ablation
removed after the v2fix9 INS/complex regression.

## Gap-analysis conclusions (ranked)

1. **Breakpoint merge precedence** — original keeps clip/split coordinates when
   PE intervals touch; naive union collapsed true DEL endpoints. Fixed by
   precedence merge (always on).
2. **Pair-only false calls** — not fixed by precedence alone; need scoring /
   filters (separate from missing TPs).
3. **`max_paths=1000` vs original 2000** — paper mode; can skip dense loci.
4. **Runtime/memory** — in-memory full-chr pairs vs original streaming (see
   `bench/RESULTS.md`).

## Other known differences

1. **Insert KDE** — Gaussian histogram convolution vs `sklearn` KDE `rtol=1e-6`.
2. **Path neighbor order** — sorted by node id (deterministic); igraph order may differ.
3. **Affected length** — every co-optimal traceback both string directions; original unions one direction.
4. **`bp_uncertainty`** — uncovered half-open gaps; original closed `hi - lo - 2`.
5. **Locus merge** — expand only when one starts strictly inside the previous (`pyinter.open` match).
6. **`simplify_blocks_diploid`** — ported with sequence + flank-preserving guard (caller-v1 / v2fix11).
7. **SV blacklist BED omitted (ML-0)** — the 10x Genomics GRCh38 `svblacklist.bed`
   (`cf.10xgenomics.com/.../svblacklist.bed`) returned **HTTP 403** as of 2026-10-05
   (`bench/logs/ml0/s3_filter_beds.log`). ML-0 configs (`label.ml0_s3.yaml`,
   `label.ml0_s4.yaml`) set `blacklist: ""` rather than substituting another track.
   Merge-filter still applies `simple_repeats`, `segdups`, and `gaps`. Dataset context
   column `frac_blacklist` is therefore 0 for all rows.
8. **Blacklist fallback** — `scripts/fetch_filter_beds.py` still requests the 10x
   SV blacklist first. When that download fails it uses ENCODE blacklist v2
   (Boyle-Lab/Blacklist commit `d0527cdea356bf57d30535b48783804ca8d21031`,
   `lists/hg38-blacklist.v2.bed.gz`, sha256
   `c92e763af17271446194991e71917ac220593a5a3d40a06667be24178ef08cf2`).
   ENCODE v2 is a genomic-assay blacklist (low mappability and high-signal
   regions). It is a substitute for the 10x SV blacklist, not that track.
   `blacklist.source.txt` next to the BED records which file was written.
   ML-0 `frac_blacklist` stays 0 until that BED is applied in dataset v1.
   Recomputing it needs only call coordinates (`minbp`, `maxbp`) and the BED.

## Oracle

Target (§10.2): ≥99% of original calls on the same BAM/region with the same
rearrangement and genotype, relative score difference &lt; 1e-6.

`tests/oracle/setup_oracle.sh` builds the original env;
`tests/oracle/run_oracle.py` writes the oracle table. Full tree has not yet been
run against `example/expected_output.tab`.
