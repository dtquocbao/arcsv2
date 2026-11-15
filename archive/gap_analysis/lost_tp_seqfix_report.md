# Lost TPs: v2fix10 A_current → v2fix11_seqfix A_current

## Aggregate (pre-fix seqfix)

| arm | validated | n01 | precision |
|:--|--:|--:|--:|
| v2fix10 A_current | 38 | 90 | 42.2% |
| seqfix A_current (broken) | 35 | 90 | 38.9% |
| delta | −3 | 0 | −3.3 pp |

## Per window (pre-fix)

- **HG01884**: fix10 23/50 → seqfix 20/50 (−3)
- **NA19700**: fix10 15/40 → seqfix 15/40 (2 lost, 2 gained)

Script: `scratch/gap_analysis/compare_lost_tp_seqfix.py`

## The 3 net lost calls (HG01884; explain aggregate −3)

Also 2 NA19700 losses offset by 2 gains (see below). Full lost set = 5.

### 1. HG01884 DEL 62270773–62270915 (`chr20_62269058-62273373`)

| field | v2fix10 | seqfix (broken) |
|:--|:--|:--|
| svtype / rearr | DEL `ABCD`/`ACD` | DEL `ABC`/`AC` |
| label / reason | 1 validated | 0 low_identity |
| identity | 0.9755 (hap2) | 0.8237 (hap1) |
| label_bp_shift | 0 | 0 |
| bp_point_source | clip,clip,clip | clip,clip |

- full allele coords identical: **True**
- flank-trimmed alt identical: **False** (647 vs 300 bp)

### 2. HG01884 BND 63442864–63443024 (`chr20_63441145-63445413`)

| field | v2fix10 | seqfix (broken) |
|:--|:--|:--|
| svtype / rearr | BND complex multi-block | DUP `ABC`/`ABBC` |
| label / reason | 1 validated | 0 clipped |
| identity / mapq | 0.9433 / 60 | 0.9030 / 51 |
| label_bp_shift | 0 | 0 |

- flank-trimmed alt identical: **False** (massive outer fusion)

### 3. HG01884 DEL 63553544–63553675 (`chr20_63546658-63554646`)

| field | v2fix10 | seqfix (broken) |
|:--|:--|:--|
| svtype / rearr | DEL `ABCDEFG`/`ABCFG` | DEL `ABC`/`AC` |
| label / reason | 1 validated | 0 low_identity |
| identity | 0.9675 (hap2) | 0.8108 (hap1) |
| bp_point_source | …,mid | clip,clip (mid fused away) |

- full allele coords identical: **True**
- flank-trimmed alt identical: **False**

### Also lost (NA19700; offset by gains)

4. NA19700 DEL 63559447–63559720 — validated → low_identity (flank break)
5. NA19700 DEL 63836236–63836302 — validated → low_mapq (flank break)

### Gains in broken seqfix (NA19700)

- G1. DEL 61641710–61642014 — clipped → validated (shift=−85; flank-changing simplify helped)
- G2. DEL 62741271–62743423 — large_gap → validated

## Root cause

`simplify_blocks_diploid` preserves **full** `blocks_point` allele bases (caller RuntimeError guard), but the labeler trims each **outer** block to 150 bp (`label._prepared_blocks`). Fusing a kept internal block into an outer flank (e.g. DEL `ACD`→`AC`) moves the flank cut, so the query sequence changes and identity collapses. Mid junctions for `label_bp_shift` can also disappear.

Hypothesis mid-shift-only: **partial** — call 3 lost a mid junction, but all five losses share flank-trim breakage; calls 1–2 had shift=0 even in v2fix10.

Original also simplifies before writing blocks (`orig:sv_validate.simplify_blocks_diploid` via `sv_inference`) but uses a different validation scheme — not ARC-SV2's flank-trimmed labeler.

## Fix applied

In `caller._row_for_genotype`: after the full-allele identity check, also compare **150 bp flank-trimmed** allele coords. If flanks changed, **revert simplify** and emit pre-fusion blocks/paths. ABCBCD→ABBC still simplifies (outer A/D unchanged).

Tests: `test_simplify_blocks_diploid_preserves_allele_base_coords` (flank OK for ABBC), `test_simplify_del_acd_changes_flank_trimmed_allele`.

## Gate after fix (force re-call)

| aggregate | validated | n01 | precision | orig TP |
|:--|--:|--:|--:|--:|
| A_current (fixed) | **38** | 90 | **42.2%** | 33/48 |
| B_v2fix8 | 34 | 66 | 51.5% | 28/48 |

**GATE PASS** (`validated ≥ 38` and `precision ≥ 42%`). ABCBC spot-check: 60792115 still **DUP ABBC**. Post-fix vs v2fix10: lost=[] / gained=[] on validated matches (same 38 TPs by coords; 60792115 stays label=1 as DUP vs former BND).

## Files

- Compare script: `scratch/gap_analysis/compare_lost_tp_seqfix.py`
- Gate outputs: `scratch/gap_analysis/v2fix11_seqfix_windows/compare_summary.json`
