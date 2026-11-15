# Caller history (v2fix → caller-v1)

Benchmark: chr20, HG01884 + NA19700, haplotype-labeled validated/precision.
Golden freeze: `export/v2fix11/`, tag `caller-v1`.

| Step | Change | Validated | Precision | Commit / notes |
|:--|:--|--:|--:|:--|
| early M2 | union merge; first dual-sample label.smk | — | — | uncommitted; 445+428 calls |
| v2fix | `merge_mode=precedence`, absorb/min_block off, Del mean insert | 157 | — | `134bc8e` |
| v2fix2 | DRP typing parity (no hanging→Inv; Dup before Ins; Ins slop 25) | 151 | — | `c6b9dce` |
| v2fix3 | signed `process_insert_len` + concordant early return + same-block split Dup | 160 | — | `693fade` |
| v2fix4 | skip-edge ref-concordance gate + `pe_support` adjacency-only | 151 | — | `6092b21` |
| v2fix5 | same-block pe junction match + `COPY_FAMILY_BAND` + `skip_edge_ref_gate` flag | — | — | `46adb80`; DUP 16/159 labeled |
| v2fix6 | split-mate join + `is_adj_satisfied` + split-path sr | 182 | ~65% | `a44611b`; DUP 31/55; build_graph slow |
| v2fix7 | indexed `overlapping_vertices`, early label length/timeout, no emit `min_affected_length` | 182 | 65% | labels match v2fix6; orig TPs 169/243 |
| v2fix8 | midpoint `len_affected` rewrite (no re-call) | 207 | 65% | orig TPs 192/243; export `export/v2fix8/` |
| v2fix9 | `skip_edge_ref_gate=False`, within-read PE, point BP cuts, flip_parity between | 198 | 56.6% | `ef23352`; regression vs v2fix8 |
| v2fix10 | gapped blocks default + 1-bp gap Dup credit (experimental) | 222 | 59% | `d4c2265`; −5.5 pp vs v2fix8 |
| **v2fix11 / caller-v1** | `simplify_blocks_diploid` with sequence + flank-preserving guard | **222** | **59.0%** | tag `caller-v1`; orig TPs 221/243; `export/v2fix11/` |

Prose detail for each step lives in git history and `archive/docs/CURSOR_TASK_*.md`.
